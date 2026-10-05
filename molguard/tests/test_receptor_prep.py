import subprocess
from pathlib import Path
import pytest

from molguard.io import receptor_prep
from molguard.io.receptor_prep import (
    _parse_excess_bond_residues,
    prepare_receptor_pdbqt,
    prepare_receptors_in_directory,
    run_prepare_command,
    sanitize_pdb_for_meeko,
)


def _atom(serial=1, atom_type="N", name="N", x=0.0):
    return (f"ATOM  {serial:5d} {name:>4} ALA A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00    {0.1:6.3f} {atom_type:>2}")


def _ready_receptor():
    return _atom() + "\n" + _atom(2, "HD", "H", 1.0) + "\n"


def test_sanitize_pdb_for_meeko_assigns_blank_chain_segments(tmp_path: Path) -> None:
    input_pdb = tmp_path / "receptor.pdb"
    output_pdb = tmp_path / "receptor.sanitized.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA     1       0.000   0.000   0.000",
                "ATOM      2  CA  ALA     1       1.000   0.000   0.000",
                "TER",
                "ATOM      3  N   LYS   443       1.200   0.000   0.000",
                "ATOM      4  CA  LYS   443       2.200   0.000   0.000",
                "TER",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    sanitize_pdb_for_meeko(input_pdb, output_pdb)
    lines = output_pdb.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    ter_lines = [line for line in lines if line.startswith("TER")]

    assert atom_lines[0][21] == "A"
    assert atom_lines[1][21] == "A"
    assert ter_lines[0][21] == "A"
    assert atom_lines[2][21] == "B"
    assert atom_lines[3][21] == "B"
    assert ter_lines[1][21] == "B"


def test_sanitize_pdb_for_meeko_resolves_altlocs_deterministically(tmp_path: Path) -> None:
    def atom(serial: int, name: str, altloc: str, resname: str, x: float) -> str:
        return (
            f"ATOM  {serial:5d} {name:<4}{altloc}{resname:>3} E{8:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}"
        )

    input_pdb = tmp_path / "altloc_receptor.pdb"
    output_pdb = tmp_path / "altloc_receptor.sanitized.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                atom(1, "N", " ", "ALA", 0.0),
                atom(2, "CA", "A", "ALA", 1.0),
                atom(3, "CA", "B", "ALA", 9.0),
                "TER",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    sanitize_pdb_for_meeko(input_pdb, output_pdb)
    atom_lines = [
        line for line in output_pdb.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]

    assert len(atom_lines) == 2
    assert [line[12:16].strip() for line in atom_lines] == ["N", "CA"]
    assert all(line[16] == " " for line in atom_lines)
    assert float(atom_lines[1][30:38]) == 1.0


def test_sanitize_pdb_for_meeko_converts_selenomethionine_atoms(tmp_path: Path) -> None:
    def atom(serial: int, name: str, resname: str, resseq: int, element: str) -> str:
        return (
            f"ATOM  {serial:5d} {name:<4} {resname:>3} A{resseq:4d}    "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2}"
        )

    input_pdb = tmp_path / "selenium_receptor.pdb"
    output_pdb = tmp_path / "selenium_receptor.sanitized.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                atom(1, "SE", "MET", 1, "SE"),
                atom(2, "SE", "MSE", 2, "SE"),
                "TER",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    sanitize_pdb_for_meeko(input_pdb, output_pdb)
    atom_lines = [
        line for line in output_pdb.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]

    assert len(atom_lines) == 2
    assert [line[17:20].strip() for line in atom_lines] == ["MET", "MET"]
    assert [line[12:16].strip() for line in atom_lines] == ["SD", "SD"]
    assert [line[76:78].strip() for line in atom_lines] == ["S", "S"]


def test_parse_excess_bond_residues_from_meeko_padding_error() -> None:
    stderr = (
        "matched with excess inter-residue bond(s): A:23\n"
        "matched with excess inter-residue bond(s): d:443\n"
        "RuntimeError: Expected 2 paddings for (A:23, d:443) "
        "with bonds [(8, 18)], but got 0"
    )

    assert _parse_excess_bond_residues(stderr) == ["A:23", "d:443"]


def test_meeko_excess_bond_retry_allows_single_residue(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_pdb = tmp_path / "receptor.sanitized.pdb"
    output_pdbqt = tmp_path / "receptor.pdbqt"
    work_dir = tmp_path / "work"
    input_pdb.write_text("ATOM\n", encoding="ascii")
    calls: list[list[str]] = []
    run_cwds: list[str | None] = []

    def fake_run(argv, **kwargs):
        calls.append(list(argv))
        run_cwds.append(kwargs.get("cwd"))
        if len(calls) == 1:
            raise subprocess.CalledProcessError(
                1,
                argv,
                stderr="matched with excess inter-residue bond(s): L:688",
            )
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(receptor_prep.subprocess, "run", fake_run)

    run_prepare_command(
        command_template="mk_prepare_receptor.py -i {input} -p {output} --allow_bad_res",
        input_path=input_pdb,
        output_path=output_pdbqt,
        seed=42,
        cwd=work_dir,
    )

    assert len(calls) == 2
    assert calls[1][-2:] == ["--delete_residues", "L:688"]
    assert run_cwds == [str(work_dir.resolve()), str(work_dir.resolve())]


def test_pdbqt_input_is_canonicalized_even_with_prepare_command(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_pdbqt = tmp_path / "ready_receptor.pdbqt"
    output_pdbqt = tmp_path / "canonical_ready_receptor.pdbqt"
    input_pdbqt.write_text(_ready_receptor(), encoding="ascii")

    def fake_canonicalize(infile: Path, outfile: Path, *, timestamp: str) -> str:
        assert infile.read_text() == input_pdbqt.read_text()
        assert timestamp == "TEST"
        outfile.write_text(_ready_receptor(), encoding="ascii")
        return "digest-pdbqt"

    def fail_prepare_command(**kwargs) -> None:
        raise AssertionError("PDBQT receptors must not run external preparation")

    monkeypatch.setattr(receptor_prep, "canonicalize_receptor", fake_canonicalize)
    monkeypatch.setattr(receptor_prep, "run_prepare_command", fail_prepare_command)

    digest = prepare_receptor_pdbqt(
        input_path=input_pdbqt,
        output_path=output_pdbqt,
        prepare_command="external-prep {input} {output}",
        seed=123,
        timestamp="TEST",
    )

    assert digest == "digest-pdbqt"
    assert output_pdbqt.read_text(encoding="ascii") == _ready_receptor()


def test_pdb_receptor_prep_uses_output_local_scratch(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_pdb = tmp_path / "raw_receptor.pdb"
    output_pdbqt = tmp_path / "prepared" / "raw_receptor.pdbqt"
    input_pdb.write_text("ATOM\n", encoding="ascii")
    prep_cwds: list[Path | None] = []

    def fake_sanitize(infile: Path, outfile: Path) -> Path:
        assert infile == input_pdb
        outfile.write_text("SANITIZED\n", encoding="ascii")
        return outfile

    def fake_prepare_command(
        *,
        command_template: str,
        input_path: Path,
        output_path: Path,
        seed: int,
        cwd: Path | None,
    ) -> None:
        prep_cwds.append(cwd)
        assert input_path.parent == cwd
        output_path.write_text(_ready_receptor(), encoding="ascii")

    def fake_canonicalize(infile: Path, outfile: Path, *, timestamp: str) -> str:
        assert infile.read_text(encoding="ascii") == _ready_receptor()
        outfile.write_text(_ready_receptor(), encoding="ascii")
        return "digest-pdb"

    monkeypatch.setattr(receptor_prep, "sanitize_pdb_for_meeko", fake_sanitize)
    monkeypatch.setattr(receptor_prep, "run_prepare_command", fake_prepare_command)
    monkeypatch.setattr(receptor_prep, "canonicalize_receptor", fake_canonicalize)

    digest = prepare_receptor_pdbqt(
        input_path=input_pdb,
        output_path=output_pdbqt,
        prepare_command="external-prep {input} {output}",
        seed=123,
        timestamp="TEST",
    )

    assert digest == "digest-pdb"
    assert output_pdbqt.read_text(encoding="ascii") == _ready_receptor()
    assert len(prep_cwds) == 1
    assert prep_cwds[0] is not None
    assert prep_cwds[0].parent == output_pdbqt.parent
    assert not prep_cwds[0].exists()


def test_prepare_receptors_in_directory_uses_input_stems_for_mixed_extensions(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_pdb = tmp_path / "source_receptor.pdb"
    ready_pdbqt = tmp_path / "prepared_receptor.pdbqt"
    ignored = tmp_path / "notes.txt"
    raw_pdb.write_text("ATOM\n", encoding="ascii")
    ready_pdbqt.write_text("PDBQT\n", encoding="ascii")
    ignored.write_text("ignore me\n", encoding="ascii")

    calls: list[tuple[str, str, str | None]] = []

    def fake_prepare_receptor_pdbqt(
        *,
        input_path: Path,
        output_path: Path,
        prepare_command: str | None,
        seed: int,
        timestamp: str,
    ) -> str:
        calls.append((input_path.name, output_path.name, prepare_command))
        output_path.write_text(f"{input_path.name} -> {output_path.name}\n", encoding="ascii")
        return f"digest-{input_path.stem}"

    monkeypatch.setattr(
        receptor_prep,
        "prepare_receptor_pdbqt",
        fake_prepare_receptor_pdbqt,
    )

    records = prepare_receptors_in_directory(
        tmp_path,
        prepare_command="external-prep {input} {output}",
        seed=7,
        timestamp="TEST",
    )

    by_source = {record.source_path.name: record for record in records}
    assert sorted(by_source) == ["prepared_receptor.pdbqt", "source_receptor.pdb"]
    assert by_source["prepared_receptor.pdbqt"].output_path.name == "prepared_receptor.pdbqt"
    assert by_source["prepared_receptor.pdbqt"].action == "validated/canonicalized"
    assert by_source["source_receptor.pdb"].output_path.name == "source_receptor.pdbqt"
    assert by_source["source_receptor.pdb"].action == "prepared"
    assert [call[0] for call in calls] == ["prepared_receptor.pdbqt", "source_receptor.pdb"]


def test_hydrogen_recovery_preserves_original_and_records_diagnostics(tmp_path, monkeypatch):
    import json
    from molguard.io import receptor_checks
    source = tmp_path / "receptor with spaces.pdbqt"
    original = _atom() + "\n"
    source.write_text(original)
    monkeypatch.setattr(receptor_checks.shutil, "which", lambda _: "/tools/obabel")

    def convert(argv, **kwargs):
        pdb = Path(argv[argv.index("-ipdb") + 1])
        assert pdb.read_text().splitlines()[0][76:78].strip() == "N"
        Path(argv[argv.index("-O") + 1]).write_text(_ready_receptor())
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="converter warning")

    monkeypatch.setattr(receptor_checks.subprocess, "run", convert)
    digest = prepare_receptor_pdbqt(input_path=source, output_path=source,
                                   prepare_command=None, seed=42, timestamp="TEST")
    report = json.loads(source.with_suffix(".prep.json").read_text())
    assert report["output_sha256"] == digest
    assert report["diagnostics"] == "converter warning"
    assert Path(report["backup"]).read_text() == original
    assert report["donor_hydrogens"] == 1
    assert report["heavy_atoms_preserved"] == 1
    # A second preparation never re-runs chemistry or overwrites provenance.
    monkeypatch.setattr(receptor_checks.subprocess, "run", lambda *a, **k: pytest.fail("reconverted"))
    assert prepare_receptor_pdbqt(input_path=source, output_path=source,
                                 prepare_command=None, seed=42, timestamp="TEST") == digest


@pytest.mark.parametrize("failure", ["missing_tool", "lost_atom", "moved_atom", "no_hydrogens", "timeout"])
def test_failed_hydrogen_recovery_never_overwrites_input(tmp_path, monkeypatch, failure):
    from molguard.io import receptor_checks
    source = tmp_path / "receptor.pdbqt"
    original = _atom() + "\n"
    source.write_text(original)
    monkeypatch.setattr(receptor_checks.shutil, "which",
                        lambda _: None if failure == "missing_tool" else "/tools/obabel")

    def convert(argv, **kwargs):
        if failure == "timeout":
            raise subprocess.TimeoutExpired(argv, 120)
        content = {
            "lost_atom": _atom(2, "HD", "H", 1.0) + "\n",
            "moved_atom": _atom(x=8.0) + "\n" + _atom(2, "HD", "H", 1.0) + "\n",
            "no_hydrogens": original,
        }[failure]
        Path(argv[argv.index("-O") + 1]).write_text(content)
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(receptor_checks.subprocess, "run", convert)
    with pytest.raises((ValueError, RuntimeError)):
        prepare_receptor_pdbqt(input_path=source, output_path=source, prepare_command=None, seed=42)
    assert source.read_text() == original
    assert not source.with_suffix(".prep.json").exists()


@pytest.mark.parametrize("content, message", [
    ("REMARK empty\n", "no ATOM"),
    ("MODEL        1\n" + _ready_receptor() + "ENDMDL\nMODEL        2\n" + _ready_receptor(), "multiple MODEL"),
    ("ROOT\n" + _ready_receptor() + "ENDROOT\nTORSDOF 0\n", "torsion"),
    (_atom(atom_type="XX") + "\n", "unsupported receptor atom type"),
    (_atom().replace("   0.000", "     nan", 1) + "\n", "non-finite"),
])
def test_unsafe_receptors_are_rejected_before_overwrite(tmp_path, content, message):
    source = tmp_path / "unsafe.pdbqt"
    source.write_text(content)
    with pytest.raises((ValueError, RuntimeError), match=message):
        prepare_receptor_pdbqt(input_path=source, output_path=source, prepare_command=None, seed=42)
    assert source.read_text() == content


def test_gzipped_pdb_keeps_format_suffix(tmp_path):
    import gzip
    source = tmp_path / "receptor.pdb.gz"
    with gzip.open(source, "wt") as handle:
        handle.write("ATOM\n")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    unpacked = receptor_prep._materialize_if_gz(source, scratch)
    assert unpacked.name == "receptor.pdb"
    assert unpacked.read_text() == "ATOM\n"


def test_command_template_keeps_unquoted_space_paths_in_one_argument():
    argv = receptor_prep.format_command_template(
        "prep -i {input} -o '{output}' --seed {seed}",
        Path("/tmp/input with 'quotes'.pdb"), Path("/tmp/output with spaces.pdbqt"), 7)
    assert argv == ["prep", "-i", "/tmp/input with 'quotes'.pdb", "-o",
                    "/tmp/output with spaces.pdbqt", "--seed", "7"]


def test_directory_reports_all_failures_instead_of_silently_continuing(tmp_path):
    for name in ("a", "b"):
        (tmp_path / f"{name}.pdbqt").write_text("REMARK empty\n")
    with pytest.raises(RuntimeError) as exc:
        prepare_receptors_in_directory(tmp_path)
    assert "a.pdbqt" in str(exc.value) and "b.pdbqt" in str(exc.value)


def test_prepare_repairs_numeric_formatting_without_changing_values(tmp_path):
    source = tmp_path / "numeric.pdbqt"
    first = _atom()
    source.write_text(first[:30] + " 1.00e+1" + first[38:] + "\n" + _atom(2, "HD", "H", 11.0) + "\n")
    prepare_receptor_pdbqt(input_path=source, output_path=source, prepare_command=None, seed=42)
    from molguard.io.pdbqt import pdbqt_check
    assert pdbqt_check(source).ok
    from molguard.io.receptor_checks import receptor_atom_lines
    assert sorted(float(line[30:38]) for line in receptor_atom_lines(source)) == [10.0, 11.0]


@pytest.mark.parametrize("force, expected", [(False, "target.pdbqt"), (True, "target.pdb")])
def test_sibling_formats_are_processed_once(tmp_path, monkeypatch, force, expected):
    for name in ("target.pdb", "target.mol2", "target.pdbqt", "target.pdbqt.gz"):
        (tmp_path / name).write_text("placeholder")
    calls = []

    def prepare(**kwargs):
        calls.append(kwargs["input_path"].name)
        return "digest"

    monkeypatch.setattr(receptor_prep, "prepare_receptor_pdbqt", prepare)
    records = prepare_receptors_in_directory(tmp_path, force=force)
    assert calls == [expected]
    assert len(records) == 1


def test_uppercase_receptor_extension_is_visible_to_docking(tmp_path):
    (tmp_path / "target.PDBQT").write_text(_ready_receptor())
    records = prepare_receptors_in_directory(tmp_path)
    assert records[0].output_path == tmp_path / "target.pdbqt"
    assert len(list(tmp_path.glob("*.pdbqt"))) == 1


@pytest.mark.parametrize("field", [(30, 38), (70, 76)])
def test_missing_required_numbers_are_not_invented(tmp_path, field):
    start, end = field
    line = _atom()
    content = line[:start] + " " * (end - start) + line[end:] + "\n"
    source = tmp_path / "missing.pdbqt"
    source.write_text(content)
    with pytest.raises(ValueError, match="missing coordinate or partial charge"):
        prepare_receptor_pdbqt(input_path=source, output_path=source, prepare_command=None, seed=42)
    assert source.read_text() == content


def test_blank_metadata_fields_get_defaults(tmp_path):
    line = _atom()
    source = tmp_path / "metadata.pdbqt"
    source.write_text(line[:54] + " " * 12 + line[66:] + "\n" + _atom(2, "HD", "H", 1.0) + "\n")
    prepare_receptor_pdbqt(input_path=source, output_path=source, prepare_command=None, seed=42)
    from molguard.io.receptor_checks import receptor_atom_lines
    assert all(line[54:66] == "  1.00  0.00" for line in receptor_atom_lines(source))


def test_mol2_uses_a_backend_that_reads_mol2(monkeypatch):
    monkeypatch.setattr(receptor_prep.shutil, "which", lambda name: f"/tools/{name}")
    assert receptor_prep._auto_receptor_prepare_command(".mol2").startswith("obabel ")


def test_gzipped_ready_pdbqt_can_be_prepared_repeatedly(tmp_path):
    import gzip
    source = tmp_path / "target.pdbqt.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(_ready_receptor())
    first = prepare_receptors_in_directory(tmp_path)
    second = prepare_receptors_in_directory(tmp_path)
    assert len(first) == len(second) == 1
    assert first[0].digest == second[0].digest
    assert second[0].source_path == tmp_path / "target.pdbqt"
