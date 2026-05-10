import subprocess
from pathlib import Path

from molguard.io import receptor_prep
from molguard.io.receptor_prep import (
    _parse_excess_bond_residues,
    prepare_receptor_pdbqt,
    prepare_receptors_in_directory,
    run_prepare_command,
    sanitize_pdb_for_meeko,
)


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
    input_pdbqt.write_text("RECEPTOR\n", encoding="ascii")

    def fake_canonicalize(infile: Path, outfile: Path, *, timestamp: str) -> str:
        assert infile == input_pdbqt
        assert timestamp == "TEST"
        outfile.write_text("CANONICAL\n", encoding="ascii")
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
    assert output_pdbqt.read_text(encoding="ascii") == "CANONICAL\n"


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
        output_path.write_text("PDBQT\n", encoding="ascii")

    def fake_canonicalize(infile: Path, outfile: Path, *, timestamp: str) -> str:
        assert infile.read_text(encoding="ascii") == "PDBQT\n"
        outfile.write_text("CANONICAL\n", encoding="ascii")
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
    assert output_pdbqt.read_text(encoding="ascii") == "CANONICAL\n"
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
    assert by_source["prepared_receptor.pdbqt"].action == "canonicalized"
    assert by_source["source_receptor.pdb"].output_path.name == "source_receptor.pdbqt"
    assert by_source["source_receptor.pdb"].action == "prepared"
    assert [call[0] for call in calls] == ["prepared_receptor.pdbqt", "source_receptor.pdb"]
