#!/usr/bin/env python3
"""Deterministic prep helpers for DUD-E style docking benchmarks.

Receptor preparation is owned by ``molguard.io.receptor_prep`` and re-exported
from this module for old benchmark imports. Ligand preparation remains here.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from molguard.io.pdbqt import canonicalize_receptor, pdbqt_normalize


SUPPORTED_LIGAND_SUFFIXES = {".pdbqt", ".mol2"}
SUPPORTED_RECEPTOR_SUFFIXES = {".pdbqt", ".pdb", ".mol2"}


def _effective_suffix(path: Path) -> str:
    """Return the payload suffix, treating .gz as transparent compression."""
    suffixes = [part.lower() for part in path.suffixes]
    if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2:
        return suffixes[-2]
    return path.suffix.lower()


def _stem_without_gz(path: Path) -> str:
    """Return a stable stem that ignores a trailing .gz suffix."""
    if path.suffix.lower() == ".gz":
        return Path(path.stem).stem
    return path.stem


def _materialize_if_gz(input_path: Path, tmp_dir: Path) -> Path:
    """Unpack a gzipped file into tmp_dir, otherwise return the original path."""
    if input_path.suffix.lower() != ".gz":
        return input_path
    unpacked = tmp_dir / _stem_without_gz(input_path)
    with gzip.open(input_path, "rb") as src, unpacked.open("wb") as dst:
        dst.write(src.read())
    return unpacked


def _infer_element_from_atom_name(atom_name_field: str) -> str:
    """Infer a plausible element symbol from the PDB atom-name field."""
    field = atom_name_field[:4].ljust(4)
    token = "".join(ch for ch in field.strip() if ch.isalpha())
    if not token:
        return ""

    # PDB alignment rule: one-letter elements are usually right-justified in
    # the 4-char atom-name field (first char blank), while two-letter elements
    # typically start in column 13. This prevents " CA " (alpha carbon) from
    # being mistaken for calcium.
    if field[0].isspace():
        return token[0].upper()

    token = token.upper()
    if len(token) >= 2 and token[:2] in {
        "CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CU", "NI", "CO", "CD",
    }:
        return token[:2].title()
    if token[0] == "H":
        return "H"
    return token[0]


def _infer_element_for_record(atom_name_field: str, residue_name: str) -> str:
    """Infer an element symbol, preserving common monoatomic ions."""
    residue_name = residue_name.strip().upper()
    if residue_name in _ION_RESIDUES:
        token = "".join(ch for ch in atom_name_field[:4].strip() if ch.isalpha()).upper()
        if token == residue_name:
            return residue_name.title()
    return _infer_element_from_atom_name(atom_name_field)


# Standard amino acid and nucleotide residue names that Meeko/ProDy recognise.
# Any ATOM record with a residue name NOT in this set is silently dropped to
# avoid "unknown residue" errors in mk_prepare_receptor.
_STANDARD_RESIDUES: set[str] = {
    # 20 standard amino acids (+ common protonation variants)
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # common caps / terminals
    "ACE", "NME", "NH2",
    # nucleotides
    "DA", "DT", "DC", "DG", "A", "U", "C", "G",
    # water / ions (often present in DUD-E PDBs)
    "HOH", "WAT", "TIP", "TIP3",
}

# Common DUD-E receptor aliases that are chemically close enough to remap onto
# standard residue templates instead of deleting them outright before Meeko.
_RESIDUE_RENAMES: dict[str, str] = {
    "ASH": "ASP",
    "ASQ": "ASP",
    "DIC": "ASP",
    "DID": "ASP",
    "GLO": "GLU",
    "GLZ": "GLY",
    "HIZ": "HIS",
    "HIY": "HIS",
    "LEV": "LEU",
    "MEU": "MET",
}

_SOLVENT_RESIDUES: set[str] = {"HOH", "WAT", "TIP", "TIP3", "WAM"}

_ION_RESIDUES: set[str] = {
    "CA",
    "CD",
    "CO",
    "CU",
    "FE",
    "K",
    "MG",
    "MN",
    "NA",
    "NI",
    "ZN",
}

_SYNTHETIC_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _synthetic_chain_id(segment_index: int) -> str:
    """Return a deterministic one-character chain ID for blank-chain PDB segments."""
    return _SYNTHETIC_CHAIN_IDS[segment_index % len(_SYNTHETIC_CHAIN_IDS)]


def _print_receptor_integrity_warning(title: str, details: Iterable[str]) -> None:
    """Print a loud warning when receptor prep removes receptor atoms/residues."""
    border = "!" * 78
    print(f"\n{border}")
    print(f"!!! RECEPTOR PREP WARNING: {title}")
    for detail in details:
        print(f"!!! {detail}")
    print("!!! This changes the receptor model and must be reported with results.")
    print(f"{border}\n")


def sanitize_pdb_for_meeko(input_path: Path, output_path: Path) -> Path:
    """
    Rewrite short PDB records into a PDB variant Meeko/ProDy can parse.

    DUD-E receptor PDBs commonly omit occupancy, B-factor, and element fields.
    This fills those columns deterministically while preserving coordinates.

    Non-standard residues (e.g. co-crystal ligand fragments left as ATOM
    records) are silently dropped so that Meeko does not choke on unknown
    residue templates.
    """
    lines_out: list[str] = []
    dropped_residues: set[str] = set()
    dropped_solvents: set[str] = set()
    preserved_ions: set[str] = set()
    remapped_residues: dict[tuple[str, str], int] = {}
    assigned_blank_chain_segments: set[str] = set()
    blank_segment_index = 0
    blank_segment_has_atoms = False

    for raw in input_path.read_text(encoding="utf-8").splitlines():
        if raw.startswith("TER"):
            line = raw.rstrip("\n").ljust(80)
            if blank_segment_has_atoms and not line[21].strip():
                chain_id = _synthetic_chain_id(blank_segment_index)
                line = f"{line[:21]}{chain_id}{line[22:]}"
            lines_out.append(line.rstrip())
            if blank_segment_has_atoms:
                blank_segment_index += 1
                blank_segment_has_atoms = False
            continue

        if not raw.startswith(("ATOM", "HETATM")):
            lines_out.append(raw)
            continue

        line = raw.rstrip("\n").ljust(80)

        # Some DUD-E PDBs use blank chain IDs plus many TER records. ProDy can
        # still infer covalent bonds across those breaks when fragments are
        # spatially close, which makes Meeko fail while padding the polymer.
        # Assign deterministic chain IDs to blank-chain segments so the breaks
        # survive receptor preparation.
        if not line[21].strip():
            chain_id = _synthetic_chain_id(blank_segment_index)
            line = f"{line[:21]}{chain_id}{line[22:]}"
            assigned_blank_chain_segments.add(chain_id)
            blank_segment_has_atoms = True

        # Extract residue name (columns 17-20 in PDB format)
        resname = line[17:20].strip()

        if resname in _SOLVENT_RESIDUES:
            dropped_solvents.add(resname)
            continue

        if resname in _RESIDUE_RENAMES:
            mapped = _RESIDUE_RENAMES[resname]
            line = f"{line[:17]}{mapped:>3}{line[20:]}"
            remapped_residues[(resname, mapped)] = remapped_residues.get((resname, mapped), 0) + 1
            resname = mapped

        if line.startswith("ATOM") and resname in _ION_RESIDUES:
            line = "HETATM" + line[6:]
            preserved_ions.add(resname)

        # Drop non-standard residues from ATOM records
        if line.startswith("ATOM") and resname and resname not in _STANDARD_RESIDUES:
            dropped_residues.add(resname)
            continue

        occ = line[54:60].strip()
        bfac = line[60:66].strip()
        element = line[76:78].strip()

        if not occ:
            line = f"{line[:54]}{1.00:6.2f}{line[60:]}"
        else:
            try:
                float(occ)
            except ValueError:
                line = f"{line[:54]}{1.00:6.2f}{line[60:]}"

        if not bfac:
            line = f"{line[:60]}{0.00:6.2f}{line[66:]}"
        else:
            try:
                float(bfac)
            except ValueError:
                line = f"{line[:60]}{0.00:6.2f}{line[66:]}"

        if not element:
            inferred = _infer_element_for_record(line[12:16], resname)
            if inferred:
                line = f"{line[:76]}{inferred:>2}{line[78:]}"
                element = inferred

        # Strip hydrogens so Meeko can safely rebuild the protonation states
        # without encountering disulfide padding errors or valency conflicts.
        if element.upper() == "H" or line[12:16].strip().startswith("H"):
            continue

        lines_out.append(line.rstrip())

    if dropped_residues:
        removed = ", ".join(sorted(dropped_residues))
        print(f"  [sanitize] dropped non-standard residues: {removed}")
        _print_receptor_integrity_warning(
            "non-standard ATOM residues were removed during sanitization",
            [
                f"Removed residue names: {removed}",
                "This is intended only for residues Meeko cannot template safely.",
            ],
        )
    if dropped_solvents:
        print(f"  [sanitize] dropped solvent residues: {', '.join(sorted(dropped_solvents))}")
    if remapped_residues:
        details = ", ".join(
            f"{src}->{dst} ({count})"
            for (src, dst), count in sorted(remapped_residues.items())
        )
        print(f"  [sanitize] remapped residues: {details}")
    if preserved_ions:
        print(f"  [sanitize] preserved ions as HETATM: {', '.join(sorted(preserved_ions))}")
    if assigned_blank_chain_segments:
        print(
            "  [sanitize] assigned synthetic chain IDs to "
            f"{len(assigned_blank_chain_segments)} blank-chain segment(s)"
        )

    output_path.write_text("\n".join(lines_out) + "\n", encoding="utf-8")
    return output_path


def _auto_receptor_prepare_command() -> str | None:
    """Auto-detect a receptor .pdb --> .pdbqt conversion tool."""
    if shutil.which("mk_prepare_receptor.py"):
        return "mk_prepare_receptor.py -i {input} -p {output} --allow_bad_res"
    if shutil.which("mk_prepare_receptor"):
        return "mk_prepare_receptor -i {input} -p {output} --allow_bad_res"

    try:
        subprocess.run([sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return f"{sys.executable} -m meeko.cli.mk_prepare_receptor -i {{input}} -p {{output}} --allow_bad_res"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    if shutil.which("obabel"):
        return "obabel {input} -O {output} -xr"
    return None


def _auto_ligand_prepare_command() -> str | None:
    """Auto-detect a ligand .mol2 --> .pdbqt conversion tool."""
    if shutil.which("mk_prepare_ligand.py"):
        return "mk_prepare_ligand.py -i {input} -o {output}"
    if shutil.which("mk_prepare_ligand"):
        return "mk_prepare_ligand -i {input} -o {output}"

    try:
        subprocess.run([sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return f"{sys.executable} -m meeko.cli.mk_prepare_ligand -i {{input}} -o {{output}}"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    if shutil.which("obabel"):
        return "obabel {input} -O {output}"
    return None


def deterministic_env(seed: int) -> dict[str, str]:
    """Environment variables that reduce avoidable run-to-run variation."""
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", str(seed))
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    return env


def format_command_template(template: str, input_path: Path, output_path: Path, seed: int) -> list[str]:
    """Expand a user-provided command template into argv."""
    rendered = template.format(
        input=str(input_path),
        output=str(output_path),
        seed=seed,
        input_stem=input_path.stem,
        output_stem=output_path.stem,
    )
    return shlex.split(rendered)


def _format_prepare_failure(
    argv: list[str],
    returncode: int,
    stdout: str,
    stderr: str,
) -> str:
    details = [f"prepare command failed: {' '.join(argv)}", f"exit code: {returncode}"]
    if stdout:
        details.append(f"stdout:\n{stdout}")
    if stderr:
        details.append(f"stderr:\n{stderr}")
    return "\n\n".join(details)


def _is_meeko_receptor_prepare_argv(argv: list[str]) -> bool:
    rendered = " ".join(argv).lower()
    return "mk_prepare_receptor" in rendered or "meeko.cli.mk_prepare_receptor" in rendered


def _has_delete_residues_arg(argv: list[str]) -> bool:
    return any(
        arg in {"-d", "--delete_residues"} or arg.startswith("--delete_residues=")
        for arg in argv
    )


def _parse_excess_bond_residues(stderr: str) -> list[str]:
    residues: list[str] = []
    padding_pairs = re.findall(
        r"Expected\s+\d+\s+paddings\s+for\s+\(([^()]+)\)\s+with bonds",
        stderr,
    )
    for pair in padding_pairs:
        for token in pair.split(","):
            residue = token.strip()
            if residue and residue not in residues:
                residues.append(residue)

    if residues:
        return residues

    for residue in re.findall(r"matched with excess inter-residue bond\(s\):\s*(\S+)", stderr):
        residue = residue.strip()
        if residue and residue not in residues:
            residues.append(residue)
    return residues


def run_prepare_command(
    *,
    command_template: str,
    input_path: Path,
    output_path: Path,
    seed: int,
) -> None:
    argv = format_command_template(command_template, input_path, output_path, seed)
    try:
        subprocess.run(
            argv,
            check=True,
            env=deterministic_env(seed),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        stdout = (exc.stdout or "").strip()
        stderr = (exc.stderr or "").strip()
        excess_residues = _parse_excess_bond_residues(stderr)
        if (
            len(excess_residues) >= 2
            and _is_meeko_receptor_prepare_argv(argv)
            and not _has_delete_residues_arg(argv)
        ):
            retry_argv = [*argv, "--delete_residues", ",".join(excess_residues)]
            removed = ", ".join(excess_residues)
            _print_receptor_integrity_warning(
                "Meeko retry will delete receptor residues",
                [
                    "The initial Meeko receptor prep failed because of excess inter-residue bonds.",
                    f"Deleting residues for retry: {removed}",
                    "Use this receptor only as a documented prep rescue, not as an untouched structure.",
                ],
            )
            try:
                subprocess.run(
                    retry_argv,
                    check=True,
                    env=deterministic_env(seed),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                return
            except subprocess.CalledProcessError as retry_exc:
                retry_stdout = (retry_exc.stdout or "").strip()
                retry_stderr = (retry_exc.stderr or "").strip()
                message = "\n\n".join(
                    [
                        _format_prepare_failure(argv, exc.returncode, stdout, stderr),
                        "retry after deleting excess-bond residues also failed",
                        _format_prepare_failure(
                            retry_argv,
                            retry_exc.returncode,
                            retry_stdout,
                            retry_stderr,
                        ),
                    ]
                )
                raise RuntimeError(message) from retry_exc

        raise RuntimeError(_format_prepare_failure(argv, exc.returncode, stdout, stderr)) from exc


def prepare_receptor_pdbqt(
    *,
    input_path: Path,
    output_path: Path,
    prepare_command: str | None,
    seed: int,
    timestamp: str = "BENCHMARK",
) -> str:
    """
    Convert a receptor into a canonical PDBQT.

    If ``input_path`` is already ``.pdbqt``, only the in-repo canonicalization
    step is applied, even when a conversion command is configured.
    """
    input_path = input_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dude-receptor-") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        prepared_input = _materialize_if_gz(input_path, tmp_dir)
        if _effective_suffix(input_path) == ".pdb":
            prepared_input = sanitize_pdb_for_meeko(
                prepared_input,
                tmp_dir / f"{_stem_without_gz(input_path)}.sanitized.pdb",
            )

        if _effective_suffix(input_path) == ".pdbqt":
            return canonicalize_receptor(prepared_input, output_path, timestamp=timestamp)

        if not prepare_command:
            prepare_command = _auto_receptor_prepare_command()
            if not prepare_command:
                raise ValueError(
                    "Receptor input is not a .pdbqt file and no conversion tool found.\n"
                    "Install Meeko (pip install meeko) or Open Babel (apt install openbabel),\n"
                    "or provide --receptor-prepare-command with {input}/{output} placeholders."
                )

        tmp_pdbqt = tmp_dir / f"{_stem_without_gz(input_path)}.pdbqt"
        run_prepare_command(
            command_template=prepare_command,
            input_path=prepared_input,
            output_path=tmp_pdbqt,
            seed=seed,
        )
        return canonicalize_receptor(tmp_pdbqt, output_path, timestamp=timestamp)


# Receptor preparation is shared with the production pipeline through molguard.
# Keep these names here as benchmark compatibility re-exports so older benchmark
# scripts still import the same symbols, but do not maintain a second active
# receptor-prep implementation in benchmarks/.
from molguard.io.receptor_prep import (  # noqa: E402
    SUPPORTED_RECEPTOR_SUFFIXES,
    _parse_excess_bond_residues,
    prepare_receptor_pdbqt,
    sanitize_pdb_for_meeko,
)


def prepare_ligand_pdbqt(
    *,
    input_path: Path,
    output_path: Path,
    prepare_command: str | None,
    seed: int,
) -> None:
    """
    Convert a ligand into a normalized PDBQT.

    If ``input_path`` is already ``.pdbqt`` and no external prepare command is
    provided, only the in-repo numeric normalization step is applied.
    """
    input_path = input_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dude-ligand-") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        prepared_input = _materialize_if_gz(input_path, tmp_dir)

        if _effective_suffix(input_path) == ".pdbqt" and prepare_command is None:
            pdbqt_normalize(prepared_input, output_path)
            return

        if not prepare_command:
            prepare_command = _auto_ligand_prepare_command()
            if not prepare_command:
                raise ValueError(
                    "Ligand input is not a .pdbqt file and no conversion tool found.\n"
                    "Install Meeko (pip install meeko) or Open Babel (apt install openbabel),\n"
                    "or provide --ligand-prepare-command with {input}/{output} placeholders."
                )

        tmp_pdbqt = tmp_dir / f"{_stem_without_gz(input_path)}.pdbqt"
        run_prepare_command(
            command_template=prepare_command,
            input_path=prepared_input,
            output_path=tmp_pdbqt,
            seed=seed,
        )
        pdbqt_normalize(tmp_pdbqt, output_path)


@dataclass(frozen=True)
class LigandInput:
    """A single logical ligand, which may come from a file or a MOL2 chunk."""

    source_path: Path
    stem: str
    mol2_block: str | None = None


def _safe_stem(name: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", (name or "").strip())
    cleaned = cleaned.strip("._-")
    return cleaned or fallback


def split_mol2_records(path: Path) -> list[LigandInput]:
    """Split a possibly multi-record MOL2 file into logical ligand entries."""
    if path.suffix.lower() == ".gz":
        text = gzip.open(path, "rt", encoding="utf-8").read()
    else:
        text = path.read_text(encoding="utf-8")
    marker = "@<TRIPOS>MOLECULE"
    starts = [match.start() for match in re.finditer(re.escape(marker), text)]
    if not starts:
        raise ValueError(f"{path} does not contain any {marker} records")

    starts.append(len(text))
    records: list[LigandInput] = []
    for index in range(len(starts) - 1):
        chunk = text[starts[index] : starts[index + 1]]
        lines = chunk.splitlines()
        title = ""
        if len(lines) > 1:
            title = lines[1].strip()
        records.append(
            LigandInput(
                source_path=path,
                stem=_safe_stem(title, f"{_stem_without_gz(path)}_{index + 1:05d}"),
                mol2_block=chunk,
            )
        )
    return records


def enumerate_ligand_inputs(paths: Iterable[Path]) -> list[LigandInput]:
    """Expand directories and multi-record MOL2 files into a flat ligand list."""
    entries: list[LigandInput] = []
    for raw_path in paths:
        path = Path(raw_path).resolve()
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if (
                    child.suffix.lower() in SUPPORTED_LIGAND_SUFFIXES
                    or _effective_suffix(child) in SUPPORTED_LIGAND_SUFFIXES
                ):
                    entries.extend(enumerate_ligand_inputs([child]))
            continue

        suffix = _effective_suffix(path)
        if suffix not in SUPPORTED_LIGAND_SUFFIXES:
            continue
        if suffix == ".mol2":
            entries.extend(split_mol2_records(path))
            continue
        entries.append(LigandInput(source_path=path, stem=_stem_without_gz(path)))
    return entries


def _prepare_one_ligand(
    entry_source_path: str,
    entry_stem: str,
    entry_mol2_block: str | None,
    out_path_str: str,
    prepare_command: str | None,
    seed: int,
) -> tuple[str | None, str | None]:
    """Prepare a single ligand.  Returns (prepared_path, error_msg)."""
    out_path = Path(out_path_str)
    try:
        if entry_mol2_block is None:
            prepare_ligand_pdbqt(
                input_path=Path(entry_source_path),
                output_path=out_path,
                prepare_command=prepare_command,
                seed=seed,
            )
        else:
            with tempfile.TemporaryDirectory(prefix="dude-mol2-record-") as tmp_dir:
                temp_input = Path(tmp_dir) / f"{entry_stem}.mol2"
                temp_input.write_text(entry_mol2_block, encoding="utf-8")
                prepare_ligand_pdbqt(
                    input_path=temp_input,
                    output_path=out_path,
                    prepare_command=prepare_command,
                    seed=seed,
                )
        if out_path.exists() and out_path.stat().st_size > 0:
            return (str(out_path), None)
        else:
            return (None, f"{entry_stem}: output file empty or missing")
    except Exception as exc:
        # Clean up partial output
        if out_path.exists():
            out_path.unlink(missing_ok=True)
        return (None, f"{entry_stem}: {exc}")


def prepare_ligand_inputs(
    *,
    inputs: Iterable[Path],
    output_dir: Path,
    prepare_command: str | None,
    seed: int,
    output_prefix: str = "",
    manifest_path: Path | None = None,
    max_workers: int | None = None,
) -> list[Path]:
    """Prepare all ligands in parallel and return the generated ``.pdbqt`` paths.

    Bad ligands are skipped with a warning instead of failing the entire batch.
    """
    import concurrent.futures

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = enumerate_ligand_inputs(inputs)
    if not entries:
        return []

    # Determine worker count
    if max_workers is None:
        max_workers = min(os.cpu_count() or 1, len(entries))
    max_workers = max(1, max_workers)

    print(f"  [ligand-prep] {len(entries)} ligands, {max_workers} workers")

    # Build tasks
    tasks: list[tuple[str, str, str | None, str, str | None, int]] = []
    entry_map: list[tuple[str, str]] = []  # (stem, source_path) for manifest
    for entry in entries:
        out_path = output_dir / f"{output_prefix}{entry.stem}.pdbqt"
        tasks.append((
            str(entry.source_path),
            entry.stem,
            entry.mol2_block,
            str(out_path),
            prepare_command,
            seed,
        ))
        entry_map.append((f"{output_prefix}{entry.stem}", str(entry.source_path)))

    # Run in parallel
    prepared_paths: list[Path] = []
    manifest_rows: list[dict[str, str]] = []
    skipped = 0

    def _record_result(idx: int, result_path: str | None, error_msg: str | None) -> None:
        nonlocal skipped
        ligand_id, source_path = entry_map[idx]
        if result_path is not None:
            prepared_paths.append(Path(result_path))
            manifest_rows.append({
                "source_path": source_path,
                "prepared_path": result_path,
                "ligand_id": ligand_id,
            })
        else:
            skipped += 1
            if skipped <= 20:
                print(f"  [SKIP] {error_msg}")
            elif skipped == 21:
                print(f"  [SKIP] ... suppressing further warnings")

    if max_workers == 1:
        for idx, task in enumerate(tasks):
            result_path, error_msg = _prepare_one_ligand(*task)
            _record_result(idx, result_path, error_msg)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(_prepare_one_ligand, *task): idx
                for idx, task in enumerate(tasks)
            }

            for future in concurrent.futures.as_completed(futures):
                idx = futures[future]
                try:
                    result_path, error_msg = future.result()
                except Exception as exc:
                    error_msg = str(exc)
                    result_path = None

                _record_result(idx, result_path, error_msg)

    if skipped:
        print(f"  [ligand-prep] skipped {skipped}/{len(entries)} bad ligands")
    print(f"  [ligand-prep] prepared {len(prepared_paths)} ligands")

    # Sort for deterministic ordering
    prepared_paths.sort(key=lambda p: p.name)

    if manifest_path is not None:
        manifest_rows.sort(key=lambda r: r["ligand_id"])
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest_rows, indent=2), encoding="utf-8")

    return prepared_paths
