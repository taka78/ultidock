"""Shared deterministic receptor preparation for Ultidock.

This module is the single receptor-prep path used by both the normal docking
pipeline and the benchmark scripts. It discovers receptor inputs by extension,
canonicalizes existing PDBQT receptors, and converts raw PDB/MOL2 receptors into
canonical receptor PDBQT files with the same sanitizer and Meeko fallback rules
everywhere unless a caller explicitly provides different options.
"""

from __future__ import annotations

import gzip
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

from molguard.io.pdbqt import canonicalize_receptor

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


def receptor_output_stem(path: Path) -> str:
    """Return the receptor stem used for generated PDBQT files and site folders."""
    return _stem_without_gz(Path(path))


def receptor_pdbqt_output_name(path: Path) -> str:
    """Return the canonical PDBQT filename for a receptor input path."""
    return f"{receptor_output_stem(path)}.pdbqt"


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

    if field[0].isspace():
        return token[0].upper()

    token = token.upper()
    if len(token) >= 2 and token[:2] in {
        "CL",
        "BR",
        "NA",
        "MG",
        "ZN",
        "FE",
        "MN",
        "CU",
        "NI",
        "CO",
        "CD",
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


_STANDARD_RESIDUES: set[str] = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "HSD",
    "HSE",
    "HSP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ACE",
    "NME",
    "NH2",
    "DA",
    "DT",
    "DC",
    "DG",
    "A",
    "U",
    "C",
    "G",
    "HOH",
    "WAT",
    "TIP",
    "TIP3",
}

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


@dataclass(frozen=True)
class ReceptorPrepRecord:
    """A receptor-prep action performed for a single input file."""

    source_path: Path
    output_path: Path
    digest: str
    action: str


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

    The sanitizer fills missing occupancy, B-factor, and element columns,
    removes waters/hydrogens, remaps common DUD-E residue aliases, preserves
    monoatomic ions as HETATM, and assigns synthetic chain IDs to blank-chain
    TER-delimited fragments so preparation respects chain breaks.
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

        if not line[21].strip():
            chain_id = _synthetic_chain_id(blank_segment_index)
            line = f"{line[:21]}{chain_id}{line[22:]}"
            assigned_blank_chain_segments.add(chain_id)
            blank_segment_has_atoms = True

        resname = line[17:20].strip()

        if resname in _SOLVENT_RESIDUES:
            dropped_solvents.add(resname)
            continue

        if resname in _RESIDUE_RENAMES:
            mapped = _RESIDUE_RENAMES[resname]
            line = f"{line[:17]}{mapped:>3}{line[20:]}"
            remapped_residues[(resname, mapped)] = (
                remapped_residues.get((resname, mapped), 0) + 1
            )
            resname = mapped

        if line.startswith("ATOM") and resname in _ION_RESIDUES:
            line = "HETATM" + line[6:]
            preserved_ions.add(resname)

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
    """Auto-detect a receptor .pdb -> .pdbqt conversion tool."""
    if shutil.which("mk_prepare_receptor.py"):
        return "mk_prepare_receptor.py -i {input} -p {output} --allow_bad_res"
    if shutil.which("mk_prepare_receptor"):
        return "mk_prepare_receptor -i {input} -p {output} --allow_bad_res"

    try:
        subprocess.run(
            [sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "-h"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return (
            f"{sys.executable} -m meeko.cli.mk_prepare_receptor "
            "-i {input} -p {output} --allow_bad_res"
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    if shutil.which("obabel"):
        return "obabel {input} -O {output} -xr"
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
    """Extract residue IDs from Meeko excess-bond padding errors."""
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
    """Run a receptor conversion command with the shared Meeko retry policy."""
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
    timestamp: str = "PIPELINE",
) -> str:
    """
    Convert a receptor into a canonical PDBQT.

    If ``input_path`` is already ``.pdbqt``, only canonicalization is applied,
    even when a conversion command is configured. PDB inputs are sanitized before
    external conversion.
    """
    input_path = input_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="molguard-receptor-") as tmp_dir_name:
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
                    "or provide a receptor prepare command with {input}/{output} placeholders."
                )

        tmp_pdbqt = tmp_dir / receptor_pdbqt_output_name(input_path)
        run_prepare_command(
            command_template=prepare_command,
            input_path=prepared_input,
            output_path=tmp_pdbqt,
            seed=seed,
        )
        return canonicalize_receptor(tmp_pdbqt, output_path, timestamp=timestamp)


def receptor_inputs_in_directory(macro_dir: Path) -> list[Path]:
    """Return receptor inputs in a directory, including gzip-wrapped payloads."""
    entries: list[Path] = []
    for path in sorted(macro_dir.iterdir()):
        if not path.is_file():
            continue
        if _effective_suffix(path) in SUPPORTED_RECEPTOR_SUFFIXES:
            entries.append(path)
    return entries


def prepare_receptors_in_directory(
    macro_dir: Path,
    *,
    prepare_command: str | None = None,
    seed: int = 42,
    timestamp: str = "PIPELINE",
    force: bool = False,
) -> list[ReceptorPrepRecord]:
    """
    Prepare/canonicalize every receptor input in ``macro_dir``.

    Existing PDBQT files are canonicalized in place.  PDB/MOL2 inputs are
    converted to sibling ``.pdbqt`` files unless such an output already exists
    and ``force`` is false.
    """
    macro_dir = Path(macro_dir).resolve()
    macro_dir.mkdir(parents=True, exist_ok=True)
    inputs = receptor_inputs_in_directory(macro_dir)
    if not inputs:
        return []

    existing_pdbqt_stems = {
        receptor_output_stem(path)
        for path in inputs
        if _effective_suffix(path) == ".pdbqt" and path.suffix.lower() != ".gz"
    }

    records: list[ReceptorPrepRecord] = []
    for source in inputs:
        suffix = _effective_suffix(source)
        stem = receptor_output_stem(source)
        if suffix == ".pdbqt":
            if source.suffix.lower() == ".gz":
                output = macro_dir / receptor_pdbqt_output_name(source)
                if output.exists() and not force:
                    continue
            else:
                output = source
            action = "canonicalized"
        else:
            output = macro_dir / receptor_pdbqt_output_name(source)
            if stem in existing_pdbqt_stems and not force:
                print(
                    f"  [receptor-prep] skipping {source.name}: "
                    f"{output.name} already exists"
                )
                continue
            if output.exists() and not force:
                print(
                    f"  [receptor-prep] reusing existing {output.name}; "
                    "use force to regenerate"
                )
                continue
            action = "prepared"

        try:
            digest = prepare_receptor_pdbqt(
                input_path=source,
                output_path=output,
                prepare_command=prepare_command,
                seed=seed,
                timestamp=timestamp,
            )
        except Exception as exc:
            print(f"  [receptor-prep] FAIL {source.name}: {exc}")
            continue
        records.append(
            ReceptorPrepRecord(
                source_path=source,
                output_path=output,
                digest=digest,
                action=action,
            )
        )

    return records
