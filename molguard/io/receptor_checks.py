"""Validation and conservative hydrogen recovery for rigid AD4 receptors."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import shutil
import subprocess

from molguard.io.pdbqt import LintError, pdbqt_check, _reformat_atom_line

# Receptor types supported by the bundled AD4 parameter library. Ligand maps
# are a separate list: in particular H, HS, NS and OS can occur in receptors.
AD4_RECEPTOR_ELEMENTS = {
    "H": "H", "HD": "H", "HS": "H", "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O",
    "S": "S", "SA": "S", "P": "P", "F": "F", "Cl": "Cl", "Br": "Br",
    "I": "I", "Mg": "Mg", "Mn": "Mn", "Zn": "Zn", "Ca": "Ca", "Fe": "Fe",
    "MG": "Mg", "MN": "Mn", "ZN": "Zn", "CA": "Ca", "FE": "Fe",
    "CL": "Cl", "BR": "Br",
}


def receptor_atom_lines(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="ascii").splitlines()
            if line[:6].strip() in {"ATOM", "HETATM"}]


def has_donor_hydrogens(path: Path) -> bool:
    return any(line[77:79].strip() in {"HD", "HS"} for line in receptor_atom_lines(path))


def validate_receptor(path: Path, *, require_donors: bool = True) -> None:
    """Reject unrepresentable receptors before canonicalization loses context."""
    report = pdbqt_check(path)
    if report.errors:
        details = "; ".join(f"line {e.line_no} {e.column}: {e.message}" for e in report.errors[:5])
        raise LintError(f"{path.name}: {details}")
    lines = path.read_text(encoding="ascii").splitlines()
    if sum(line[:6].strip() == "MODEL" for line in lines) > 1:
        raise ValueError(f"{path.name}: multiple MODEL structures; select one model per receptor file")
    if any(line.split() and line.split()[0] in
           {"ROOT", "BRANCH", "TORSDOF", "BEGIN_RES"} for line in lines):
        raise ValueError(f"{path.name}: flexible/ligand torsion records in a rigid receptor; "
                         "supply the rigid receptor separately")
    atoms = receptor_atom_lines(path)
    if not atoms:
        raise ValueError(f"{path.name}: no ATOM/HETATM records")
    if len(atoms) > 99999:
        raise ValueError(f"{path.name}: more than 99999 atoms cannot fit PDBQT serial columns")
    for line in atoms:
        atom_type = line[77:79].strip()
        if atom_type not in AD4_RECEPTOR_ELEMENTS:
            raise ValueError(f"{path.name}: unsupported receptor atom type {atom_type!r}")
        # Blank numeric fields pass the generic linter, but are not usable here.
        for start, end in ((30, 38), (38, 46), (46, 54), (70, 76)):
            if not line[start:end].strip():
                raise ValueError(f"{path.name}: missing coordinate or partial charge in {line!r}")
    if require_donors and not has_donor_hydrogens(path):
        raise ValueError(
            f"{path.name}: no HD or HS donor hydrogens. Run 'molguard receptor prepare "
            f"\"{path}\" -o \"{path}\"' with Open Babel installed, or prepare a "
            "hydrogen-complete receptor from the source structure. Receptors with no "
            "chemical donors are unsupported by the bundled AutoGrid affinity-map builder."
        )


def normalize_receptor_input(source: Path, destination: Path) -> None:
    """Repair representable numeric formatting without guessing missing values."""
    report = pdbqt_check(source)
    fatal = [e for e in report.errors if e.code not in {"EXPONENT", "NO_DECIMAL"}]
    if fatal:
        details = "; ".join(f"line {e.line_no} {e.column}: {e.message}" for e in fatal[:5])
        raise LintError(f"{source.name}: {details}")
    lines = source.read_text(encoding="ascii").splitlines()
    normalized = []
    for number, line in enumerate(lines, 1):
        if line[:6].strip() in {"ATOM", "HETATM"}:
            for start, end in ((30, 38), (38, 46), (46, 54), (70, 76)):
                if not line[start:end].strip():
                    raise ValueError(f"{source.name}: line {number}: missing coordinate or partial charge")
            # These metadata fields have harmless defaults; coordinates/charges do not.
            if not line[54:60].strip():
                line = line[:54] + "  1.00" + line[60:]
            if not line[60:66].strip():
                line = line[:60] + "  0.00" + line[66:]
            line = _reformat_atom_line(line)
        normalized.append(line)
    destination.write_text("\n".join(normalized) + "\n", encoding="ascii")
    validate_receptor(destination, require_donors=False)


def _heavy_atoms(path: Path) -> Counter:
    # Preserve identity, multiplicity, elements and coordinates (PDB precision).
    return Counter(
        (line[12:27], tuple(float(line[s:s + 8]) for s in (30, 38, 46)),
         AD4_RECEPTOR_ELEMENTS[line[77:79].strip()])
        for line in receptor_atom_lines(path)
        if AD4_RECEPTOR_ELEMENTS[line[77:79].strip()] != "H"
    )


def recover_donor_hydrogens(source: Path, destination: Path) -> dict:
    """Re-perceive bonds, add H and recalculate charges; never delete heavy atoms.

    PDBQT does not retain bond orders. Open Babel's geometry-based recovery is
    recorded explicitly and is not a guarantee of correct protonation chemistry.
    """
    executable = shutil.which("obabel")
    if not executable:
        raise ValueError(f"{source.name}: missing donor hydrogens; install Open Babel "
                         "(obabel) for recovery, or supply a prepared receptor/source PDB")
    pdb = destination.with_suffix(".pdb")
    lines = []
    for line in receptor_atom_lines(source):
        element = AD4_RECEPTOR_ELEMENTS[line[77:79].strip()]
        if element != "H":
            lines.append(line[:66] + " " * 10 + f"{element:>2}")
    pdb.write_text("\n".join(lines) + "\nEND\n", encoding="ascii")
    argv = [executable, "-ipdb", str(pdb), "-opdbqt", "-O", str(destination),
            "-xr", "-xn", "-h", "--partialcharge", "gasteiger"]
    try:
        result = subprocess.run(argv, check=True, capture_output=True, text=True, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"Hydrogen recovery failed: {exc}\n{exc.stderr or ''}") from exc
    if not destination.is_file():
        raise RuntimeError("Hydrogen recovery produced no receptor PDBQT")
    normalized = destination.with_suffix(".normalized.pdbqt")
    normalize_receptor_input(destination, normalized)
    validate_receptor(normalized)
    if _heavy_atoms(source) != _heavy_atoms(normalized):
        raise ValueError("Hydrogen recovery changed heavy-atom identities or coordinates; "
                         "original receptor preserved. Prepare from the source structure.")
    normalized.replace(destination)
    diagnostics = (result.stdout + "\n" + result.stderr).strip()
    print("  [receptor-prep] recovered donor hydrogens with Open Babel; "
          "recalculated atom types and Gasteiger charges; preserved all heavy atoms")
    if diagnostics:
        print(f"  [receptor-prep] converter diagnostics:\n{diagnostics}")
    return {"method": "Open Babel geometry-based hydrogen recovery",
            "command": argv, "diagnostics": diagnostics,
            "heavy_atoms_preserved": sum(_heavy_atoms(source).values()),
            "donor_hydrogens": sum(line[77:79].strip() in {"HD", "HS"}
                                   for line in receptor_atom_lines(destination)),
            "limitation": "Review protonation and bond orders, especially converter warnings; "
                          "missing heavy atoms are not reconstructed."}
