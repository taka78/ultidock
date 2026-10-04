"""Restore ligand chemistry from an explicit atom map; never infer bonds from PDBQT."""

from pathlib import Path

import numpy as np


AD_ELEMENTS = {"A": "C", "OA": "O", "OS": "O", "NA": "N", "NS": "N",
               "SA": "S", "HD": "H", "HS": "H"}


def atom_xyz(line: str) -> np.ndarray:
    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not np.isfinite(xyz).all():
        raise ValueError("Non-finite molecular coordinates")
    return xyz


def dock_atoms(lines: list[str]) -> dict[int, tuple[str, np.ndarray]]:
    atoms = {}
    for line in lines:
        serial = int(line[6:11])
        element = AD_ELEMENTS.get(line.split()[-1], line.split()[-1])
        if serial in atoms:
            raise ValueError(f"Duplicate docking atom serial {serial}")
        atoms[serial] = (element, atom_xyz(line))
    return atoms


def posed_molecule(sdf: Path, atom_map: dict, lines: list[str], charge: int, transform):
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
    if len(supplier) != 1 or supplier[0] is None:
        raise ValueError(f"Expected exactly one valid chemical structure in {sdf}")
    original = supplier[0]
    if any(label == "?" for _, label in Chem.FindMolChiralCenters(original, includeUnassigned=True)):
        raise ValueError("Specify ligand stereochemistry explicitly in the chemical source")
    if len(Chem.GetMolFrags(original)) != 1:
        raise ValueError("Ligand source must be one connected molecule, without counterions")
    if Chem.GetFormalCharge(original) != charge:
        raise ValueError("Requested ligand net charge disagrees with the SDF formal charge")
    if any(atom.GetNumRadicalElectrons() or atom.GetSymbol() not in
           {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
           for atom in original.GetAtoms()):
        raise ValueError("Automatic GAFF2 route supports closed-shell organic ligands only")
    declared = {a.GetIdx(): a.GetProp("_CIPCode") for a in original.GetAtoms()
                if a.HasProp("_CIPCode")}
    for atom in original.GetAtoms():
        atom.SetIntProp("source_index", atom.GetIdx())
    mol = Chem.RemoveHs(original)
    source_indices = {a.GetIntProp("source_index"): a.GetIdx() for a in mol.GetAtoms()}
    all_docked = dock_atoms(lines)
    atoms = {serial: value for serial, value in all_docked.items() if value[0] != "H"}
    mapping = {int(serial): int(index) for serial, index in atom_map.items()}
    if set(mapping) != set(atoms) or set(mapping.values()) != set(source_indices):
        raise ValueError("atom_map must cover each docking heavy-atom serial and each SDF "
                         "heavy-atom index exactly once (SDF indices are zero-based)")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Duplicate atom-map target")
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    matrix = np.asarray(transform)
    for serial, old_index in mapping.items():
        index = source_indices[old_index]
        element, xyz = atoms[serial]
        if mol.GetAtomWithIdx(index).GetSymbol() != element:
            raise ValueError(f"Atom map element mismatch at docking serial {serial}")
        placed = (matrix @ np.r_[xyz, 1])[:3]
        conf.SetAtomPosition(index, Point3D(*map(float, placed)))
    mol.AddConformer(conf)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = np.linalg.norm(np.array(conf.GetAtomPosition(a)) - np.array(conf.GetAtomPosition(b)))
        radii = Chem.GetPeriodicTable()
        expected = sum(radii.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in (a, b))
        if not 0.65 * expected < distance < 1.45 * expected:
            raise ValueError("Mapped docking atoms have implausible bond lengths; check atom_map")
    stereo = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(stereo, replaceExistingTags=True)
    for old_index, label in declared.items():
        atom = stereo.GetAtomWithIdx(source_indices[old_index])
        if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != label:
            raise ValueError("Docking coordinates disagree with the SDF stereochemistry")
    mol = Chem.AddHs(mol, addCoords=True)
    # Retain explicitly docked hydrogens (usually polar H in PDBQT) where the
    # reviewed chemical graph supplies the corresponding bonded hydrogen.
    used = set()
    conf = mol.GetConformer()
    for element, xyz in all_docked.values():
        if element != "H":
            continue
        distances = sorted((np.linalg.norm(xyz - position), serial)
                           for serial, (_, position) in atoms.items())
        if not distances or not 0.6 < distances[0][0] < 1.4:
            raise ValueError("Docked hydrogen has no plausible heavy-atom bond")
        parent = source_indices[mapping[distances[0][1]]]
        candidates = [a.GetIdx() for a in mol.GetAtomWithIdx(parent).GetNeighbors()
                      if a.GetAtomicNum() == 1 and a.GetIdx() not in used]
        if not candidates:
            raise ValueError("Docked hydrogen count disagrees with the chemical protonation state")
        placed = (matrix @ np.r_[xyz, 1])[:3]
        hydrogen = min(candidates, key=lambda i: np.linalg.norm(
            np.array(conf.GetAtomPosition(i)) - placed))
        conf.SetAtomPosition(hydrogen, Point3D(*map(float, placed)))
        used.add(hydrogen)
    return mol


def convert_mol2(mol2: Path, mol) -> dict[str, np.ndarray]:
    """Name atoms explicitly so ACPYPE's topology reordering cannot scramble the pose."""
    lines = mol2.read_text().splitlines()
    section, found, positions = "", [], mol.GetConformer().GetPositions()
    names = {}
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            section = line
        elif section == "@<TRIPOS>ATOM" and line.strip():
            fields = line.split()
            index = int(fields[0]) - 1
            if index != len(found) or index >= mol.GetNumAtoms():
                raise ValueError("Open Babel changed the chemical atom ordering")
            xyz = np.array(list(map(float, fields[2:5])))
            if np.linalg.norm(xyz - positions[index]) > 0.002:
                raise ValueError("Open Babel moved ligand coordinates")
            name = f"{mol.GetAtomWithIdx(index).GetSymbol()}{index + 1}"
            if len(name) > 5:
                raise ValueError("Ligand atom name exceeds GRO format limit")
            fields[1] = name
            fields[6:8] = ["1", "LIG"]
            lines[i] = " ".join(fields)
            names[name] = positions[index] / 10  # Angstrom -> nm
            found.append(index)
    if len(found) != mol.GetNumAtoms():
        raise ValueError("Incomplete MOL2 conversion")
    mol2.write_text("\n".join(lines) + "\n")
    return names


def validate_receptor(pdb: Path, pdbqt: Path, tolerance: float = 0.02) -> None:
    """Require the reviewed protein to be in the docking coordinate frame."""
    def collect(path, docking):
        atoms = {}
        for line in path.read_text().splitlines():
            if line.startswith("MODEL") or (not docking and line.startswith("HETATM")):
                raise ValueError("Use a single-model, protein-only reviewed receptor PDB; "
                                 "cofactors/structural waters require a custom topology")
            if not line.startswith("ATOM  "):
                continue
            element = line.split()[-1] if docking else line[76:78].strip()
            if element in {"H", "HD", "HS"} or line[12:16].strip().lstrip("0123456789").startswith("H"):
                continue
            if line[16:17] not in {" ", ""}:
                raise ValueError("Resolve alternate receptor conformations before MD")
            key = (line[21:27], line[12:16].strip())
            if key in atoms:
                raise ValueError("Duplicate receptor atom identity")
            atoms[key] = atom_xyz(line)
        return atoms
    source, docked = collect(pdb, False), collect(pdbqt, True)
    if not docked or set(source) != set(docked):
        raise ValueError("Reviewed receptor heavy atoms differ from the docking receptor; "
                         "repair the receptor and redock rather than silently changing the complex")
    if any(np.linalg.norm(source[key] - xyz) > tolerance for key, xyz in docked.items()):
        raise ValueError("Reviewed receptor is outside the docking coordinate frame")


def transform_pdb(source: Path, target: Path, transform) -> None:
    matrix, rows = np.asarray(transform), []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = (matrix @ np.r_[atom_xyz(line), 1])[:3]
            line = line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]
        rows.append(line)
    target.write_text("\n".join(rows) + "\n")
