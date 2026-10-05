"""Coordinate/topology assembly with explicit atom-order and membrane checks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass
class GroAtom:
    residue: int
    resname: str
    name: str
    xyz: np.ndarray


def read_gro(path: Path) -> tuple[list[GroAtom], np.ndarray]:
    lines = path.read_text().splitlines()
    count = int(lines[1])
    if len(lines) != count + 3:
        raise ValueError(f"Invalid GRO atom count: {path}")
    atoms = [GroAtom(int(s[:5]), s[5:10].strip(), s[10:15].strip(),
                     np.array([float(s[20:28]), float(s[28:36]), float(s[36:44])]))
             for s in lines[2:2 + count]]
    values = list(map(float, lines[-1].split()))
    if len(values) not in {3, 9} or not np.isfinite(values).all():
        raise ValueError(f"Invalid GRO box: {path}")
    box = np.diag(values[:3])
    if len(values) == 9:
        box[0, 1], box[0, 2], box[1, 0], box[1, 2], box[2, 0], box[2, 1] = values[3:]
    if any(not np.isfinite(a.xyz).all() for a in atoms):
        raise ValueError(f"Non-finite GRO coordinates: {path}")
    return atoms, box


def write_gro(path: Path, atoms: list[GroAtom], box) -> None:
    box = np.asarray(box)
    values = [box[0, 0], box[1, 1], box[2, 2], box[0, 1], box[0, 2],
              box[1, 0], box[1, 2], box[2, 0], box[2, 1]]
    rows = ["Ultidock protein/ligand assembly", str(len(atoms))]
    for index, atom in enumerate(atoms, 1):
        if len(atom.resname) > 5 or len(atom.name) > 5 or np.max(np.abs(atom.xyz)) >= 999:
            raise ValueError("Coordinates or atom names exceed GRO format limits")
        rows.append(f"{atom.residue % 100000:5d}{atom.resname:>5s}{atom.name:>5s}"
                    f"{index % 100000:5d}" + "".join(f"{v:8.3f}" for v in atom.xyz))
    rows.append(" ".join(f"{v:.7f}" for v in values))
    path.write_text("\n".join(rows) + "\n")


def split_itp(text: str) -> tuple[str, str]:
    if re.search(r"\[\s*(defaults|system|molecules)\s*\]", text, re.I):
        raise ValueError("Include topology must not redefine defaults/system/molecules")
    match = re.search(r"^\s*\[\s*moleculetype\s*\]", text, re.M | re.I)
    if not match:
        return text, ""
    globals_, molecule = text[:match.start()], text[match.start():]
    if re.search(r"\[\s*(atomtypes|bondtypes|angletypes|dihedraltypes|nonbond_params|pairtypes)\s*\]",
                 molecule, re.I):
        raise ValueError("Global parameter sections must precede molecule sections")
    return globals_, molecule


def itp_atoms(text: str) -> tuple[str, list[list[str]]]:
    section, molecules, atoms = "", [], []
    for raw in text.splitlines():
        line = raw.split(";")[0].strip()
        if line.startswith("["):
            section = line.strip("[] ").lower()
        elif line and not line.startswith("#"):
            if section == "moleculetype":
                molecules.append(line.split()[0])
            elif section == "atoms":
                atoms.append(line.split())
    if len(molecules) != 1 or not atoms:
        raise ValueError("Expected one ligand moleculetype with an atoms section")
    return molecules[0], atoms


def restore_ligand(itp: Path, gro: Path, positions: dict, charge: int, output: Path) -> str:
    text = itp.read_text()
    # ACPYPE may generate its own conditional restraints. Use our heavy-atom list instead.
    text = re.sub(r"#ifdef\s+POSRES\b.*?#endif[^\n]*", "", text, flags=re.S)
    molecule, rows = itp_atoms(text)
    atoms, box = read_gro(gro)
    names = [row[4] for row in rows]
    if [int(row[0]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Non-contiguous ligand topology atom indices")
    if names != [atom.name for atom in atoms] or len(set(names)) != len(names):
        raise ValueError("ACPYPE coordinate/topology atom order mismatch")
    if set(names) != set(positions):
        raise ValueError("ACPYPE renamed, added or removed ligand atoms")
    total_charge = sum(float(row[6]) for row in rows)
    if not np.isfinite(total_charge) or abs(total_charge - charge) > 0.01:
        raise ValueError("Ligand topology charge differs from requested chemical state")
    heavy = []
    for index, atom in enumerate(atoms, 1):
        atom.xyz = np.asarray(positions[atom.name])
        atom.resname, atom.residue = "LIG", 1
        if not atom.name.startswith("H"):
            heavy.append(index)
    write_gro(output / "ligand.gro", atoms, box)
    globals_, body = split_itp(text)
    (output / "ligand_atomtypes.itp").write_text(globals_)
    (output / "ligand.itp").write_text(body + '\n#ifdef POSRES_LIG\n'
                                     '#include "posre_ligand.itp"\n#endif\n')
    (output / "posre_ligand.itp").write_text(
        "[ position_restraints ]\n; atom type fx fy fz\n" +
        "".join(f"{index} 1 1000 1000 1000\n" for index in heavy))
    return molecule


def add_topology(topology: Path, global_includes: list[str],
                 molecule_includes: list[str], molecules: list[tuple[str, int]]) -> None:
    text = topology.read_text()
    ff = re.search(r'^\s*#include\s+"[^"\n]+/forcefield.itp"[^\n]*$', text, re.M)
    system = re.search(r"^\s*\[\s*system\s*\]", text, re.M | re.I)
    count = re.search(r"^\s*\[\s*molecules\s*\]", text, re.M | re.I)
    if not ff or not system or not count or not ff.end() < system.start() < count.start():
        raise ValueError("Unexpected pdb2gmx topology structure")
    def include(paths):
        return "".join(f'\n#include "{p}"' for p in paths) + "\n"
    text = (text[:ff.end()] + include(global_includes) + text[ff.end():system.start()] +
            include(molecule_includes) + text[system.start():])
    text += "\n" + "".join(f"{name} {number}\n" for name, number in molecules)
    topology.write_text(text)


def embed_bilayer(complex_atoms: list[GroAtom], config: dict):
    """Remove whole overlapping lipid molecules from a reviewed periodic bilayer."""
    membrane = config["membrane"]
    lipids, box = read_gro(Path(membrane["gro"]))
    definitions = {}
    for file in membrane["include_files"]:
        text = (Path(membrane["topology_dir"]) / file).read_text()
        section, current = "", None
        for raw in text.splitlines():
            line = raw.split(";")[0].strip()
            if line.startswith("["):
                section = line.strip("[] ").lower()
            elif line and not line.startswith("#"):
                if section == "moleculetype":
                    current = line.split()[0]
                    if current in definitions:
                        raise ValueError(f"Duplicate lipid moleculetype {current}")
                    definitions[current] = []
                elif section == "atoms" and current:
                    definitions[current].append(line.split()[4])
    lengths = np.diag(box)
    # A rectangular triclinic cell is supported. Tilted bilayers require a separate builder.
    cutoff = max(membrane["mdp_nonbonded"][k] for k in ("rvdw", "rcoulomb"))
    if not np.allclose(box, np.diag(lengths)) or (lengths <= 2 * cutoff).any():
        raise ValueError("Membrane seed must have an orthorhombic periodic box larger than twice the cutoff")
    low, high = membrane["hydrophobic_z_nm"]
    if not (0 < low < high < lengths[2]):
        raise ValueError("Invalid membrane hydrophobic_z_nm bounds")
    if config["box_shape"] == "cubic" and not np.allclose(lengths, lengths[0]):
        raise ValueError("Requested cubic box does not match the membrane seed box")
    xyz = np.array([a.xyz for a in complex_atoms])
    padding = config["padding_nm"]
    if (xyz.min(axis=0) < padding).any() or (xyz.max(axis=0) > lengths - padding).any():
        raise ValueError("Oriented complex lacks the requested periodic padding in the bilayer box")
    heavy = np.array([a.xyz for a in complex_atoms if not a.name.lstrip("0123456789").startswith("H")])
    from scipy.spatial import cKDTree
    tree = cKDTree(heavy % lengths, boxsize=lengths)
    kept, counts, offset = [], [], 0
    report = []
    for block in membrane["molecules"]:
        name, count, size = block["name"], block["count"], block["atom_count"]
        if any(not isinstance(v, int) or isinstance(v, bool) or v < 1 for v in (count, size)):
            raise ValueError("Membrane molecule count/atom_count must be positive integers")
        if name in {"SOL", "NA", "CL"} or len(definitions.get(name, [])) != size:
            raise ValueError(f"Lipid atom count differs from the {name} topology")
        retained = 0
        upper, lower = 0, 0
        for _ in range(count):
            molecule = lipids[offset:offset + size]
            if len(molecule) != size:
                raise ValueError("Membrane molecule blocks exceed coordinate atom count")
            offset += size
            if [a.name for a in molecule] != definitions[name]:
                raise ValueError(f"Seed coordinate/topology atom order mismatch for {name}")
            distances, _ = tree.query(np.array([a.xyz for a in molecule]) % lengths)
            if distances.min() >= membrane["clash_distance_nm"]:
                kept.extend(molecule)
                retained += 1
                midpoint = (low + high) / 2
                # Keep a molecule whole for leaflet assignment even when its
                # coordinates straddle the periodic boundary.
                z = np.array([a.xyz[2] for a in molecule])
                dz = z - z[0]
                center_z = (z[0] + np.mean(dz - lengths[2] * np.rint(dz / lengths[2]))) % lengths[2]
                if center_z > midpoint:
                    upper += 1
                else:
                    lower += 1
        if retained == 0:
            raise ValueError(f"No {name} lipids remain; check orientation and box size")
        counts.append((name, retained))
        report.append({"molecule": name, "initial": count, "retained": retained,
                       "upper_leaflet": upper, "lower_leaflet": lower})
    if offset != len(lipids):
        raise ValueError("Membrane must contain exactly the specified lipid molecules, no water/ions")
    if not all(sum(row[key] for row in report) for key in ("upper_leaflet", "lower_leaflet")):
        raise ValueError("Both membrane leaflets must retain lipids; check orientation and clashes")
    return kept, box, counts, report


def remove_membrane_core_water(path: Path, membrane: dict) -> int:
    """Solvate aqueous regions, keeping explicit user-declared pore waters possible."""
    atoms, box = read_gro(path)
    lengths = np.diag(box)
    if not np.allclose(box, np.diag(lengths)) or (lengths <= 0).any():
        raise ValueError("Membrane water filtering requires an orthorhombic periodic box")
    low, high = membrane["hydrophobic_z_nm"]
    if not (0 < low < high < box[2, 2]):
        raise ValueError("Invalid membrane hydrophobic_z_nm bounds")
    kept, waters, i = [], 0, 0
    sites = membrane["water_sites"]
    if sites not in (3, 4):
        raise ValueError("Membrane water filtering supports three- or four-site water")
    pores = membrane.get("water_pores", [])
    while i < len(atoms):
        atom = atoms[i]
        if atom.resname != "SOL":
            kept.append(atom)
            i += 1
            continue
        water = atoms[i:i + sites]
        if (len(water) != sites or any(a.resname != "SOL" or a.residue != atom.residue for a in water)
                or not atom.name.startswith("O") or len({a.name for a in water}) != sites):
            raise ValueError("Solvent atom grouping does not match the selected water model")
        oxygen = water[0].xyz % lengths
        pore = False
        for region in pores:
            delta = oxygen[:2] - np.asarray(region["xy_nm"])
            delta -= lengths[:2] * np.rint(delta / lengths[:2])
            if np.linalg.norm(delta) <= region["radius_nm"]:
                pore = True
                break
        if not low < oxygen[2] < high or pore:
            kept.extend(water)
            waters += 1
        i += sites
    if not waters:
        raise ValueError("No aqueous solvent remains")
    write_gro(path, kept, box)
    return waters


def set_water_count(topology: Path, count: int) -> None:
    text = topology.read_text()
    before, marker, body = text.partition("[ molecules ]")
    if not marker:
        match = re.search(r"\[\s*molecules\s*\]", text)
        if not match:
            raise ValueError("Missing molecules section")
        before, marker, body = text[:match.start()], text[match.start():match.end()], text[match.end():]
    body, changes = re.subn(r"^\s*SOL\s+\d+[^\n]*$", f"SOL {count}", body, flags=re.M)
    if changes != 1:
        raise ValueError("Expected one solvent molecule count")
    topology.write_text(before + marker + body)


def topology_charge(text: str) -> float:
    """Read the fully preprocessed topology, including water/ion definitions."""
    section, current, charges, molecules = "", None, {}, []
    checked_defaults = False
    for raw in text.splitlines():
        line = raw.split(";")[0].strip()
        if line.startswith("["):
            section = line.strip("[] ").lower()
        elif line and not line.startswith("#"):
            fields = line.split()
            if section == "defaults":
                if (fields[:3] != ["1", "2", "yes"] or
                        abs(float(fields[3]) - 0.5) > 1e-6 or
                        abs(float(fields[4]) - 1 / 1.2) > 5e-5):
                    raise ValueError("Protein/lipid defaults are incompatible with Amber GAFF2")
                checked_defaults = True
            elif section == "moleculetype":
                current = fields[0]
                charges[current] = 0
            elif section == "atoms" and current:
                charges[current] += float(fields[6])
            elif section == "molecules":
                molecules.append((fields[0], int(fields[1])))
    if not checked_defaults or not molecules:
        raise ValueError("Incomplete preprocessed topology")
    total = sum(charges[name] * count for name, count in molecules)
    if not np.isfinite(total) or abs(total - round(total)) > 0.01:
        raise ValueError("System topology has a non-integer or non-finite net charge")
    return round(total)
