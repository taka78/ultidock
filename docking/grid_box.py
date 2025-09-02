# grid_box.py
import os, hashlib, math
import numpy as np

ODD = lambda n: n if n % 2 == 1 else n + 1

def _hash_file(path, block=1<<20):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()[:12]

def _load_coords_from_pdb_like(path, atom_selector=lambda name: True):
    coords = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            if not atom_selector(name):  # e.g., ignore hydrogens
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
            except:  # robust against malformed lines
                continue
    if not coords:
        raise ValueError(f"No atom coordinates parsed from: {path}")
    return np.asarray(coords, dtype=float)

def _bbox(coords):
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    return mins, maxs, (maxs + mins) / 2.0

def _extents(mins, maxs):
    return maxs - mins

def _npts_for_size(side_len_angstrom, spacing=0.375, clamp=(25, 126)):
    # Convert physical Å length to grid points (odd), clamp to safe range.
    n = max(ODD(int(math.ceil(side_len_angstrom / spacing))), clamp[0])
    return min(n, clamp[1])

def ligand_mode_box(ref_ligand_path, margin=5.0, spacing=0.375, ignore_h=True):
    atoms = _load_coords_from_pdb_like(
        ref_ligand_path,
        atom_selector=(lambda n: (not n.startswith("H"))) if ignore_h else (lambda n: True),
    )
    mins, maxs, center = _bbox(atoms)
    size = _extents(mins, maxs) + margin*2.0  # Å per axis
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def residues_mode_box(receptor_path, residue_predicate, margin=5.0, spacing=0.375):
    # residue_predicate: function(line)->bool (select lines belonging to the site residues)
    coords = []
    with open(receptor_path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
            if not residue_predicate(line): continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
            except:
                pass
    if not coords:
        raise ValueError("No residue coords selected for binding site.")
    coords = np.asarray(coords, dtype=float)
    mins, maxs, center = _bbox(coords)
    size = _extents(mins, maxs) + margin*2.0
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def centers_file_mode(centers_tsv_path, key=None, spacing=0.375, npts=(81,81,81)):
    """
    centers.tsv format:
    key    cx    cy    cz    [nx ny nz spacing]
    If key is None and file has a single row, use it. Otherwise match by key (e.g., receptor name).
    """
    rows = []
    with open(centers_tsv_path, "r", errors="ignore") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            parts = ln.split()
            rows.append(parts)
    if not rows:
        raise ValueError("centers file empty.")
    if key is None and len(rows) == 1:
        parts = rows[0]
    else:
        found = None
        for parts in rows:
            if parts[0] == key:
                found = parts; break
        if not found:
            raise ValueError(f"Key {key} not found in centers file.")
        parts = found
    # parse
    if len(parts) < 4:
        raise ValueError("centers row must have at least: key cx cy cz")
    cx, cy, cz = map(float, parts[1:4])
    if len(parts) >= 7:
        nx, ny, nz = map(int, parts[4:7])
    else:
        nx, ny, nz = npts
    sp = spacing if len(parts) < 8 else float(parts[7])
    return np.array([cx, cy, cz], float), np.array([nx, ny, nz], int), sp

def blind_box(receptor_path, cap=50.0, spacing=0.375):
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, center = _bbox(rec)
    size = _extents(mins, maxs)
    # Cap each axis length to avoid massive grids
    size = np.minimum(size, np.array([cap, cap, cap]))
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)# grid_box.py
import os, hashlib, math
import numpy as np

ODD = lambda n: n if n % 2 == 1 else n + 1

def _hash_file(path, block=1<<20):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()[:12]

def _load_coords_from_pdb_like(path, atom_selector=lambda name: True):
    coords = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            if not atom_selector(name):  # e.g., ignore hydrogens
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
            except:  # robust against malformed lines
                continue
    if not coords:
        raise ValueError(f"No atom coordinates parsed from: {path}")
    return np.asarray(coords, dtype=float)

def _bbox(coords):
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    return mins, maxs, (maxs + mins) / 2.0

def _extents(mins, maxs):
    return maxs - mins

def _npts_for_size(side_len_angstrom, spacing=0.375, clamp=(25, 126)):
    # Convert physical Å length to grid points (odd), clamp to safe range.
    n = max(ODD(int(math.ceil(side_len_angstrom / spacing))), clamp[0])
    return min(n, clamp[1])

def ligand_mode_box(ref_ligand_path, margin=5.0, spacing=0.375, ignore_h=True):
    atoms = _load_coords_from_pdb_like(
        ref_ligand_path,
        atom_selector=(lambda n: (not n.startswith("H"))) if ignore_h else (lambda n: True),
    )
    mins, maxs, center = _bbox(atoms)
    size = _extents(mins, maxs) + margin*2.0  # Å per axis
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def residues_mode_box(receptor_path, residue_predicate, margin=5.0, spacing=0.375):
    # residue_predicate: function(line)->bool (select lines belonging to the site residues)
    coords = []
    with open(receptor_path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
            if not residue_predicate(line): continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
            except:
                pass
    if not coords:
        raise ValueError("No residue coords selected for binding site.")
    coords = np.asarray(coords, dtype=float)
    mins, maxs, center = _bbox(coords)
    size = _extents(mins, maxs) + margin*2.0
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def centers_file_mode(centers_tsv_path, key=None, spacing=0.375, npts=(81,81,81)):
    """
    centers.tsv format:
    key    cx    cy    cz    [nx ny nz spacing]
    If key is None and file has a single row, use it. Otherwise match by key (e.g., receptor name).
    """
    rows = []
    with open(centers_tsv_path, "r", errors="ignore") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            parts = ln.split()
            rows.append(parts)
    if not rows:
        raise ValueError("centers file empty.")
    if key is None and len(rows) == 1:
        parts = rows[0]
    else:
        found = None
        for parts in rows:
            if parts[0] == key:
                found = parts; break
        if not found:
            raise ValueError(f"Key {key} not found in centers file.")
        parts = found
    # parse
    if len(parts) < 4:
        raise ValueError("centers row must have at least: key cx cy cz")
    cx, cy, cz = map(float, parts[1:4])
    if len(parts) >= 7:
        nx, ny, nz = map(int, parts[4:7])
    else:
        nx, ny, nz = npts
    sp = spacing if len(parts) < 8 else float(parts[7])
    return np.array([cx, cy, cz], float), np.array([nx, ny, nz], int), sp

def blind_box(receptor_path, cap=50.0, spacing=0.375):
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, center = _bbox(rec)
    size = _extents(mins, maxs)
    # Cap each axis length to avoid massive grids
    size = np.minimum(size, np.array([cap, cap, cap]))
    npts = np.array([_npts_for_size(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def validate_center_inside_receptor(center, receptor_path):
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, _ = _bbox(rec)
    inside = np.all(center >= mins) and np.all(center <= maxs)
    return bool(inside), mins, maxs

def grid_signature(receptor_path, mode_name, params_dict):
    key = f"{_hash_file(receptor_path)}|{mode_name}|{sorted(params_dict.items())}"
    return hashlib.sha1(key.encode()).hexdigest()[:16]

    return center, npts, spacing

def validate_center_inside_receptor(center, receptor_path):
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, _ = _bbox(rec)
    inside = np.all(center >= mins) and np.all(center <= maxs)
    return bool(inside), mins, maxs

def grid_signature(receptor_path, mode_name, params_dict):
    key = f"{_hash_file(receptor_path)}|{mode_name}|{sorted(params_dict.items())}"
    return hashlib.sha1(key.encode()).hexdigest()[:16]