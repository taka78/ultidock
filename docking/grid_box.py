# grid_box.py
import os, hashlib, math
import numpy as np

# --- box policy (Å) ---
MIN_BOX_A = 21.0
MAX_BOX_A = 50.0

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
            except:
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

# --- NEW: clamp by Å first, then convert to odd npts ---
def _npts_from_size_clamped(size_A, spacing, min_box_A=MIN_BOX_A, max_box_A=MAX_BOX_A):
    """
    Clamp a physical side length (Å) into [min_box_A, max_box_A], then
    convert to an odd grid count at the given spacing.
    """
    side = float(np.clip(size_A, float(min_box_A), float(max_box_A)))
    n = max(1, int(math.ceil(side / float(spacing))))
    if n % 2 == 0:
        n += 1
    return n

def ligand_mode_box(ref_ligand_path, margin=5.0, spacing=0.375, ignore_h=True):
    atoms = _load_coords_from_pdb_like(
        ref_ligand_path,
        atom_selector=(lambda n: (not n.startswith("H"))) if ignore_h else (lambda n: True),
    )
    mins, maxs, center = _bbox(atoms)
    size = _extents(mins, maxs) + margin*2.0  # Å per axis
    npts = np.array([_npts_from_size_clamped(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def residues_mode_box(receptor_path, residue_predicate, margin=5.0, spacing=0.375):
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
    npts = np.array([_npts_from_size_clamped(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def centers_file_mode(centers_tsv_path, key=None, spacing=0.375, npts=(81,81,81)):
    """
    centers.tsv format:
    key  cx  cy  cz  [nx ny nz spacing]
    If nx/ny/nz are missing, compute them from policy (Å→npts).
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

    if len(parts) < 4:
        raise ValueError("centers row must have at least: key cx cy cz")
    cx, cy, cz = map(float, parts[1:4])

    if len(parts) >= 8:
        # full row with nx ny nz spacing
        nx, ny, nz = map(int, parts[4:7])
        sp = float(parts[7])
    elif len(parts) >= 7:
        nx, ny, nz = map(int, parts[4:7])
        sp = float(spacing)
    else:
        # no npts given -> enforce policy using current spacing
        sp = float(spacing)
        # choose a default within policy (e.g., middle of the range) or MIN
        target_side = max(MIN_BOX_A, min(MAX_BOX_A, 32.0))
        n = _npts_from_size_clamped(target_side, sp)
        nx = ny = nz = n

    return np.array([cx, cy, cz], float), np.array([nx, ny, nz], int), float(sp)

def blind_box(receptor_path, cap=30.0, spacing=0.375):
    """
    Blind box over the whole receptor, but clamp each axis to [MIN_BOX_A, MAX_BOX_A] Å.
    Also respect the 'cap' (Å) ceiling to avoid massive boxes.
    """
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, center = _bbox(rec)
    size = _extents(mins, maxs)
    size = np.minimum(size, np.array([cap, cap, cap], dtype=float))
    # Enforce policy range in Å
    size = np.clip(size, MIN_BOX_A, MAX_BOX_A)
    npts = np.array([_npts_from_size_clamped(s, spacing) for s in size], dtype=int)
    return center, npts, spacing

def validate_center_inside_receptor(center, receptor_path):
    rec = _load_coords_from_pdb_like(receptor_path, atom_selector=lambda n: True)
    mins, maxs, _ = _bbox(rec)
    inside = np.all(center >= mins) and np.all(center <= maxs)
    return bool(inside), mins, maxs

def grid_signature(receptor_path, mode_name, params_dict):
    key = f"{_hash_file(receptor_path)}|{mode_name}|{sorted(params_dict.items())}"
    return hashlib.sha1(key.encode()).hexdigest()[:16]
