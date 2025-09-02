# make_grids.py
import os, json, subprocess
import numpy as np
from grid_box import (
    ligand_mode_box, residues_mode_box, centers_file_mode, blind_box,
    validate_center_inside_receptor, grid_signature
)

def write_gpf(gpf_path, receptor_pdbqt, center, npts, spacing, atom_types):
    # AutoGrid4 GPF template: change as needed (maps for each atom type)
    with open(gpf_path, "w") as g:
        g.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
        g.write(f"spacing {spacing:.3f}\n")
        g.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
        g.write(f"receptor_types { ' '.join(sorted(set(atom_types))) }\n")
        g.write(f"receptor {os.path.basename(receptor_pdbqt)}\n")
        g.write("smooth 0.5\n")
        g.write("map grid\n")  # optional; your pipeline may set maps per type explicitly

def autodetect_receptor_types(receptor_pdbqt):
    types = set()
    with open(receptor_pdbqt, "r", errors="ignore") as f:
        for ln in f:
            if ln.startswith("ATOM") or ln.startswith("HETATM"):
                t = ln[77:79].strip()  # AD4 keeps atom type in columns 78-79 usually
                if t:
                    types.add(t)
    # fallback minimal set if parsing fails:
    return types or {"C", "A", "HD", "N", "NA", "OA", "SA", "Cl", "F", "P", "S"}

def run_autogrid4(gpf_path, workdir, autogrid4_bin="autogrid4"):
    subprocess.run([autogrid4_bin, "-p", os.path.basename(gpf_path)], cwd=workdir, check=True)

def ensure_grids(
    receptor_pdbqt, out_dir, mode,
    ref_ligand=None, centers_file=None, residues_predicate=None,
    spacing=0.375, margin=5.0, cap=50.0, autogrid4_bin="autogrid4"
):
    os.makedirs(out_dir, exist_ok=True)
    params = {"spacing": spacing, "margin": margin, "cap": cap, "mode": mode}
    # pick mode
    if mode == "ligand":
        center, npts, sp = ligand_mode_box(ref_ligand, margin=margin, spacing=spacing)
    elif mode == "residues":
        center, npts, sp = residues_mode_box(receptor_pdbqt, residues_predicate, margin=margin, spacing=spacing)
    elif mode == "centers":
        key = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
        center, npts, sp = centers_file_mode(centers_file, key=key, spacing=spacing)
    elif mode == "blind":
        center, npts, sp = blind_box(receptor_path=receptor_pdbqt, cap=cap, spacing=spacing)
    else:
        raise ValueError("mode must be ligand|residues|centers|blind")

    inside, mins, maxs = validate_center_inside_receptor(center, receptor_pdbqt)
    if not inside:
        # Clamp the center to receptor bounds (robustness)
        center = center.copy()
        center = np.maximum(center, mins)
        center = np.minimum(center, maxs)

    atom_types = autodetect_receptor_types(receptor_pdbqt)
    sig = grid_signature(receptor_pdbqt, mode, {
        "center": tuple(round(float(x),3) for x in center),
        "npts": tuple(int(x) for x in npts),
        "spacing": float(sp),
        "types": tuple(sorted(atom_types)),
    })
    stamp = os.path.join(out_dir, f".grid_{sig}.json")
    gpf_path = os.path.join(out_dir, "grid.gpf")

    if os.path.exists(stamp):
        # Maps already up-to-date for these exact params
        return {"center": center.tolist(), "npts": npts.tolist(), "spacing": sp, "types": sorted(atom_types), "signature": sig, "regenerated": False}

    # Clean stale stamps so cache is consistent
    for f in os.listdir(out_dir):
        if f.startswith(".grid_") and f.endswith(".json"):
            try: os.remove(os.path.join(out_dir, f))
            except: pass

    write_gpf(gpf_path, receptor_pdbqt, center, npts, sp, atom_types)
    run_autogrid4(gpf_path, workdir=out_dir, autogrid4_bin=autogrid4_bin)

    with open(stamp, "w") as s:
        json.dump({"center": center.tolist(), "npts": npts.tolist(), "spacing": sp, "types": sorted(atom_types)}, s, indent=2)

    return {"center": center.tolist(), "npts": npts.tolist(), "spacing": sp, "types": sorted(atom_types), "signature": sig, "regenerated": True}