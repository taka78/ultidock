#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualize pockets/cavities and the composite hotspot score as heatmaps.

Inputs
------
--seed_dir   dir with grid.gpf and *.map
--centers    centers.tsv   (receptor_key site_id cx cy cz nx ny nz spacing ...)
--receptor   receptor key/stem to filter rows (e.g., 5i6x-edited)
--axis       x|y|z (projection axis; default z)
--out        output PNG

Notes
-----
- Map ASCII is Fortran-ordered (AutoDock/AutoGrid).
- We plot in *world coordinates* using (gridcenter, npts, spacing).
"""

import argparse, os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

# ---------------------- IO helpers ----------------------

def parse_gpf(gpf_path):
    nx = ny = nz = None
    sp = None
    gc = None
    receptor_path = None
    with open(gpf_path, "r", errors="ignore") as f:
        for ln in f:
            t = ln.strip().split()
            if not t:
                continue
            k = t[0].lower()
            if k == "npts" and len(t) >= 4:
                nx, ny, nz = map(int, t[1:4])
            elif k == "spacing" and len(t) >= 2:
                sp = float(t[1])
            elif k == "gridcenter" and len(t) >= 4:
                gc = tuple(map(float, t[1:4]))
            elif k == "receptor" and len(t) >= 2:
                receptor_path = t[1]
    if None in (nx, ny, nz) or sp is None or gc is None:
        raise ValueError(f"GPF missing fields: {gpf_path}")
    return (nx, ny, nz), sp, gc, receptor_path

def read_map_ascii(map_path, gpf_path):
    npts, sp, gc, _ = parse_gpf(gpf_path)
    nx, ny, nz = npts
    vals = []
    with open(map_path, "r", errors="ignore") as f:
        for ln in f:
            for tok in ln.strip().split():
                try: vals.append(float(tok))
                except: pass
    arr = np.asarray(vals, dtype=np.float32)
    expect = nx*ny*nz
    if arr.size < expect:
        raise ValueError(f"{map_path}: not enough values ({arr.size} < {expect})")
    if arr.size > expect:
        arr = arr[-expect:]  # tolerate text headers
    vol = arr.reshape((nx, ny, nz), order="F")
    return vol, npts, sp, gc

def load_maps(seed_dir, stem):
    gpf = os.path.join(seed_dir, "grid.gpf")
    e_map = os.path.join(seed_dir, f"{stem}.e.map")
    d_map = os.path.join(seed_dir, f"{stem}.d.map")
    c_map = os.path.join(seed_dir, f"{stem}.C.map")
    for p in (gpf, e_map, d_map, c_map):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing: {p}")
    vol_e, npts, sp, gc = read_map_ascii(e_map, gpf)
    vol_d, _,   _,  _   = read_map_ascii(d_map, gpf)
    vol_c, _,   _,  _   = read_map_ascii(c_map, gpf)
    return (vol_e, vol_d, vol_c), npts, sp, gc, gpf

def parse_centers(centers_path, receptor_key):
    rows = []
    with open(centers_path, "r", errors="ignore") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            s = s.replace(",", " ")
            parts = s.split()
            if len(parts) >= 9 and parts[0] == receptor_key:
                _, site_id, cx, cy, cz, *_ = parts
                rows.append((site_id, float(cx), float(cy), float(cz)))
    return rows

# ---------------------- geometry / scoring ----------------------

def zscore(a):
    a = np.asarray(a, dtype=np.float32)
    m = np.nanmean(a); s = np.nanstd(a) + 1e-8
    return (a - m) / s

def build_heavy_coords_from_pdbqt(pdbqt_path):
    xyz = []
    if not os.path.exists(pdbqt_path):
        return np.empty((0,3), dtype=np.float32)
    with open(pdbqt_path, "r", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(("ATOM", "HETATM")): continue
            el = ln[76:78].strip().upper()
            if not el:
                name = ln[12:16].strip()
                el = ''.join(ch for ch in name if ch.isalpha())[:2].upper() or 'C'
            if el == 'H': 
                continue
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            xyz.append((x,y,z))
    if not xyz:
        raise ValueError(f"No heavy atoms parsed from {pdbqt_path}")
    return np.asarray(xyz, dtype=np.float32)

def edt_to_heavy(npts, sp, gc, heavy_xyz):
    nx, ny, nz = npts
    xs = gc[0] + (np.arange(nx) - (nx - 1)/2.0) * sp
    ys = gc[1] + (np.arange(ny) - (ny - 1)/2.0) * sp
    zs = gc[2] + (np.arange(nz) - (nz - 1)/2.0) * sp
    kdt = cKDTree(heavy_xyz)
    dist = np.empty((nx, ny, nz), dtype=np.float32)
    for i, x in enumerate(xs):
        P = np.stack(np.meshgrid([x], ys, zs, indexing="ij"), axis=-1).reshape(-1, 3)
        d, _ = kdt.query(P, k=1, workers=-1)
        dist[i, :, :] = d.reshape((ny, nz))
    return dist

def composite_score(vol_e, vol_d, vol_c, sp, sigma_A=1.2, wE=0.4, wD=0.2, wC=1.0):
    sigma_vox = max(0.1, float(sigma_A) / float(sp))
    e = gaussian_filter(vol_e, sigma=sigma_vox, mode="nearest")
    d = gaussian_filter(vol_d, sigma=sigma_vox, mode="nearest")
    c = gaussian_filter(vol_c, sigma=sigma_vox, mode="nearest")
    e = np.minimum(e, 0.0); d = np.minimum(d, 0.0); c = np.minimum(c, 0.0)
    return (wC * (-zscore(c))) + (wE * (-zscore(e))) + (wD * (-zscore(d)))

def rescale01(a):
    a = np.asarray(a, dtype=np.float32)
    lo = np.nanpercentile(a, 1.0); hi = np.nanpercentile(a, 99.0)
    return np.clip((a - lo) / (hi - lo + 1e-8), 0.0, 1.0)

# ---------------------- world extents + projection ----------------------

def world_bounds(npts, sp, gc):
    size = (np.array(npts) - 1) * sp
    xmin, xmax = gc[0] - size[0]/2, gc[0] + size[0]/2
    ymin, ymax = gc[1] - size[1]/2, gc[1] + size[1]/2
    zmin, zmax = gc[2] - size[2]/2, gc[2] + size[2]/2
    return (xmin, xmax, ymin, ymax, zmin, zmax)

def plane_extent(npts, sp, gc, axis='z'):
    xmin, xmax, ymin, ymax, zmin, zmax = world_bounds(npts, sp, gc)
    if axis == 'z': return [xmin, xmax, ymin, ymax]   # X vs Y
    if axis == 'y': return [xmin, xmax, zmin, zmax]   # X vs Z
    if axis == 'x': return [ymin, ymax, zmin, zmax]   # Y vs Z
    raise ValueError("axis must be x|y|z")

def sample_sphere(vol, center_xyz, npts, sp, gc, r=2.0):
    """mean value in a ~r Å sphere (grid coords)"""
    cx, cy, cz = center_xyz
    nx, ny, nz = npts
    xs = gc[0] + (np.arange(nx) - (nx - 1)/2.0) * sp
    ys = gc[1] + (np.arange(ny) - (ny - 1)/2.0) * sp
    zs = gc[2] + (np.arange(nz) - (nz - 1)/2.0) * sp
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    m = (X-cx)**2 + (Y-cy)**2 + (Z-cz)**2 <= r*r
    vals = vol[m]
    return float(np.nanmean(vals)) if vals.size else np.nan

def print_site_scores(rows, mE, mD, mC, npts, sp, gc):
    print("\n[sites]  mean within 2.0 Å sphere (AD maps)")
    print("site   E(map)   D(map)   C(map)   composite")
    for sid, cx, cy, cz in rows:
        e = sample_sphere(mE,(cx,cy,cz),npts,sp,gc)
        d = sample_sphere(mD,(cx,cy,cz),npts,sp,gc)
        c = sample_sphere(mC,(cx,cy,cz),npts,sp,gc)
        # same weights as the heatmap, negative is favorable before z-score
        comp = 1.0*(-c) + 0.4*(-e) + 0.2*(-d)
        print(f"{sid:>3s}  {e:7.3f}  {d:7.3f}  {c:7.3f}   {comp:8.3f}")

def mip(vol, axis='z', how='max'):
    ax = {'x':0,'y':1,'z':2}[axis.lower()]
    return (np.nanmax if how == 'max' else np.nanmean)(vol, axis=ax)

def project_points(pts_xyz, axis='z'):
    out = []
    for (x,y,z) in pts_xyz:
        if axis == 'z': out.append((x, y))
        elif axis == 'y': out.append((x, z))
        else: out.append((y, z))  # axis == 'x'
    return np.array(out, float) if pts_xyz else np.zeros((0,2), float)

# ---------------------- main ----------------------

def main():
    ap = argparse.ArgumentParser(description="Cavity & hotspot heatmap visualizer")
    ap.add_argument("--seed_dir", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--centers", required=True)
    ap.add_argument("--axis", default="z", choices=list("xyz"))
    ap.add_argument("--out", default="hotspot_heatmaps.png")
    ap.add_argument("--sigmaA", type=float, default=1.2)
    ap.add_argument("--wE", type=float, default=0.4)
    ap.add_argument("--wD", type=float, default=0.2)
    ap.add_argument("--wC", type=float, default=1.0)
    ap.add_argument("--shell_min", type=float, default=1.8)
    ap.add_argument("--shell_max", type=float, default=6.0)
    ap.add_argument("--receptor_pdbqt", default=None)
    args = ap.parse_args()

    stem = args.receptor
    (mE, mD, mC), npts, sp, gc, gpf = load_maps(args.seed_dir, stem)

    # receptor for EDT: CLI > receptor path in GPF > seed_dir/<stem>.pdbqt
    _, _, _, rec_from_gpf = parse_gpf(gpf)
    rec_candidates = [args.receptor_pdbqt, rec_from_gpf, os.path.join(args.seed_dir, f"{stem}.pdbqt")]
    rec_pdbqt = next((p for p in rec_candidates if p and os.path.exists(p)), None)
    if not rec_pdbqt:
        raise FileNotFoundError("Receptor PDBQT not found (pass --receptor_pdbqt or ensure it’s in GPF/seed_dir).")

    heavy = build_heavy_coords_from_pdbqt(rec_pdbqt)
    dist  = edt_to_heavy(npts, sp, gc, heavy)

    shell = (dist >= float(args.shell_min)) & (dist <= float(args.shell_max))
    if not np.any(shell):
        print("[warn] cavity shell is empty with current shell_min/max; consider relaxing thresholds.")
    cavity = np.where(shell, dist, np.nan)

    score     = composite_score(mE, mD, mC, sp, sigma_A=args.sigmaA, wE=args.wE, wD=args.wD, wC=args.wC)
    heat_cmp  = rescale01(score)
    heat_e    = rescale01(-zscore(np.minimum(gaussian_filter(mE, max(0.1, args.sigmaA/sp)), 0.0)))
    heat_c    = rescale01(-zscore(np.minimum(gaussian_filter(mC, max(0.1, args.sigmaA/sp)), 0.0)))
    cav_heat  = rescale01(np.nan_to_num(cavity, nan=0.0))

    axis = args.axis.lower()
    Hcav = mip(cav_heat,  axis=axis, how="max")
    Hcmp = mip(heat_cmp,  axis=axis, how="max")
    He   = mip(heat_e,    axis=axis, how="max")
    Hc   = mip(heat_c,    axis=axis, how="max")
    ext  = plane_extent(npts, sp, gc, axis=axis)

    # hotspots
    rows  = parse_centers(args.centers, stem)
    print_site_scores(rows, mE, mD, mC, npts, sp, gc)
    pts3d = [(cx,cy,cz) for (_,cx,cy,cz) in rows]
    pts2d = project_points(pts3d, axis=axis)

    # plot
    fig, axs = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    panels = [
        (axs[0,0], Hcav, "Cavityness (EDT within shell)"),
        (axs[0,1], Hcmp, "Composite heat (e+d+C)"),
        (axs[1,0], He,   "E-only (favorable)"),
        (axs[1,1], Hc,   "C-only (favorable)"),
    ]
    for ax, H, title in panels:
        im = ax.imshow(H.T, origin="lower", extent=ext, aspect="equal")
        if pts2d.size:
            ax.scatter(pts2d[:,0], pts2d[:,1], s=30, c="white",
                       edgecolors="black", linewidths=0.8, zorder=3)
            for (site_id, x, y, z), (px, py) in zip(rows, pts2d):
                ax.text(px, py, site_id, color="white", fontsize=8, weight="bold",
                        ha="left", va="bottom", zorder=4)
        ax.set_title(f"{title}  (axis={axis})")
        ax.set_xlabel("Å"); ax.set_ylabel("Å")
        fig.colorbar(im, ax=ax, shrink=0.85)

    fig.suptitle(f"{stem} — pocket/cavity heatmaps (spacing={sp:.3f} Å)", y=1.02, fontsize=14)
    out = Path(args.out)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"[viz] wrote {out.resolve()}")

if __name__ == "__main__":
    main()
