import numpy as np
from scipy.ndimage import binary_closing, distance_transform_edt, label, generate_binary_structure
from collections import namedtuple
import re, os, argparse
from pathlib import Path
from scipy.ndimage import gaussian_filter, maximum_filter
import mmap
from grid_box import blind_box
from config import AUTODOCK_GPU_DIR, DOCKING_DIR, RESULTS_DIR, NUMWI, LIGANDS_DIR, CENTERS_TSV, GRID_MODE, GRID_SPACING, GRID_MARGIN, GRID_CAP, R_MIN_CAVITY_A, HOTSPOT_BOX_ANGLE

spacing = float(GRID_SPACING)

# ===== Multi-center wiring for dock_v02.py =====
from pathlib import Path
import subprocess


def _is_float_token(tok: str) -> bool:
    try:
        float(tok)
        return True
    except Exception:
        return False


def _parse_centers_tsv(centers_tsv_path, receptor_key=None, default_spacing=GRID_SPACING, default_npts=(81, 81, 81)):
    """
    Parses centers.tsv supporting:
      receptor site_id cx cy cz nx ny nz spacing [r_peak F ...]
    or
      legacy: key cx cy cz [nx ny nz spacing]

    Returns list of dicts, one per row (optionally filtered by receptor_key).
    """
    rows = []
    with open(centers_tsv_path, "r", errors="ignore") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if receptor_key and parts[0] != receptor_key:
                continue
            has_site = len(parts) >= 2 and (not _is_float_token(parts[1]))
            i0 = 2 if has_site else 1
            if len(parts) < i0 + 3:
                continue
            cx, cy, cz = map(float, parts[i0:i0 + 3])
            if len(parts) >= i0 + 6:
                try:
                    nx, ny, nz = map(int, parts[i0 + 3:i0 + 6])
                except Exception:
                    nx, ny, nz = default_npts
            else:
                nx, ny, nz = default_npts
            if len(parts) >= i0 + 7 and _is_float_token(parts[i0 + 6]):
                sp = float(parts[i0 + 6])
            else:
                sp = float(default_spacing)
            sd = {
                "receptor": parts[0],
                "site_id": parts[1] if has_site else f"S{len(rows) + 1}",
                "center": (cx, cy, cz),
                "npts": (nx, ny, nz),
                "spacing": sp,
            }
            rows.append(sd)
    return rows


def _ensure_odd_clamped(nxyz, clamp=(60, 255)):
    import numpy as _np
    n = _np.array(nxyz, int)
    n = n + (n % 2 == 0)
    n = _np.clip(n, clamp[0], clamp[1])
    return tuple(int(x) for x in n.tolist())


def _write_gpf(gpf_path, receptor_pdbqt, center, npts, spacing):
    """Write a per-site GPF; **do not** modify npts here."""
    cx, cy, cz = map(float, center)  # correct: center → (cx,cy,cz)
    nx, ny, nz = (int(npts[0]), int(npts[1]), int(npts[2]))  # correct: npts → (nx,ny,nz)
    rec_path = Path(receptor_pdbqt).resolve()
    rec_stem = rec_path.stem
    fld_name = f"{rec_stem}.maps.fld"
    with open(gpf_path, "w") as f:
        # some builds are picky: put gridfld before gridcenter
        f.write(f"gridfld {fld_name}\n")
        f.write(f"npts {nx} {ny} {nz}\n")
        f.write(f"spacing {float(spacing):.3f}\n")
        f.write(f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}\n")
        f.write(f"receptor {rec_path}\n")
        f.write("receptor_types " + " ".join(_AD4_TYPES) + "\n")
        f.write("ligand_types " + " ".join(_AD4_TYPES) + "\n")
        for t in _AD4_TYPES:
            f.write(f"map {rec_stem}.{t}.map\n")
        f.write(f"elecmap {rec_stem}.e.map\n")
        f.write(f"dsolvmap {rec_stem}.d.map\n")
    return str(gpf_path)


def autogenerate_centers_tsv(
    receptor_pdbqt: str,
    out_root: str,
    centers_tsv_path: str,
    n_sites: int = 6,
    default_npts=(81, 81, 81),
    default_spacing: float = GRID_SPACING,
    blind_cap: float = 50.0,
    autogrid4_bin: str = "autogrid4",
    hotspot_box_ang: float = 35.0,
    hotspot_sigma_A: float = 1.0,  # (kept for API compatibility; not used directly below)
    hotspot_bury_z: float = 0.35,  # (ditto)
    mode: str = "hybrid",
):
    """
    Generate centers.tsv under out_root using already-present maps or by making whole-protein maps.
    Reuse existing centers.tsv if it already has rows for this receptor.
    Fallback order: internal cavities → maps hotspots → looser maps → blind box.
    """
    from pathlib import Path

    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    centers_tsv_path = Path(centers_tsv_path)
    rec_stem = Path(receptor_pdbqt).stem

    # 0) Reuse if present
    if centers_tsv_path.exists():
        rows = _parse_centers_tsv(centers_tsv_path, receptor_key=rec_stem)
        if rows:
            print(f"[centers] reusing {centers_tsv_path} with {len(rows)} rows")
            return str(centers_tsv_path)

        # 1) Find (or build) whole-protein maps (.fld)
    fld_candidates = list(out_root.glob("**/*.fld"))
    if fld_candidates:
        def grid_volume(f):
            try:
                m = load_fld_or_map_meta(str(f))
                sh = m["shape"]
                return int(sh[0] * sh[1] * sh[2])
            except Exception:
                return -1
        fld_candidates.sort(key=grid_volume, reverse=True)
        fld = fld_candidates[0]
    else:
        # bootstrap: generate whole-protein maps once
        fld_path = ensure_whole_protein_maps(
            receptor_pdbqt=receptor_pdbqt,
            out_root=out_root,
            spacing=default_spacing,
            cap_ang=blind_cap,
            autogrid4_bin=autogrid4_bin,
        )
        fld_candidates = [Path(fld_path)]

    fld = fld_candidates[0]
    meta = load_fld_or_map_meta(str(fld))
    origin, spacing, shape = meta["origin"], meta["spacing"], meta["shape"]
    mp = meta["map_paths"]
    C = load_map_ascii(mp["C"], shape)
    E = load_map_ascii(mp["E"], shape)
    D = load_map_ascii(mp["D"], shape)

    sites = []

    # 2) Internal cavities (EDT) first if requested
    if mode in ("internal", "hybrid"):
        try:
            sites = pick_centers(
                receptor_pdbqt,
                {"C": C, "E": E, "D": D},
                map_origin=origin,
                map_spacing=spacing,
                voxel_spacing=0.375,
                r_min=R_MIN_CAVITY_A,  # more permissive than 1.8
                min_sep_A=5.0,
                k_box=4.0,
                max_sites=n_sites,
                inflate_A=0.0,
            )
        except Exception as e:
            print(f"[centers/internal] skipped due to error: {e}")

    # 3) Maps hotspots if hybrid or maps-only, when internal found none
    if mode in ("maps", "hybrid") and len(sites) == 0:
        print("[centers] using maps-mode hotspots")
        sites = detect_maps_hotspots(
            C, E, D, origin, spacing, tau_rel=0.52, min_sep_A=5.0, max_sites=n_sites, half_size_A=hotspot_box_ang
        )

    # 4) Loosen and retry maps if still empty
    if len(sites) == 0:
        print("[centers] no sites found; retrying maps-mode with looser params")
        # We already ensured whole-protein maps above; reuse meta/maps
        sites = detect_maps_hotspots(
            C, E, D, origin, spacing, tau_rel=0.52,  # looser
            min_sep_A=5.0,
            max_sites=n_sites,
            half_size_A=hotspot_box_ang,
        )

    # 5) Final fallback: blind box around the protein, so pipeline never dies
    if len(sites) == 0:
        print("[centers] maps still empty; falling back to single blind box")
        center, npts, sp = blind_box(receptor_pdbqt, cap=blind_cap, spacing=default_spacing)
        sites = [
            dict(
                site_id="S1",
                cx=center[0],
                cy=center[1],
                cz=center[2],
                nx=npts[0],
                ny=npts[1],
                nz=npts[2],
                spacing=sp,
                r_peak=float(blind_cap) / 4.0,
                F=0.0,
            )
        ]

    # 6) Write TSV
    with open(centers_tsv_path, "w") as f:
        f.write("# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tr_peak\tF\n")
        for s in sites:
            f.write(
                f"{rec_stem}\t{s['site_id']}\t{s['cx']:.3f}\t{s['cy']:.3f}\t{s['cz']:.3f}\t"
                f"{s['nx']}\t{s['ny']}\t{s['nz']}\t{s['spacing']:.3f}\t{s.get('r_peak', 2.5):.2f}\t{s.get('F', 1.0):.3f}\n"
            )
    print(f"[centers] wrote {len(sites)} → {centers_tsv_path}")
    return str(centers_tsv_path)


def ensure_grids_multi_centers(
    receptor_pdbqt: str,
    out_root: str,
    centers_tsv: str,
    autogrid4_bin: str = "autogrid4",
):
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rec_stem = Path(receptor_pdbqt).stem

    rows = _parse_centers_tsv(
        centers_tsv,
        receptor_key=rec_stem,
        default_spacing=GRID_SPACING,
        default_npts=(81, 81, 81),
    )
    if not rows:
        raise ValueError(f"No centers in {centers_tsv} for receptor {rec_stem}")

    sites = []
    for r in rows:
        site_id  = r["site_id"]
        center   = r["center"]
        npts     = r["npts"]
        spacing  = r["spacing"]

        site_dir = out_root / site_id
        site_dir.mkdir(parents=True, exist_ok=True)

        gpf_path = site_dir / "grid.gpf"
        _write_gpf(gpf_path, receptor_pdbqt, center, npts, spacing)

        log_path = site_dir / "grid.glg"
        cmd = [autogrid4_bin, "-p", gpf_path.name, "-l", log_path.name]
        print(f"[autogrid] {site_id}: running in {site_dir}")

        res = subprocess.run(cmd, cwd=site_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print(f"[autogrid/{site_id}] FAILED\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
            raise RuntimeError(f"autogrid failed for {site_id}")

        # Prefer the exact expected fld name first (deterministic)
        expected_fld = site_dir / f"{rec_stem}.maps.fld"
        if expected_fld.exists():
            fld_path = str(expected_fld.resolve())
        else:
            fld_candidates = list(site_dir.glob("*.fld"))
            if not fld_candidates:
                fld_candidates = list(site_dir.glob(f"{rec_stem}*.fld"))
            if not fld_candidates:
                raise FileNotFoundError(f"No .fld produced in {site_dir}; see {log_path}")
            fld_path = str(fld_candidates[0].resolve())

        # NOW meta is safe
        meta = load_fld_or_map_meta(fld_path)
        origin = meta["origin"]
        sp = meta["spacing"]
        shape = meta["shape"]

        # if you want your clash distance grid:
        atoms = pdbqt_atoms(receptor_pdbqt)
        occ, distA, origin, sp = make_grid_from_map(
            atoms,
            map_origin=origin,
            map_spacing=sp,
            map_shape=shape,
            inflate_A=0.25,
        )
        # occ must be a boolean numpy array here
        if isinstance(occ, tuple):
            # common bug: occ = make_grid_from_map(...) without unpacking
            # or occ = (occ_grid, origin, spacing, ...)
            # pick the first ndarray inside the tuple
            for x in occ:
                if isinstance(x, np.ndarray):
                    occ = x
                    break
            else:
                raise TypeError(f"occ is tuple with no ndarray inside: {[type(x) for x in occ]}")

        occ = np.asarray(occ)
        if occ.dtype != np.bool_:
            # IMPORTANT: distance_transform_edt(~occ) only makes sense if occ is boolean occupancy
            occ = occ.astype(bool, copy=False)

        distA = distance_transform_edt(~occ) * float(sp)

        distA = distance_transform_edt(~occ) * float(sp)
        save_clash_dist_grid(site_dir / "clash_dist.npz", distA, origin, sp)


        sites.append({
            "site_id": site_id,
            "center": center,
            "npts": tuple(map(int, npts)),
            "spacing": float(spacing),
            "fld_path": fld_path,
            "out_dir": str(site_dir.resolve()),
        })

    print(f"[autogrid] prepared {len(sites)} site grids")
    return sites


def save_clash_dist_grid(out_path, distA, origin, spacing):
    """
    Save distance-to-receptor grid aligned with AutoGrid maps.
    distA: (nx,ny,nz) float32 in Å
    """
    out_path = Path(out_path)
    np.save(out_path, distA.astype(np.float32))
    meta_path = out_path.with_suffix(".meta.txt")
    with open(meta_path, "w") as f:
        f.write(f"origin {origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}\n")
        f.write(f"spacing {float(spacing):.6f}\n")
        f.write(f"shape {distA.shape[0]} {distA.shape[1]} {distA.shape[2]}\n")


def _robust_z(x):
    x = np.asarray(x, dtype=np.float32)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med)) + 1e-6
    return (x - med) / mad


def _sigmoid_stable(x):
    # clip to avoid overflow; equivalent to logistic for our purposes
    x = np.clip(x, -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-x))


def detect_maps_hotspots(C, E, D, origin, spacing, tau_rel=0.52, min_sep_A=5.0, max_sites=6, half_size_A=18.0):
    """Surface-pocket finder from C/E/D maps with stable math and matched masks."""
    C = np.asarray(C, dtype=np.float32)
    E = np.asarray(E, dtype=np.float32)
    D = np.asarray(D, dtype=np.float32)

    # robust-normalize C and E
    Cn, En = _robust_z(C), _robust_z(E)

    # scale D into 0..1 using headroom
    dmin = float(np.nanmin(D))
    d95 = float(np.nanpercentile(D, 95))
    Dn = np.clip((D - dmin) / (d95 - dmin + 1e-6), 0.0, 1.0)

    # squashes → scores (stable)
    sC = _sigmoid_stable(Cn)
    sE = _sigmoid_stable(+0.7 * (-En))  # more negative E → higher
    sD = Dn
    F = 0.45 * sC + 0.25 * sE + 0.30 * sD
    F = gaussian_filter(F, sigma=1.0)

    # pocket-ish mask
    mask = (Dn > 0.35) & (Cn > -0.5)

    # apply mask by filling -inf outside; this keeps shapes aligned
    Fm = np.full_like(F, -np.inf, dtype=np.float32)
    Fm[mask] = F[mask]

    # local maxima on masked field
    Fmax = maximum_filter(Fm, size=5, mode="nearest")

    # threshold relative to best finite value INSIDE the mask
    finite = Fm[np.isfinite(Fm)]
    if finite.size == 0:
        return []
    thr = float(finite.max()) * float(tau_rel)

    # IMPORTANT: use the SAME mask for peaks and scores
    peakmask = (Fm == Fmax) & (Fm > thr)
    if not np.any(peakmask):
        return []

    peaks = np.argwhere(peakmask)  # N x 3
    scores = Fm[peakmask].astype(np.float32)  # length N
    order = np.argsort(-scores, kind="stable")
    peaks = peaks[order]
    scores = scores[order]

    # NMS in Å
    ox, oy, oz = origin
    sp = float(spacing)

    def vox2world(ijk):
        i, j, k = ijk
        return np.array([ox + i * sp, oy + j * sp, oz + k * sp], dtype=np.float32)

    kept_xyz, sites = [], []
    for ijk, sc in zip(peaks, scores):
        wp = vox2world(ijk)
        if all(np.linalg.norm(wp - q) >= min_sep_A for q in kept_xyz):
            kept_xyz.append(wp)
            half = float(half_size_A)
            n = int(np.ceil(2 * half / sp))
            if n % 2 == 0:
                n += 1
            n = max(60, min(255, n))
            sites.append(dict(
                site_id=f"S{len(kept_xyz)}",
                cx=float(wp[0]),
                cy=float(wp[1]),
                cz=float(wp[2]),
                nx=n, ny=n, nz=n,
                spacing=sp,
                r_peak=half / 4.0,  # placeholder for maps-mode
                F=float(sc)
            ))
            if len(sites) == int(max_sites):
                break
    return sites


def pdbqt_atoms(path):
    atoms = []
    with open(path, "r") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        for line in iter(mm.readline, b""):
            if line.startswith((b"ATOM", b"HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                t = line[77:79].strip().decode()
                atoms.append((x, y, z, t))
    return atoms


AD4_RADII = {
    # minimal table; extend as needed
    "C": 2.0, "A": 2.0, "HD": 1.0, "N": 1.9, "NA": 1.9, "OA": 1.7, "SA": 2.0, "S": 2.0, "P": 2.1,
    "Cl": 2.0, "F": 1.8, "Br": 2.2, "I": 2.3, "Zn": 1.2
}


def make_grid_from_map(atoms, map_origin, map_spacing, map_shape, inflate_A=0.25):
    """
    Build occupancy aligned to the map grid. 'inflate_A' thickens each atom radius by that many Å
    to make the envelope more conservative (use 0.0–0.5 typically).
    """
    ox, oy, oz = map_origin
    sp = float(map_spacing)
    nx, ny, nz = map_shape
    occ = np.zeros((nx, ny, nz), dtype=bool)
    gx = ox + np.arange(nx) * sp
    gy = oy + np.arange(ny) * sp
    gz = oz + np.arange(nz) * sp

    for (x, y, z, t) in atoms:
        r = AD4_RADII.get(t, 2.0) + float(inflate_A)
        ix0 = max(0, int(np.floor((x - r - ox) / sp)))
        ix1 = min(nx - 1, int(np.ceil((x + r - ox) / sp)))
        iy0 = max(0, int(np.floor((y - r - oy) / sp)))
        iy1 = min(ny - 1, int(np.ceil((y + r - oy) / sp)))
        iz0 = max(0, int(np.floor((z - r - oz) / sp)))
        iz1 = min(nz - 1, int(np.ceil((z + r - oz) / sp)))
        if ix0 > ix1 or iy0 > iy1 or iz0 > iz1:
            continue
        dx2 = (gx[ix0:ix1 + 1] - x) ** 2
        dy2 = (gy[iy0:iy1 + 1] - y) ** 2
        dz2 = (gz[iz0:iz1 + 1] - z) ** 2
        mask = (dx2[:, None, None] + dy2[None, :, None] + dz2[None, None, :]) <= (r * r)
        occ[ix0:ix1 + 1, iy0:iy1 + 1, iz0:iz1 + 1] |= mask

    occ = binary_closing(occ, structure=np.ones((3, 3, 3), bool))
    dist_vox = distance_transform_edt(~occ)
    dist_A = dist_vox * sp
    return occ, dist_A, (ox, oy, oz), sp


def internal_cavities(occ):
    free = ~occ

    # exterior flood from borders
    ext = np.zeros_like(free, bool)
    ext[[0, -1], :, :] = True
    ext[:, [0, -1], :] = True
    ext[:, :, [0, -1]] = True

    from scipy.ndimage import binary_dilation
    prev = np.zeros_like(free, bool)
    cur = ext & free
    while True:
        nxt = (binary_dilation(cur, structure=np.ones((3, 3, 3))) & free) | cur | ext
        if np.array_equal(nxt, cur):
            break
        cur = nxt
    exterior = cur
    internal = free & (~exterior)

    # distance (in voxels)
    dist = distance_transform_edt(internal)

    # components
    cc, ncc = label(internal, structure=generate_binary_structure(3, 2))
    return dist, cc, ncc


#### MAP SAMPLING AND SCORING ####
def trilinear_sample(vol, xyz, origin, sp):
    x, y, z = xyz
    ox, oy, oz = origin
    fx = (x - ox) / sp
    fy = (y - oy) / sp
    fz = (z - oz) / sp
    i, j, k = np.floor([fx, fy, fz]).astype(int)
    dx, dy, dz = fx - i, fy - j, fz - k
    nx, ny, nz = vol.shape
    if not (0 <= i < nx - 1 and 0 <= j < ny - 1 and 0 <= k < nz - 1):
        return np.nan
    c000 = vol[i, j, k]
    c100 = vol[i + 1, j, k]
    c010 = vol[i, j + 1, k]
    c110 = vol[i + 1, j + 1, k]
    c001 = vol[i, j, k + 1]
    c101 = vol[i + 1, j, k + 1]
    c011 = vol[i, j + 1, k + 1]
    c111 = vol[i + 1, j + 1, k + 1]
    c00 = c000 * (1 - dx) + c100 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy
    return c0 * (1 - dz) + c1 * dz


def robust_norm(x):
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med)) + 1e-6
    return (x - med) / mad


def site_score(E, C, D, origin, sp, center):
    # sample small ball statistics
    valsE, valsC, valsD = [], [], []
    # neighborhood: 2-voxel radius
    for dx in range(-2, 3):
        for dy in range(-2, 3):
            for dz in range(-2, 3):
                pt = (center[0] + dx * sp, center[1] + dy * sp, center[2] + dz * sp)
                valsE.append(trilinear_sample(E, pt, origin, sp))
                valsC.append(trilinear_sample(C, pt, origin, sp))
                valsD.append(trilinear_sample(D, pt, origin, sp))
    vE = np.array(valsE)
    vC = np.array(valsC)
    vD = np.array(valsD)

    # robust z
    Ez = robust_norm(vE)
    Cz = robust_norm(vC)

    # scale D to 0..1
    dmin, d95 = np.nanmin(vD), np.nanpercentile(vD, 95)
    Dn = np.clip((vD - dmin) / (d95 - dmin + 1e-6), 0, 1)

    # squashes
    sC = 1.0 / (1.0 + np.exp(-Cz))
    sE = 1.0 / (1.0 + np.exp(+0.7 * Ez * -1.0))  # more negative E -> higher
    sD = Dn

    # combine
    F = 0.45 * np.nanmedian(sC) + 0.30 * np.nanmedian(sD) + 0.25 * np.nanmedian(sE)
    return float(F), float(np.nanmedian(sC)), float(np.nanmedian(sD)), float(np.nanmedian(sE))


### CENTERING AND BOXING ###
def voxel_to_world(ijk, origin, sp):
    return origin[0] + ijk[0] * sp, origin[1] + ijk[1] * sp, origin[2] + ijk[2] * sp


def pick_centers(
    pdbqt,
    maps,
    map_origin,
    map_spacing,
    voxel_spacing=0.375,
    r_min=R_MIN_CAVITY_A,
    min_sep_A=5.0,
    k_box=4.0,
    max_sites=6,
    inflate_A=0.0
):
    atoms = pdbqt_atoms(pdbqt)

    # align geometry grid to the map grid; use inflation
    occ, distA, (ox, oy, oz), sp = make_grid_from_map(
        atoms, map_origin, map_spacing, maps["C"].shape, inflate_A=inflate_A
    )
    print(f"[grid] map shape={maps['C'].shape}, spacing={map_spacing:.3f} Å")
    print(f"[grid] occ voxels={int(occ.sum())}, free voxels={int((~occ).sum())}")

    dist, cc, ncc = internal_cavities(occ)
    print(f"[cav] internal components={ncc}, max_r_peak_vox={float(dist.max()):.2f} (Å)={float(dist.max() * sp):.2f}")

    candidates = []
    for lab in range(1, ncc + 1):
        mask = (cc == lab)
        if not mask.any():
            continue
        # peak voxel
        i, j, k = np.unravel_index(np.argmax(dist * mask), dist.shape)
        r_peak = dist[i, j, k] * sp
        if r_peak < r_min:
            continue
        cx, cy, cz = voxel_to_world((i, j, k), (ox, oy, oz), sp)
        F, sC, sD, sE = site_score(maps["E"], maps["C"], maps["D"], map_origin, map_spacing, (cx, cy, cz))
        # simple gates
        # if sD < 0.45 or sC < 0.55: continue ##################fix these thresholds
        candidates.append((np.array([cx, cy, cz]), r_peak, F))

    # NMS
    candidates.sort(key=lambda x: -x[2])
    kept = []
    for c in candidates:
        if all(np.linalg.norm(c[0] - k[0]) >= min_sep_A for k in kept):
            kept.append(c)
        if len(kept) == max_sites:
            break

    # boxes
    sites = []
    for idx, (ctr, r_peak, F) in enumerate(kept, 1):
        half = max(map_spacing, k_box * r_peak)
        n = int(np.ceil(2 * half / map_spacing))
        if n % 2 == 0:
            n += 1
        n = max(60, min(255, n))
        sites.append(dict(
            site_id=f"S{idx}",
            cx=ctr[0], cy=ctr[1], cz=ctr[2],
            nx=n, ny=n, nz=n,
            spacing=map_spacing,
            r_peak=r_peak,
            F=F
        ))
    return sites


# ---------- GLUE: FLD/MAP loaders, TSV writer, CLI ----------
_FLOAT_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")


def load_fld_or_map_meta(fld_or_dir):
    """
    Parse AutoGrid AVS .fld:
    - dims: dim1, dim2, dim3 (these are npts+1)
    - #NELEMENTS nx ny nz (true npts)
    - #SPACING s
    - #CENTER cx cy cz
    - variable i file=...map + label lines

    Returns:
    {
      'origin': (ox,oy,oz),  # lower-left-back corner in Å
      'spacing': spacing,    # Å
      'shape': (nx,ny,nz),   # for reshaping .map arrays (Fortran order)
      'dir': Path,           # directory of the fld
      'map_paths': {'C': Path, 'E': Path, 'D': Path}
    }
    """
    import re
    from pathlib import Path

    p = Path(fld_or_dir)
    if not p.is_file():
        # try to find a .fld inside the directory
        cands = list(Path(fld_or_dir).glob("*.fld"))
        if not cands:
            raise FileNotFoundError(f"No .fld in {fld_or_dir}")
        p = cands[0]
    d = p.parent
    txt = open(p, "r", encoding="utf-8", errors="ignore").read()

    # pull commented metadata first (preferred)
    def floats_after(tag, n=None):
        m = re.search(rf"^{re.escape(tag)}\s+(.+)$", txt, flags=re.M | re.I)
        if not m:
            return None
        vals = [float(t) for t in m.group(1).split()]
        if n and len(vals) < n:
            return None
        return vals

    sp_list = floats_after("#SPACING", 1)
    cen_list = floats_after("#CENTER", 3)
    nelem_lst = floats_after("#NELEMENTS", 3)
    spacing = sp_list[0] if sp_list else None
    center = tuple(cen_list) if cen_list else None
    nelems = tuple(int(v) for v in nelem_lst) if nelem_lst else None

    # parse dim1/2/3 (= npts+1)
    def int_after(k):
        m = re.search(rf"^{k}\s*=\s*([0-9]+)", txt, flags=re.M | re.I)
        return int(m.group(1)) if m else None

    dim1 = int_after("dim1")
    dim2 = int_after("dim2")
    dim3 = int_after("dim3")
    if not (dim1 and dim2 and dim3):
        raise RuntimeError("Could not parse dim1/dim2/dim3 from .fld")

    shape = (dim1, dim2, dim3)

    # compute origin from center + spacing + npts (true npts, NOT dims)
    if spacing is None or center is None:
        raise RuntimeError("Missing #SPACING or #CENTER in .fld header comments.")
    if nelems is None:
        # fallback: infer npts = dims - 1
        nelems = (dim1 - 1, dim2 - 1, dim3 - 1)

    # after parsing 'nelems' (nx,ny,nz), 'center', and 'spacing sp'
    # ...
    # shape must be NELEMENTS (true npts), not dim*
    shape = nelems  # (nx, ny, nz)

    # correct origin: center ± (npts-1)*spacing/2
    nx, ny, nz = nelems
    cx, cy, cz = center
    ox = cx - ((nx - 1) * spacing) / 2.0
    oy = cy - ((ny - 1) * spacing) / 2.0
    oz = cz - ((nz - 1) * spacing) / 2.0
    origin = (ox, oy, oz)

    # map variable -> label + file
    # example lines:
    # label=Electrostatics
    # variable 18 file=5i6x-edited.e.map filetype=ascii skip=6
    # #
    # We’ll link files by looking for label lines near variable index.
    var_file = {}
    for m in re.finditer(r"variable\s+(\d+)\s+file=(\S+)", txt, flags=re.I):
        idx = int(m.group(1))
        path = (d / m.group(2)).resolve()
        var_file[idx] = path

    # collect labels (in order of appearance)
    labels = []
    for m in re.finditer(r"label\s*=\s*(.+)", txt, flags=re.I):
        labels.append(m.group(1).strip())  # e.g., "C-affinity", "Electrostatics", "Desolvation"

    # Build name->file mapping using typical AutoGrid ordering:
    # veclen=19; variables 1..19 correspond to labels 1..19
    # safeguard: only take those indices that exist in var_file
    name_to_path = {}
    for idx, lab in enumerate(labels, start=1):
        if idx in var_file:
            name_to_path[lab] = var_file[idx]

    # Resolve C/E/D paths
    # - Contact/“C-affinity” → label contains "C-affinity" (but not "Cl-affinity")
    # - Electrostatics → label == "Electrostatics"
    # - Desolvation → label == "Desolvation"
    C_path = None  # pick exact "C-affinity" (not Cl, Ca…); use word boundary trick
    for k, v in name_to_path.items():
        if k.lower() == "c-affinity":
            C_path = v
            break
    if C_path is None:
        # fallback: first label ending with "-affinity" starting with "C-"
        for k, v in name_to_path.items():
            if k.lower().startswith("c-affinity"):
                C_path = v
                break
    E_path = name_to_path.get("Electrostatics") or name_to_path.get("electrostatics")
    D_path = name_to_path.get("Desolvation") or name_to_path.get("desolvation")

    if not C_path:
        # fallback glob (avoid Cl/Br etc.)
        cands = sorted(d.glob("*C.map"))
        if cands:
            C_path = cands[0]
    if not E_path:
        # often lowercase 'e.map'
        cands = sorted(list(d.glob("*E.map")) + list(d.glob("*.e.map")))
        if cands:
            E_path = cands[0]
    if not D_path:
        cands = sorted(list(d.glob("*D.map")) + list(d.glob("*.d.map")))
        if cands:
            D_path = cands[0]

    missing = [n for n, p in {"C": C_path, "E": E_path, "D": D_path}.items() if p is None]
    if missing:
        raise FileNotFoundError(f"Could not resolve map files for: {', '.join(missing)}")

    return {
        "origin": origin,
        "spacing": spacing,
        "shape": shape,  # use (dim1,dim2,dim3) when reshaping maps
        "dir": d,
        "map_paths": {"C": C_path, "E": E_path, "D": D_path},
    }


def load_map_ascii(map_path, shape):
    """
    Robust ASCII loader for AutoGrid .map (x fastest).
    Skips non-numeric header lines.
    Reshapes with Fortran order so that index [i,j,k] corresponds to x,y,z (i fastest).
    """
    vals = []
    with open(map_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            # skip header lines (start with letters)
            c0 = s[0]
            if c0 not in "0123456789-+.":
                continue
            # append floats on this line
            for tok in s.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
    arr = np.asarray(vals, dtype=np.float32)
    need = int(np.prod(shape))
    if arr.size > need:
        arr = arr[-need:]  # keep trailing data block
    arr = arr[:need].reshape(shape, order="F")  # x fastest
    return arr


def write_centers_tsv(path, receptor_stem, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tr_peak\tF\n")
        for r in rows:
            f.write(
                f"{receptor_stem}\t{r['site_id']}\t{r['cx']:.3f}\t{r['cy']:.3f}\t{r['cz']:.3f}\t"
                f"{r['nx']}\t{r['ny']}\t{r['nz']}\t{r['spacing']:.3f}\t{r['r_peak']:.2f}\t{r['F']:.3f}\n"
            )
    print(f"[OK] wrote {len(rows)} sites → {path}")


# ===== Single-site grid builder (ligand / residues / blind) =====
import math


def _odd(n: int) -> int:
    return n if (n % 2 == 1) else (n + 1)


def _npts_for_size(side_len_A: float, spacing: float, clamp=(60, 255)) -> int:
    n = int(math.ceil(side_len_A - 1 / float(spacing)))
    n = _odd(n)
    return max(clamp[0], min(clamp[1], n))


def _load_coords_pdb_like(path, atom_filter=lambda name: True, line_filter=None):
    coords = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line_filter and (not line_filter(line)):
                continue
            name = line[12:16].strip()
            if not atom_filter(name):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except Exception:
                pass
    if not coords:
        raise ValueError(f"No atom coordinates parsed from: {path}")
    arr = np.asarray(coords, dtype=float)
    mins = arr.min(axis=0)
    maxs = arr.max(axis=0)
    center = (maxs + mins) / 2.0
    size = (maxs - mins)
    return center, size


def ensure_grids(
    receptor_pdbqt: str,
    out_dir: str,
    mode: str,
    ref_ligand: str = None,
    centers_file: str = None,  # unused here; kept for signature compatibility
    residues_predicate=None,  # function(line)->bool for residues mode
    spacing: float = GRID_SPACING,
    margin: float = 5.0,
    cap: float = 50.0,
    autogrid4_bin: str = "autogrid4",
):
    """
    Build a single grid box and run AutoGrid for modes:
    - 'ligand'   : box around reference ligand (+margin)
    - 'residues' : box around selected receptor residues (+margin)
    - 'blind'    : box over receptor AABB, per-axis capped by 'cap'

    Returns:
    {
      'center': (cx,cy,cz),
      'npts': (nx,ny,nz),
      'spacing': float,
      'fld_path': str,
      'out_dir': str
    }
    """
    mode = (mode or "").lower()
    if mode not in {"ligand", "residues", "blind"}:
        raise ValueError(f"ensure_grids: unsupported mode '{mode}'")

    # Determine center & box size
    if mode == "ligand":
        if not ref_ligand:
            raise ValueError("ligand mode requires ref_ligand path")
        center, size = _load_coords_pdb_like(ref_ligand, atom_filter=lambda n: (not n.startswith("H")))
        size = size + float(margin) * 2.0
    elif mode == "residues":
        if residues_predicate is None:
            raise ValueError("residues mode requires residues_predicate(line)->bool")
        center, size = _load_coords_pdb_like(receptor_pdbqt, atom_filter=lambda n: True, line_filter=residues_predicate)
        size = size + float(margin) * 2.0
    else:
        # blind
        center, size = _load_coords_pdb_like(receptor_pdbqt, atom_filter=lambda n: True)
        size = np.minimum(size, np.array([cap, cap, cap], dtype=float))

    # Convert physical Å to grid points (odd, clamped)
    npts = np.array([_npts_for_size(size[i], spacing) for i in range(3)], dtype=int)
    npts = np.array(_ensure_odd_clamped(npts), dtype=int)

    # Write GPF & run AutoGrid
    out_dir_p = Path(out_dir)
    out_dir_p.mkdir(parents=True, exist_ok=True)
    _write_gpf(out_dir_p / "grid.gpf", receptor_pdbqt, tuple(map(float, center)), tuple(map(int, npts)), float(spacing))
    log = out_dir_p / "grid.glg"
    cmd = [autogrid4_bin, "-p", "grid.gpf", "-l", str(log.name)]
    print(f"[autogrid] single-site ({mode}) in {out_dir_p}")
    res = subprocess.run(cmd, cwd=out_dir_p, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"[autogrid/{mode}] FAILED\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        raise RuntimeError("autogrid failed")

    # Pick up resulting FLD
    fld_candidates = list(out_dir_p.glob("*.fld"))
    if not fld_candidates:
        raise FileNotFoundError(f"No .fld produced in {out_dir_p}; see {log}")
    fld_path = str(fld_candidates[0].resolve())
    return {
        "center": (float(center[0]), float(center[1]), float(center[2])),
        "npts": (int(npts[0]), int(npts[1]), int(npts[2])),
        "spacing": float(spacing),
        "fld_path": fld_path,
        "out_dir": str(out_dir_p.resolve()),
    }


# --- Whole-protein bootstrap (GPF + AutoGrid) ---
import math, subprocess


def _odd(n):
    # odd grid points
    return n if n % 2 == 1 else n + 1


def _load_coords_from_pdb_like(path):
    xs, ys, zs = [], [], []
    with open(path, "r", errors="ignore") as f:
        for ln in f:
            if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                continue
            try:
                xs.append(float(ln[30:38]))
                ys.append(float(ln[38:46]))
                zs.append(float(ln[46:54]))
            except:
                pass
    if not xs:
        raise ValueError(f"No atoms parsed from {path}")
    import numpy as np
    xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
    mins = np.array([xs.min(), ys.min(), zs.min()], float)
    maxs = np.array([xs.max(), ys.max(), zs.max()], float)
    ctr = (mins + maxs) / 2.0
    ext = (maxs - mins)
    return ctr, ext


def _npts_for_extent(extent_A, spacing, cap):
    # cap each axis length to <= cap Å, convert to odd grid points, clamp to [25, 129]
    side = float(min(extent_A, cap))
    n = int(math.ceil((side - 1) / float(spacing)))
    n = _odd(max(n, 25))
    return min(n, 129)


# ========= HotspotGPFGenerator: whole maps + centers.tsv + per-site grids =========
from pathlib import Path
import os
import subprocess

_AD4_TYPES = [
    "A", "C", "HD", "N", "NA", "OA", "SA", "S", "P", "Cl", "F", "Br", "I", "Zn", "Fe", "Mg", "Mn", "Ca"
]


def _ensure_odd_clamped(nxyz, clamp=(60, 255)):
    import numpy as _np
    n = _np.array(nxyz, int)
    n = n + (n % 2 == 0)
    n = _np.clip(n, clamp[0], clamp[1])
    return tuple(int(x) for x in n.tolist())


def _write_site_gpf(gpf_path, receptor_pdbqt, center, npts, spacing):
    """GPF writer that matches ligand_types ⇔ map lines 1:1 (+ elec/dsolv maps)."""
    rec_stem = Path(receptor_pdbqt).stem
    cx, cy, cz = [float(v) for v in center]
    nx, ny, nz = _ensure_odd_clamped(npts)
    with open(gpf_path, "w") as f:
        f.write(f"gridfld {rec_stem}.maps.fld\n")
        f.write(f"npts {nx} {ny} {nz}\n")
        f.write(f"spacing {float(spacing):.3f}\n")
        f.write(f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}\n")
        f.write(f"receptor {Path(receptor_pdbqt).resolve()}\n")
        f.write(f"receptor_types {' '.join(_AD4_TYPES)}\n")
        f.write(f"ligand_types {' '.join(_AD4_TYPES)}\n")
        for t in _AD4_TYPES:
            f.write(f"map {rec_stem}.{t}.map\n")
        f.write(f"elecmap {rec_stem}.e.map\n")
        f.write(f"dsolvmap {rec_stem}.d.map\n")
        ##### f.write("dielectric -0.1465\n")
    return str(gpf_path)

def receptor_aabb_A(pdbqt_path: str):
    mins = np.array([+np.inf, +np.inf, +np.inf], dtype=float)
    maxs = np.array([-np.inf, -np.inf, -np.inf], dtype=float)
    with open(pdbqt_path, "r", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                except Exception:
                    continue
                mins = np.minimum(mins, [x, y, z])
                maxs = np.maximum(maxs, [x, y, z])
    if not np.isfinite(mins).all():
        raise ValueError(f"No atom coords parsed from {pdbqt_path}")
    ctr = (mins + maxs) / 2.0
    ext = (maxs - mins)
    return ctr, ext, mins, maxs

def _odd_int(n: int) -> int:
    return n if (n % 2 == 1) else (n + 1)

def compute_whole_box_auto(
    receptor_pdbqt: str,
    base_spacing: float,
    margin_A: float = 8.0,
    npts_max: int = 255,
):
    """
    Returns center, (nx,ny,nz) for GPF npts, and spacing.
    Uses larger spacing if needed so the box fully covers the receptor AABB (+margin)
    while keeping npts <= npts_max.
    """
    ctr, ext, mins, maxs = receptor_aabb_A(receptor_pdbqt)
    side = ext + 2.0 * float(margin_A)

    # choose spacing so max(axis side)/spacing <= npts_max
    need_sp = float(side.max()) / float(npts_max)
    sp = max(float(base_spacing), need_sp)

    npts = []
    for s in side:
        n = int(math.ceil(float(s) / sp))
        n = _odd_int(n)
        n = max(25, min(npts_max, n))
        npts.append(n)

    return (float(ctr[0]), float(ctr[1]), float(ctr[2])), tuple(npts), float(sp)

def validate_box_covers_receptor(
    receptor_pdbqt: str,
    center,
    npts,
    spacing,
    margin_tol_A: float = 0.5,
):
    """
    Validate that receptor AABB is inside [center - (npts*sp)/2, center + (npts*sp)/2]
    (within a small tolerance).
    """
    ctr, ext, mins, maxs = receptor_aabb_A(receptor_pdbqt)
    center = np.array(center, dtype=float)
    npts = np.array(npts, dtype=float)
    sp = float(spacing)

    half = 0.5 * (npts * sp)  # because side length = npts * spacing
    gmin = center - half
    gmax = center + half

    ok = np.all(mins >= (gmin - margin_tol_A)) and np.all(maxs <= (gmax + margin_tol_A))
    if not ok:
        msg = (
            f"[whole-map] grid does NOT cover receptor!\n"
            f"  receptor mins: {mins}\n"
            f"  receptor maxs: {maxs}\n"
            f"  grid mins:     {gmin}\n"
            f"  grid maxs:     {gmax}\n"
            f"  npts={tuple(map(int,npts))} spacing={sp:.3f}\n"
            f"Fix: increase WHOLE margin/cap or allow larger npts or increase spacing."
        )
        raise ValueError(msg)


def ensure_whole_protein_maps(
    receptor_pdbqt: str,
    out_root: str,
    spacing: float = GRID_SPACING,
    cap_ang: float = None,            # keep arg for API compatibility
    margin_A: float = 8.0,
    npts_max: int = 255,
    autogrid4_bin: str = "autogrid4",
):
    """
    Build whole-receptor maps robustly:
    - auto box from receptor AABB (+margin)
    - auto spacing if needed to keep npts <= npts_max
    - validate coverage before/after
    """
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rec_stem = Path(receptor_pdbqt).stem

    center, npts, sp = compute_whole_box_auto(
        receptor_pdbqt=receptor_pdbqt,
        base_spacing=float(spacing),
        margin_A=float(margin_A),
        npts_max=int(npts_max),
    )

    # Optional: if you REALLY want a hard cap, apply it as an upper bound *only if it still covers receptor*
    # (I recommend leaving cap_ang=None for whole maps, and controlling size via npts_max/auto spacing)
    validate_box_covers_receptor(receptor_pdbqt, center, npts, sp)

    gpf_path = out_root / "grid.gpf"
    _write_site_gpf(gpf_path, receptor_pdbqt, center, npts, sp)

    log_path = out_root / "grid.glg"
    cmd = [autogrid4_bin, "-p", gpf_path.name, "-l", log_path.name]
    res = subprocess.run(cmd, cwd=out_root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"autogrid whole-protein failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )

    fld_candidates = list(out_root.glob(f"{rec_stem}*.fld"))
    if not fld_candidates:
        raise FileNotFoundError(f"No .fld produced in {out_root}. Check {log_path}")
    fld_path = str(fld_candidates[0].resolve())

    # (Optional) post-check: make sure the produced fld bounds cover receptor too (uses your meta loader)
    # meta = load_fld_or_map_meta(fld_path)
    # validate_grid_bounds_vs_receptor(receptor_pdbqt, meta)

    return fld_path


class HotspotGPFGenerator:
    """
    Orchestrates:
    (A) whole-protein maps → .fld
    (B) centers.tsv (hybrid or chosen mode)
    (C) per-site GPF + AutoGrid → per-site .fld
    """

    def __init__(self, autogrid4_bin="autogrid4"):
        self.autogrid4_bin = autogrid4_bin

    def prepare_centers_and_grids(
        self,
        receptor_pdbqt: str,
        out_root: str,
        centers_tsv_path: str,
        mode: str = "hybrid",  # "internal" | "maps" | "hybrid"
        n_sites: int = 6,
        whole_spacing: float = GRID_SPACING,
        whole_cap_ang: float = 100.0,
        hotspot_box_ang: float = 35.0,
        tau_rel: float = 0.60,
        min_sep_A: float = 7.0,
    ):
        """
        Returns: a list of site dicts (site_id, center, npts, spacing, fld_path, out_dir)
        """
        out_root = Path(out_root)
        out_root.mkdir(parents=True, exist_ok=True)

        # (A) one-time whole-protein maps (so we have consistent C/E/D for hotspot finding)
        fld_whole = ensure_whole_protein_maps(
            receptor_pdbqt=receptor_pdbqt,
            out_root=str(out_root),
            spacing=whole_spacing,
            cap_ang=whole_cap_ang,
            autogrid4_bin=self.autogrid4_bin,
        )

        # (B) detect hotspots and write centers.tsv (internal/maps/hybrid)
        autogenerate_centers_tsv(
            receptor_pdbqt=receptor_pdbqt,
            out_root=str(out_root),
            centers_tsv_path=str(centers_tsv_path),
            n_sites=n_sites,
            default_npts=(81, 81, 81),
            default_spacing=whole_spacing,
            blind_cap=whole_cap_ang,
            autogrid4_bin=self.autogrid4_bin,
            hotspot_box_ang=hotspot_box_ang,
            hotspot_sigma_A=1.0,
            hotspot_bury_z=0.35,
            mode=mode,
        )

        # (C) per-site GPF + AutoGrid into <out_root>/<site_id>/*
        sites = ensure_grids_multi_centers(
            receptor_pdbqt=receptor_pdbqt,
            out_root=str(out_root),
            centers_tsv=str(centers_tsv_path),
            autogrid4_bin=self.autogrid4_bin,
        )
        return sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fld_or_dir", help=".fld file or directory containing whole-protein maps")
    ap.add_argument("--receptor-pdbqt", required=True, help="macromolecule .pdbqt for geometry EDT")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--voxel-spacing", type=float, default=0.375, help="EDT voxel size (Å) for geometry stage")
    ap.add_argument("--r-min", type=float, default=R_MIN_CAVITY_A, help="min inscribed-sphere radius (Å)")
    ap.add_argument("--min-sep", type=float, default=3.0, help="NMS min separation between hotspots (Å)")
    ap.add_argument("--k-box", type=float, default=10.0, help="box half-size multiplier × r_peak")
    ap.add_argument("--max-sites", type=int, default=8)
    ap.add_argument("--inflate", type=float, default=0.25, help="inflate atom radii by this many Å when voxelizing (0.0–0.5 typical)")
    ap.add_argument("--mode", choices=["internal", "maps", "hybrid"], default="hybrid", help="internal=EDT cavities; maps=C/E/D peaks; hybrid=internal then fallback to maps")
    ap.add_argument("--tau-rel", type=float, default=0.52, help="relative peak threshold for maps mode (0.58–0.62 typical)")
    ap.add_argument("--half-size", type=float, default=12.0, help="box half-size (Å) for maps mode")
    args = ap.parse_args()

    meta = load_fld_or_map_meta(args.fld_or_dir)
    origin, spacing, shape = meta["origin"], meta["spacing"], meta["shape"]
    mp = meta["map_paths"]

    # load maps
    C = load_map_ascii(mp["C"], shape)
    E = load_map_ascii(mp["E"], shape)
    D = load_map_ascii(mp["D"], shape)
    maps = {"C": C, "E": E, "D": D}

    # pick centers according to mode
    if args.mode == "internal":
        sites = pick_centers(
            args.receptor_pdbqt, maps,
            map_origin=origin, map_spacing=spacing,
            voxel_spacing=args.voxel_spacing, r_min=args.r_min,
            min_sep_A=args.min_sep, k_box=args.k_box, max_sites=args.max_sites,
            inflate_A=args.inflate
        )
    elif args.mode == "maps":
        print("[maps] detecting surface-pocket hotspots from C/E/D maps")
        sites = detect_maps_hotspots(maps["C"], maps["E"], maps["D"], origin, spacing,
                                     tau_rel=args.tau_rel, min_sep_A=args.min_sep,
                                     max_sites=args.max_sites, half_size_A=args.half_size)
    else:
        # hybrid
        print("[hybrid] internal cavities first…")
        sites = pick_centers(
            args.receptor_pdbqt, maps,
            map_origin=origin, map_spacing=spacing,
            voxel_spacing=0.5, r_min=R_MIN_CAVITY_A,  # <-- enforce 20 Å minimum cavity radius
            min_sep_A=7.0, k_box=4.0, max_sites=args.max_sites, inflate_A=0.0,
        )
        if len(sites) == 0:
            print("[hybrid] no internal sites → falling back to maps hotspots")
            sites = detect_maps_hotspots(maps["C"], maps["E"], maps["D"], origin, spacing,
                                         tau_rel=args.tau_rel, min_sep_A=args.min_sep,
                                         max_sites=args.max_sites, half_size_A=args.half_size)

    # write TSV
    rec_stem = Path(args.receptor_pdbqt).stem
    out_tsv = Path(args.out) / "centers.tsv"
    write_centers_tsv(out_tsv, rec_stem, sites)


if __name__ == "__main__":
    main()