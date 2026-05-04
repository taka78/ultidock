#!/usr/bin/env python3
"""
profile_receptors.py  —  Ultidock per-receptor config generator

Scans every .pdbqt in MACRO_MOL_DIR, runs a lightweight geometry analysis
(bounding box, residue count, EDT-based cavity depth), and writes a
receptor_name.config.toml sidecar next to each receptor.

Already-existing sidecars are skipped by default. Use --force to overwrite.

Usage:
    python3 docking/profile_receptors.py
    python3 docking/profile_receptors.py --force
    python3 docking/profile_receptors.py --receptor path/to/specific.pdbqt
    python3 docking/profile_receptors.py --dry-run
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

# ── path bootstrap so we can import config ────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import config as _config
from config import (
    MACRO_MOL_DIR,
    GRID_SPACING,
    HOTSPOT_BOX_ANGLE,
    HOTSPOT_NMS_MINSEP_A,
    AUTOSITES,
    GRID_MARGIN,
    GRID_CAP,
)

ADAPTIVE_R_MIN_PERCENTILE = float(getattr(_config, "ADAPTIVE_R_MIN_PERCENTILE", 85.0))
ADAPTIVE_R_MIN_PEAK_PERCENTILE = float(getattr(_config, "ADAPTIVE_R_MIN_PEAK_PERCENTILE", 10.0))
ADAPTIVE_R_MIN_FLOOR_A = float(getattr(_config, "ADAPTIVE_R_MIN_FLOOR_A", 2.0))
ADAPTIVE_R_MIN_CEIL_A = float(getattr(_config, "ADAPTIVE_R_MIN_CEIL_A", 5.0))
ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A = float(getattr(_config, "ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A", 8.0))
ADAPTIVE_R_MIN_PEAK_WINDOW_A = float(getattr(_config, "ADAPTIVE_R_MIN_PEAK_WINDOW_A", 4.0))
MAPS_POCKET_MAX_A = float(getattr(_config, "MAPS_POCKET_MAX_A", 15.0))
HOTSPOT_NMS_BOX_FRACTION = float(getattr(_config, "HOTSPOT_NMS_BOX_FRACTION", 0.40))
HOTSPOT_NMS_MIN_A = float(getattr(_config, "HOTSPOT_NMS_MIN_A", 4.0))
HOTSPOT_NMS_MAX_A = float(getattr(_config, "HOTSPOT_NMS_MAX_A", 25.0))

# ── geometry helpers ───────────────────────────────────────────────────────

def _parse_pdbqt_atoms(path: Path) -> tuple[np.ndarray, int]:
    """
    Return (coords, n_residues).
    coords     : (N, 3) float32 array of all heavy-atom positions
    n_residues : number of unique (chain, res_seq) pairs
    """
    coords   = []
    residues: set[tuple[str, str]] = set()

    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append((x, y, z))
            chain   = line[21:22].strip()
            res_seq = line[22:26].strip()
            residues.add((chain, res_seq))

    if not coords:
        raise ValueError(f"No heavy-atom coordinates found in {path}")

    return np.array(coords, dtype=np.float32), len(residues)


def _bounding_box(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mins   = coords.min(axis=0)
    maxs   = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    return mins, maxs, center


def _build_edt(
    coords: np.ndarray,
    spacing: float = 0.8,
    *,
    max_dim: int = 200,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Voxelise the receptor at `spacing` Å and compute EDT (distance from each
    empty voxel to the nearest occupied voxel, in Å).

    Uses a coarse grid (0.8 Å default) purely for speed — this is a profiling
    pass, not a docking grid.

    Returns (edt, origin, actual_spacing) where edt is a float32 ndarray in Å.
    """
    from scipy.ndimage import distance_transform_edt as _edt

    pad    = 6.0
    origin = coords.min(axis=0) - pad
    maxs   = coords.max(axis=0) + pad

    extents = maxs - origin
    actual_spacing = max(float(spacing), float(extents.max()) / float(max(2, max_dim - 1)))
    shape = np.ceil(extents / actual_spacing).astype(int) + 1
    shape = np.maximum(shape, 10)

    occ = np.zeros(shape, dtype=bool)

    atom_radius = 1.8
    gx = origin[0] + np.arange(shape[0]) * actual_spacing
    gy = origin[1] + np.arange(shape[1]) * actual_spacing
    gz = origin[2] + np.arange(shape[2]) * actual_spacing

    for x, y, z in coords:
        ix0 = max(0, int(math.floor((x - atom_radius - origin[0]) / actual_spacing)))
        ix1 = min(shape[0] - 1, int(math.ceil((x + atom_radius - origin[0]) / actual_spacing)))
        iy0 = max(0, int(math.floor((y - atom_radius - origin[1]) / actual_spacing)))
        iy1 = min(shape[1] - 1, int(math.ceil((y + atom_radius - origin[1]) / actual_spacing)))
        iz0 = max(0, int(math.floor((z - atom_radius - origin[2]) / actual_spacing)))
        iz1 = min(shape[2] - 1, int(math.ceil((z + atom_radius - origin[2]) / actual_spacing)))
        dx2 = (gx[ix0:ix1 + 1] - x) ** 2
        dy2 = (gy[iy0:iy1 + 1] - y) ** 2
        dz2 = (gz[iz0:iz1 + 1] - z) ** 2
        mask = (dx2[:, None, None] + dy2[None, :, None] + dz2[None, None, :]) <= atom_radius ** 2
        occ[ix0:ix1 + 1, iy0:iy1 + 1, iz0:iz1 + 1] |= mask

    edt = _edt(~occ).astype(np.float32) * actual_spacing
    return edt, origin, actual_spacing


def _cavity_stats(edt: np.ndarray, spacing: float = 0.8) -> dict:
    """
    Extract pocket-relevant statistics from an EDT grid.

    We focus on voxels in the "pocket zone": EDT value between 1.4 Å (too
    tight for any drug-like atom) and 8.0 Å (beyond this you're in open
    solvent, not a pocket).  This window captures real binding sites and
    ignores both steric clashes and bulk solvent.

    Returns a dict with:
        p25, median, p75  — percentiles of pocket-zone EDT values (Å)
        n_peaks           — number of distinct local EDT maxima (cavity count)
        median_peak_r     — median EDT value at detected peaks (Å)
    """
    from scipy.ndimage import maximum_filter, label

    POCKET_MIN = 1.4   # Å — smaller than this is a clash, not a pocket
    POCKET_MAX = 8.0   # Å — larger than this is open solvent

    pocket_mask = (edt >= POCKET_MIN) & (edt <= POCKET_MAX)

    if not np.any(pocket_mask):
        # fallback: use whatever's available
        pocket_vals = edt[edt > 0].ravel()
    else:
        pocket_vals = edt[pocket_mask].ravel()

    p25    = float(np.percentile(pocket_vals, 25))
    median = float(np.percentile(pocket_vals, 50))
    p75    = float(np.percentile(pocket_vals, 75))

    # ── count distinct cavity peaks ───────────────────────────────────────
    # A "peak" is a local EDT maximum inside the pocket zone.
    # We run maximum_filter with a window that corresponds to ~4 Å separation
    # (so two peaks must be at least 4 Å apart to be counted separately).
    # Then we label connected regions of peak voxels.
    win = max(3, int(round(4.0 / spacing)))   # 4 Å separation — resolves adjacent sub-pockets
    local_max = maximum_filter(edt, size=win, mode="constant", cval=0.0)
    peak_mask = (edt == local_max) & pocket_mask

    # erode to single-voxel representatives
    labeled, n_peaks = label(peak_mask)

    # median EDT at peak centres
    if n_peaks > 0:
        peak_vals = [float(edt[labeled == i].max()) for i in range(1, n_peaks + 1)]
        # Only peaks with EDT >= p25 are "significant" (real drug-binding
        # cavities, not shallow surface crevices).  The raw n_peaks includes
        # surface noise; n_significant drives AUTOSITES.
        significant_peak_vals = [v for v in peak_vals if v >= p25]
        n_significant = len(significant_peak_vals)
        median_peak_r = float(np.median(peak_vals))
        significant_peak_p10 = float(np.percentile(significant_peak_vals, 10)) if significant_peak_vals else p25
    else:
        n_significant = 0
        median_peak_r = median
        significant_peak_p10 = p25

    return {
        "p25":            round(p25,    2),
        "median":         round(median, 2),
        "p75":            round(p75,    2),
        "n_peaks":        int(n_peaks),
        "n_significant":  int(n_significant),
        "median_peak_r":  round(median_peak_r, 2),
        "significant_peak_p10": round(significant_peak_p10, 2),
    }


def _select_r_min_cavity(stats: dict) -> float:
    """
    Pick the minimum accepted cavity-centre radius from peak statistics, not
    from the raw distribution of all pocket voxels.

    The old p25(all pocket voxels) rule was dominated by shallow surface-shell
    geometry and systematically underestimated the depth of real cavity centres.
    The 10th percentile of significant peak radii is a better proxy for the
    shallowest *actual* pockets that should still be retained.
    """
    peak_floor = float(stats.get("significant_peak_p10", stats["p25"]))
    return round(float(np.clip(peak_floor, ADAPTIVE_R_MIN_FLOOR_A, ADAPTIVE_R_MIN_CEIL_A)), 2)


# ── per-receptor parameter derivation ─────────────────────────────────────

def derive_params(path: Path) -> dict:
    """
    Analyse a receptor PDBQT and return a dict of suggested config overrides.

    Key design decisions
    --------------------
    R_MIN_CAVITY_A
        The main pipeline now profiles this adaptively at center-generation
        time, using the same map-aligned receptor grid as the site finder.
        This script still controls the clamp by writing receptor-specific
        adaptive floor/ceiling bounds and keeps the fixed r_min value as a
        diagnostic suggestion only.

    AUTOSITES
        Count *significant* EDT peaks (those with EDT ≥ p25) rather than
        all local maxima.  Raw peak counts are dominated by shallow surface
        noise; filtering by p25 keeps only genuine drug-binding cavities.
        Clamped to [2, 8].

    HOTSPOT_BOX_ANGLE
        A docking box must enclose the pocket plus approach vectors for the
        ligand.  We use the larger of median_peak_r and p75 as the cavity
        radius proxy (median_peak_r alone is dragged down by surface noise),
        then scale by 5× to give ample sampling room.  Clamped to [30, 55] Å
        — below 30 truncates moderate pockets, above 55 wastes GPU time.

    ADAPTIVE_R_MIN_*
        Bounds around the runtime adaptive r_min estimate. The profiler uses
        coarse receptor geometry to keep the runtime value in a receptor-local
        range, while make_grids still computes the final threshold on the
        actual map-aligned grid.

    HOTSPOT_NMS_MINSEP_A
        The non-maximum-suppression radius should avoid merging distinct
        pockets without rejecting valid nearby sub-pockets.  It is tied to
        the generated box side (default 0.4 × box side) so multiple accepted
        boxes do not collapse onto the same region.  Clamped by the configured
        HOTSPOT_NMS_MIN_A / HOTSPOT_NMS_MAX_A bounds.

    GRID_SPACING
        Coarser for large proteins (saves AutoGrid time), finer for small
        ones (better pose resolution).  Breakpoints at 200 and 400 residues.
    """
    coords, n_residues = _parse_pdbqt_atoms(path)
    mins, maxs, center = _bounding_box(coords)
    extents             = maxs - mins
    longest_axis        = float(extents.max())
    volume_proxy        = float(extents.prod())

    edt, _origin, edt_spacing = _build_edt(coords, spacing=0.8)
    stats        = _cavity_stats(edt, spacing=edt_spacing)

    p25           = stats["p25"]
    p75           = stats["p75"]
    median_peak_r = stats["median_peak_r"]
    n_peaks       = stats["n_peaks"]
    n_significant = stats.get("n_significant", n_peaks)

    # ── R_MIN_CAVITY_A ────────────────────────────────────────────────────
    # Use the lower tail of significant cavity peaks rather than all pocket
    # voxels; this tracks real cavity-centre depths instead of shallow walls.
    r_min = _select_r_min_cavity(stats)

    # Keep runtime adaptive r_min receptor-aware without freezing it.  The
    # sidecar supplies bounds; make_grids computes the final map-aligned value.
    adaptive_floor = round(float(np.clip(r_min - 0.75, 1.5, 2.5)), 2)
    adaptive_ceil = round(float(np.clip(max(r_min + 1.75, p75 + 1.0), 4.0, 7.0)), 2)
    adaptive_zone_max = round(float(max(ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A, adaptive_ceil + 2.0)), 1)

    # ── AUTOSITES ─────────────────────────────────────────────────────────
    # Use significant peaks (EDT ≥ p25) — raw peak count includes noise.
    autosites = int(np.clip(n_significant, 2, 5))

    # ── HOTSPOT_BOX_ANGLE ─────────────────────────────────────────────────
    # Use the larger of median_peak_r and p75 as cavity radius proxy;
    # median_peak_r alone is dragged down by surface pockets in large proteins.
    cavity_r = max(median_peak_r, p75)
    box_angle = round(float(np.clip(5.0 * cavity_r, 30.0, 55.0)), 1)

    # ── HOTSPOT_NMS_MINSEP_A ──────────────────────────────────────────────
    # Non-Maximum Suppression (NMS) minimum separation between accepted
    # binding-site centres.
    #
    # WHY scale with box size?
    #   The purpose of NMS is to ensure that each accepted site covers a
    #   *distinct* region of the protein.  Two sites whose docking boxes
    #   overlap > 80% are essentially sampling the same conformational space
    #   and waste GPU time.  Box overlap is determined by the ratio of
    #   inter-site distance to box size:
    #     overlap ≈ 1 - (distance / box_diameter)
    #   For meaningful coverage diversity, we need overlap < 60%, which
    #   requires distance > 0.4 × box_diameter ≈ 0.4 × box_angle.
    #
    # WHY NOT scale with pocket depth (the old formula)?
    #   The previous formula (0.8 × p25, clamped to [1.5, 5.0] Å) produced
    #   NMS = 2.3 Å for a 63 Å box — meaning two site centres only 2.3 Å
    #   apart (one grid voxel!) were accepted as "distinct".  This led to
    #   4-5 nearly identical sites clustered within 4-8 Å of each other,
    #   all sampling the same protein surface region.
    #
    # Bounds: [4, 25] Å ensures at minimum two sites are 4 Å apart (allowing
    # nearby sub-pockets), and at most 25 Å apart (for the largest 55 Å boxes,
    # preventing excessive site rejection on very large proteins).
    nms_minsep = round(
        float(
            np.clip(
                HOTSPOT_NMS_BOX_FRACTION * box_angle,
                HOTSPOT_NMS_MIN_A,
                HOTSPOT_NMS_MAX_A,
            )
        ),
        1,
    )

    # Deep transmembrane/channel-like targets benefit from the 15 Å maps shell;
    # compact soluble proteins can stay tighter to reduce surface noise.
    maps_pocket_max = 15.0 if (n_residues > 300 or longest_axis > 55.0 or p75 >= 4.5) else 12.0

    # ── GRID_SPACING ──────────────────────────────────────────────────────
    if n_residues < 200:
        grid_spacing = 0.25
    elif n_residues > 400:
        grid_spacing = 0.375
    else:
        grid_spacing = float(GRID_SPACING)

    return {
        # ── derived ──────────────────────────────────────────────────────
        "AUTOSITES":            autosites,
        "HOTSPOT_BOX_ANGLE":    box_angle,
        "HOTSPOT_NMS_MINSEP_A": nms_minsep,
        "ADAPTIVE_R_MIN_PERCENTILE": ADAPTIVE_R_MIN_PERCENTILE,
        "ADAPTIVE_R_MIN_PEAK_PERCENTILE": ADAPTIVE_R_MIN_PEAK_PERCENTILE,
        "ADAPTIVE_R_MIN_FLOOR_A": adaptive_floor,
        "ADAPTIVE_R_MIN_CEIL_A": adaptive_ceil,
        "ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A": adaptive_zone_max,
        "ADAPTIVE_R_MIN_PEAK_WINDOW_A": ADAPTIVE_R_MIN_PEAK_WINDOW_A,
        "MAPS_POCKET_MAX_A": maps_pocket_max,
        "HOTSPOT_NMS_BOX_FRACTION": HOTSPOT_NMS_BOX_FRACTION,
        "HOTSPOT_NMS_MIN_A": HOTSPOT_NMS_MIN_A,
        "HOTSPOT_NMS_MAX_A": HOTSPOT_NMS_MAX_A,
        "GRID_SPACING":         grid_spacing,
        # ── passthrough (hand-tunable) ────────────────────────────────────
        "GRID_MARGIN":          float(GRID_MARGIN),
        "GRID_CAP":             float(GRID_CAP),
        # ── diagnostics (written as comments, not loaded by pipeline) ─────
        "_n_residues":          n_residues,
        "_longest_axis_A":      round(float(longest_axis), 1),
        "_volume_proxy_A3":     round(volume_proxy, 0),
        "_edt_p25_A":           stats["p25"],
        "_edt_median_A":        stats["median"],
        "_edt_p75_A":           stats["p75"],
        "_edt_spacing_A":       round(float(edt_spacing), 3),
        "_n_cavity_peaks":      n_peaks,
        "_n_significant_peaks": n_significant,
        "_median_peak_r_A":     stats["median_peak_r"],
        "_significant_peak_p10_A": stats["significant_peak_p10"],
        "_suggested_R_MIN_CAVITY_A": r_min,
    }


# ── TOML writer (stdlib tomllib only reads; we write manually) ─────────────

def _write_toml(path: Path, params: dict, receptor_name: str) -> None:
    diag_keys = {k for k in params if k.startswith("_")}
    cfg_keys  = [k for k in params if k not in diag_keys]

    lines = [
        f"# {receptor_name}.config.toml",
        "# Auto-generated by profile_receptors.py — review and adjust as needed.",
        "# Delete this file to fall back to global config.py defaults.",
        "#",
        "# Diagnostic info (read-only, not loaded by the pipeline):",
    ]
    for k in sorted(diag_keys):
        lines.append(f"#   {k[1:]:30s} = {params[k]}")

    lines += [
        "#",
        "# Per-receptor overrides (all values in Ångström unless noted):",
        "#   AUTOSITES         — number of binding sites to detect",
        "#   R_MIN_CAVITY_A    — optional fixed cavity radius override; omit for runtime adaptive",
        "#   ADAPTIVE_R_MIN_*  — receptor-specific bounds for runtime adaptive r_min",
        "#   MAPS_POCKET_MAX_A — maximum EDT depth included in maps/surface pocket shell",
        "#   HOTSPOT_BOX_ANGLE — minimum docking box side length (Å)",
        "#   HOTSPOT_NMS_MINSEP_A — minimum distance between two hotspots (Å)",
        "#   GRID_SPACING      — AutoGrid point spacing (Å); smaller = finer but slower",
        "",
    ]
    for k in cfg_keys:
        v = params[k]
        lines.append(f"{k} = {v}")

    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


# ── main ──────────────────────────────────────────────────────────────────

def profile_one(pdbqt: Path, *, force: bool = False, dry_run: bool = False) -> bool:
    """Profile a single receptor. Returns True if a sidecar was written."""
    sidecar = pdbqt.with_suffix(".config.toml")

    if sidecar.exists() and not force:
        print(f"  [SKIP] {pdbqt.name}  (sidecar exists; use --force to overwrite)")
        return False

    print(f"  [PROFILE] {pdbqt.name} ...", end=" ", flush=True)
    try:
        params = derive_params(pdbqt)
    except Exception as exc:
        print(f"FAILED — {exc}")
        return False

    print(
        f"residues={params['_n_residues']}  "
        f"peaks={params['_n_cavity_peaks']}({params['_n_significant_peaks']} sig)  "
        f"p25={params['_edt_p25_A']:.1f}Å  "
        f"median_r={params['_median_peak_r_A']:.1f}Å  "
        f"--> AUTOSITES={params['AUTOSITES']}  "
        f"R_MIN=auto[{params['ADAPTIVE_R_MIN_FLOOR_A']}-{params['ADAPTIVE_R_MIN_CEIL_A']}] "
        f"(suggest {params['_suggested_R_MIN_CAVITY_A']})  "
        f"BOX={params['HOTSPOT_BOX_ANGLE']}  "
        f"NMS={params['HOTSPOT_NMS_MINSEP_A']}"
    )

    if not dry_run:
        _write_toml(sidecar, params, pdbqt.stem)
        print(f"           wrote --> {sidecar.name}")

    return True


def profile_all(macro_dir: Path, *, force: bool = False, dry_run: bool = False) -> int:
    """Profile all receptors in macro_dir. Returns the number of sidecars written."""
    pdbqt_files = sorted(macro_dir.glob("*.pdbqt"))
    if not pdbqt_files:
        print(f"[profile] No .pdbqt files found in {macro_dir}")
        return 0

    print(f"[profile] Found {len(pdbqt_files)} receptor(s) in {macro_dir}")
    written = 0
    for pdbqt in pdbqt_files:
        if profile_one(pdbqt, force=force, dry_run=dry_run):
            written += 1

    action = "would write" if dry_run else "wrote"
    print(f"[profile] Done — {action} {written} sidecar(s).")
    return written


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Generate per-receptor .config.toml sidecars for Ultidock"
    )
    ap.add_argument(
        "--macro-mol-dir",
        default=MACRO_MOL_DIR,
        help="Directory containing receptor PDBQT files (default: from config.py)",
    )
    ap.add_argument(
        "--receptor",
        metavar="PATH",
        help="Profile a single receptor instead of the whole directory",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing sidecars",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print derived params without writing any files",
    )
    args = ap.parse_args(argv)

    if args.receptor:
        pdbqt = Path(args.receptor).resolve()
        if not pdbqt.exists():
            print(f"[error] Receptor not found: {pdbqt}")
            return 1
        profile_one(pdbqt, force=args.force, dry_run=args.dry_run)
    else:
        macro_dir = Path(args.macro_mol_dir).resolve()
        if not macro_dir.exists():
            print(f"[error] MACRO_MOL_DIR not found: {macro_dir}")
            return 1
        profile_all(macro_dir, force=args.force, dry_run=args.dry_run)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
