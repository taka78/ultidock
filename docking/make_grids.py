import numpy as np
from scipy.ndimage import binary_closing, distance_transform_edt, label, generate_binary_structure
from collections import namedtuple
import hashlib
import re, os, argparse
from pathlib import Path
from scipy.ndimage import gaussian_filter, maximum_filter
import mmap
from grid_box import blind_box
from config import (
    AUTODOCK_GPU_DIR,
    DOCKING_DIR,
    RESULTS_DIR,
    NUMWI,
    LIGANDS_DIR,
    CENTERS_TSV,
    GRID_MODE,
    GRID_SPACING,
    GRID_MARGIN,
    GRID_CAP,
    R_MIN_CAVITY_A,
    HOTSPOT_BOX_ANGLE,
    HOTSPOT_NMS_MINSEP_A,
    SURFACE_SHELL__MIN_A,
    SURFACE_SHELL__MAX_A,
    SURFACE_NMS_MINSEP_A,
    MAX_CENTER_DIST_A,
    CONTACT_SHELL_A,
    MIN_SURFACE_FRAC,
    AUTOSITES,
)
import config as _config

spacing = float(GRID_SPACING)

ADAPTIVE_R_MIN_PERCENTILE = float(getattr(_config, "ADAPTIVE_R_MIN_PERCENTILE", 85.0))
ADAPTIVE_R_MIN_PEAK_PERCENTILE = float(getattr(_config, "ADAPTIVE_R_MIN_PEAK_PERCENTILE", 10.0))
ADAPTIVE_R_MIN_FLOOR_A = float(getattr(_config, "ADAPTIVE_R_MIN_FLOOR_A", 2.0))
ADAPTIVE_R_MIN_CEIL_A = float(getattr(_config, "ADAPTIVE_R_MIN_CEIL_A", 5.0))
ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A = float(getattr(_config, "ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A", 8.0))
ADAPTIVE_R_MIN_PEAK_WINDOW_A = float(getattr(_config, "ADAPTIVE_R_MIN_PEAK_WINDOW_A", 4.0))
MAPS_POCKET_MAX_A = float(getattr(_config, "MAPS_POCKET_MAX_A", 15.0))
MAPS_CONVOLUTION_RADIUS_A = float(getattr(_config, "MAPS_CONVOLUTION_RADIUS_A", 8.0))
MAPS_GAUSSIAN_SIGMA_A = float(
    getattr(_config, "MAPS_GAUSSIAN_SIGMA_A", max(1.0, 0.5 * MAPS_CONVOLUTION_RADIUS_A))
)
MAPS_C_WEIGHT = float(getattr(_config, "MAPS_C_WEIGHT", 0.42))
MAPS_E_WEIGHT = float(getattr(_config, "MAPS_E_WEIGHT", 0.28))
MAPS_D_WEIGHT = float(getattr(_config, "MAPS_D_WEIGHT", 0.30))
MAPS_EDT_WEIGHT_FLOOR = float(getattr(_config, "MAPS_EDT_WEIGHT_FLOOR", 0.15))
MAPS_EDT_WEIGHT_POWER = float(getattr(_config, "MAPS_EDT_WEIGHT_POWER", 1.25))
HOTSPOT_NMS_BOX_FRACTION = float(getattr(_config, "HOTSPOT_NMS_BOX_FRACTION", 0.40))
HOTSPOT_NMS_MIN_A = float(getattr(_config, "HOTSPOT_NMS_MIN_A", 4.0))
HOTSPOT_NMS_MAX_A = float(getattr(_config, "HOTSPOT_NMS_MAX_A", 25.0))
INTERNAL_MEDOID_MIN_DEPTH_FRACTION = float(
    getattr(_config, "INTERNAL_MEDOID_MIN_DEPTH_FRACTION", 0.50)
)

# AutoGrid4 accepts at most 1024 intervals (1025 grid points) per axis.  Keep
# docking maps within AutoDock-GPU's 256-point limit, but let CaV-EMPS retain
# the configured base resolution over much larger whole-receptor boxes.
AUTOGRID4_NPTS_MAX = 1024
DOCKING_GRID_NPTS_MAX = 255
CAV_EMPS_WHOLE_NPTS_MAX = int(getattr(_config, "CAV_EMPS_WHOLE_NPTS_MAX", AUTOGRID4_NPTS_MAX))
if not 25 <= CAV_EMPS_WHOLE_NPTS_MAX <= AUTOGRID4_NPTS_MAX:
    raise ValueError(
        "CAV_EMPS_WHOLE_NPTS_MAX must be between 25 and "
        f"{AUTOGRID4_NPTS_MAX}, found {CAV_EMPS_WHOLE_NPTS_MAX}"
    )


# ── AD4 parameter-file locator ────────────────────────────────────────────
_AD4_PARAM_FILE_CACHE: str | None = None

def _find_ad4_parameter_file(autogrid4_bin: str = "autogrid4") -> str | None:
    """
    Locate the AD4_parameters.dat shipped with AutoGrid.

    Search order:
    1. <autogrid4 binary dir>/../ad4_shared/AD4_parameters.dat
    2. <autogrid4 binary dir>/ad4_shared/AD4_parameters.dat
    3. <AUTODOCK_GPU_DIR>/autogrid/ad4_shared/AD4_parameters.dat
    4. Fall back to AD4.1_bound.dat in the same locations
    """
    global _AD4_PARAM_FILE_CACHE
    if _AD4_PARAM_FILE_CACHE is not None:
        return _AD4_PARAM_FILE_CACHE

    candidates = []
    # Resolve the binary path
    ag_path = Path(autogrid4_bin).resolve()
    ag_dir = ag_path.parent

    for base_dir in [ag_dir.parent, ag_dir, Path(AUTODOCK_GPU_DIR) / "autogrid"]:
        for name in ["AD4_parameters.dat", "AD4.1_bound.dat"]:
            p = base_dir / "ad4_shared" / name
            candidates.append(p)
            # Also check Tests/ dir
            candidates.append(base_dir / "Tests" / name)

    for p in candidates:
        if p.is_file():
            _AD4_PARAM_FILE_CACHE = str(p.resolve())
            return _AD4_PARAM_FILE_CACHE

    print("[WARNING] Could not locate AD4_parameters.dat — GPF will lack parameter_file directive")
    return None


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


def _parse_centers_metadata(centers_tsv_path) -> dict[str, str]:
    meta: dict[str, str] = {}
    try:
        with open(centers_tsv_path, "r", errors="ignore") as handle:
            for line in handle:
                line = line.strip()
                if not line.startswith("# meta "):
                    continue
                for token in line[len("# meta "):].split():
                    if "=" not in token:
                        continue
                    key, value = token.split("=", 1)
                    meta[key.strip()] = value.strip()
    except OSError:
        return {}
    return meta


def _ensure_odd_clamped(nxyz, clamp=(60, 255)):
    import numpy as _np
    n = _np.array(nxyz, int)
    n = n + (n % 2 == 0)
    n = _np.clip(n, clamp[0], clamp[1])
    return tuple(int(x) for x in n.tolist())


def receptor_atom_types(pdbqt_path):
    found = []
    valid = set(_AD4_TYPES)
    with open(pdbqt_path, "r", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_type = line[77:79].strip()
            if not atom_type:
                raise ValueError(f"{pdbqt_path}: missing PDBQT atom type in line: {line.rstrip()}")
            if atom_type not in valid:
                raise ValueError(
                    f"{pdbqt_path}: unsupported receptor atom type {atom_type!r}; "
                    f"expected one of {', '.join(_AD4_TYPES)}"
                )
            if atom_type not in found:
                found.append(atom_type)
    if not found:
        raise ValueError(f"{pdbqt_path}: no ATOM/HETATM records with PDBQT atom types found")
    return tuple(found)


def _write_gpf(gpf_path, receptor_pdbqt, center, npts, spacing, autogrid4_bin="autogrid4", map_types=None):
    """Write a per-site GPF; **do not** modify npts here."""
    cx, cy, cz = map(float, center)  # correct: center → (cx,cy,cz)
    nx, ny, nz = (int(npts[0]), int(npts[1]), int(npts[2]))  # correct: npts → (nx,ny,nz)
    rec_path = Path(receptor_pdbqt).resolve()
    rec_stem = rec_path.stem
    receptor_types = receptor_atom_types(rec_path)
    ligand_types = _normalize_autogrid_map_types(map_types)
    fld_name = f"{rec_stem}.maps.fld"
    param_file = _find_ad4_parameter_file(autogrid4_bin)
    with open(gpf_path, "w") as f:
        # parameter_file MUST come first — AutoGrid reads parameters before
        # processing any atom types.  Without it, all maps are zero.
        if param_file:
            f.write(f"parameter_file {param_file}\n")
        # AutoGrid reads GPFs sequentially. The receptor atom types must be
        # declared before the receptor is parsed, and the receptor must be
        # parsed before gridcenter so coordinates are translated into the grid
        # frame. If gridcenter comes first, far-from-origin receptors produce
        # zero C/D maps because atoms never enter the nonbonded cutoff.
        f.write(f"gridfld {fld_name}\n")
        f.write(f"npts {nx} {ny} {nz}\n")
        f.write(f"spacing {float(spacing):.3f}\n")
        f.write("receptor_types " + " ".join(receptor_types) + "\n")
        f.write(f"receptor {rec_path}\n")
        f.write(f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}\n")
        f.write("ligand_types " + " ".join(ligand_types) + "\n")
        f.write("smooth 0.500\n")
        for t in ligand_types:
            f.write(f"map {rec_stem}.{t}.map\n")
        f.write(f"elecmap {rec_stem}.e.map\n")
        f.write(f"dsolvmap {rec_stem}.d.map\n")
        f.write("dielectric -0.1465\n")
    return str(gpf_path)


def autogenerate_centers_tsv(
    receptor_pdbqt: str,
    out_root: str,
    centers_tsv_path: str,
    n_sites: int = AUTOSITES,
    default_npts=(81, 81, 81),
    default_spacing: float = GRID_SPACING,
    blind_cap: float = GRID_CAP,
    autogrid4_bin: str = "autogrid4",
    hotspot_box_ang: float = HOTSPOT_BOX_ANGLE,
    hotspot_sigma_A: float = 1.0,  # (kept for API compatibility; not used directly below)
    hotspot_bury_z: float = 0.35,  # (ditto)
    tau_rel: float = 0.52,
    min_sep_A: float = HOTSPOT_NMS_MINSEP_A,
    k_box: float = 4.0,
    r_min: float | None = R_MIN_CAVITY_A,
    mode: str = "receptor_search",
    adaptive_r_min_params: dict | None = None,
    maps_pocket_max_A: float | None = None,
    map_types=None,
    nms_box_fraction: float = HOTSPOT_NMS_BOX_FRACTION,
    nms_min_A: float = HOTSPOT_NMS_MIN_A,
    nms_max_A: float = HOTSPOT_NMS_MAX_A,
    whole_map_npts_max: int = CAV_EMPS_WHOLE_NPTS_MAX,
):
    """
    Generate centers.tsv under out_root using already-present maps or by making whole-protein maps.
    Reuse existing centers.tsv only when its metadata matches the requested
    search policy. Policy options:
      - receptor_search  : compact portfolio of complementary receptor-derived sites
      - exhaustive_search: broader multi-site search across all families
      - internal         : internal-cavity family only
      - surface / maps   : surface-cleft family only
      - hybrid           : consensus family only
    """
    from pathlib import Path

    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    centers_tsv_path = Path(centers_tsv_path)
    rec_stem = Path(receptor_pdbqt).stem
    policy = (mode or "receptor_search").lower()
    if policy == "maps":
        policy = "surface"
    expected_site_count = _expected_site_count(policy, n_sites)
    requested_r_min = _coerce_optional_float(r_min)
    adaptive_params = _adaptive_r_min_params(adaptive_r_min_params)
    maps_pocket_max_value = MAPS_POCKET_MAX_A if maps_pocket_max_A is None else float(maps_pocket_max_A)
    whole_map_npts_max = int(whole_map_npts_max)
    if not 25 <= whole_map_npts_max <= AUTOGRID4_NPTS_MAX:
        raise ValueError(
            "CaV-EMPS whole_map_npts_max must be between 25 and "
            f"{AUTOGRID4_NPTS_MAX}, found {whole_map_npts_max}"
        )
    autogrid_map_types = _normalize_autogrid_map_types(map_types)
    cav_emps_policies = {
        "internal",
        "surface",
        "hybrid",
        "receptor_search",
        "exhaustive_search",
    }
    if policy in cav_emps_policies:
        full_types = tuple(_AD4_TYPES)
        if autogrid_map_types != full_types:
            raise ValueError(
                "CaV-EMPS requires full AD4 ligand_types because "
                "AutoGrid dsolvmap depends on the requested type set."
            )
        autogrid_map_types = full_types
    effective_min_sep_A = _effective_hotspot_min_sep_A(
        min_sep_A,
        hotspot_box_ang,
        box_fraction=nms_box_fraction,
        min_A=nms_min_A,
        max_A=nms_max_A,
    )
    candidate_min_sep_A = _candidate_harvest_min_sep_A(effective_min_sep_A, expected_site_count)
    requested_meta = {
        "policy": policy,
        "site_count": str(expected_site_count),
        "box_side_A": f"{float(hotspot_box_ang):.3f}",
        "min_sep_A": f"{effective_min_sep_A:.3f}",
        "candidate_min_sep_A": f"{candidate_min_sep_A:.3f}",
        "surface_relaxation": "progressive_tau_budget_v2",
        "surface_center_refinement": "edt_component_core_centroid_v3",
        "surface_region_centroid": "cluster_lobes_v1",
        "surface_contact_gate": "refined_fixed_shell_v1",
        "core_reserve": "replace_last_v1",
        "axis_reserve": "replace_penultimate_v1",
        "underfilled_rescue": "edt_tail_fill_v1",
        "geometry_reserve_fill": "underfilled_tail_v1",
        "surface_scoring": "isotropic_gaussian_CED_edt_v1",
        "surface_kernel": f"gaussian_sigma_A={MAPS_GAUSSIAN_SIGMA_A:.3f}",
        "surface_weights": f"C={MAPS_C_WEIGHT:.3f},E={MAPS_E_WEIGHT:.3f},D={MAPS_D_WEIGHT:.3f}",
        "surface_edt_weight": f"floor={MAPS_EDT_WEIGHT_FLOOR:.3f},power={MAPS_EDT_WEIGHT_POWER:.3f}",
        "ranking_score": "common_regional_physics_v1",
        "ranking_normalization": "receptor_field_max_v1",
        "legacy_ranking_score": "portfolio_role_prior_v3",
        "map_loader": "autogrid_dim_plus_one_exact_dsolv_v2",
        "dsolv_context": "full_ad4_ligand_types_v1",
        "autogrid_ligand_types": ",".join(autogrid_map_types),
        "autogrid_map_types": ",".join(autogrid_map_types),
        "whole_map_policy": "adaptive_spacing_autogrid_cap_v1",
        "whole_map_npts_max": str(whole_map_npts_max),
        "r_min_A": "auto" if requested_r_min is None else f"{requested_r_min:.3f}",
        "rmin_floor_A": f"{adaptive_params['floor_A']:.3f}",
        "rmin_ceil_A": f"{adaptive_params['ceil_A']:.3f}",
        "rmin_peak_pct": f"{adaptive_params['peak_percentile']:.3f}",
        "rmin_zone_max_A": f"{adaptive_params['pocket_zone_max_A']:.3f}",
        "maps_pocket_max_A": f"{maps_pocket_max_value:.3f}",
    }

    def _metadata_matches(meta: dict[str, str]) -> bool:
        return all(str(meta.get(key, "")) == value for key, value in requested_meta.items())

    # 0) Reuse if present
    if centers_tsv_path.exists():
        rows = _parse_centers_tsv(centers_tsv_path, receptor_key=rec_stem)
        meta = _parse_centers_metadata(centers_tsv_path)
        if rows and _metadata_matches(meta):
            print(f"[centers] reusing {centers_tsv_path} with {len(rows)} rows")
            return str(centers_tsv_path)
        if rows:
            print(
                f"[centers] regenerating {centers_tsv_path}: "
                f"stored policy={meta.get('policy', 'unknown')} site_count={meta.get('site_count', 'unknown')} "
                f"requested policy={policy} site_count={expected_site_count}"
            )

    # 1) Find (or build) whole-protein maps (.fld)
    fld_candidates = sorted(out_root.glob("**/*.fld"))
    if fld_candidates:
        if len(fld_candidates) != 1:
            raise RuntimeError(
                f"Expected one compatible FLD, found {len(fld_candidates)}: "
                f"{fld_candidates}"
            )
        fld = fld_candidates[0]
        existing_meta = load_fld_or_map_meta(str(fld))
        if tuple(existing_meta["affinity_types"]) != tuple(autogrid_map_types):
            raise RuntimeError(
                f"Incompatible cached FLD {fld}: expected affinity types "
                f"{autogrid_map_types}, found "
                f"{existing_meta['affinity_types']}. Use a clean work directory."
            )
        _, _, expected_spacing = compute_whole_box_auto(
            receptor_pdbqt=receptor_pdbqt,
            base_spacing=float(default_spacing),
            margin_A=8.0,
            npts_max=whole_map_npts_max,
        )
        expected_spacing = float(f"{expected_spacing:.3f}")
        if not np.isclose(float(existing_meta["spacing"]), expected_spacing, rtol=0.0, atol=1e-9):
            raise RuntimeError(
                f"Incompatible cached FLD {fld}: CaV-EMPS expects spacing "
                f"{expected_spacing:.3f} A with npts_max={whole_map_npts_max}, found "
                f"{float(existing_meta['spacing']):.3f} A. Use a clean work directory."
            )
    else:
        # bootstrap: generate whole-protein maps once
        fld_path = ensure_whole_protein_maps(
            receptor_pdbqt=receptor_pdbqt,
            out_root=out_root,
            spacing=default_spacing,
            cap_ang=blind_cap,
            autogrid4_bin=autogrid4_bin,
            map_types=autogrid_map_types,
            npts_max=whole_map_npts_max,
        )
        fld_candidates = [Path(fld_path)]

    fld = fld_candidates[0]
    meta = load_fld_or_map_meta(str(fld))
    origin, spacing, shape = meta["origin"], meta["spacing"], meta["shape"]
    mp = meta["map_paths"]
    C = load_map_ascii(mp["C"], shape)
    E = load_map_ascii(mp["E"], shape)
    D = load_map_ascii(mp["D"], shape)
    map_diagnostics = {}
    map_diagnostics.update(_map_diagnostics("C", mp["C"], C))
    map_diagnostics.update(_map_diagnostics("E", mp["E"], E))
    map_diagnostics.update(_map_diagnostics("D", mp["D"], D))

    profiled_r_min = requested_r_min
    r_min_is_profiled = False
    if profiled_r_min is None:
        profiled_r_min = _profile_adaptive_r_min_from_maps(
            receptor_pdbqt,
            origin,
            spacing,
            shape,
            adaptive_params=adaptive_params,
        )
        r_min_is_profiled = profiled_r_min is not None

    if effective_min_sep_A > float(min_sep_A) + 1e-6:
        print(
            f"[nms] clamped site separation from {float(min_sep_A):.2f} Å "
            f"to {effective_min_sep_A:.2f} Å for {float(hotspot_box_ang):.1f} Å boxes"
        )
    if candidate_min_sep_A + 1e-6 < effective_min_sep_A:
        print(
            f"[nms] harvesting candidates at {candidate_min_sep_A:.2f} Å; "
            f"final site separation target is {effective_min_sep_A:.2f} Å"
        )

    if policy in {"receptor_search", "hybrid", "exhaustive_search"}:
        candidate_budget = max(int(n_sites), min(96, max(40, int(n_sites) * 8)))
    else:
        candidate_budget = max(3, int(n_sites))
    surface_target_count = min(
        candidate_budget,
        max(expected_site_count, int(expected_site_count) * 2),
    )
    internal_sites: list[dict] = []
    surface_sites: list[dict] = []
    hybrid_sites: list[dict] = []
    ranking_score_context: dict | None = None

    if policy in {"internal", "hybrid", "receptor_search", "exhaustive_search"}:
        try:
            internal_sites = pick_centers(
                receptor_pdbqt,
                {"C": C, "E": E, "D": D},
                map_origin=origin,
                map_spacing=spacing,
                voxel_spacing=0.375,
                r_min=profiled_r_min,
                r_min_is_profiled=r_min_is_profiled,
                adaptive_r_min_params=adaptive_params,
                min_sep_A=candidate_min_sep_A,
                k_box=float(k_box),
                max_sites=candidate_budget,
                inflate_A=0.0,
            )
        except Exception as exc:
            print(f"[centers/internal] skipped due to error: {exc}")

    def _collect_surface_sites(threshold: float) -> tuple[list[dict], dict | None]:
        return detect_maps_hotspots(
            C,
            E,
            D,
            origin,
            spacing,
            tau_rel=float(threshold),
            min_sep_A=max(candidate_min_sep_A, float(SURFACE_NMS_MINSEP_A)),
            max_sites=candidate_budget,
            box_side_A=float(hotspot_box_ang),
            r_min_A=profiled_r_min,
            adaptive_r_min_params=adaptive_params,
            pocket_max_A=maps_pocket_max_value,
            receptor_pdbqt=receptor_pdbqt,
            min_peak_candidates=candidate_budget,
            return_score_context=True,
        )

    if policy in {"surface", "hybrid", "receptor_search", "exhaustive_search"}:
        print("[centers] collecting surface-pocket candidates")
        surface_sites, ranking_score_context = _collect_surface_sites(float(tau_rel))
        if len(surface_sites) < surface_target_count:
            print(
                f"[centers] surface search underfilled after adaptive thresholding "
                f"({len(surface_sites)}/{surface_target_count})"
            )

    surface_region_sites: list[dict] = []
    if surface_sites and policy in {"hybrid", "receptor_search", "exhaustive_search"}:
        surface_region_sites = _build_surface_region_sites(
            surface_sites,
            cluster_radius_A=float(hotspot_box_ang) * 0.60,
            min_sep_A=candidate_min_sep_A,
            max_sites=candidate_budget,
        )
        if surface_region_sites:
            print(f"[centers] built {len(surface_region_sites)} surface-region candidates")

    if policy in {"hybrid", "receptor_search", "exhaustive_search"}:
        hybrid_sites = _build_hybrid_sites(
            internal_sites,
            surface_sites,
            min_sep_A=max(candidate_min_sep_A, float(MAX_CENTER_DIST_A) * 0.5),
            max_sites=candidate_budget,
        )
        if surface_region_sites:
            hybrid_sites = _rank_union_sites(surface_region_sites, hybrid_sites)[:candidate_budget]

    if policy == "internal":
        sites = _relabel_sites(internal_sites[: int(n_sites)])
    elif policy == "surface":
        sites = _relabel_sites(surface_sites[: int(n_sites)])
    elif policy == "hybrid":
        sites = _relabel_sites(hybrid_sites[: int(n_sites)])
    elif policy == "receptor_search":
        sites = _assemble_receptor_search_sites(
            _normalize_ranked_sites(internal_sites, family="internal"),
            _normalize_ranked_sites(surface_sites, family="surface"),
            _normalize_ranked_sites(hybrid_sites, family="hybrid"),
            min_sep_A=effective_min_sep_A,
            target_count=expected_site_count,
        )
    elif policy == "exhaustive_search":
        sites = _assemble_exhaustive_search_sites(
            _normalize_ranked_sites(internal_sites, family="internal"),
            _normalize_ranked_sites(surface_sites, family="surface"),
            _normalize_ranked_sites(hybrid_sites, family="hybrid"),
            min_sep_A=candidate_min_sep_A,
            target_count=int(n_sites),
        )
    else:
        raise ValueError(f"Unsupported center-generation policy: {policy}")

    if policy == "receptor_search" and expected_site_count >= 6 and len(sites) >= expected_site_count:
        core_site = _receptor_core_reserve_site(
            receptor_pdbqt,
            box_side_A=float(hotspot_box_ang),
            spacing=float(spacing),
        )
        axis_site = _receptor_axis_reserve_site(
            receptor_pdbqt,
            preserved_sites=sites[: max(0, expected_site_count - 2)],
            box_side_A=float(hotspot_box_ang),
            spacing=float(spacing),
        )
        sites = _replace_tail_with_reserve_sites(
            sites,
            axis_site=axis_site,
            core_site=core_site,
            min_recenter_A=2.0,
        )

    # 5) If normal receptor-search families found some hypotheses but fewer
    # than requested, fill the tail with geometry-only EDT surface-cleft sites.
    # This preserves the existing ranked hypotheses and improves recall for
    # receptors where map scoring underfills the requested portfolio.
    if policy == "receptor_search" and 0 < len(sites) < expected_site_count:
        try:
            pocket_min_A, shell_max_A = _maps_pocket_shell_bounds(
                r_min_A=profiled_r_min,
                pocket_max_A=maps_pocket_max_value,
            )
            rescue_sites = _edt_surface_sites(
                receptor_pdbqt,
                origin,
                spacing,
                shape,
                min_sep_A=candidate_min_sep_A,
                max_sites=max(expected_site_count * 4, 24),
                box_side_A=float(hotspot_box_ang),
                pocket_min_A=pocket_min_A,
                pocket_max_A=shell_max_A,
            )
            added = _append_distinct_tail_sites(
                sites,
                rescue_sites,
                target_count=expected_site_count,
                portfolio_role="edt_tail_rescue",
                min_sep_A=_portfolio_min_sep_A(effective_min_sep_A),
            )
            if added:
                sites = _relabel_sites(sites)
                print(f"[centers/rescue] added {added} EDT tail site(s)")
        except Exception as exc:
            print(f"[centers/rescue] EDT tail fill skipped due to error: {exc}")

    # 6) Final fallback: keep the blind box, then add EDT-surface rescue sites.
    # This path is intentionally narrow: it only runs when all normal
    # receptor-derived families failed. Keeping blind as S1 preserves the old
    # fallback behavior while allowing extra hypotheses for difficult cases
    # such as transmembrane channels.
    if len(sites) == 0:
        print("[centers] search yielded no sites; falling back to blind box plus EDT rescue")
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
                family="blind_fallback",
                portfolio_role="blind_fallback",
            )
        ]
        if policy in {"surface", "hybrid", "receptor_search", "exhaustive_search"}:
            try:
                pocket_min_A, shell_max_A = _maps_pocket_shell_bounds(
                    r_min_A=profiled_r_min,
                    pocket_max_A=maps_pocket_max_value,
                )
                rescue_sites = _edt_surface_sites(
                    receptor_pdbqt,
                    origin,
                    spacing,
                    shape,
                    min_sep_A=candidate_min_sep_A,
                    max_sites=max(expected_site_count * 4, 24),
                    box_side_A=float(hotspot_box_ang),
                    pocket_min_A=pocket_min_A,
                    pocket_max_A=shell_max_A,
                )
                _append_distinct_tail_sites(
                    sites,
                    rescue_sites,
                    target_count=expected_site_count,
                    portfolio_role="edt_tail_rescue",
                    min_sep_A=_portfolio_min_sep_A(effective_min_sep_A),
                )
                sites = _relabel_sites(sites)
            except Exception as exc:
                print(f"[centers/rescue] EDT rescue skipped due to error: {exc}")

    # 7) Last-resort geometry reserve fill. This is deliberately lower priority
    # than map, hybrid, internal, and EDT rescue hypotheses; it exists to honor
    # the requested portfolio size without promoting reserve boxes as confident
    # binding-site predictions.
    if policy == "receptor_search" and 0 < len(sites) < expected_site_count:
        reserve_candidates = []
        axis_sites = _receptor_axis_reserve_sites(
            receptor_pdbqt,
            preserved_sites=sites,
            box_side_A=float(hotspot_box_ang),
            spacing=float(spacing),
        )
        core_site = _receptor_core_reserve_site(
            receptor_pdbqt,
            box_side_A=float(hotspot_box_ang),
            spacing=float(spacing),
        )
        if core_site is not None:
            reserve_candidates.append(core_site)
        reserve_candidates.extend(axis_sites)
        added = _append_distinct_tail_sites(
            sites,
            reserve_candidates,
            target_count=expected_site_count,
            portfolio_role="reserve",
            min_sep_A=_portfolio_min_sep_A(effective_min_sep_A),
        )
        if added:
            sites = _relabel_sites(sites)
            print(f"[centers/reserve] added {added} geometry reserve site(s)")

    if policy == "receptor_search":
        sites = _assign_receptor_search_fitness_scores(sites)
        sites = _assign_common_physics_ranking_scores(
            sites,
            score_context=ranking_score_context,
        )

    # 8) Write TSV
    diagnostic_columns = [
        "original_peak_x",
        "original_peak_y",
        "original_peak_z",
        "refined_center_x",
        "refined_center_y",
        "refined_center_z",
        "refinement_shift_A",
        "refinement_component_size",
        "refinement_component_count",
        "internal_anchor_x",
        "internal_anchor_y",
        "internal_anchor_z",
        "surface_anchor_x",
        "surface_anchor_y",
        "surface_anchor_z",
        "anchor_separation_A",
        "hybrid_center_x",
        "hybrid_center_y",
        "hybrid_center_z",
        "cluster_size",
        "cluster_spread_A",
        "unprojected_centroid_x",
        "unprojected_centroid_y",
        "unprojected_centroid_z",
    ]

    def _format_diagnostic_value(value) -> str:
        if value is None:
            return ""
        if isinstance(value, str):
            return value
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        if isinstance(value, (float, np.floating)):
            if not np.isfinite(float(value)):
                return ""
            return f"{float(value):.3f}"
        return str(value)

    with open(centers_tsv_path, "w") as f:
        f.write(
            "# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tr_peak\tF\t"
            "raw_F\tlegacy_F\tcommon_physics_score\tcommon_physics_raw\tcommon_edt_A\t"
            "ranking_basis\tfamily\tportfolio_role\tselection_score\tcenter_closeness\t"
            + "\t".join(diagnostic_columns)
            + "\n"
        )
        write_meta = dict(requested_meta)
        write_meta.update(map_diagnostics)
        write_meta["site_count"] = str(len(sites))
        f.write("# meta " + " ".join(f"{key}={value}" for key, value in write_meta.items()) + "\n")
        for s in sites:
            f.write(
                f"{rec_stem}\t{s['site_id']}\t{s['cx']:.3f}\t{s['cy']:.3f}\t{s['cz']:.3f}\t"
                f"{s['nx']}\t{s['ny']}\t{s['nz']}\t{s['spacing']:.3f}\t"
                f"{s.get('r_peak', 2.5):.2f}\t{s.get('F', 1.0):.6f}\t"
                f"{s.get('raw_F', s.get('F', 1.0)):.6f}\t"
                f"{s.get('legacy_F', s.get('F', 1.0)):.6f}\t"
                f"{s.get('common_physics_score', s.get('F', 1.0)):.6f}\t"
                f"{s.get('common_physics_raw', 0.0):.8f}\t"
                f"{s.get('common_edt_A', 0.0):.4f}\t{s.get('ranking_basis', '')}\t"
                f"{s.get('family', '')}\t{s.get('portfolio_role', '')}\t"
                f"{s.get('selection_score', '')}\t{s.get('center_closeness', '')}\t"
                + "\t".join(_format_diagnostic_value(s.get(key, "")) for key in diagnostic_columns)
                + "\n"
            )
    print(f"[centers] wrote {len(sites)} --> {centers_tsv_path}")
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

        expected_fld = site_dir / f"{rec_stem}.maps.fld"
        if not expected_fld.exists():
            fld_candidates = sorted(site_dir.glob("*.fld"))
            raise FileNotFoundError(
                f"Expected AutoGrid FLD {expected_fld}, found {fld_candidates}; see {log_path}"
            )
        fld_path = str(expected_fld.resolve())

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


def _positive_signal_unit_interval(values: np.ndarray, *, floor: float = 1e-6) -> np.ndarray:
    """Robustly scale nonnegative map support to [0, 1]."""
    values = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(values)
    signal = np.zeros_like(values, dtype=np.float32)
    positive = values[finite & (values > floor)]
    if positive.size < 8:
        return signal
    lo = float(np.percentile(positive, 5.0))
    hi = float(np.percentile(positive, 95.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo + 1e-6:
        lo = 0.0
        hi = float(np.max(positive))
    if hi <= lo + 1e-6:
        return signal
    signal[finite] = np.clip((values[finite] - lo) / (hi - lo + 1e-6), 0.0, 1.0)
    return signal.astype(np.float32)


def _normalize_positive_weights(weights: dict[str, float]) -> dict[str, float]:
    cleaned = {key: max(0.0, float(value)) for key, value in weights.items()}
    total = sum(cleaned.values())
    if total <= 0.0:
        return {key: 1.0 / len(cleaned) for key in cleaned}
    return {key: value / total for key, value in cleaned.items()}


def _continuous_edt_weight(
    edt_grid: np.ndarray,
    *,
    pocket_min_A: float,
    pocket_max_A: float,
    floor: float = MAPS_EDT_WEIGHT_FLOOR,
    power: float = MAPS_EDT_WEIGHT_POWER,
) -> np.ndarray:
    """Continuous enclosure/depth weight for pocket voxels."""
    span = max(1e-6, float(pocket_max_A) - float(pocket_min_A))
    depth = np.clip((edt_grid.astype(np.float32) - float(pocket_min_A)) / span, 0.0, 1.0)
    floor = float(np.clip(floor, 0.0, 1.0))
    return floor + (1.0 - floor) * np.power(depth, max(0.1, float(power)))


def _maps_are_degenerate(C, E, D) -> bool:
    """
    Return True if the AutoGrid maps lack enough signal for meaningful
    scoring.  Common causes: AutoGrid ran without a GPF parameter_file
    directive, or the receptor has unusual atom types.

    We check:
      - E-map all zero / constant → no electrostatic signal
      - C/D maps have effectively no favorable or desolvation signal anywhere

    Important: whole-protein affinity maps are expected to be sparse.  A low
    nonzero fraction by itself is not evidence of a broken grid; large positive
    steric walls also make raw C-map ranges look extreme even when the map is
    valid.  We therefore look for the presence of a meaningful favorable tail,
    not density of nonzero voxels.
    """
    e_range = float(np.nanmax(E)) - float(np.nanmin(E))
    c_favorable_frac = float(np.mean(np.asarray(C) < -0.01))
    d_signal_frac = float(np.mean(np.abs(np.asarray(D)) > 0.01))

    if e_range < 1e-6:
        print("[maps/degenerate] E-map is constant (all zeros) — maps scoring is unreliable")
        return True
    if c_favorable_frac < 1e-5 and d_signal_frac < 1e-5:
        print(
            "[maps/degenerate] C/D maps have no meaningful favorable signal "
            f"(C< -0.01 frac={c_favorable_frac:.2e}, |D|>0.01 frac={d_signal_frac:.2e})"
        )
        return True
    return False


def _maps_pocket_shell_bounds(
    r_min_A: float | None = R_MIN_CAVITY_A,
    pocket_max_A: float | None = None,
) -> tuple[float, float]:
    """
    Derive the maps-mode EDT shell from the receptor-specific cavity threshold.

    `R_MIN_CAVITY_A` is the minimum radius for accepting an internal cavity
    centre. Surface-contact voxels for real ligand poses sit closer to the
    protein wall than that, so maps-mode should search a shallower shell rather
    than reuse the full cavity-centre cutoff. We therefore anchor the lower EDT
    bound to half of the profiled cavity threshold, with conservative clamps to
    avoid collapsing into pure surface noise.
    """
    if r_min_A is None and R_MIN_CAVITY_A is not None:
        r_min_A = float(R_MIN_CAVITY_A)
    elif r_min_A is None:
        r_min_A = 3.0
    if pocket_max_A is None:
        pocket_max_A = MAPS_POCKET_MAX_A
    lower = float(np.clip(0.5 * float(r_min_A), 0.75, 2.5))
    upper = min(float(SURFACE_SHELL__MAX_A), max(float(pocket_max_A), lower + 2.0))
    return lower, upper


def _effective_hotspot_min_sep_A(
    min_sep_A: float,
    box_side_A: float,
    *,
    box_fraction: float = HOTSPOT_NMS_BOX_FRACTION,
    min_A: float = HOTSPOT_NMS_MIN_A,
    max_A: float = HOTSPOT_NMS_MAX_A,
) -> float:
    """
    Clamp site NMS to the docking-box scale so automatic sites do not collapse
    into several nearly identical boxes around the same pocket.
    """
    requested = max(0.0, float(min_sep_A))
    box_scaled = float(np.clip(float(box_fraction) * float(box_side_A), float(min_A), float(max_A)))
    return round(max(requested, box_scaled), 2)


def _adaptive_r_min_params(overrides: dict | None = None) -> dict:
    params = {
        "percentile": ADAPTIVE_R_MIN_PERCENTILE,
        "peak_percentile": ADAPTIVE_R_MIN_PEAK_PERCENTILE,
        "floor_A": ADAPTIVE_R_MIN_FLOOR_A,
        "ceil_A": ADAPTIVE_R_MIN_CEIL_A,
        "pocket_zone_max_A": ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A,
        "peak_window_A": ADAPTIVE_R_MIN_PEAK_WINDOW_A,
    }
    if overrides:
        for key in params:
            if key in overrides and overrides[key] is not None:
                params[key] = float(overrides[key])
    if params["ceil_A"] < params["floor_A"]:
        params["ceil_A"] = params["floor_A"]
    if params["pocket_zone_max_A"] < params["floor_A"]:
        params["pocket_zone_max_A"] = params["floor_A"]
    return params


def _coerce_optional_float(value) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and value.strip().lower() in {"", "none", "auto", "adaptive"}:
        return None
    return float(value)


def _profile_adaptive_r_min_from_maps(
    receptor_pdbqt: str | Path | None,
    origin,
    spacing,
    shape,
    *,
    adaptive_params: dict | None = None,
) -> float | None:
    if receptor_pdbqt is None:
        return None
    try:
        atoms = pdbqt_atoms(receptor_pdbqt)
        occ, _, _, sp = make_grid_from_map(atoms, origin, spacing, shape, inflate_A=0.0)
        return adaptive_r_min_cavity(occ, sp, **_adaptive_r_min_params(adaptive_params))
    except Exception as exc:
        print(f"[r_min] adaptive profiling failed: {exc}")
        return None


def _site_center(site: dict) -> np.ndarray:
    return np.array([site["cx"], site["cy"], site["cz"]], dtype=np.float32)


def _site_score_value(site: dict) -> float:
    return float(site.get("selection_score", site.get("F", 0.0)))


def _site_surface_focus_score(site: dict) -> float:
    centrality = float(site.get("center_closeness", 0.0))
    affinity = _site_score_value(site)
    return 0.85 * centrality + 0.15 * affinity


def _site_hybrid_focus_score(site: dict) -> float:
    centrality = float(site.get("center_closeness", 0.0))
    affinity = _site_score_value(site)
    return 0.70 * centrality + 0.30 * affinity


def _clone_site(site: dict, **updates) -> dict:
    cloned = dict(site)
    cloned.update(updates)
    return cloned


def _relabel_sites(sites: list[dict]) -> list[dict]:
    relabeled = []
    for idx, site in enumerate(sites, start=1):
        relabeled.append(_clone_site(site, site_id=f"S{idx}"))
    return relabeled


def _normalize_ranked_sites(sites: list[dict], *, family: str) -> list[dict]:
    if not sites:
        return []
    scores = np.array([float(site.get("F", 0.0)) for site in sites], dtype=np.float32)
    if np.allclose(scores.max(), scores.min()):
        norm = np.ones_like(scores, dtype=np.float32)
    else:
        norm = (scores - scores.min()) / (scores.max() - scores.min() + 1e-6)
    ranked = []
    for idx, (site, norm_score) in enumerate(zip(sites, norm), start=1):
        ranked.append(
            _clone_site(
                site,
                family=family,
                family_rank=idx,
                selection_score=float(norm_score),
            )
        )
    return ranked


def _take_best_distinct(
    candidates: list[dict],
    selected: list[dict],
    *,
    min_sep_A: float,
) -> dict | None:
    for site in candidates:
        ctr = _site_center(site)
        if all(np.linalg.norm(ctr - _site_center(existing)) >= min_sep_A for existing in selected):
            return site
    return None


def _append_distinct_tail_sites(
    selected: list[dict],
    candidates: list[dict],
    *,
    target_count: int,
    portfolio_role: str,
    min_sep_A: float,
) -> int:
    """
    Fill only the tail of an underfilled portfolio.

    The normal portfolio separation is intentionally conservative for docking
    boxes.  When a difficult receptor underfills, ask for many rescue candidates
    and relax the tail separation gradually rather than returning too few sites.
    """
    if len(selected) >= int(target_count) or not candidates:
        return 0
    primary = max(0.0, float(min_sep_A))
    schedule = [
        primary,
        max(5.0, 0.75 * primary),
        max(3.0, 0.50 * primary),
        2.0,
    ]
    added = 0
    for site in candidates:
        if len(selected) >= int(target_count):
            return added
        for sep in schedule:
            if _take_best_distinct([site], selected, min_sep_A=sep) is not None:
                selected.append(_clone_site(site, portfolio_role=portfolio_role))
                added += 1
                break
    return added


def _rank_union_sites(*site_groups: list[dict]) -> list[dict]:
    union = [site for group in site_groups for site in group]
    return sorted(union, key=lambda site: _site_score_value(site), reverse=True)


def _rank_surface_focus_sites(surface_sites: list[dict]) -> list[dict]:
    return sorted(surface_sites, key=_site_surface_focus_score, reverse=True)


def _rank_hybrid_focus_sites(hybrid_sites: list[dict]) -> list[dict]:
    return sorted(hybrid_sites, key=_site_hybrid_focus_score, reverse=True)


def _rank_hybrid_low_score_focus_sites(hybrid_sites: list[dict]) -> list[dict]:
    if len(hybrid_sites) <= 2:
        return _rank_hybrid_focus_sites(hybrid_sites)
    split_idx = max(1, len(hybrid_sites) // 2)
    low_score_sites = hybrid_sites[split_idx:]
    high_score_sites = hybrid_sites[:split_idx]
    return (
        sorted(low_score_sites, key=_site_hybrid_focus_score, reverse=True)
        + sorted(high_score_sites, key=_site_hybrid_focus_score, reverse=True)
    )


def _normalized_site_values(sites: list[dict], field: str) -> dict[int, float]:
    values = np.array([float(site.get(field, 0.0) or 0.0) for site in sites], dtype=np.float32)
    if values.size == 0:
        return {}
    if np.allclose(float(values.max()), float(values.min())):
        return {index: 0.5 for index in range(len(sites))}
    norm = (values - values.min()) / (values.max() - values.min() + 1e-6)
    return {index: float(value) for index, value in enumerate(norm)}


def _assign_receptor_search_fitness_scores(sites: list[dict]) -> list[dict]:
    """
    Preserve the historical v3 role-prior score for benchmark ablations.

    Raw map/internal scores are useful inside a family, but not comparable
    across surface peaks, hybrid consensus sites, and reserve/core boxes. This
    score is retained in ``legacy_F`` by the common-physics ranker and is no
    longer the publication ranking score.
    """
    if not sites:
        return sites

    role_prior = {
        "surface_focus": 0.92,
        "surface_primary": 0.88,
        "hybrid_primary": 0.82,
        "internal_primary": 0.76,
        "hybrid_low_score_focus": 0.68,
        "surface": 0.64,
        "edt_tail_rescue": 0.42,
        "reserve": 0.22,
        "axis_reserve": 0.14,
        "core_reserve": 0.10,
        "blind_fallback": 0.05,
    }
    raw_norm = _normalized_site_values(sites, "F")
    r_norm = _normalized_site_values(sites, "r_peak")
    calibrated: list[dict] = []
    for index, site in enumerate(sites):
        raw_f = float(site.get("F", 0.0) or 0.0)
        role = str(site.get("portfolio_role") or site.get("family") or "reserve")
        centrality = float(site.get("center_closeness", 0.0) or 0.0)
        contact = float(site.get("contact_frac", 0.0) or 0.0)
        cluster_support = min(1.0, float(site.get("cluster_size", 1) or 1) / 4.0)
        score = (
            role_prior.get(role, 0.20)
            + 0.10 * raw_norm.get(index, 0.0)
            + 0.04 * r_norm.get(index, 0.0)
            + 0.04 * centrality
            + 0.03 * contact
            + 0.02 * cluster_support
        )
        calibrated.append(
            _clone_site(
                site,
                raw_F=raw_f,
                ranking_score=float(score),
                F=float(score),
            )
        )
    return calibrated


def _sample_site_grid_value(site: dict, grid: np.ndarray | None, origin, spacing: float) -> float:
    if grid is None:
        return 0.0
    value = trilinear_sample(
        grid,
        (float(site["cx"]), float(site["cy"]), float(site["cz"])),
        origin,
        float(spacing),
    )
    return float(value) if np.isfinite(value) else 0.0


def _assign_common_physics_ranking_scores(
    sites: list[dict],
    *,
    score_context: dict | None,
) -> list[dict]:
    """
    Rank every emitted center in one shared physical scoring space.

    Candidate-family scores remain useful while harvesting surface, internal,
    and hybrid hypotheses, but their numerical scales are not comparable.  At
    the portfolio boundary, sample the same regional C/E/D/EDT field at every
    final center and normalize only by the strongest field value in that
    receptor.  Portfolio roles are retained as diagnostics and do not
    contribute to the new ranking score.

    If map scoring is unavailable or wholly degenerate, use the same EDT shell
    support at every center as a deterministic geometry-only fallback.
    """
    if not sites:
        return sites

    context = score_context or {}
    origin = context.get("origin")
    spacing = float(context.get("spacing", 1.0))
    score_field = context.get("score_field")
    edt_grid = context.get("edt_grid")
    pocket_min_A = float(context.get("pocket_min_A", 0.0))
    pocket_max_A = float(context.get("pocket_max_A", np.inf))

    if origin is None:
        origin = (0.0, 0.0, 0.0)

    raw_physics = [
        max(0.0, _sample_site_grid_value(site, score_field, origin, spacing))
        for site in sites
    ]
    edt_values = [
        max(0.0, _sample_site_grid_value(site, edt_grid, origin, spacing))
        for site in sites
    ]

    field_max = float(context.get("score_max", 0.0) or 0.0)
    if not np.isfinite(field_max) or field_max <= 1e-12:
        field_max = max(raw_physics, default=0.0)

    if field_max > 1e-12:
        scores = [float(np.clip(value / field_max, 0.0, 1.0)) for value in raw_physics]
        ranking_basis = "regional_CED_EDT"
    else:
        geometry_support = [
            value if pocket_min_A <= value <= pocket_max_A else 0.0
            for value in edt_values
        ]
        geometry_max = max(geometry_support, default=0.0)
        if geometry_max > 1e-12:
            scores = [float(value / geometry_max) for value in geometry_support]
        else:
            scores = [0.0 for _ in sites]
        ranking_basis = "EDT_fallback"

    ranked = []
    for site, score, raw_score, edt_A in zip(sites, scores, raw_physics, edt_values):
        legacy_score = float(site.get("F", 0.0) or 0.0)
        ranked.append(
            _clone_site(
                site,
                legacy_F=legacy_score,
                legacy_ranking_score=float(site.get("ranking_score", legacy_score) or 0.0),
                common_physics_score=float(score),
                common_physics_raw=float(raw_score),
                common_edt_A=float(edt_A),
                ranking_basis=ranking_basis,
                ranking_score=float(score),
                F=float(score),
            )
        )
    return ranked


def _build_surface_region_sites(
    surface_sites: list[dict],
    *,
    cluster_radius_A: float,
    min_sep_A: float,
    max_sites: int,
) -> list[dict]:
    """
    Build docking-box centers from clusters of nearby surface hotspots.

    Surface hotspots are wall/contact points. A practical docking box often
    belongs between several such contact lobes, especially in elongated ligand
    pockets, so create additional region-centroid hypotheses from local hotspot
    clusters.
    """
    if len(surface_sites) < 2:
        return []

    ordered = sorted(surface_sites, key=_site_score_value, reverse=True)
    regions = []
    for anchor in ordered:
        anchor_center = _site_center(anchor)
        neighbors = [
            site for site in ordered
            if np.linalg.norm(_site_center(site) - anchor_center) <= float(cluster_radius_A)
        ]
        if len(neighbors) < 2:
            continue

        centers = np.array([_site_center(site) for site in neighbors], dtype=np.float32)
        centroid = centers.mean(axis=0)
        cluster_spread_A = float(np.max(np.linalg.norm(centers - centroid, axis=1)))
        f_values = [float(site.get("F", 0.0)) for site in neighbors]
        r_values = [float(site.get("r_peak", 0.0)) for site in neighbors]
        template = max(neighbors, key=_site_score_value)
        support_score = float(np.mean(f_values) * np.sqrt(len(neighbors)))
        regions.append(
            _clone_site(
                template,
                cx=float(centroid[0]),
                cy=float(centroid[1]),
                cz=float(centroid[2]),
                family="surface_region",
                F=support_score,
                selection_score=support_score,
                r_peak=float(np.mean(r_values)) if r_values else float(template.get("r_peak", 1.0)),
                cluster_size=len(neighbors),
                cluster_spread_A=cluster_spread_A,
                unprojected_centroid_x=float(centroid[0]),
                unprojected_centroid_y=float(centroid[1]),
                unprojected_centroid_z=float(centroid[2]),
                center_closeness=max(
                    float(site.get("center_closeness", 0.0)) for site in neighbors
                ),
            )
        )

    regions.sort(key=lambda site: (int(site.get("cluster_size", 0)), _site_score_value(site)), reverse=True)
    selected = []
    for site in regions:
        if _take_best_distinct([site], selected, min_sep_A=float(min_sep_A)) is not None:
            selected.append(site)
        if len(selected) >= int(max_sites):
            break
    return selected


def _receptor_core_reserve_site(
    receptor_pdbqt,
    *,
    box_side_A: float,
    spacing: float,
) -> dict | None:
    try:
        center, _, sp = blind_box(receptor_pdbqt, cap=float(box_side_A), spacing=float(spacing))
    except Exception as exc:
        print(f"[centers/core] skipped due to error: {exc}")
        return None
    n = _npts_for_box_side(float(box_side_A), float(spacing))
    return dict(
        site_id="S_core",
        cx=float(center[0]),
        cy=float(center[1]),
        cz=float(center[2]),
        nx=n,
        ny=n,
        nz=n,
        spacing=float(sp),
        r_peak=max(1.0, float(box_side_A) * 0.15),
        F=0.0,
        family="core_reserve",
    )


def _receptor_axis_reserve_site(
    receptor_pdbqt,
    *,
    preserved_sites: list[dict],
    box_side_A: float,
    spacing: float,
    offset_fraction: float = 0.22,
    offset_min_A: float = 8.0,
    offset_max_A: float = 14.0,
) -> dict | None:
    sites = _receptor_axis_reserve_sites(
        receptor_pdbqt,
        preserved_sites=preserved_sites,
        box_side_A=box_side_A,
        spacing=spacing,
        offset_fraction=offset_fraction,
        offset_min_A=offset_min_A,
        offset_max_A=offset_max_A,
    )
    return sites[0] if sites else None


def _receptor_axis_reserve_sites(
    receptor_pdbqt,
    *,
    preserved_sites: list[dict],
    box_side_A: float,
    spacing: float,
    offset_fraction: float = 0.22,
    offset_min_A: float = 8.0,
    offset_max_A: float = 14.0,
) -> list[dict]:
    try:
        atoms = pdbqt_atoms(receptor_pdbqt)
        if not atoms:
            return []
        coords = np.array([(a[0], a[1], a[2]) for a in atoms], dtype=np.float32)
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)
        extents = maxs - mins
        core, _, sp = blind_box(receptor_pdbqt, cap=float(box_side_A), spacing=float(spacing))
        core = np.asarray(core, dtype=np.float32)
    except Exception as exc:
        print(f"[centers/axis] skipped due to error: {exc}")
        return []

    preserved_centers = [_site_center(site) for site in preserved_sites]
    candidates: list[tuple[float, np.ndarray, int, int]] = []
    for axis in range(3):
        offset = float(np.clip(float(offset_fraction) * float(extents[axis]), offset_min_A, offset_max_A))
        for sign in (-1, 1):
            center = core.copy()
            center[axis] += float(sign) * offset
            center = np.maximum(mins, np.minimum(maxs, center))
            if preserved_centers:
                coverage_gap = min(float(np.linalg.norm(center - old)) for old in preserved_centers)
            else:
                coverage_gap = float("inf")
            candidates.append((coverage_gap, center, axis, sign))

    if not candidates:
        return []

    n = _npts_for_box_side(float(box_side_A), float(spacing))
    reserve_sites = []
    for coverage_gap, center, axis, sign in sorted(candidates, key=lambda item: item[0], reverse=True):
        reserve_sites.append(
            dict(
                site_id=f"S_axis_{int(axis)}_{int(sign)}",
                cx=float(center[0]),
                cy=float(center[1]),
                cz=float(center[2]),
                nx=n,
                ny=n,
                nz=n,
                spacing=float(sp),
                r_peak=max(1.0, float(box_side_A) * 0.15),
                F=0.0,
                family="axis_reserve",
                reserve_axis=int(axis),
                reserve_sign=int(sign),
                reserve_coverage_gap_A=float(coverage_gap),
            )
        )
    return reserve_sites


def _replace_tail_with_reserve_sites(
    sites: list[dict],
    *,
    axis_site: dict | None,
    core_site: dict | None,
    min_recenter_A: float,
) -> list[dict]:
    if not sites:
        return sites
    updated = list(sites)

    def _far_enough(site: dict, others: list[dict]) -> bool:
        center = _site_center(site)
        return all(np.linalg.norm(center - _site_center(old)) >= float(min_recenter_A) for old in others)

    changed = False
    if axis_site is not None and len(updated) >= 2 and _far_enough(axis_site, updated[:-2]):
        updated[-2] = _clone_site(axis_site, portfolio_role="axis_reserve")
        changed = True
    if core_site is not None and _far_enough(core_site, updated[:-1]):
        updated[-1] = _clone_site(core_site, portfolio_role="core_reserve")
        changed = True

    return _relabel_sites(updated) if changed else sites


def _portfolio_min_sep_A(min_sep_A: float) -> float:
    return float(np.clip(0.55 * float(min_sep_A), 4.0, 8.0))


def _candidate_harvest_min_sep_A(min_sep_A: float, target_count: int) -> float:
    """
    Use a looser NMS while collecting raw candidates.

    The box-scaled clamp is meant to keep final docking boxes from collapsing
    onto the same site. Applying that full separation during candidate harvest
    can starve multi-site searches before the portfolio assembler gets a chance
    to choose diverse hypotheses.
    """
    if int(target_count) <= 3:
        return float(min_sep_A)
    return min(float(min_sep_A), _portfolio_min_sep_A(min_sep_A))


def _surface_contact_fraction(
    edt_grid: np.ndarray,
    ijk: tuple[int, int, int],
    spacing: float,
    *,
    shell_min_A: float,
    contact_shell_A: float = CONTACT_SHELL_A,
) -> float:
    radius_vox = max(2, int(round(float(contact_shell_A) / float(spacing))))
    slices = []
    for axis, center in enumerate(ijk):
        lo = max(0, center - radius_vox)
        hi = min(int(edt_grid.shape[axis]), center + radius_vox + 1)
        slices.append(slice(lo, hi))
    block = edt_grid[tuple(slices)]
    if block.size == 0:
        return 0.0
    contact_mask = (block >= shell_min_A) & (block <= float(contact_shell_A))
    return float(contact_mask.mean())


def _refine_surface_hotspot_ijk(
    ijk: np.ndarray,
    *,
    score: np.ndarray,
    edt_grid: np.ndarray | None,
    pocket_mask: np.ndarray,
    border_mask: np.ndarray,
    spacing: float,
    raw_score: float,
    radius_A: float = 8.0,
    score_keep_fraction: float = 0.55,
    edt_core_percentile: float = 85.0,
    return_diagnostics: bool = False,
) -> tuple[int, int, int] | tuple[tuple[int, int, int], dict]:
    """
    Move a surface-score peak toward the geometric center of the same pocket.

    AutoGrid energy maxima sit near pocket walls. Docking boxes, however, should
    be centered nearer the ligand-sized free volume between those walls. Keep
    the refinement local and require retained map support so it does not drift
    into unrelated solvent space.
    """
    seed_ijk = (int(ijk[0]), int(ijk[1]), int(ijk[2]))

    def _return(result: tuple[int, int, int], **diagnostics):
        if return_diagnostics:
            return result, diagnostics
        return result

    if edt_grid is None:
        return _return(
            seed_ijk,
            refinement_component_count=0,
            refinement_component_size=0,
        )

    sp = float(spacing)
    radius_vox = max(1, int(round(float(radius_A) / sp)))
    slices = []
    for axis, center in enumerate(ijk):
        lo = max(0, int(center) - radius_vox)
        hi = min(int(score.shape[axis]), int(center) + radius_vox + 1)
        slices.append(slice(lo, hi))
    sl = tuple(slices)

    local_score = score[sl]
    local_edt = edt_grid[sl]
    local_mask = pocket_mask[sl] & border_mask[sl]

    grids = np.ogrid[
        slices[0].start:slices[0].stop,
        slices[1].start:slices[1].stop,
        slices[2].start:slices[2].stop,
    ]
    local_dist_A = np.sqrt(
        (grids[0] - int(ijk[0])) ** 2
        + (grids[1] - int(ijk[1])) ** 2
        + (grids[2] - int(ijk[2])) ** 2
    ) * sp

    score_floor = max(0.01, float(raw_score) * float(score_keep_fraction))
    local_mask &= (local_dist_A <= float(radius_A)) & (local_score >= score_floor)
    if not np.any(local_mask):
        return _return(
            seed_ijk,
            refinement_component_count=0,
            refinement_component_size=0,
        )

    component_labels, component_count = label(
        local_mask,
        structure=generate_binary_structure(3, 2),
    )
    seed_local = tuple(
        int(ijk[axis]) - int(slices[axis].start)
        for axis in range(3)
    )
    seed_component = int(component_labels[seed_local])
    if seed_component <= 0:
        valid_coords = np.argwhere(local_mask)
        if valid_coords.size:
            seed_array = np.asarray(seed_local, dtype=np.int32)
            nearest_index = int(np.argmin(np.sum((valid_coords - seed_array) ** 2, axis=1)))
            nearest_coord = tuple(int(value) for value in valid_coords[nearest_index])
            seed_component = int(component_labels[nearest_coord])
    if seed_component > 0:
        local_mask &= component_labels == seed_component
    if not np.any(local_mask):
        return _return(
            seed_ijk,
            refinement_component_count=int(component_count),
            refinement_component_size=0,
        )
    component_size = int(np.count_nonzero(local_mask))

    core_floor = float(np.percentile(local_edt[local_mask], float(edt_core_percentile)))
    core_mask = local_mask & (local_edt >= core_floor)
    if not np.any(core_mask):
        core_mask = local_mask

    weights = local_edt[core_mask] * (local_score[core_mask] + 1e-3)
    if weights.size == 0 or not np.isfinite(weights).all() or float(weights.sum()) <= 0.0:
        return _return(
            seed_ijk,
            refinement_component_count=int(component_count),
            refinement_component_size=component_size,
        )

    local_coords = np.where(core_mask)
    centroid = np.array(
        [np.average(local_coords[axis], weights=weights) for axis in range(3)],
        dtype=np.float32,
    )
    dist_to_centroid = (
        (local_coords[0] - centroid[0]) ** 2
        + (local_coords[1] - centroid[1]) ** 2
        + (local_coords[2] - centroid[2]) ** 2
    )
    closest = int(np.argmin(dist_to_centroid))
    refined_ijk = tuple(
        int(local_coords[axis][closest] + slices[axis].start)
        for axis in range(3)
    )
    return _return(
        refined_ijk,
        refinement_component_count=int(component_count),
        refinement_component_size=component_size,
    )


def _build_hybrid_sites(
    internal_sites: list[dict],
    surface_sites: list[dict],
    *,
    min_sep_A: float,
    max_sites: int,
) -> list[dict]:
    if not internal_sites and not surface_sites:
        return []

    internal_ranked = _normalize_ranked_sites(internal_sites, family="internal")
    surface_ranked = _normalize_ranked_sites(surface_sites, family="surface")
    candidate_pool = internal_ranked + surface_ranked

    hybrid_candidates = []
    decay = max(2.0, float(MAX_CENTER_DIST_A))

    for candidate in candidate_pool:
        ctr = _site_center(candidate)

        def _support(pool: list[dict]) -> tuple[float, dict | None]:
            best_score = 0.0
            best_site = None
            for site in pool:
                dist = float(np.linalg.norm(ctr - _site_center(site)))
                score = float(site["selection_score"]) * float(np.exp(-dist / decay))
                if score > best_score:
                    best_score = score
                    best_site = site
            return best_score, best_site

        internal_support, internal_anchor = _support(internal_ranked)
        surface_support, surface_anchor = _support(surface_ranked)
        if internal_support <= 0.0 and surface_support <= 0.0:
            continue

        agreement = float(np.sqrt(max(0.0, internal_support * surface_support)))
        coverage = 0.5 * (internal_support + surface_support)
        hybrid_score = 0.65 * agreement + 0.35 * coverage

        anchors = [site for site in (internal_anchor, surface_anchor) if site is not None]
        if anchors:
            weights = np.array([max(1e-3, float(site["selection_score"])) for site in anchors], dtype=np.float32)
            centers = np.array([_site_center(site) for site in anchors], dtype=np.float32)
            ctr = np.average(centers, axis=0, weights=weights)
        internal_anchor_center = _site_center(internal_anchor) if internal_anchor is not None else None
        surface_anchor_center = _site_center(surface_anchor) if surface_anchor is not None else None
        if internal_anchor_center is not None and surface_anchor_center is not None:
            anchor_separation_A = float(np.linalg.norm(internal_anchor_center - surface_anchor_center))
        else:
            anchor_separation_A = ""

        template = max(anchors + [candidate], key=lambda site: _site_score_value(site))
        center_closeness = max(
            [float(site.get("center_closeness", 0.0)) for site in anchors + [candidate]],
            default=0.0,
        )
        hybrid_candidates.append(
            _clone_site(
                template,
                cx=float(ctr[0]),
                cy=float(ctr[1]),
                cz=float(ctr[2]),
                family="hybrid",
                selection_score=float(hybrid_score),
                F=float(hybrid_score),
                center_closeness=float(center_closeness),
                internal_anchor_x="" if internal_anchor_center is None else float(internal_anchor_center[0]),
                internal_anchor_y="" if internal_anchor_center is None else float(internal_anchor_center[1]),
                internal_anchor_z="" if internal_anchor_center is None else float(internal_anchor_center[2]),
                surface_anchor_x="" if surface_anchor_center is None else float(surface_anchor_center[0]),
                surface_anchor_y="" if surface_anchor_center is None else float(surface_anchor_center[1]),
                surface_anchor_z="" if surface_anchor_center is None else float(surface_anchor_center[2]),
                anchor_separation_A=anchor_separation_A,
                hybrid_center_x=float(ctr[0]),
                hybrid_center_y=float(ctr[1]),
                hybrid_center_z=float(ctr[2]),
            )
        )

    hybrid_candidates.sort(key=lambda site: _site_score_value(site), reverse=True)
    selected = []
    for site in hybrid_candidates:
        if _take_best_distinct([site], selected, min_sep_A=min_sep_A) is not None:
            selected.append(site)
        if len(selected) == int(max_sites):
            break
    return selected


def _assemble_receptor_search_sites(
    internal_sites: list[dict],
    surface_sites: list[dict],
    hybrid_sites: list[dict],
    *,
    min_sep_A: float,
    target_count: int = 3,
) -> list[dict]:
    selected: list[dict] = []
    selection_min_sep_A = float(min_sep_A) if int(target_count) <= 3 else _portfolio_min_sep_A(min_sep_A)
    surface_focus_sites = _rank_surface_focus_sites(surface_sites)
    hybrid_low_score_focus_sites = _rank_hybrid_low_score_focus_sites(hybrid_sites)
    if int(target_count) <= 3:
        portfolio_pools = (
            ("hybrid_primary", hybrid_sites),
            ("internal_primary", internal_sites),
            ("surface_primary", surface_sites),
        )
    else:
        portfolio_pools = (
            ("hybrid_primary", hybrid_sites),
            ("internal_primary", internal_sites),
            ("surface_primary", surface_sites),
            ("surface_focus", surface_focus_sites),
            ("hybrid_low_score_focus", hybrid_low_score_focus_sites),
        )
    for role, pool in portfolio_pools:
        site = _take_best_distinct(pool, selected, min_sep_A=selection_min_sep_A)
        if site is not None:
            selected.append(_clone_site(site, portfolio_role=role))
        if len(selected) == target_count:
            break

    if len(selected) < target_count:
        reserves = _rank_union_sites(internal_sites, surface_sites, hybrid_sites)
        for site in reserves:
            if _take_best_distinct([site], selected, min_sep_A=selection_min_sep_A) is not None:
                selected.append(_clone_site(site, portfolio_role="reserve"))
            if len(selected) == target_count:
                break

    return _relabel_sites(selected[:target_count])


def _assemble_exhaustive_search_sites(
    internal_sites: list[dict],
    surface_sites: list[dict],
    hybrid_sites: list[dict],
    *,
    min_sep_A: float,
    target_count: int,
) -> list[dict]:
    selected: list[dict] = []
    for site in _rank_union_sites(hybrid_sites, internal_sites, surface_sites):
        if _take_best_distinct([site], selected, min_sep_A=min_sep_A) is not None:
            selected.append(site)
        if len(selected) == int(target_count):
            break
    return _relabel_sites(selected)


def _expected_site_count(mode: str, n_sites: int) -> int:
    policy = (mode or "").lower()
    if policy == "receptor_search":
        return max(3, int(n_sites))
    return int(n_sites)


def _select_internal_search_r_min(
    requested_r_min_A: float | None,
    peak_radii_A: list[float],
    *,
    fallback_floor_A: float = 1.5,
    percentile: float = 90.0,
) -> float | None:
    """
    Pick the internal-family search threshold from the observed enclosed-cavity
    peak distribution.

    `R_MIN_CAVITY_A` describes the shallowest cavity centres we consider
    *credible* for a receptor profile. Reusing it unchanged as the internal
    candidate-generation gate is too strict: if the enclosed-cavity geometry is
    slightly shallower than the profiled threshold, the internal family drops to
    zero sites and contributes no search hypothesis at all.

    For internal candidate generation we therefore use a softer, recall-oriented
    gate:
      - never stricter than the requested/profiled threshold
      - if that threshold would eliminate the whole internal peak distribution,
        relax to the upper tail of the observed enclosed peaks
      - never relax below a small physical floor, to avoid flooding with
        single-voxel noise
    """
    if requested_r_min_A is None:
        return None
    if not peak_radii_A:
        return float(requested_r_min_A)

    observed_tail_A = float(np.percentile(np.asarray(peak_radii_A, dtype=np.float32), percentile))
    relaxed_r_min_A = max(float(fallback_floor_A), observed_tail_A)
    return round(float(min(float(requested_r_min_A), relaxed_r_min_A)), 2)


def _edt_pocket_centers(
    receptor_pdbqt, origin, spacing, shape,
    *,
    min_sep_A, max_sites, box_side_A,
    # pocket_min_A: minimum EDT to accept as a pocket centre.
    #   2.5 Å ≈ radius of a water molecule (1.4 Å) + one C–C bond (1.5 Å).
    #   Cavities smaller than this cannot accommodate a drug-like fragment.
    pocket_min_A=2.5,
    # pocket_max_A: maximum EDT included in the pocket zone.
    #   Set to 15 Å (was 10 Å) to include deep transmembrane channels that are
    #   common in GPCR, transporter (SERT/DAT/NET), and ion channel targets.
    #   The SERT vestibule, for example, extends ~13 Å from the nearest
    #   protein atom to the channel axis.  Capping at 10 Å excluded these
    #   pharmacologically critical sites.
    pocket_max_A=15.0,
):
    """
    Pure EDT-based pocket detection fallback for degenerate maps.

    Strategy (in priority order):
    1. **Internal cavities**: flood-fill the exterior, keep only enclosed
       free-space voxels (true pockets/channels), pick the deepest local
       EDT maxima.  This is the most reliable for buried binding sites.
    2. **Surface buriedness**: if no internal cavities pass the minimum
       depth threshold, score each pocket-zone voxel by *enclosure* —
       what fraction of 26 neighbour directions are blocked by protein
       within a probe radius.  This finds clefts and grooves on the
       surface that are partially buried and likely druggable.
    """
    atoms = pdbqt_atoms(receptor_pdbqt)
    occ, _, (ox, oy, oz), sp_actual = make_grid_from_map(
        atoms, origin, spacing, shape, inflate_A=0.0
    )
    sp = float(spacing)
    edt = distance_transform_edt(~occ) * sp

    # ── Strategy 1: internal cavities ──────────────────────────────────────
    try:
        dist_vox, cc, ncc = internal_cavities(occ)
        dist_A = dist_vox * sp

        candidates = []
        for lab in range(1, ncc + 1):
            mask = (cc == lab)
            if not mask.any():
                continue
            i, j, k = np.unravel_index(np.argmax(dist_A * mask), dist_A.shape)
            r_peak = float(dist_A[i, j, k])
            if r_peak < pocket_min_A:
                continue
            cx = ox + i * sp
            cy = oy + j * sp
            cz = oz + k * sp
            candidates.append((np.array([cx, cy, cz]), r_peak))

        if candidates:
            # Sort by depth (deepest first), then NMS
            candidates.sort(key=lambda x: -x[1])
            kept_xyz, sites = [], []
            for ctr, r_peak in candidates:
                if all(np.linalg.norm(ctr - q) >= min_sep_A for q in kept_xyz):
                    kept_xyz.append(ctr)
                    n = _npts_for_box_side(box_side_A, sp)
                    sites.append(dict(
                        site_id=f"S{len(kept_xyz)}",
                        cx=float(ctr[0]), cy=float(ctr[1]), cz=float(ctr[2]),
                        nx=n, ny=n, nz=n,
                        spacing=sp,
                        r_peak=r_peak,
                        F=r_peak,   # score by cavity depth
                    ))
                    if len(sites) == int(max_sites):
                        break
            if sites:
                print(f"[EDT-fallback/internal] found {len(sites)} enclosed pockets from {ncc} components")
                return sites
    except Exception as exc:
        print(f"[EDT-fallback/internal] failed: {exc}")

    # ── Strategy 2: surface buriedness ────────────────────────────────────
    # For each pocket-zone voxel, cast rays in 26 directions and count
    # how many hit protein within a probe radius.  High enclosure =
    # buried cleft = likely binding site.
    print("[EDT-fallback] no internal cavities found, using surface buriedness")
    pocket = (edt >= pocket_min_A) & (edt <= pocket_max_A)

    # Compute buriedness: fraction of 26 directions blocked by protein
    # within a shell of ~12 Å.  Uses binary dilation as a fast proxy.
    from scipy.ndimage import binary_dilation, generate_binary_structure

    # Probe shell: 12 Å / spacing voxels
    probe_vox = max(3, int(round(12.0 / sp)))
    struct = generate_binary_structure(3, 2)  # 26-connected

    # Dilate the occupancy and count overlap with pocket zone
    # "buriedness" ≈ how much protein is nearby
    occ_dilated = binary_dilation(occ, iterations=probe_vox, structure=struct)
    # Voxels that are in the pocket zone AND surrounded by protein
    # are the ones where dilation of protein reaches them from many sides

    # Better: score by EDT value × inverse distance to protein center-of-mass
    # (prefer pockets closer to the protein core)
    atom_coords = np.array([(a[0], a[1], a[2]) for a in atoms], dtype=np.float32)
    com = atom_coords.mean(axis=0)

    # Build distance-from-COM grid
    gi = ox + np.arange(shape[0]) * sp
    gj = oy + np.arange(shape[1]) * sp
    gk = oz + np.arange(shape[2]) * sp
    dist_from_com = np.sqrt(
        (gi[:, None, None] - com[0])**2 +
        (gj[None, :, None] - com[1])**2 +
        (gk[None, None, :] - com[2])**2
    ).astype(np.float32)

    # Score: EDT depth × closeness to protein core
    # Deeper pocket + closer to core = higher score
    max_com_dist = float(dist_from_com[pocket].max()) + 1e-6
    closeness = 1.0 - (dist_from_com / max_com_dist)
    score = edt * closeness
    score[~pocket] = 0.0
    score = gaussian_filter(score, sigma=2.0)
    score[~pocket] = 0.0

    # Border exclusion
    BORDER = 5
    border_mask = np.zeros(shape, dtype=bool)
    border_mask[BORDER:-BORDER, BORDER:-BORDER, BORDER:-BORDER] = True

    local_max = maximum_filter(score, size=7, mode="constant", cval=0.0)
    peak_mask = (score == local_max) & (score > 0.5) & border_mask & pocket

    peaks = np.argwhere(peak_mask)
    scores_arr = score[peak_mask].astype(np.float32)
    order = np.argsort(-scores_arr, kind="stable")
    peaks = peaks[order]
    scores_arr = scores_arr[order]

    kept_xyz, sites = [], []
    for ijk, sc in zip(peaks, scores_arr):
        i, j, k = ijk
        wp = np.array([ox + i * sp, oy + j * sp, oz + k * sp], dtype=np.float32)
        if all(np.linalg.norm(wp - q) >= min_sep_A for q in kept_xyz):
            kept_xyz.append(wp)
            n = _npts_for_box_side(box_side_A, sp)
            r_peak = float(edt[i, j, k])
            sites.append(dict(
                site_id=f"S{len(kept_xyz)}",
                cx=float(wp[0]),
                cy=float(wp[1]),
                cz=float(wp[2]),
                nx=n, ny=n, nz=n,
                spacing=sp,
                r_peak=r_peak,
                F=float(sc),
            ))
            if len(sites) == int(max_sites):
                break

    print(f"[EDT-fallback/surface] found {len(sites)} pocket centres from {len(peaks)} surface peaks")
    return sites


def _edt_surface_sites(
    receptor_pdbqt,
    origin,
    spacing,
    shape,
    *,
    min_sep_A,
    max_sites,
    box_side_A,
    pocket_min_A,
    pocket_max_A,
):
    """Surface-cleft fallback based purely on EDT buriedness near the protein wall."""
    atoms = pdbqt_atoms(receptor_pdbqt)
    occ, _, (ox, oy, oz), _ = make_grid_from_map(
        atoms, origin, spacing, shape, inflate_A=0.0
    )
    sp = float(spacing)
    edt = distance_transform_edt(~occ) * sp
    pocket = (edt >= pocket_min_A) & (edt <= pocket_max_A)

    atom_coords = np.array([(a[0], a[1], a[2]) for a in atoms], dtype=np.float32)
    com = atom_coords.mean(axis=0)

    gi = ox + np.arange(shape[0]) * sp
    gj = oy + np.arange(shape[1]) * sp
    gk = oz + np.arange(shape[2]) * sp
    dist_from_com = np.sqrt(
        (gi[:, None, None] - com[0]) ** 2 +
        (gj[None, :, None] - com[1]) ** 2 +
        (gk[None, None, :] - com[2]) ** 2
    ).astype(np.float32)

    if not np.any(pocket):
        return []

    max_com_dist = float(dist_from_com[pocket].max()) + 1e-6
    closeness = 1.0 - (dist_from_com / max_com_dist)
    score = gaussian_filter(edt * closeness, sigma=2.0)
    score[~pocket] = 0.0

    border_mask = np.zeros(shape, dtype=bool)
    BORDER = 5
    border_mask[BORDER:-BORDER, BORDER:-BORDER, BORDER:-BORDER] = True

    local_max = maximum_filter(score, size=7, mode="constant", cval=0.0)
    peak_mask = (score == local_max) & (score > 0.5) & border_mask & pocket
    peaks = np.argwhere(peak_mask)
    scores_arr = score[peak_mask].astype(np.float32)
    order = np.argsort(-scores_arr, kind="stable")
    peaks = peaks[order]
    scores_arr = scores_arr[order]

    kept_xyz, sites = [], []
    for ijk, sc in zip(peaks, scores_arr):
        i, j, k = ijk
        wp = np.array([ox + i * sp, oy + j * sp, oz + k * sp], dtype=np.float32)
        if all(np.linalg.norm(wp - q) >= min_sep_A for q in kept_xyz):
            contact_frac = _surface_contact_fraction(
                edt,
                (int(i), int(j), int(k)),
                sp,
                shell_min_A=pocket_min_A,
            )
            if contact_frac < float(MIN_SURFACE_FRAC):
                continue
            kept_xyz.append(wp)
            n = _npts_for_box_side(box_side_A, sp)
            sites.append(
                dict(
                    site_id=f"S{len(kept_xyz)}",
                    cx=float(wp[0]),
                    cy=float(wp[1]),
                    cz=float(wp[2]),
                    nx=n,
                    ny=n,
                    nz=n,
                    spacing=sp,
                    r_peak=float(edt[i, j, k]),
                    F=float(sc),
                    family="surface",
                    contact_frac=float(contact_frac),
                    center_closeness=float(np.clip(closeness[i, j, k], 0.0, 1.0)),
                )
            )
            if len(sites) == int(max_sites):
                break
    print(f"[EDT-surface] found {len(sites)} pocket centres from {len(peaks)} surface peaks")
    return sites


def detect_maps_hotspots(
    C, E, D, origin, spacing,
    tau_rel=0.52,
    min_sep_A=HOTSPOT_NMS_MINSEP_A,
    max_sites=AUTOSITES,
    box_side_A=HOTSPOT_BOX_ANGLE,
    r_min_A=R_MIN_CAVITY_A,
    adaptive_r_min_params: dict | None = None,
    pocket_max_A: float | None = None,
    receptor_pdbqt=None,        # if supplied, EDT r_peak is computed from true geometry
    min_peak_candidates: int | None = None,
    return_score_context: bool = False,
):
    """
    Binding-site finder using pocket-averaged interaction energy from AutoGrid maps.

    Scientific rationale
    --------------------
    AutoGrid computes the Lennard-Jones (C-map) and Coulombic (E-map) interaction
    energy that a probe atom would experience at each grid point.  These energies
    are most negative (most favorable) at the van-der-Waals contact distance
    (~3.4–3.8 Å) from protein heavy atoms — i.e., at the **protein surface**,
    not at the geometric center of a binding pocket.

    Naïvely scoring individual voxels therefore picks surface atoms rather than
    pocket interiors.  A real binding pocket is characterised by *many* favorable
    contacts surrounding a ligand-sized cavity from multiple sides.  We therefore:

    1. **Clip/scale** C and E to favorable contact support, discarding steric
       clashes (large positives from VDW overlap) that dominate the raw range.
    2. **Add** D-map desolvation support so hydrophobic/dehydrating crypts are
       not invisible to the regional score.
    3. **Convolve** the combined field with an isotropic Gaussian over a
       ligand-sized length scale.  This is rotationally invariant, unlike a
       cubic box filter.
    4. **Weight/mask** by the receptor-specific EDT pocket shell, using a
       continuous enclosure/depth weight before excluding bulk solvent and
       steric-clash voxels.

    This regional scoring replaces the old per-voxel sigmoid scoring, which was
    too sensitive to isolated surface contacts and too weak at distinguishing
    enclosed pockets from exposed grooves.

    If the maps are degenerate (E all-zero, C mostly-zero with extreme outliers),
    falls back to pure EDT-based pocket detection from the receptor geometry.

    Parameters
    ----------
    receptor_pdbqt : str | Path | None
        If provided, the receptor occupancy grid is built (aligned to the C/E/D
        maps) and distance_transform_edt is run once.  Each accepted peak is then
        assigned the true geometric r_peak (EDT value = distance to nearest
        protein atom surface in Å).
    """
    C = np.asarray(C, dtype=np.float32)
    E = np.asarray(E, dtype=np.float32)
    D = np.asarray(D, dtype=np.float32)

    sp = float(spacing)

    def _result(found_sites: list[dict], context: dict | None = None):
        if return_score_context:
            return found_sites, context
        return found_sites

    # ── Build EDT for pocket masking ─────────────────────────────────────────
    edt_grid = None
    occ_grid = None
    protein_center = None
    protein_radius_A = None
    if receptor_pdbqt is not None:
        try:
            atoms = pdbqt_atoms(receptor_pdbqt)
            atom_coords = np.array([(a[0], a[1], a[2]) for a in atoms], dtype=np.float32)
            if atom_coords.size:
                protein_center = atom_coords.mean(axis=0)
                protein_radius_A = float(
                    np.max(np.linalg.norm(atom_coords - protein_center, axis=1))
                ) + 1e-6
            occ, _, _, _sp = make_grid_from_map(
                atoms, origin, spacing, C.shape, inflate_A=0.0
            )
            occ_grid = occ
            edt_grid = distance_transform_edt(~occ_grid) * sp
        except Exception as exc:
            print(f"[maps/EDT] failed: {exc}")

    if r_min_A is None and occ_grid is not None:
        r_min_A = adaptive_r_min_cavity(occ_grid, sp, **_adaptive_r_min_params(adaptive_r_min_params))

    # ── Degenerate-map guard ─────────────────────────────────────────────────
    if receptor_pdbqt is not None and _maps_are_degenerate(C, E, D):
        pocket_min_A, shell_max_A = _maps_pocket_shell_bounds(r_min_A=r_min_A, pocket_max_A=pocket_max_A)
        fallback_sites = _edt_surface_sites(
            receptor_pdbqt,
            origin,
            spacing,
            C.shape,
            min_sep_A=min_sep_A,
            max_sites=max_sites,
            box_side_A=box_side_A,
            pocket_min_A=pocket_min_A,
            pocket_max_A=shell_max_A,
        )
        return _result(
            fallback_sites,
            {
                "origin": tuple(float(value) for value in origin),
                "spacing": sp,
                "score_field": None,
                "score_max": 0.0,
                "edt_grid": edt_grid,
                "pocket_min_A": float(pocket_min_A),
                "pocket_max_A": float(shell_max_A),
            },
        )

    # ── Pocket-averaged favorability scoring ─────────────────────────────────
    #
    # WHY NOT per-voxel scoring?
    #   AutoGrid's Lennard-Jones potential has its energy minimum at the VDW
    #   contact distance (~3.4 Å for carbon).  The most negative C-map values
    #   are therefore RIGHT AT the protein surface, not at pocket centres.
    #   Per-voxel scoring (the old robust_z → sigmoid approach) picks surface
    #   atoms as "hotspots" rather than the enclosed cavities that ligands
    #   actually bind in.  It also tends to compress scores into a narrow,
    #   poorly discriminating range across unrelated surface grooves.
    #
    # WHY isotropic regional convolution?
    #   A druggable binding pocket (Kd < 10 µM) is defined by *complementarity*:
    #   the ligand must form simultaneous favorable contacts with protein residues
    #   on multiple sides.  Convolving the interaction support over a ligand-sized
    #   length scale (~8 Å support, typical for MW 300–500 Da drug-like molecules)
    #   effectively measures this complementarity without introducing grid-axis
    #   bias.  Isolated favorable voxels at exposed surfaces wash out, while
    #   enclosed pockets with walls on 3+ sides accumulate high scores.
    #
    # WHY clip to negative values only?
    #   Positive C/E values represent steric clashes (VDW repulsion) or
    #   unfavorable electrostatics.  The raw C-map ranges from -0.94 to
    #   +302,000 kcal/mol.  Without clipping, the enormous positive outliers
    #   dominate any normalization (robust_z, sigmoid) and collapse the useful
    #   signal into numeric noise.  Clipping to [-5, 0] for C and [-50, 0] for
    #   E keeps only the physically meaningful favorable region of the energy
    #   landscape.  The asymmetric clip limits reflect that electrostatic
    #   interactions (H-bonds, salt bridges) span a wider energy range than
    #   van-der-Waals contacts.
    #
    # Step 1: clip to favorable-only interaction support
    # C_fav: van-der-Waals favorable contacts.  Typical VDW well depths for
    #   drug-like atoms are -0.1 to -0.5 kcal/mol; -5.0 accommodates H-bonds.
    # E_fav: electrostatic favorable contacts.  Salt bridges can reach -40
    #   kcal/mol in vacuum; -50.0 is a generous upper bound.
    C_fav = np.clip(-np.clip(C, -5.0, 0.0) / 5.0, 0.0, 1.0).astype(np.float32)
    E_fav = np.clip(-np.clip(E, -50.0, 0.0) / 50.0, 0.0, 1.0).astype(np.float32)
    D_fav = _positive_signal_unit_interval(np.abs(D))

    weights = _normalize_positive_weights(
        {"C": MAPS_C_WEIGHT, "E": MAPS_E_WEIGHT, "D": MAPS_D_WEIGHT}
    )
    fav = (
        weights["C"] * C_fav
        + weights["E"] * E_fav
        + weights["D"] * D_fav
    ).astype(np.float32)

    # Step 2: restrict to a receptor-specific "pocket shell" — voxels near the
    # protein surface, but not inside steric-clash volume.  The lower EDT bound
    # is derived from the profiled cavity threshold instead of being a fixed
    # 3 Å cutoff, because shallow surface pockets often place ligand atoms at
    # ~1–2.5 Å from the nearest receptor surface.
    #   < shell_min: protein interior or steric clash zone
    #   shell_min–shell_max: pocket-contact zone where ligands make wall contacts
    #   > shell_max: bulk solvent — no enclosure, no binding complementarity
    pocket_min_A, shell_max_A = _maps_pocket_shell_bounds(r_min_A=r_min_A, pocket_max_A=pocket_max_A)
    if edt_grid is not None:
        pocket = (edt_grid >= pocket_min_A) & (edt_grid <= shell_max_A)
    else:
        # Fallback when EDT is unavailable: use desolvation map as pocket proxy
        dmin = float(np.nanmin(D))
        d95 = float(np.nanpercentile(D, 95))
        Dn = np.clip((D - dmin) / (d95 - dmin + 1e-6), 0.0, 1.0)
        pocket = Dn > 0.35

    # Step 3: isotropic ligand-scale convolution.  The Gaussian is separable and
    # rotationally symmetric, so the score does not depend on how the receptor
    # sits relative to the grid axes.  The truncate value keeps the effective
    # support near MAPS_CONVOLUTION_RADIUS_A.
    sigma_vox = max(0.5, float(MAPS_GAUSSIAN_SIGMA_A) / sp)
    truncate = max(
        2.0,
        float(MAPS_CONVOLUTION_RADIUS_A) / max(float(MAPS_GAUSSIAN_SIGMA_A), 1e-6),
    )
    fav_smooth = gaussian_filter(
        fav,
        sigma=sigma_vox,
        mode="constant",
        cval=0.0,
        truncate=truncate,
    ).astype(np.float32)

    # Step 4: continuously reward enclosure/depth inside the allowed pocket
    # shell.  This keeps the hard shell as a safety mask while avoiding a binary
    # all-or-none treatment of shallow grooves versus deep vestibules.
    score = fav_smooth
    if edt_grid is not None:
        score = score * _continuous_edt_weight(
            edt_grid,
            pocket_min_A=pocket_min_A,
            pocket_max_A=shell_max_A,
        )

    # Step 5: exclude outermost 5 voxels on every face.  AutoGrid's finite-
    # difference solver uses boundary conditions that produce artefactual
    # energy values at grid edges (often maximal, causing false peaks).
    BORDER = 5
    border_mask = np.zeros(C.shape, dtype=bool)
    border_mask[BORDER:-BORDER, BORDER:-BORDER, BORDER:-BORDER] = True

    score[~(pocket & border_mask)] = 0.0
    score_context = {
        "origin": tuple(float(value) for value in origin),
        "spacing": sp,
        "score_field": score,
        "score_max": float(np.nanmax(score)),
        "edt_grid": edt_grid,
        "pocket_min_A": float(pocket_min_A),
        "pocket_max_A": float(shell_max_A),
    }

    # Local maxima. Start with the requested relative threshold, but if that
    # yields too few raw peaks, relax it before NMS. This keeps the default
    # strict for strong pockets while preserving recall for weaker, buried
    # pockets whose absolute map score is lower than exposed surface patches.
    score_max = float(np.nanmax(score))
    local_max = maximum_filter(score, size=7, mode="constant", cval=0.0)
    target_peak_count = max(
        1,
        min(int(max_sites), int(min_peak_candidates) if min_peak_candidates else int(max_sites)),
    )
    threshold_schedule = []
    for threshold in (
        float(tau_rel),
        min(float(tau_rel) - 0.12, 0.40),
        min(float(tau_rel) - 0.22, 0.30),
        min(float(tau_rel) - 0.32, 0.22),
    ):
        threshold = max(0.15, float(threshold))
        if threshold_schedule and threshold >= threshold_schedule[-1] - 1e-6:
            continue
        if all(abs(threshold - seen) > 1e-6 for seen in threshold_schedule):
            threshold_schedule.append(threshold)

    peak_mask = np.zeros(score.shape, dtype=bool)
    peak_floor = 0.0
    selected_tau_rel = threshold_schedule[0]
    for threshold in threshold_schedule:
        rel_floor = float(np.clip(threshold, 0.0, 1.0)) * score_max
        peak_floor = max(0.01, rel_floor)
        peak_mask = (score == local_max) & (score >= peak_floor)
        selected_tau_rel = threshold
        if int(np.count_nonzero(peak_mask)) >= target_peak_count:
            break

    if selected_tau_rel + 1e-6 < float(tau_rel):
        print(
            f"[maps/hotspots] relaxed tau_rel from {float(tau_rel):.2f} "
            f"to {selected_tau_rel:.2f} to collect candidate peaks"
        )

    if not np.any(peak_mask):
        print("[maps/hotspots] no peaks found in favorability landscape")
        if receptor_pdbqt is not None:
            fallback_sites = _edt_surface_sites(
                receptor_pdbqt,
                origin,
                spacing,
                C.shape,
                min_sep_A=min_sep_A,
                max_sites=max_sites,
                box_side_A=box_side_A,
                pocket_min_A=pocket_min_A,
                pocket_max_A=shell_max_A,
            )
            return _result(fallback_sites, score_context)
        return _result([], score_context)

    peaks = np.argwhere(peak_mask)
    scores = score[peak_mask].astype(np.float32)
    order = np.argsort(-scores, kind="stable")
    peaks = peaks[order]
    scores = scores[order]

    print(f"[maps/hotspots] {len(peaks)} peaks, score range {scores[-1]:.4f}–{scores[0]:.4f}")

    # NMS in Å
    ox, oy, oz = origin

    def vox2world(ijk):
        i, j, k = ijk
        return np.array([ox + i * sp, oy + j * sp, oz + k * sp], dtype=np.float32)

    def center_closeness(wp):
        if protein_center is None or protein_radius_A is None:
            return 0.0
        dist = float(np.linalg.norm(wp - protein_center))
        return float(np.clip(1.0 - dist / protein_radius_A, 0.0, 1.0))

    kept_xyz, sites = [], []
    for ijk, sc in zip(peaks, scores):
        refined_ijk = (int(ijk[0]), int(ijk[1]), int(ijk[2]))
        refinement_diagnostics = {
            "refinement_component_count": "",
            "refinement_component_size": "",
        }
        if edt_grid is not None:
            refined_ijk, refinement_diagnostics = _refine_surface_hotspot_ijk(
                ijk,
                score=score,
                edt_grid=edt_grid,
                pocket_mask=pocket,
                border_mask=border_mask,
                spacing=sp,
                raw_score=float(sc),
                return_diagnostics=True,
            )
        original_wp = vox2world(ijk)
        wp = vox2world(refined_ijk)
        if all(np.linalg.norm(wp - q) >= min_sep_A for q in kept_xyz):
            contact_frac = 1.0
            if edt_grid is not None:
                contact_frac = _surface_contact_fraction(
                    edt_grid,
                    refined_ijk,
                    sp,
                    shell_min_A=pocket_min_A,
                )
                if contact_frac < float(MIN_SURFACE_FRAC):
                    continue
            kept_xyz.append(wp)
            n = _npts_for_box_side(box_side_A, sp)
            # r_peak from EDT if available
            i_p, j_p, k_p = refined_ijk
            if edt_grid is not None:
                r_peak = float(edt_grid[i_p, j_p, k_p])
            else:
                r_peak = max(1.0, float(box_side_A) * 0.15)
            refined_score = float(score[i_p, j_p, k_p])
            sites.append(dict(
                site_id=f"S{len(kept_xyz)}",
                cx=float(wp[0]),
                cy=float(wp[1]),
                cz=float(wp[2]),
                nx=n, ny=n, nz=n,
                spacing=sp,
                r_peak=r_peak,
                F=refined_score,
                family="surface",
                contact_frac=float(contact_frac),
                center_closeness=center_closeness(wp),
                original_peak_x=float(original_wp[0]),
                original_peak_y=float(original_wp[1]),
                original_peak_z=float(original_wp[2]),
                refined_center_x=float(wp[0]),
                refined_center_y=float(wp[1]),
                refined_center_z=float(wp[2]),
                refinement_shift_A=float(np.linalg.norm(wp - original_wp)),
                refinement_component_size=refinement_diagnostics.get("refinement_component_size", ""),
                refinement_component_count=refinement_diagnostics.get("refinement_component_count", ""),
            ))
            if len(sites) == int(max_sites):
                break
    return _result(sites, score_context)


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


def _component_edt_weighted_medoid(mask, dist):
    """Choose a real cavity voxel near the EDT-weighted center of a component."""
    indices = np.argwhere(mask)
    if indices.size == 0:
        return None

    weights = dist[mask].astype(np.float64, copy=False)
    peak_pos = int(np.argmax(weights))
    peak_i, peak_j, peak_k = (int(v) for v in indices[peak_pos])
    peak_r_vox = float(weights[peak_pos])
    total_weight = float(weights.sum())
    if total_weight <= 0.0 or peak_r_vox <= 0.0:
        return peak_i, peak_j, peak_k, peak_r_vox

    weighted_center = np.average(indices.astype(np.float64), axis=0, weights=weights)
    depth_fraction = float(np.clip(INTERNAL_MEDOID_MIN_DEPTH_FRACTION, 0.0, 1.0))
    deep = weights >= (depth_fraction * peak_r_vox)
    if not np.any(deep):
        return peak_i, peak_j, peak_k, peak_r_vox

    deep_indices = indices[deep]
    deltas = deep_indices.astype(np.float64) - weighted_center
    medoid_pos = int(np.argmin(np.einsum("ij,ij->i", deltas, deltas)))
    i, j, k = (int(v) for v in deep_indices[medoid_pos])
    return i, j, k, peak_r_vox


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
    sC = _sigmoid_stable(Cz)
    sE = _sigmoid_stable(-0.7 * Ez)  # more negative E -> higher
    sD = Dn

    # combine
    F = 0.45 * np.nanmedian(sC) + 0.30 * np.nanmedian(sD) + 0.25 * np.nanmedian(sE)
    return float(F), float(np.nanmedian(sC)), float(np.nanmedian(sD)), float(np.nanmedian(sE))


### CENTERING AND BOXING ###
def voxel_to_world(ijk, origin, sp):
    return origin[0] + ijk[0] * sp, origin[1] + ijk[1] * sp, origin[2] + ijk[2] * sp


def adaptive_r_min_cavity(
    occ: np.ndarray,
    spacing: float,
    *,
    percentile: float = ADAPTIVE_R_MIN_PERCENTILE,
    peak_percentile: float = ADAPTIVE_R_MIN_PEAK_PERCENTILE,
    floor_A: float = ADAPTIVE_R_MIN_FLOOR_A,
    ceil_A: float = ADAPTIVE_R_MIN_CEIL_A,
    pocket_zone_max_A: float = ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A,
    peak_window_A: float = ADAPTIVE_R_MIN_PEAK_WINDOW_A,
) -> float:
    """
    Compute a pocket-size threshold from the protein's own EDT distribution.

    Instead of using a global fixed cutoff, derive the threshold from
    pocket-like EDT peaks. The earlier all-free-voxel percentile was biased by
    open solvent in whole-protein grids, especially for large receptors. This
    version first looks for local EDT maxima in the drug-sized pocket zone and
    uses the lower tail of significant peak radii as the shallowest cavity
    centre worth keeping.

    Parameters
    ----------
    occ       : boolean occupancy grid (True = protein atom present)
    spacing   : voxel spacing in Å
    percentile: fallback percentile when no usable local peaks are found
    peak_percentile: percentile of significant local EDT peak radii
    floor_A   : minimum r_min in Å
    ceil_A    : maximum r_min in Å

    Returns
    -------
    r_min_A : float, in Å
    """
    edt_A = distance_transform_edt(~occ).astype(np.float32) * float(spacing)
    free_edt_A = edt_A[~occ].ravel()
    if free_edt_A.size == 0:
        return floor_A

    pocket_zone = (~occ) & (edt_A >= float(floor_A)) & (edt_A <= float(pocket_zone_max_A))
    peak_radii_A = np.array([], dtype=np.float32)
    if np.any(pocket_zone):
        win = max(3, int(round(float(peak_window_A) / float(spacing))))
        local_max = maximum_filter(edt_A, size=win, mode="constant", cval=0.0)
        peak_mask = (edt_A == local_max) & pocket_zone
        peak_radii_A = edt_A[peak_mask].astype(np.float32)

    if peak_radii_A.size:
        zone_vals = edt_A[pocket_zone].astype(np.float32)
        significant_floor = max(float(floor_A), float(np.percentile(zone_vals, 25)))
        significant_peaks = peak_radii_A[peak_radii_A >= significant_floor]
        if significant_peaks.size:
            raw_r_min_A = float(np.percentile(significant_peaks, peak_percentile))
            source = f"p{peak_percentile:.0f} significant EDT peaks"
        else:
            raw_r_min_A = float(np.percentile(peak_radii_A, peak_percentile))
            source = f"p{peak_percentile:.0f} EDT peaks"
    else:
        clipped_free = free_edt_A[
            (free_edt_A >= float(floor_A)) & (free_edt_A <= float(pocket_zone_max_A))
        ]
        if clipped_free.size == 0:
            clipped_free = free_edt_A
        raw_r_min_A = float(np.percentile(clipped_free, percentile))
        source = f"p{percentile:.0f} pocket-zone EDT"

    r_min_A = raw_r_min_A
    r_min_A = float(np.clip(r_min_A, floor_A, ceil_A))
    print(f"[r_min] adaptive: {source} = {r_min_A:.2f} Å  "
          f"(free EDT range {float(free_edt_A.min()):.2f}–{float(free_edt_A.max()):.2f} Å, "
          f"clamped to [{floor_A}, {ceil_A}])")
    return r_min_A


def pick_centers(
    pdbqt,
    maps,
    map_origin,
    map_spacing,
    voxel_spacing=0.375,
    r_min=None,              # None → computed adaptively from the EDT distribution
    r_min_is_profiled=False,
    adaptive_r_min_params: dict | None = None,
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

    # Adaptive r_min: derive from the protein's own EDT distribution so we
    # never need to tune a global constant across different receptor sizes.
    if r_min is None:
        r_min = adaptive_r_min_cavity(occ, sp, **_adaptive_r_min_params(adaptive_r_min_params))
    elif r_min_is_profiled:
        print(f"[r_min] using adaptive r_min={float(r_min):.2f} Å")
    else:
        print(f"[r_min] using explicit r_min={float(r_min):.2f} Å (override)")

    dist, cc, ncc = internal_cavities(occ)
    print(f"[cav] internal components={ncc}, max_r_peak_vox={float(dist.max()):.2f} (Å)={float(dist.max() * sp):.2f}, r_min={r_min:.2f} Å")

    peak_records = []
    for lab in range(1, ncc + 1):
        mask = (cc == lab)
        if not mask.any():
            continue
        medoid = _component_edt_weighted_medoid(mask, dist)
        if medoid is None:
            continue
        i, j, k, r_peak_vox = medoid
        r_peak = r_peak_vox * sp
        peak_records.append((int(lab), int(i), int(j), int(k), float(r_peak)))

    search_r_min = _select_internal_search_r_min(
        float(r_min) if r_min is not None else None,
        [record[4] for record in peak_records],
    )
    if search_r_min is not None and search_r_min + 1e-6 < float(r_min):
        print(
            f"[r_min/internal] relaxing internal search threshold from {float(r_min):.2f} Å "
            f"to {float(search_r_min):.2f} Å based on enclosed-cavity peaks"
        )

    candidates = []
    for lab, i, j, k, r_peak in peak_records:
        if search_r_min is not None and r_peak < search_r_min:
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
        n = _npts_for_box_side(2.0 * half, map_spacing)
        if n == 255:
            actual_half = (255 - 1) * map_spacing / 2.0
            print(f"[box/S{idx}] WARNING: npts hit 255 cap — box half={actual_half:.1f} Å "
                  f"but cavity r_peak={r_peak:.1f} Å wants half={half:.1f} Å")
        sites.append(dict(
            site_id=f"S{idx}",
            cx=ctr[0], cy=ctr[1], cz=ctr[2],
            nx=n, ny=n, nz=n,
            spacing=map_spacing,
            r_peak=r_peak,
            F=F,
            family="internal",
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
      'shape': (dim1,dim2,dim3),  # for reshaping .map arrays (Fortran order)
      'dir': Path,           # directory of the fld
      'map_paths': {'C': Path, 'E': Path, 'D': Path},
      'affinity_types': tuple[str, ...],
    }
    """
    import re
    from pathlib import Path

    p = Path(fld_or_dir)
    if not p.is_file():
        # try to find a .fld inside the directory
        cands = sorted(Path(fld_or_dir).glob("*.fld"))
        if len(cands) != 1:
            raise RuntimeError(f"Expected one FLD in {fld_or_dir}, found {len(cands)}: {cands}")
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

    dims = (dim1, dim2, dim3)

    # compute origin from center + spacing + npts (true npts, NOT dims)
    if spacing is None or center is None:
        raise RuntimeError("Missing #SPACING or #CENTER in .fld header comments.")
    if nelems is None:
        # fallback: infer npts = dims - 1
        nelems = (dim1 - 1, dim2 - 1, dim3 - 1)
    expected_dims = tuple(n + 1 for n in nelems)
    if dims != expected_dims:
        raise ValueError(
            f"FLD dimensions {dims} do not equal "
            f"NELEMENTS+1 {expected_dims}"
        )

    shape = dims
    origin = tuple(
        float(center_value) - 0.5 * float(nelement) * float(spacing)
        for center_value, nelement in zip(center, nelems)
    )

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
    affinity_types = []
    suffix = "-affinity"
    for label in labels:
        label_clean = label.strip()
        if label_clean.lower().endswith(suffix):
            affinity_types.append(label_clean[: -len(suffix)])

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
    def _path_for_label(label_name):
        wanted = str(label_name).strip().lower()
        for label, path in name_to_path.items():
            if str(label).strip().lower() == wanted:
                return path
        return None

    def _map_candidates_by_suffixes(*suffixes):
        matches = []
        for map_path in d.glob("*.map"):
            if any(map_path.name.endswith(suffix) for suffix in suffixes):
                matches.append(map_path)
        return sorted(matches)

    E_path = _path_for_label("Electrostatics")
    D_path = _path_for_label("Desolvation")

    if not C_path:
        # fallback glob (avoid Cl/Br etc.)
        cands = _map_candidates_by_suffixes(".C.map")
        if cands:
            C_path = cands[0]
    if not E_path:
        # often lowercase 'e.map'
        cands = _map_candidates_by_suffixes(".E.map", ".e.map")
        if cands:
            E_path = cands[0]
    if not D_path:
        cands = _map_candidates_by_suffixes(".D.map", ".d.map")
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
        "affinity_types": tuple(affinity_types),
    }


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _map_diagnostics(label, map_path, arr):
    values = np.asarray(arr)
    finite_mask = np.isfinite(values)
    finite = values[finite_mask]
    prefix = f"map_{label}"
    shape_value = "x".join(str(int(x)) for x in values.shape)
    diagnostics = {
        f"{prefix}_filename": Path(map_path).name,
        f"{prefix}_sha256": _sha256_file(map_path),
        f"{prefix}_shape": shape_value,
        f"{prefix}_value_count": str(int(values.size)),
        f"{prefix}_finite_count": str(int(finite.size)),
        f"{prefix}_nonfinite_count": str(int(values.size - finite.size)),
    }
    if finite.size:
        p01, median, p99 = np.percentile(finite, [1.0, 50.0, 99.0])
        diagnostics.update(
            {
                f"{prefix}_min": f"{float(np.min(finite)):.9g}",
                f"{prefix}_max": f"{float(np.max(finite)):.9g}",
                f"{prefix}_p01": f"{float(p01):.9g}",
                f"{prefix}_median": f"{float(median):.9g}",
                f"{prefix}_p99": f"{float(p99):.9g}",
            }
        )
    else:
        diagnostics.update(
            {
                f"{prefix}_min": "nan",
                f"{prefix}_max": "nan",
                f"{prefix}_p01": "nan",
                f"{prefix}_median": "nan",
                f"{prefix}_p99": "nan",
            }
        )
    return diagnostics


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
    need = int(np.prod(shape))
    arr = np.asarray(vals, dtype=np.float32)
    if arr.size != need:
        raise ValueError(
            f"{map_path}: expected {need} values for shape {shape}, "
            f"found {arr.size}"
        )
    arr = arr.reshape(shape, order="F")  # x fastest
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
    print(f"[OK] wrote {len(rows)} sites --> {path}")


# ===== Single-site grid builder (ligand / residues / blind) =====
import math


def _odd(n: int) -> int:
    return n if (n % 2 == 1) else (n + 1)


def _npts_for_box_side(side_len_A: float, spacing: float, clamp=(60, 255)) -> int:
    """Convert a physical box side length in Angstrom to odd AutoGrid npts."""
    n = int(math.ceil(float(side_len_A) / float(spacing)))
    n = _odd(n)
    return max(clamp[0], min(clamp[1], n))


def _npts_for_size(side_len_A: float, spacing: float, clamp=(60, 255)) -> int:
    """Backward-compatible alias for side-length based box sizing."""
    return _npts_for_box_side(side_len_A, spacing, clamp=clamp)


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
    margin: float = GRID_MARGIN,
    cap: float = GRID_CAP,
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

    expected_fld = out_dir_p / f"{Path(receptor_pdbqt).stem}.maps.fld"
    if not expected_fld.exists():
        fld_candidates = sorted(out_dir_p.glob("*.fld"))
        raise FileNotFoundError(f"Expected AutoGrid FLD {expected_fld}, found {fld_candidates}; see {log}")
    fld_path = str(expected_fld.resolve())
    return {
        "site_id": "S1",
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
    n = int(math.ceil(side / float(spacing)))
    n = _odd(max(n, 25))
    return min(n, 129)


# ========= HotspotGPFGenerator: whole maps + centers.tsv + per-site grids =========
from pathlib import Path
import os
import subprocess

_AD4_TYPES = [
    "A", "C", "HD", "N", "NA", "OA", "SA", "S", "P", "Cl", "F", "Br", "I", "Zn", "Fe", "Mg", "Mn", "Ca"
]


def _normalize_autogrid_map_types(map_types=None) -> tuple[str, ...]:
    """
    Return ligand map types to ask AutoGrid to compute.

    None preserves historical behavior and generates every AD4 affinity map.
    CaV-EMPS policies reject selective maps because AutoGrid dsolvmap depends
    on the requested ligand type set.
    """
    if map_types is None:
        return tuple(_AD4_TYPES)
    if isinstance(map_types, str):
        raw = [part.strip() for part in re.split(r"[,\s]+", map_types) if part.strip()]
    else:
        raw = [str(part).strip() for part in map_types if str(part).strip()]
    if not raw:
        return tuple(_AD4_TYPES)
    valid = set(_AD4_TYPES)
    normalized = []
    for atom_type in raw:
        if atom_type not in valid:
            raise ValueError(f"unsupported AutoGrid map type {atom_type!r}; expected one of {', '.join(_AD4_TYPES)}")
        if atom_type not in normalized:
            normalized.append(atom_type)
    return tuple(normalized)


def _ensure_odd_clamped(nxyz, clamp=(60, 255)):
    import numpy as _np
    n = _np.array(nxyz, int)
    n = n + (n % 2 == 0)
    n = _np.clip(n, clamp[0], clamp[1])
    return tuple(int(x) for x in n.tolist())


def _write_site_gpf(
    gpf_path,
    receptor_pdbqt,
    center,
    npts,
    spacing,
    autogrid4_bin="autogrid4",
    map_types=None,
    npts_max: int = DOCKING_GRID_NPTS_MAX,
):
    """GPF writer that matches ligand_types ⇔ map lines 1:1 (+ elec/dsolv maps)."""
    rec_path = Path(receptor_pdbqt).resolve()
    rec_stem = rec_path.stem
    cx, cy, cz = [float(v) for v in center]
    npts_max = int(npts_max)
    if not 25 <= npts_max <= AUTOGRID4_NPTS_MAX:
        raise ValueError(
            f"npts_max must be between 25 and {AUTOGRID4_NPTS_MAX}, found {npts_max}"
        )
    nx, ny, nz = _ensure_odd_clamped(npts, clamp=(60, npts_max))
    receptor_types = receptor_atom_types(rec_path)
    ligand_types = _normalize_autogrid_map_types(map_types)
    param_file = _find_ad4_parameter_file(autogrid4_bin)
    with open(gpf_path, "w") as f:
        # parameter_file MUST come first — AutoGrid reads parameters before
        # processing any atom types.  Without it, all maps are zero.
        if param_file:
            f.write(f"parameter_file {param_file}\n")
        f.write(f"gridfld {rec_stem}.maps.fld\n")
        f.write(f"npts {nx} {ny} {nz}\n")
        f.write(f"spacing {float(spacing):.3f}\n")
        f.write(f"receptor_types {' '.join(receptor_types)}\n")
        f.write(f"receptor {rec_path}\n")
        f.write(f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}\n")
        f.write(f"ligand_types {' '.join(ligand_types)}\n")
        f.write("smooth 0.500\n")
        for t in ligand_types:
            f.write(f"map {rec_stem}.{t}.map\n")
        f.write(f"elecmap {rec_stem}.e.map\n")
        f.write(f"dsolvmap {rec_stem}.d.map\n")
        f.write("dielectric -0.1465\n")
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
    npts_max: int = DOCKING_GRID_NPTS_MAX,
):
    """
    Returns center, (nx,ny,nz) for GPF npts, and spacing.
    Uses larger spacing if needed so the box fully covers the receptor AABB (+margin)
    while keeping npts <= npts_max.
    """
    npts_max = int(npts_max)
    if not 25 <= npts_max <= AUTOGRID4_NPTS_MAX:
        raise ValueError(
            f"npts_max must be between 25 and {AUTOGRID4_NPTS_MAX}, found {npts_max}"
        )

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
    npts_max: int = DOCKING_GRID_NPTS_MAX,
    autogrid4_bin: str = "autogrid4",
    map_types=None,
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
    _write_site_gpf(
        gpf_path,
        receptor_pdbqt,
        center,
        npts,
        sp,
        autogrid4_bin=autogrid4_bin,
        map_types=map_types,
        npts_max=npts_max,
    )

    log_path = out_root / "grid.glg"
    cmd = [autogrid4_bin, "-p", gpf_path.name, "-l", log_path.name]
    res = subprocess.run(cmd, cwd=out_root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"autogrid whole-protein failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )

    expected_fld = out_root / f"{rec_stem}.maps.fld"
    if not expected_fld.exists():
        fld_candidates = sorted(out_root.glob("*.fld"))
        raise FileNotFoundError(f"Expected AutoGrid FLD {expected_fld}, found {fld_candidates}. Check {log_path}")
    fld_path = str(expected_fld.resolve())

    # (Optional) post-check: make sure the produced fld bounds cover receptor too (uses your meta loader)
    # meta = load_fld_or_map_meta(fld_path)
    # validate_grid_bounds_vs_receptor(receptor_pdbqt, meta)

    return fld_path


class HotspotGPFGenerator:
    """
    Orchestrates:
    (A) whole-protein maps → .fld
    (B) centers.tsv (receptor_search, exhaustive_search, or chosen family)
    (C) per-site GPF + AutoGrid → per-site .fld
    """

    def __init__(self, autogrid4_bin="autogrid4"):
        self.autogrid4_bin = autogrid4_bin

    def prepare_centers_and_grids(
        self,
        receptor_pdbqt: str,
        out_root: str,
        centers_tsv_path: str,
        mode: str = "receptor_search",  # receptor_search | exhaustive_search | internal | surface | hybrid
        n_sites: int = AUTOSITES,
        whole_spacing: float = GRID_SPACING,
        whole_cap_ang: float = GRID_CAP,
        hotspot_box_ang: float = HOTSPOT_BOX_ANGLE,
        tau_rel: float = 0.60,
        min_sep_A: float = HOTSPOT_NMS_MINSEP_A,
        r_min: float | None = R_MIN_CAVITY_A,
        k_box: float = 4.0,
        adaptive_r_min_params: dict | None = None,
        maps_pocket_max_A: float | None = None,
        nms_box_fraction: float = HOTSPOT_NMS_BOX_FRACTION,
        nms_min_A: float = HOTSPOT_NMS_MIN_A,
        nms_max_A: float = HOTSPOT_NMS_MAX_A,
        whole_npts_max: int = CAV_EMPS_WHOLE_NPTS_MAX,
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
            npts_max=whole_npts_max,
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
            tau_rel=tau_rel,
            min_sep_A=min_sep_A,
            k_box=k_box,
            r_min=r_min,
            mode=mode,
            adaptive_r_min_params=adaptive_r_min_params,
            maps_pocket_max_A=maps_pocket_max_A,
            nms_box_fraction=nms_box_fraction,
            nms_min_A=nms_min_A,
            nms_max_A=nms_max_A,
            whole_map_npts_max=whole_npts_max,
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
    ap.add_argument("--r-min", type=float, default=R_MIN_CAVITY_A, help="min inscribed-sphere radius (Å); None uses adaptive profiling")
    ap.add_argument("--min-sep", type=float, default=HOTSPOT_NMS_MINSEP_A, help="NMS min separation between hotspots (Å)")
    ap.add_argument("--k-box", type=float, default=10.0, help="box half-size multiplier × r_peak")
    ap.add_argument("--max-sites", type=int, default=8)
    ap.add_argument("--inflate", type=float, default=0.25, help="inflate atom radii by this many Å when voxelizing (0.0–0.5 typical)")
    ap.add_argument(
        "--mode",
        choices=["receptor_search", "exhaustive_search", "internal", "surface", "maps", "hybrid"],
        default="receptor_search",
        help="receptor_search=internal+surface+hybrid; exhaustive_search=multi-site union; internal/surface/hybrid run one family",
    )
    ap.add_argument("--tau-rel", type=float, default=0.52, help="relative peak threshold for maps mode (0.58–0.62 typical)")
    ap.add_argument("--half-size", type=float, default=12.0, help="box half-size (Å) for maps mode")
    args = ap.parse_args()
    effective_min_sep_A = _effective_hotspot_min_sep_A(args.min_sep, 2.0 * args.half_size)
    if effective_min_sep_A > float(args.min_sep) + 1e-6:
        print(
            f"[nms] clamped site separation from {float(args.min_sep):.2f} Å "
            f"to {effective_min_sep_A:.2f} Å for {2.0 * float(args.half_size):.1f} Å boxes"
        )
    candidate_min_sep_A = _candidate_harvest_min_sep_A(effective_min_sep_A, int(args.max_sites))
    if candidate_min_sep_A + 1e-6 < effective_min_sep_A:
        print(
            f"[nms] harvesting candidates at {candidate_min_sep_A:.2f} Å; "
            f"final site separation target is {effective_min_sep_A:.2f} Å"
        )

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
            min_sep_A=candidate_min_sep_A, k_box=args.k_box, max_sites=args.max_sites,
            inflate_A=args.inflate
        )
    elif args.mode in {"maps", "surface"}:
        print("[maps] detecting surface-pocket hotspots from C/E/D maps")
        sites = detect_maps_hotspots(maps["C"], maps["E"], maps["D"], origin, spacing,
                                     tau_rel=args.tau_rel, min_sep_A=max(candidate_min_sep_A, float(SURFACE_NMS_MINSEP_A)),
                                     max_sites=args.max_sites, box_side_A=2.0 * args.half_size,
                                     r_min_A=args.r_min)
    elif args.mode == "hybrid":
        print("[hybrid] building consensus sites from internal + surface candidates")
        internal_sites = pick_centers(
            args.receptor_pdbqt, maps,
            map_origin=origin, map_spacing=spacing,
            voxel_spacing=0.5, r_min=None,  # adaptive: derived from each receptor's own EDT
            min_sep_A=max(candidate_min_sep_A, 7.0), k_box=4.0, max_sites=args.max_sites, inflate_A=0.0,
        )
        surface_sites = detect_maps_hotspots(
            maps["C"],
            maps["E"],
            maps["D"],
            origin,
            spacing,
            tau_rel=args.tau_rel,
            min_sep_A=max(candidate_min_sep_A, float(SURFACE_NMS_MINSEP_A)),
            max_sites=args.max_sites,
            box_side_A=2.0 * args.half_size,
            r_min_A=args.r_min,
        )
        sites = _build_hybrid_sites(
            internal_sites,
            surface_sites,
            min_sep_A=max(candidate_min_sep_A, float(MAX_CENTER_DIST_A) * 0.5),
            max_sites=args.max_sites,
        )
    else:
        print(f"[{args.mode}] generating centers via policy")
        autogenerate_centers_tsv(
            receptor_pdbqt=args.receptor_pdbqt,
            out_root=args.out,
            centers_tsv_path=str(Path(args.out) / "centers.tsv"),
            n_sites=args.max_sites,
            default_spacing=spacing,
            hotspot_box_ang=2.0 * args.half_size,
            tau_rel=args.tau_rel,
            min_sep_A=effective_min_sep_A,
            r_min=args.r_min,
            mode=args.mode,
        )
        return

    # write TSV
    rec_stem = Path(args.receptor_pdbqt).stem
    out_tsv = Path(args.out) / "centers.tsv"
    write_centers_tsv(out_tsv, rec_stem, sites)


if __name__ == "__main__":
    main()
