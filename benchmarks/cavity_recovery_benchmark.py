#!/usr/bin/env python3
"""Batch benchmark Ultidock auto-sites against DUD-E crystal-ligand centers."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import importlib
import json
import math
import os
import subprocess
import sys
import tempfile
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_ROOT = REPO_ROOT / "docking"
DEFAULT_DATASET_ROOT = REPO_ROOT / "benchmarks" / "datasets"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "benchmarks" / "results" / "cavity_recovery"
DEFAULT_AUTOGRID4 = DOCKING_ROOT / "AUTODOCK_GPU_DIR" / "autogrid" / "autogrid4"
DEFAULT_THRESHOLDS = (4.0, 5.0, 10.0)


def _default_benchmark_jobs() -> int:
    for env_name in ("ULTIDOCK_BENCHMARK_WORKERS", "ULTIDOCK_WORKERS"):
        raw_value = os.environ.get(env_name)
        if raw_value:
            try:
                return max(1, int(raw_value))
            except ValueError:
                pass
    cpu_count = os.cpu_count() or 1
    return max(1, min(4, cpu_count // 2 or 1))

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from benchmarks.evaluate_cavities import euclidean_distance, parse_centers_tsv
except ModuleNotFoundError:
    from evaluate_cavities import euclidean_distance, parse_centers_tsv  # type: ignore

from molguard.io.receptor_prep import (
    SUPPORTED_RECEPTOR_SUFFIXES,
    prepare_receptor_pdbqt,
    receptor_pdbqt_output_name,
)


SITE_POLICIES = ("receptor_search", "exhaustive_search", "internal", "surface", "hybrid")


def _site_sort_key(site_id: str) -> tuple[int, str]:
    label = site_id.strip()
    if label.upper().startswith("S") and label[1:].isdigit():
        return int(label[1:]), label
    if label.isdigit():
        return int(label), label
    return 10**9, label


def _first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def _load_metadata(target_dir: Path) -> dict[str, Any]:
    metadata_path = target_dir / "benchmark.json"
    if not metadata_path.exists():
        return {}
    return json.loads(metadata_path.read_text(encoding="utf-8"))


def _reference_candidates(target_dir: Path, metadata: dict[str, Any]) -> list[Path]:
    candidates: list[Path] = []
    if metadata.get("reference_ligand"):
        candidates.append(target_dir / str(metadata["reference_ligand"]))
    candidates.extend(
        target_dir / name
        for name in (
            "crystal_ligand.mol2",
            "crystal_ligand.mol2.gz",
            "crystal_ligand.pdbqt",
            "crystal_ligand.pdbqt.gz",
            "crystal_ligand.pdb",
            "crystal_ligand.pdb.gz",
            "ligand.mol2",
            "ligand.mol2.gz",
            "ref_ligand.mol2",
            "ref_ligand.mol2.gz",
        )
    )
    return candidates


def _receptor_candidates(target_dir: Path, metadata: dict[str, Any]) -> list[Path]:
    candidates: list[Path] = []
    if metadata.get("receptor"):
        candidates.append(target_dir / str(metadata["receptor"]))
    candidates.extend(
        target_dir / f"receptor{suffix}" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    candidates.extend(
        target_dir / f"receptor{suffix}.gz" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    candidates.extend(
        target_dir / f"protein{suffix}" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    candidates.extend(
        target_dir / f"protein{suffix}.gz" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    return candidates


def load_target_files(dataset_root: Path, target_name: str) -> tuple[Path, Path]:
    target_dir = (dataset_root / target_name).resolve()
    if not target_dir.exists():
        raise FileNotFoundError(f"target directory not found: {target_dir}")
    metadata = _load_metadata(target_dir)
    receptor = _first_existing(_receptor_candidates(target_dir, metadata))
    reference = _first_existing(_reference_candidates(target_dir, metadata))
    if receptor is None:
        raise FileNotFoundError(f"{target_name}: receptor file not found")
    if reference is None:
        raise FileNotFoundError(f"{target_name}: crystal/reference ligand file not found")
    return receptor, reference


def discover_targets(dataset_root: Path) -> list[str]:
    targets: list[str] = []
    if not dataset_root.exists():
        return targets
    for target_dir in sorted(path for path in dataset_root.iterdir() if path.is_dir()):
        metadata = _load_metadata(target_dir)
        if _first_existing(_receptor_candidates(target_dir, metadata)) and _first_existing(
            _reference_candidates(target_dir, metadata)
        ):
            targets.append(target_dir.name.lower())
    return targets


def parse_targets(value: str | None, dataset_root: Path) -> list[str]:
    discovered = discover_targets(dataset_root)
    if value is None or value.strip().lower() == "all":
        return discovered
    requested = [part.strip().lower() for part in value.split(",") if part.strip()]
    return requested


def parse_thresholds(value: str) -> list[float]:
    thresholds: list[float] = []
    for part in value.split(","):
        stripped = part.strip()
        if not stripped:
            continue
        threshold = float(stripped)
        if threshold <= 0:
            raise ValueError("thresholds must be positive")
        thresholds.append(threshold)
    if not thresholds:
        raise ValueError("at least one threshold is required")
    return sorted(set(thresholds))


def threshold_label(value: float) -> str:
    label = f"{value:g}".replace("-", "m").replace(".", "p")
    return f"{label}a"


def _effective_suffix(path: Path) -> str:
    suffixes = [part.lower() for part in path.suffixes]
    if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2:
        return suffixes[-2]
    return path.suffix.lower()


def calculate_reference_center(path: Path) -> tuple[float, float, float]:
    """Calculate a geometric center from MOL2, PDB, or PDBQT coordinates."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    suffix = _effective_suffix(path)
    coords: list[tuple[float, float, float]] = []

    if suffix == ".mol2":
        in_atom_block = False
        with opener(path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.rstrip("\n")
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atom_block = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atom_block:
                    break
                if not in_atom_block:
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                except ValueError:
                    continue
    elif suffix in {".pdb", ".pdbqt"}:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
    else:
        raise ValueError(f"unsupported reference ligand format: {path}")

    if not coords:
        raise ValueError(f"no coordinates found in {path}")
    count = float(len(coords))
    return (
        sum(coord[0] for coord in coords) / count,
        sum(coord[1] for coord in coords) / count,
        sum(coord[2] for coord in coords) / count,
    )


def git_value(*args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def git_status_short() -> str | None:
    return git_value("status", "--short")


def _setup_equivalent_config_defaults() -> dict[str, Any]:
    """Fallback defaults matching the generated config written by setup.py."""
    return {
        "AUTODOCK_GPU_DIR": str(DOCKING_ROOT / "AUTODOCK_GPU_DIR"),
        "DOCKING_DIR": str(DOCKING_ROOT / "DOCKING_DIR"),
        "RESULTS_DIR": str(DOCKING_ROOT / "RESULTS_DIR"),
        "NUMWI": "128",
        "LIGANDS_DIR": str(DOCKING_ROOT / "LIGANDS_DIR"),
        "CENTERS_TSV": None,
        "GRID_MODE": "centers",
        "GRID_SPACING": 0.375,
        "GRID_MARGIN": 5.0,
        "GRID_CAP": 150.0,
        "R_MIN_CAVITY_A": None,
        "HOTSPOT_BOX_ANGLE": 35,
        "HOTSPOT_NMS_MINSEP_A": 14.0,
        "SURFACE_SHELL__MIN_A": 2.0,
        "SURFACE_SHELL__MAX_A": 20.0,
        "SURFACE_NMS_MINSEP_A": 5.0,
        "MAX_CENTER_DIST_A": 10.0,
        "CONTACT_SHELL_A": 4.0,
        "MIN_SURFACE_FRAC": 0.01,
        "AUTOSITES": 6,
        "ADAPTIVE_R_MIN_PERCENTILE": 85.0,
        "ADAPTIVE_R_MIN_PEAK_PERCENTILE": 10.0,
        "ADAPTIVE_R_MIN_FLOOR_A": 2.0,
        "ADAPTIVE_R_MIN_CEIL_A": 5.0,
        "ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A": 8.0,
        "ADAPTIVE_R_MIN_PEAK_WINDOW_A": 4.0,
        "MAPS_POCKET_MAX_A": 15.0,
        "HOTSPOT_NMS_BOX_FRACTION": 0.40,
        "HOTSPOT_NMS_MIN_A": 4.0,
        "HOTSPOT_NMS_MAX_A": 25.0,
        "SITE_POLICY": "receptor_search",
        "AUTO_GRID_BIN": str(DEFAULT_AUTOGRID4),
        "REF_LIGAND_PDB": None,
        "VINA_CPU": 2,
        "VINA_SEED": None,
        "VINA_EXHAUSTIVENESS": 8,
        "VINA_NUM_MODES": 9,
        "BENCHMARK_MODE": True,
    }


def _load_docking_config_defaults() -> dict[str, Any]:
    """
    Load algorithm knob defaults from the real ``docking/config.py``.

    This keeps the benchmark in sync with the mainline pipeline: any parameter
    tuned in ``docking/config.py`` is automatically picked up here without
    needing a matching CLI override in the benchmark script. If the generated
    config was removed by ``ultidock clean --all``, fall back to the same
    defaults that ``setup.py`` writes.
    """
    import importlib.util

    values = _setup_equivalent_config_defaults()
    real_config_path = DOCKING_ROOT / "config.py"
    if not real_config_path.exists():
        print(
            f"[config] {real_config_path} not found; using setup-equivalent "
            "fallback defaults"
        )
        return values
    spec = importlib.util.spec_from_file_location("_real_docking_config", real_config_path)
    mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    values.update({
        key: getattr(mod, key)
        for key in dir(mod)
        if key.isupper() and not key.startswith("_")
    })
    return values


def write_runtime_config(config_dir: Path, args: argparse.Namespace) -> Path:
    """
    Write a scratch ``config.py`` for the benchmark run.

    Algorithm knobs (``R_MIN_CAVITY_A``, ``HOTSPOT_BOX_ANGLE``, etc.) are
    inherited from the real ``docking/config.py`` so the benchmark uses the
    **same tuned defaults as the mainline pipeline**.  Only path/directory
    variables and ``AUTOSITES`` are overridden for the sandbox.

    If ``--r-min`` was explicitly supplied on the CLI it still takes priority.
    """
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "config.py"
    scratch_root = Path(args.output_dir).resolve() / "_runtime"

    # Start from real docking/config.py so algorithm knobs are inherited.
    values = _load_docking_config_defaults()

    # Override only path/runtime variables for the benchmark sandbox.
    values.update(
        {
            "AUTODOCK_GPU_DIR": str(Path(args.autodock_gpu_dir).resolve()),
            "DOCKING_DIR": str(scratch_root / "DOCKING_DIR"),
            "RESULTS_DIR": str(scratch_root / "RESULTS_DIR"),
            "NUMWI": "128",
            "LIGANDS_DIR": str(scratch_root / "LIGANDS_DIR"),
            "CENTERS_TSV": None,
            "GRID_MODE": "centers",
            "GRID_SPACING": float(args.grid_spacing),
            "GRID_MARGIN": 5.0,
            "GRID_CAP": float(args.blind_cap),
            "AUTO_GRID_BIN": str(Path(args.autogrid4_bin).resolve()),
            "AUTOSITES": int(args.autosites),
        }
    )

    # Honour an explicit --r-min override; None means use adaptive per receptor.
    if args.r_min is not None:
        values["R_MIN_CAVITY_A"] = float(args.r_min)

    lines = [
        "# Auto-generated by benchmarks/cavity_recovery_benchmark.py",
        "# Algorithm knobs inherited from docking/config.py unless explicitly overridden.",
    ]
    for key, value in values.items():
        lines.append(f"{key} = {value!r}")
    config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return config_path


def normalize_runtime_paths(args: argparse.Namespace) -> None:
    if args.autogrid4_bin is None:
        args.autogrid4_bin = str(
            Path(args.autodock_gpu_dir).resolve() / "autogrid" / "autogrid4"
        )


def _is_executable_file(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def ensure_generation_runtime(args: argparse.Namespace) -> None:
    if args.skip_generation:
        return

    missing: list[str] = []
    config_path = DOCKING_ROOT / "config.py"
    autogrid4_bin = Path(args.autogrid4_bin).resolve()

    if not config_path.exists():
        missing.append(str(config_path))
    if not _is_executable_file(autogrid4_bin):
        missing.append(str(autogrid4_bin))
    if not missing:
        return
    if args.no_auto_setup:
        raise RuntimeError(
            "Benchmark runtime is incomplete and --no-auto-setup was set. "
            f"Missing: {', '.join(missing)}"
        )

    setup_py = DOCKING_ROOT / "setup.py"
    command = [
        sys.executable,
        str(setup_py),
        "--mode",
        str(args.auto_setup_mode),
        "--skip-wget",
        "--skip-profile",
        "--benchmark",
        "--autodock-gpu-dir",
        str(Path(args.autodock_gpu_dir).resolve()),
        "--grid-mode",
        "centers",
        "--grid-spacing",
        str(float(args.grid_spacing)),
        "--grid-cap",
        str(float(args.blind_cap)),
        "--autosites",
        str(int(args.autosites)),
    ]
    print("[preflight] Benchmark runtime is incomplete:")
    for item in missing:
        print(f"  - missing {item}")
    print(
        "[preflight] Running minimal setup: "
        f"setup.py --mode {args.auto_setup_mode} --skip-wget --skip-profile"
    )
    try:
        subprocess.run(command, cwd=REPO_ROOT, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "Automatic benchmark setup failed. Re-run setup.py manually or pass "
            "--no-auto-setup to fail fast before generation."
        ) from exc

    if not _is_executable_file(autogrid4_bin):
        raise RuntimeError(
            "Automatic benchmark setup finished, but AutoGrid is still missing "
            f"or not executable at {autogrid4_bin}"
        )


def load_make_grids(config_dir: Path):
    config_dir_str = str(config_dir)
    docking_root_str = str(DOCKING_ROOT)
    for path in (config_dir_str, docking_root_str):
        while path in sys.path:
            sys.path.remove(path)
    sys.path.insert(0, docking_root_str)
    sys.path.insert(0, config_dir_str)
    sys.modules.pop("config", None)
    sys.modules.pop("make_grids", None)
    importlib.invalidate_caches()
    try:
        return importlib.import_module("make_grids")
    except ModuleNotFoundError as exc:
        missing = exc.name or "a required Python package"
        raise RuntimeError(
            "Could not import docking/make_grids.py because "
            f"{missing!r} is missing. Run this benchmark in the Ultidock "
            "scientific environment with numpy/scipy installed."
        ) from exc


def prepare_receptor(
    *,
    receptor_input: Path,
    output_dir: Path,
    prepare_command: str | None,
    seed: int,
    force: bool,
) -> Path:
    prepared = output_dir / "prepared" / receptor_pdbqt_output_name(receptor_input)
    if prepared.exists() and not force:
        return prepared
    prepare_receptor_pdbqt(
        input_path=receptor_input,
        output_path=prepared,
        prepare_command=prepare_command,
        seed=seed,
        timestamp="CAVITY_BENCHMARK",
    )
    return prepared


def existing_centers_candidates(
    target_name: str,
    output_dir: Path,
    centers_root: Path | None,
) -> list[Path]:
    candidates = [output_dir / target_name / "centers.tsv"]
    if centers_root is not None:
        target_root = centers_root / target_name
        candidates.append(target_root / "centers.tsv")
        candidates.extend(sorted(target_root.glob("*/centers.tsv")))
    return candidates


def generate_centers(
    *,
    make_grids,
    receptor_pdbqt: Path,
    target_output: Path,
    args: argparse.Namespace,
) -> Path:
    centers_tsv = target_output / "centers.tsv"
    if args.force and centers_tsv.exists():
        centers_tsv.unlink()
    # r_min=None triggers the adaptive_r_min_cavity() path inside make_grids
    # (same behaviour as the mainline pipeline).  A non-None value means the
    # user explicitly requested a fixed threshold via --r-min.
    r_min_arg: float | None = float(args.r_min) if args.r_min is not None else None
    make_grids.autogenerate_centers_tsv(
        receptor_pdbqt=str(receptor_pdbqt),
        out_root=str(target_output / "maps"),
        centers_tsv_path=str(centers_tsv),
        n_sites=int(args.autosites),
        default_spacing=float(args.grid_spacing),
        blind_cap=float(args.blind_cap),
        autogrid4_bin=str(Path(args.autogrid4_bin).resolve()),
        hotspot_box_ang=float(args.hotspot_box_size),
        tau_rel=float(args.tau_rel),
        min_sep_A=float(args.min_sep),
        r_min=r_min_arg,
        mode=str(args.site_policy),
    )
    if not centers_tsv.exists():
        raise FileNotFoundError(f"center generation did not produce {centers_tsv}")
    return centers_tsv


def vector_delta(
    reference: tuple[float, float, float],
    center: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (center[0] - reference[0], center[1] - reference[1], center[2] - reference[2])


def required_cube_side(delta: tuple[float, float, float]) -> float:
    return 2.0 * max(abs(delta[0]), abs(delta[1]), abs(delta[2]))


def within_cube(delta: tuple[float, float, float], box_size: float) -> bool:
    half = box_size / 2.0
    return abs(delta[0]) <= half and abs(delta[1]) <= half and abs(delta[2]) <= half


def evaluate_centers(
    *,
    target_name: str,
    reference_ligand: Path,
    centers_tsv: Path,
    thresholds: list[float],
    box_size: float,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    crystal_center = calculate_reference_center(reference_ligand)
    predicted_sites = parse_centers_tsv(centers_tsv)
    if not predicted_sites:
        raise ValueError(f"no predicted sites found in {centers_tsv}")

    site_rows: list[dict[str, Any]] = []
    for site_order, site_id in enumerate(sorted(predicted_sites, key=_site_sort_key), start=1):
        center = predicted_sites[site_id]
        distance = euclidean_distance(crystal_center, center)
        delta = vector_delta(crystal_center, center)
        row: dict[str, Any] = {
            "target": target_name,
            "site_id": site_id,
            "site_order": site_order,
            "distance_a": distance,
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "delta_x": delta[0],
            "delta_y": delta[1],
            "delta_z": delta[2],
            "required_cube_side_a": required_cube_side(delta),
            f"within_{threshold_label(box_size)}_cube": within_cube(delta, box_size),
        }
        for threshold in thresholds:
            row[f"hit_at_{threshold_label(threshold)}"] = distance <= threshold
        site_rows.append(row)

    best = min(site_rows, key=lambda row: float(row["distance_a"]))
    summary: dict[str, Any] = {
        "target": target_name,
        "status": "ok",
        "n_sites": len(site_rows),
        "crystal_center_x": crystal_center[0],
        "crystal_center_y": crystal_center[1],
        "crystal_center_z": crystal_center[2],
        "best_site_id": best["site_id"],
        "best_site_order": best["site_order"],
        "best_distance_a": best["distance_a"],
        "best_required_cube_side_a": best["required_cube_side_a"],
        f"best_within_{threshold_label(box_size)}_cube": best[
            f"within_{threshold_label(box_size)}_cube"
        ],
        "reference_ligand": str(reference_ligand.resolve()),
        "centers_tsv": str(centers_tsv.resolve()),
        "error": "",
    }
    for threshold in thresholds:
        label = threshold_label(threshold)
        summary[f"any_site_success_at_{label}"] = bool(best[f"hit_at_{label}"])

    evaluation = {
        "target": target_name,
        "reference_ligand": str(reference_ligand.resolve()),
        "centers_tsv": str(centers_tsv.resolve()),
        "crystal_center": {
            "x": crystal_center[0],
            "y": crystal_center[1],
            "z": crystal_center[2],
        },
        "thresholds_a": thresholds,
        "box_size_a": box_size,
        "summary": summary,
        "sites": site_rows,
    }
    return summary, site_rows, evaluation


def csv_value(value: Any) -> Any:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.4f}"
        return ""
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: csv_value(row.get(key, "")) for key in fieldnames})


def write_run_outputs(
    *,
    output_dir: Path,
    summary_rows: list[dict[str, Any]],
    site_rows: list[dict[str, Any]],
    thresholds: list[float],
    box_size: float,
    metadata: dict[str, Any],
) -> None:
    box_label = threshold_label(box_size)
    threshold_labels = [threshold_label(value) for value in thresholds]
    summary_fields = [
        "target",
        "status",
        "site_policy",
        "autosites",
        "n_sites",
        "crystal_center_x",
        "crystal_center_y",
        "crystal_center_z",
        "best_site_id",
        "best_site_order",
        "best_distance_a",
        "best_required_cube_side_a",
        f"best_within_{box_label}_cube",
    ]
    for label in threshold_labels:
        summary_fields.append(f"any_site_success_at_{label}")
    summary_fields.extend(
        ["reference_ligand", "receptor_input", "receptor_pdbqt", "centers_tsv", "error"]
    )

    site_fields = [
        "target",
        "site_id",
        "site_order",
        "distance_a",
        "center_x",
        "center_y",
        "center_z",
        "delta_x",
        "delta_y",
        "delta_z",
        "required_cube_side_a",
        f"within_{box_label}_cube",
    ]
    site_fields.extend(f"hit_at_{label}" for label in threshold_labels)

    write_csv(output_dir / "summary.csv", summary_rows, summary_fields)
    write_csv(output_dir / "sites.csv", site_rows, site_fields)

    summary_json = {
        "metadata": metadata,
        "n_targets": len(summary_rows),
        "n_ok": sum(1 for row in summary_rows if row.get("status") == "ok"),
        "n_failed": sum(1 for row in summary_rows if row.get("status") != "ok"),
        "targets": summary_rows,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary_json, indent=2), encoding="utf-8")


def failure_row(
    *,
    target_name: str,
    args: argparse.Namespace,
    error: str,
    reference_ligand: Path | None = None,
    receptor_input: Path | None = None,
    receptor_pdbqt: Path | None = None,
    centers_tsv: Path | None = None,
) -> dict[str, Any]:
    return {
        "target": target_name,
        "status": "failed",
        "site_policy": args.site_policy,
        "autosites": args.autosites,
        "reference_ligand": str(reference_ligand.resolve()) if reference_ligand else "",
        "receptor_input": str(receptor_input.resolve()) if receptor_input else "",
        "receptor_pdbqt": str(receptor_pdbqt.resolve()) if receptor_pdbqt else "",
        "centers_tsv": str(centers_tsv.resolve()) if centers_tsv else "",
        "error": error,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate/evaluate Ultidock auto-sites against DUD-E crystal-ligand centers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dataset-root",
        default=str(DEFAULT_DATASET_ROOT),
        help="Root directory containing DUD-E target subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Directory for generated centers and CSV/JSON outputs.",
    )
    parser.add_argument(
        "--targets",
        help="Comma-separated target names, or 'all'. Default: all discoverable targets.",
    )
    parser.add_argument("--max-targets", type=int, help="Run only the first N selected targets.")
    parser.add_argument(
        "--jobs",
        type=int,
        default=_default_benchmark_jobs(),
        help=(
            "Number of targets to process concurrently. "
            "Can also be set with ULTIDOCK_BENCHMARK_WORKERS."
        ),
    )
    parser.add_argument("--seed", type=int, default=42, help="Deterministic receptor-prep seed.")
    parser.add_argument("--autosites", type=int, default=5, help="Number of auto-sites requested.")
    parser.add_argument(
        "--site-policy",
        choices=SITE_POLICIES,
        default="receptor_search",
        help="Ultidock center-generation policy.",
    )
    parser.add_argument(
        "--thresholds",
        default=",".join(f"{value:g}" for value in DEFAULT_THRESHOLDS),
        help="Comma-separated distance thresholds in Angstrom.",
    )
    parser.add_argument(
        "--box-size",
        type=float,
        default=40.0,
        help="Cubic box side used only for post-hoc containment reporting.",
    )
    parser.add_argument("--grid-spacing", type=float, default=0.375, help="AutoGrid map spacing.")
    parser.add_argument("--blind-cap", type=float, default=150.0, help="Whole-protein map box cap.")
    parser.add_argument(
        "--hotspot-box-size",
        type=float,
        default=35.0,
        help="Generated site box side.",
    )
    parser.add_argument(
        "--tau-rel",
        type=float,
        default=0.52,
        help="Surface hotspot relative threshold.",
    )
    parser.add_argument(
        "--min-sep",
        type=float,
        default=2.0,
        help="Minimum site separation parameter.",
    )
    parser.add_argument(
        "--r-min",
        type=float,
        default=None,
        help=(
            "Minimum inscribed-sphere radius (Å) for cavity acceptance. "
            "Default: None = adaptive per receptor (mainline pipeline behaviour). "
            "Supply a positive float to force a fixed threshold."
        ),
    )
    parser.add_argument(
        "--autogrid4-bin",
        default=None,
        help="Path to autogrid4 executable. Defaults to <autodock-gpu-dir>/autogrid/autogrid4.",
    )
    parser.add_argument(
        "--autodock-gpu-dir",
        default=str(DOCKING_ROOT / "AUTODOCK_GPU_DIR"),
        help="AutoDock-GPU directory used by make_grids config.",
    )
    parser.add_argument(
        "--receptor-prepare-command",
        help="Optional receptor conversion command template with {input} and {output}.",
    )
    parser.add_argument(
        "--centers-root",
        type=Path,
        help="Optional root containing pre-existing target/centers.tsv files.",
    )
    parser.add_argument(
        "--skip-generation",
        action="store_true",
        help="Only evaluate existing centers.tsv files; do not prepare receptors or run AutoGrid.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate receptor PDBQT and centers.tsv even when cached outputs exist.",
    )
    parser.add_argument(
        "--no-auto-setup",
        action="store_true",
        help="Fail instead of running minimal setup.py when config.py or AutoGrid is missing.",
    )
    parser.add_argument(
        "--auto-setup-mode",
        choices=["cpu", "cuda", "opencl", "gpu", "auto"],
        default="cpu",
        help="setup.py mode used by automatic benchmark preflight.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print selected targets and exit.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    normalize_runtime_paths(args)
    dataset_root = Path(args.dataset_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = parse_thresholds(args.thresholds)

    target_names = parse_targets(args.targets, dataset_root)
    if args.max_targets is not None:
        target_names = target_names[: args.max_targets]
    if not target_names:
        raise SystemExit(f"No targets found in {dataset_root}")
    args.jobs = min(args.jobs, len(target_names))

    print(f"Dataset root: {dataset_root}")
    print(f"Output dir:   {output_dir}")
    print(f"Targets:      {len(target_names)}")
    print(f"Jobs:         {args.jobs}")
    for name in target_names:
        print(f"  - {name}")
    if args.dry_run:
        return

    metadata = {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "dataset_root": str(dataset_root),
        "output_dir": str(output_dir),
        "targets": target_names,
        "site_policy": args.site_policy,
        "autosites": args.autosites,
        "jobs": args.jobs,
        "thresholds_a": thresholds,
        "box_size_a": args.box_size,
        "git_commit": git_value("rev-parse", "HEAD"),
        "git_branch": git_value("branch", "--show-current"),
        "git_status_short": git_status_short(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    make_grids = None
    if not args.skip_generation:
        ensure_generation_runtime(args)
        with tempfile.TemporaryDirectory(prefix="ultidock-cavity-config-") as config_tmp:
            config_dir = Path(config_tmp)
            write_runtime_config(config_dir, args)
            make_grids = load_make_grids(config_dir)
            run_targets(
                args,
                dataset_root,
                output_dir,
                thresholds,
                target_names,
                metadata,
                make_grids,
            )
    else:
        run_targets(args, dataset_root, output_dir, thresholds, target_names, metadata, make_grids)


def run_targets(
    args: argparse.Namespace,
    dataset_root: Path,
    output_dir: Path,
    thresholds: list[float],
    target_names: list[str],
    metadata: dict[str, Any],
    make_grids,
) -> None:
    centers_root = Path(args.centers_root).resolve() if args.centers_root else None
    worker_count = min(max(1, int(args.jobs)), len(target_names))
    results_by_index: dict[int, dict[str, Any]] = {}

    if worker_count > 1:
        print(f"\n[parallel] using {worker_count} target worker(s)")

    def persist_completed_outputs() -> None:
        summary_rows, all_site_rows = ordered_benchmark_rows(results_by_index)
        write_run_outputs(
            output_dir=output_dir,
            summary_rows=summary_rows,
            site_rows=all_site_rows,
            thresholds=thresholds,
            box_size=float(args.box_size),
            metadata=metadata,
        )

    if worker_count == 1:
        for index, target_name in enumerate(target_names, start=1):
            results_by_index[index] = run_one_target(
                index=index,
                total=len(target_names),
                target_name=target_name,
                args=args,
                dataset_root=dataset_root,
                output_dir=output_dir,
                thresholds=thresholds,
                centers_root=centers_root,
                make_grids=make_grids,
            )
            persist_completed_outputs()
    else:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix="cavity-worker",
        ) as pool:
            futures = {
                pool.submit(
                    run_one_target,
                    index=index,
                    total=len(target_names),
                    target_name=target_name,
                    args=args,
                    dataset_root=dataset_root,
                    output_dir=output_dir,
                    thresholds=thresholds,
                    centers_root=centers_root,
                    make_grids=make_grids,
                ): index
                for index, target_name in enumerate(target_names, start=1)
            }
            for future in concurrent.futures.as_completed(futures):
                index = futures[future]
                try:
                    results_by_index[index] = future.result()
                except Exception as exc:
                    target_name = target_names[index - 1]
                    error = str(exc)
                    print(f"\n[{index}/{len(target_names)}] {target_name}")
                    print(f"  [FAIL] {error}")
                    results_by_index[index] = {
                        "summary": failure_row(
                            target_name=target_name,
                            args=args,
                            error=error,
                        ),
                        "site_rows": [],
                    }
                persist_completed_outputs()
                print(f"[progress] completed {len(results_by_index)}/{len(target_names)} target(s)")

    print(f"\nSummary CSV: {output_dir / 'summary.csv'}")
    print(f"Sites CSV:   {output_dir / 'sites.csv'}")


def ordered_benchmark_rows(
    results_by_index: dict[int, dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    all_site_rows: list[dict[str, Any]] = []
    for index in sorted(results_by_index):
        result = results_by_index[index]
        summary_rows.append(result["summary"])
        all_site_rows.extend(result.get("site_rows", []))
    return summary_rows, all_site_rows


def run_one_target(
    *,
    index: int,
    total: int,
    target_name: str,
    args: argparse.Namespace,
    dataset_root: Path,
    output_dir: Path,
    thresholds: list[float],
    centers_root: Path | None,
    make_grids,
) -> dict[str, Any]:
    print(f"\n[{index}/{total}] {target_name}")
    target_output = output_dir / target_name
    target_output.mkdir(parents=True, exist_ok=True)
    receptor_input: Path | None = None
    reference_ligand: Path | None = None
    receptor_pdbqt: Path | None = None
    centers_tsv: Path | None = None

    try:
        receptor_input, reference_ligand = load_target_files(dataset_root, target_name)
        if args.skip_generation:
            centers_tsv = _first_existing(
                existing_centers_candidates(target_name, output_dir, centers_root)
            )
            if centers_tsv is None:
                raise FileNotFoundError(f"{target_name}: no existing centers.tsv found")
        else:
            if make_grids is None:
                raise RuntimeError("make_grids module was not loaded")
            receptor_pdbqt = prepare_receptor(
                receptor_input=receptor_input,
                output_dir=target_output,
                prepare_command=args.receptor_prepare_command,
                seed=args.seed,
                force=args.force,
            )
            centers_tsv = generate_centers(
                make_grids=make_grids,
                receptor_pdbqt=receptor_pdbqt,
                target_output=target_output,
                args=args,
            )

        summary, site_rows, evaluation = evaluate_centers(
            target_name=target_name,
            reference_ligand=reference_ligand,
            centers_tsv=centers_tsv,
            thresholds=thresholds,
            box_size=float(args.box_size),
        )
        summary["site_policy"] = args.site_policy
        summary["autosites"] = args.autosites
        summary["receptor_input"] = str(receptor_input.resolve())
        summary["receptor_pdbqt"] = str(receptor_pdbqt.resolve()) if receptor_pdbqt else ""
        (target_output / "evaluation.json").write_text(
            json.dumps(evaluation, indent=2),
            encoding="utf-8",
        )
        print(
            "  [OK] "
            f"closest={float(summary['best_distance_a']):.2f} A "
            f"({summary['best_site_id']})"
        )
        return {"summary": summary, "site_rows": site_rows}
    except Exception as exc:
        error = str(exc)
        print(f"  [FAIL] {error}")
        (target_output / "failure.json").write_text(
            json.dumps(
                {
                    "target": target_name,
                    "error": error,
                    "traceback": traceback.format_exc(),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return {
            "summary": failure_row(
                target_name=target_name,
                args=args,
                error=error,
                reference_ligand=reference_ligand,
                receptor_input=receptor_input,
                receptor_pdbqt=receptor_pdbqt,
                centers_tsv=centers_tsv,
            ),
            "site_rows": [],
        }


if __name__ == "__main__":
    main()
