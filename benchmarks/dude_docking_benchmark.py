#!/usr/bin/env python3
"""Run a DUD-E benchmark through setup.py + dock_v02.py."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib
import json
import math
import re
import shutil
import sqlite3
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from benchmarks.dude_prep import prepare_ligand_inputs
except ModuleNotFoundError:
    from dude_prep import prepare_ligand_inputs

from molguard.io.receptor_prep import (
    SUPPORTED_RECEPTOR_SUFFIXES,
    prepare_receptor_pdbqt,
    receptor_pdbqt_output_name,
)

try:
    from benchmarks.evaluate_cavities import (
        DEFAULT_HIT_THRESHOLD_A,
        evaluate_sites,
        write_evaluation_csv,
        write_evaluation_json,
    )
except ModuleNotFoundError:
    from evaluate_cavities import (  # type: ignore
        DEFAULT_HIT_THRESHOLD_A,
        evaluate_sites,
        write_evaluation_csv,
        write_evaluation_json,
    )


DOCKING_ROOT = REPO_ROOT / "docking"
DEFAULT_DATASET_ROOT = REPO_ROOT / "benchmarks" / "datasets"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "benchmarks" / "results" / "dude"
DEFAULT_TARGETS = ("cdk2", "ache", "hivpr")
_SITE_ID_RE = re.compile(r"^S?(\d+)$", re.IGNORECASE)


@dataclass(frozen=True)
class DockingBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


@dataclass(frozen=True)
class TargetSpec:
    name: str
    receptor: Path
    reference_ligand: Path
    actives: list[Path]
    decoys: list[Path]
    box: DockingBox
    receptor_prepare_command: str | None = None
    ligand_prepare_command: str | None = None


@dataclass(frozen=True)
class LigandScore:
    ligand_id: str
    label: int
    best_affinity_kcal_mol: float
    best_binding_site: str | None
    prepared_path: Path
    dock_rows: int


def site_sort_key(site_id: str | None) -> tuple[int, str]:
    label = (site_id or "UNKNOWN").strip() or "UNKNOWN"
    match = _SITE_ID_RE.match(label)
    if match:
        return int(match.group(1)), label
    return 10**9, label


def parse_atom_records(path: Path) -> np.ndarray:
    """Read coordinates from PDB/PDBQT or MOL2 atom records."""
    suffixes = [part.lower() for part in path.suffixes]
    suffix = suffixes[-2] if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2 else path.suffix.lower()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    if suffix in {".pdb", ".pdbqt"}:
        coords: list[tuple[float, float, float]] = []
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    coords.append(
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    )
                except ValueError:
                    continue
        return np.asarray(coords, dtype=float)

    if suffix == ".mol2":
        coords = []
        in_atoms = False
        with opener(path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.rstrip("\n")
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if not in_atoms or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                except ValueError:
                    continue
        return np.asarray(coords, dtype=float)

    raise ValueError(f"Unsupported coordinate file for box inference: {path}")


def infer_box_from_reference(
    reference_ligand: Path,
    *,
    padding: float,
    min_size: float,
) -> DockingBox:
    coords = parse_atom_records(reference_ligand)
    if coords.size == 0:
        raise ValueError(f"Could not extract any coordinates from {reference_ligand}")
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = coords.mean(axis=0)
    extents = np.maximum(maxs - mins + (2.0 * padding), min_size)
    return DockingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(extents[0]),
        size_y=float(extents[1]),
        size_z=float(extents[2]),
    )


def _first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def _candidate_ligand_sources(target_dir: Path, stem: str) -> list[Path]:
    return [
        target_dir / stem,
        target_dir / f"{stem}.mol2",
        target_dir / f"{stem}.mol2.gz",
        target_dir / f"{stem}.pdbqt",
        target_dir / f"{stem}.pdbqt.gz",
        target_dir / f"{stem}.pdb",
        target_dir / f"{stem}.pdb.gz",
        target_dir / f"{stem}_final.mol2",
        target_dir / f"{stem}_final.mol2.gz",
        target_dir / f"{stem}_final.pdbqt",
        target_dir / f"{stem}_final.pdbqt.gz",
        target_dir / f"{stem}_final",
    ]


def _resolve_box(target_dir: Path, metadata: dict[str, object], args: argparse.Namespace) -> DockingBox:
    box_data = metadata.get("box")
    if isinstance(box_data, dict):
        required = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        missing = [key for key in required if key not in box_data]
        if missing:
            raise ValueError(f"{target_dir}: benchmark.json box is missing {missing}")
        return DockingBox(**{key: float(box_data[key]) for key in required})

    reference_candidates: list[Path] = []
    if metadata.get("reference_ligand"):
        reference_candidates.append(target_dir / str(metadata["reference_ligand"]))
    reference_candidates.extend(
        target_dir / name
        for name in (
            "crystal_ligand.pdbqt",
            "crystal_ligand.pdbqt.gz",
            "crystal_ligand.mol2",
            "crystal_ligand.mol2.gz",
            "crystal_ligand.pdb",
            "crystal_ligand.pdb.gz",
            "ligand.pdbqt",
            "ligand.pdbqt.gz",
            "ligand.mol2",
            "ligand.mol2.gz",
            "ligand.pdb",
            "ligand.pdb.gz",
            "ref_ligand.pdbqt",
            "ref_ligand.pdbqt.gz",
            "ref_ligand.mol2",
            "ref_ligand.mol2.gz",
            "ref_ligand.pdb",
            "ref_ligand.pdb.gz",
        )
    )
    reference = _first_existing(reference_candidates)
    if reference is None:
        raise FileNotFoundError(
            f"{target_dir}: no box provided and no reference ligand found for box inference"
        )
    return infer_box_from_reference(
        reference,
        padding=float(metadata.get("box_padding", args.box_padding)),
        min_size=float(metadata.get("box_min_size", args.box_min_size)),
    )


def load_target_spec(dataset_root: Path, target_name: str, args: argparse.Namespace) -> TargetSpec:
    target_dir = (dataset_root / target_name).resolve()
    if not target_dir.exists():
        raise FileNotFoundError(f"Target directory not found: {target_dir}")

    metadata_path = target_dir / "benchmark.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}

    receptor_candidates: list[Path] = []
    if metadata.get("receptor"):
        receptor_candidates.append(target_dir / str(metadata["receptor"]))
    receptor_candidates.extend(
        target_dir / f"receptor{suffix}" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    receptor_candidates.extend(
        target_dir / f"receptor{suffix}.gz" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    receptor_candidates.extend(
        target_dir / f"protein{suffix}" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    receptor_candidates.extend(
        target_dir / f"protein{suffix}.gz" for suffix in sorted(SUPPORTED_RECEPTOR_SUFFIXES)
    )
    receptor = _first_existing(receptor_candidates)
    if receptor is None:
        raise FileNotFoundError(f"{target_dir}: could not locate receptor input")

    reference_candidates: list[Path] = []
    if metadata.get("reference_ligand"):
        reference_candidates.append(target_dir / str(metadata["reference_ligand"]))
    reference_candidates.extend(
        target_dir / name
        for name in (
            "crystal_ligand.pdbqt",
            "crystal_ligand.pdbqt.gz",
            "crystal_ligand.mol2",
            "crystal_ligand.mol2.gz",
            "crystal_ligand.pdb",
            "crystal_ligand.pdb.gz",
            "ligand.pdbqt",
            "ligand.pdbqt.gz",
            "ligand.mol2",
            "ligand.mol2.gz",
            "ligand.pdb",
            "ligand.pdb.gz",
            "ref_ligand.pdbqt",
            "ref_ligand.pdbqt.gz",
            "ref_ligand.mol2",
            "ref_ligand.mol2.gz",
            "ref_ligand.pdb",
            "ref_ligand.pdb.gz",
        )
    )
    reference_ligand = _first_existing(reference_candidates)
    if reference_ligand is None:
        raise FileNotFoundError(f"{target_dir}: could not locate a reference ligand")

    active_candidates: list[Path] = []
    if metadata.get("actives"):
        active_candidates.append(target_dir / str(metadata["actives"]))
    active_candidates.extend(_candidate_ligand_sources(target_dir, "actives"))
    actives = _first_existing(active_candidates)

    decoy_candidates: list[Path] = []
    if metadata.get("decoys"):
        decoy_candidates.append(target_dir / str(metadata["decoys"]))
    decoy_candidates.extend(_candidate_ligand_sources(target_dir, "decoys"))
    decoys = _first_existing(decoy_candidates)
    if actives is None or decoys is None:
        raise FileNotFoundError(
            f"{target_dir}: expected actives/decoys inputs or benchmark.json overrides"
        )

    return TargetSpec(
        name=target_name,
        receptor=receptor,
        reference_ligand=reference_ligand,
        actives=[actives],
        decoys=[decoys],
        box=_resolve_box(target_dir, metadata, args),
        receptor_prepare_command=str(
            metadata.get("receptor_prepare_command") or args.receptor_prepare_command or ""
        )
        or None,
        ligand_prepare_command=str(
            metadata.get("ligand_prepare_command") or args.ligand_prepare_command or ""
        )
        or None,
    )


def average_ranks(scores: np.ndarray, *, descending: bool = False) -> np.ndarray:
    """Average ranks with stable tie handling."""
    order = np.argsort(-scores if descending else scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    sorted_scores = scores[order]
    start = 0
    while start < len(sorted_scores):
        end = start + 1
        while end < len(sorted_scores) and sorted_scores[end] == sorted_scores[start]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        ranks[order[start:end]] = avg_rank
        start = end
    return ranks


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("ROC-AUC requires at least one positive and one negative example")
    ranks = average_ranks(scores, descending=False)
    pos_ranks = ranks[labels == 1]
    u = pos_ranks.sum() - positives * (positives + 1) / 2.0
    return float(u / (positives * negatives))


def enrichment_factor(labels: np.ndarray, scores: np.ndarray, fraction: float) -> float:
    positives = int(labels.sum())
    if positives == 0:
        raise ValueError("Enrichment factor requires at least one active ligand")
    n_total = len(labels)
    top_k = max(1, math.ceil(n_total * fraction))
    order = np.argsort(-scores, kind="mergesort")
    hits = int(labels[order[:top_k]].sum())
    expected = positives * (top_k / n_total)
    if expected == 0:
        return 0.0
    return float(hits / expected)


def bedroc_like(labels: np.ndarray, scores: np.ndarray, alpha: float) -> float:
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("BEDROC requires at least one active and one decoy")
    if alpha <= 0:
        raise ValueError("BEDROC alpha must be > 0")

    order = np.argsort(-scores, kind="mergesort")
    ranks = np.nonzero(labels[order] == 1)[0] + 1
    n_total = len(labels)

    observed = float(np.exp(-alpha * ranks / n_total).mean())
    baseline = ((1.0 - math.exp(-alpha)) / n_total) / (math.exp(alpha / n_total) - 1.0)
    rie = observed / baseline

    best_ranks = np.arange(1, positives + 1, dtype=float)
    worst_ranks = np.arange(n_total - positives + 1, n_total + 1, dtype=float)
    rie_max = float(np.exp(-alpha * best_ranks / n_total).mean() / baseline)
    rie_min = float(np.exp(-alpha * worst_ranks / n_total).mean() / baseline)
    if math.isclose(rie_max, rie_min):
        return 1.0
    return float((rie - rie_min) / (rie_max - rie_min))


def log_auc(labels: np.ndarray, scores: np.ndarray, lambda_min: float = 0.001) -> float:
    """LogAUC: area under the ROC curve with a log-scaled FPR axis.

    Emphasises early enrichment (low FPR region).  A random classifier
    scores ~0.1446 when *lambda_min* = 0.001.

    Reference
    ---------
    Mysinger, M. M. & Shoichet, B. K. (2010). *Rapid context-dependent
    ligand desolvation in molecular docking*. JCIM 50(9), 1561-1573.
    """
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("LogAUC requires at least one positive and one negative example")

    order = np.argsort(-scores, kind="mergesort")
    sorted_labels = labels[order]

    tps = np.cumsum(sorted_labels)
    fps = np.cumsum(1 - sorted_labels)
    tpr = np.concatenate([[0.0], tps / positives])
    fpr = np.concatenate([[0.0], fps / negatives])

    valid = fpr >= lambda_min
    fpr_v = fpr[valid]
    tpr_v = tpr[valid]
    if len(fpr_v) < 2:
        return 0.0
    log_fpr = np.log10(fpr_v)
    integral = float(np.trapz(tpr_v, log_fpr))
    denom = -np.log10(lambda_min)  # 3.0 for lambda_min = 0.001
    return integral / denom


def write_centers_tsv(receptor_pdbqt: Path, box: DockingBox) -> Path:
    centers_dir = receptor_pdbqt.parent / receptor_pdbqt.stem
    centers_dir.mkdir(parents=True, exist_ok=True)
    centers_path = centers_dir / "centers.tsv"
    rec_stem = receptor_pdbqt.stem
    # Format: receptor_name  site_id  cx  cy  cz  (matches _parse_centers_tsv)
    centers_path.write_text(
        (
            f"{rec_stem}\tS1\t{box.center_x:.3f}\t{box.center_y:.3f}\t{box.center_z:.3f}\t"
            f"{box.size_x:.3f}\t{box.size_y:.3f}\t{box.size_z:.3f}\n"
        ),
        encoding="utf-8",
    )
    return centers_path


def clean_docking_workspace(*, dry: bool = False) -> None:
    """Clean the standard docking/ workspace using clean.py."""
    cmd = ["python3", "clean.py", "--all"]
    if not dry:
        cmd.append("-y")
    subprocess.run(
        cmd,
        cwd=DOCKING_ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )


def stage_target_inputs(
    *,
    target: TargetSpec,
    seed: int,
    center_source: str,
) -> tuple[dict[str, Path], dict[Path, int]]:
    """Prepare receptor and ligands in the standard docking/ workspace."""

    # Standard paths under docking/
    paths = {
        "root": DOCKING_ROOT,
        "ligands": DOCKING_ROOT / "LIGANDS_DIR",
        "docking": DOCKING_ROOT / "DOCKING_DIR",
        "analysis": DOCKING_ROOT / "ANALYSIS_DIR",
        "macro": DOCKING_ROOT / "MACRO_MOL_DIR",
        "results": DOCKING_ROOT / "RESULTS_DIR",
    }

    # Clean workspace between targets
    clean_docking_workspace()

    # Ensure standard directories exist
    for key in ("ligands", "docking", "analysis", "macro", "results"):
        paths[key].mkdir(parents=True, exist_ok=True)

    receptor_out = paths["macro"] / receptor_pdbqt_output_name(target.receptor)
    prepare_receptor_pdbqt(
        input_path=target.receptor,
        output_path=receptor_out,
        prepare_command=target.receptor_prepare_command,
        seed=seed,
        timestamp="BENCHMARK",
    )
    if center_source == "crystal":
        write_centers_tsv(receptor_out, target.box)

    # Prepare ligands --> LIGANDS_DIR/active__*.pdbqt + decoy__*.pdbqt
    active_paths = prepare_ligand_inputs(
        inputs=target.actives,
        output_dir=paths["ligands"],
        prepare_command=target.ligand_prepare_command,
        seed=seed,
        output_prefix="active__",
        manifest_path=DOCKING_ROOT / "actives_manifest.json",
    )
    decoy_paths = prepare_ligand_inputs(
        inputs=target.decoys,
        output_dir=paths["ligands"],
        prepare_command=target.ligand_prepare_command,
        seed=seed,
        output_prefix="decoy__",
        manifest_path=DOCKING_ROOT / "decoys_manifest.json",
    )

    label_map: dict[Path, int] = {}
    for ligand in active_paths:
        label_map[ligand.resolve()] = 1
    for ligand in decoy_paths:
        label_map[ligand.resolve()] = 0
    return paths, label_map


def _run_with_tee(
    cmd: list[str],
    log_path: Path,
    *,
    cwd: Path | None = None,
) -> None:
    """Run a command, printing output live to terminal AND saving to a log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_fh:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,  # line buffered
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_fh.write(line)
        proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


def reset_sqlite_outputs(db_path: Path) -> None:
    """Remove a SQLite database and any WAL sidecars from earlier runs."""
    for path in (
        db_path,
        Path(f"{db_path}-wal"),
        Path(f"{db_path}-shm"),
    ):
        if path.exists():
            path.unlink()


def run_setup_for_workspace(
    *,
    paths: dict[str, Path],
    args: argparse.Namespace,
    target_output: Path,
) -> Path:
    """Run setup.py using the standard docking/ workspace."""
    autosites = 1 if args.center_source == "crystal" else int(args.autosites)
    db_path = paths["results"] / "ultidock_results.db"
    reset_sqlite_outputs(db_path)

    setup_cmd = [
        "python3",
        "setup.py",
        "--benchmark",
        "--skip-wget",
        "--mode",
        args.mode,
        "--grid-mode",
        "centers",
        "--autosites",
        str(autosites),
        "--vina-cpu",
        str(args.vina_cpu),
        "--vina-seed",
        str(args.seed),
        "--vina-exhaustiveness",
        str(args.vina_exhaustiveness),
        "--vina-num-modes",
        str(args.vina_num_modes),
    ]
    setup_log = target_output / "setup.log"
    _run_with_tee(setup_cmd, setup_log, cwd=DOCKING_ROOT)
    return setup_log


def _load_dock_v02_module():
    docking_root_str = str(DOCKING_ROOT)
    if docking_root_str not in sys.path:
        sys.path.insert(0, docking_root_str)
    importlib.invalidate_caches()
    config_mod = importlib.import_module("config")
    importlib.reload(config_mod)
    dock_mod = importlib.import_module("dock_v02")
    return importlib.reload(dock_mod)


def find_prepared_receptor_pdbqt(macro_dir: Path) -> Path:
    receptors = sorted(macro_dir.glob("*.pdbqt"))
    if not receptors:
        raise FileNotFoundError(f"No prepared receptor .pdbqt found in {macro_dir}")
    if len(receptors) > 1:
        names = ", ".join(path.name for path in receptors)
        raise RuntimeError(
            "Benchmark workspace expected one prepared receptor .pdbqt, "
            f"found {len(receptors)}: {names}"
        )
    return receptors[0]


def generate_auto_sites(paths: dict[str, Path]) -> Path:
    """Force site generation before docking so it can be benchmarked."""
    dock_mod = _load_dock_v02_module()
    receptor_pdbqt = find_prepared_receptor_pdbqt(paths["macro"])
    site_root = paths["macro"] / receptor_pdbqt.stem
    dock_mod.prepare_sites_for_docking(str(receptor_pdbqt), str(site_root))
    centers_tsv = site_root / "centers.tsv"
    if not centers_tsv.exists():
        raise FileNotFoundError(f"Auto-site generation did not produce {centers_tsv}")
    return centers_tsv


def evaluate_auto_sites(
    *,
    target: TargetSpec,
    centers_tsv_path: Path,
    target_output: Path,
    hit_threshold_a: float,
) -> tuple[dict[str, object], Path, Path]:
    """Write per-target cavity-finder evaluation artifacts."""
    evaluation = evaluate_sites(
        reference_ligand_path=target.reference_ligand,
        centers_tsv_path=centers_tsv_path,
        hit_threshold_a=hit_threshold_a,
    )
    json_path = target_output / "site_summary.json"
    csv_path = target_output / "site_summary.csv"
    write_evaluation_json(json_path, evaluation)
    write_evaluation_csv(csv_path, evaluation)
    print(
        "  [sites] "
        f"best={evaluation['best_site_id']} "
        f"site_order={evaluation['best_site_order']} "
        f"dist={float(evaluation['best_distance_a']):.2f} A "
        f"any_site={'YES' if evaluation['any_site_success'] else 'NO'}"
    )
    return evaluation, json_path, csv_path


def run_docking_for_workspace(
    *,
    target: TargetSpec,
    paths: dict[str, Path],
    target_output: Path,
) -> Path:
    """Run dock_v02.py using the standard docking/ workspace."""
    db_path = paths["results"] / "ultidock_results.db"
    dock_log = target_output / "dock_v02.log"
    _run_with_tee(["python3", "dock_v02.py"], dock_log, cwd=DOCKING_ROOT)
    if not db_path.exists():
        raise FileNotFoundError(f"{target.name}: expected benchmark DB at {db_path}")
    return db_path


def load_docking_rows(
    db_path: Path,
) -> tuple[dict[Path, list[tuple[float, str | None]]], dict[str, dict[Path, list[tuple[float, str | None]]]]]:
    rows_by_file: dict[Path, list[tuple[float, str | None]]] = {}
    rows_by_site: dict[str, dict[Path, list[tuple[float, str | None]]]] = {}
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        rows = cur.execute(
            'SELECT docking_file, "binding_affinity (kcal/mol)", binding_site FROM docking_results'
        ).fetchall()

    for docking_file, affinity, binding_site in rows:
        if docking_file is None or affinity is None:
            continue
        key = Path(str(docking_file)).resolve()
        affinity_value = float(affinity)
        site_id = (binding_site or "UNKNOWN").strip() or "UNKNOWN"
        rows_by_file.setdefault(key, []).append((affinity_value, binding_site))
        rows_by_site.setdefault(site_id, {}).setdefault(key, []).append((affinity_value, binding_site))
    return rows_by_file, rows_by_site


def scores_from_grouped_rows(
    rows_by_file: dict[Path, list[tuple[float, str | None]]],
    label_map: dict[Path, int],
    *,
    warn_context: str,
) -> list[LigandScore]:
    missing = sorted(str(path) for path in label_map if path not in rows_by_file)
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "" if len(missing) <= 5 else f" ... (+{len(missing) - 5} more)"
        pct = len(missing) / len(label_map) * 100 if label_map else 0
        print(
            f"  [WARN] Missing docking results for {warn_context}: "
            f"{len(missing)}/{len(label_map)} ligand(s) ({pct:.1f}%): {preview}{suffix}"
        )

    scores: list[LigandScore] = []
    for prepared_path, label in sorted(label_map.items(), key=lambda item: item[0].name):
        if prepared_path not in rows_by_file:
            continue  # skip missing ligands
        entries = rows_by_file[prepared_path]
        best_affinity, best_site = min(entries, key=lambda item: item[0])
        scores.append(
            LigandScore(
                ligand_id=prepared_path.stem,
                label=label,
                best_affinity_kcal_mol=best_affinity,
                best_binding_site=best_site,
                prepared_path=prepared_path,
                dock_rows=len(entries),
            )
        )
    return scores


def collect_scores_from_db(db_path: Path, label_map: dict[Path, int]) -> list[LigandScore]:
    rows_by_file, _ = load_docking_rows(db_path)
    return scores_from_grouped_rows(rows_by_file, label_map, warn_context="all sites")


def collect_scores_from_db_by_site(db_path: Path, label_map: dict[Path, int]) -> dict[str, list[LigandScore]]:
    _, rows_by_site = load_docking_rows(db_path)
    site_scores: dict[str, list[LigandScore]] = {}
    for site_id, grouped_rows in sorted(rows_by_site.items(), key=lambda item: site_sort_key(item[0])):
        site_scores[site_id] = scores_from_grouped_rows(
            grouped_rows,
            label_map,
            warn_context=f"site {site_id}",
        )
    return site_scores


def write_scores_csv(path: Path, rows: Iterable[LigandScore]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["ligand_id", "label", "best_affinity_kcal_mol", "best_binding_site", "dock_rows", "prepared_path"]
        )
        for row in rows:
            writer.writerow(
                [
                    row.ligand_id,
                    row.label,
                    f"{row.best_affinity_kcal_mol:.4f}",
                    row.best_binding_site or "",
                    row.dock_rows,
                    str(row.prepared_path),
                ]
            )


def build_metrics_from_scores(
    *,
    target: TargetSpec,
    scores: list[LigandScore],
    args: argparse.Namespace,
    db_path: Path,
    scores_csv: Path,
    expected_ligands: int,
    site_aggregation: str,
    site_id: str | None = None,
) -> dict[str, object]:
    labels = np.asarray([row.label for row in scores], dtype=int)
    affinities = np.asarray([row.best_affinity_kcal_mol for row in scores], dtype=float)

    metrics: dict[str, object] = {
        "target": target.name,
        "mode": args.mode,
        "center_source": args.center_source,
        "site_aggregation": site_aggregation,
        "site_id": site_id,
        "n_expected_ligands": expected_ligands,
        "n_scored_ligands": len(scores),
        "n_missing_ligands": max(0, expected_ligands - len(scores)),
        "n_actives": int(labels.sum()) if len(labels) else 0,
        "n_decoys": int(len(labels) - labels.sum()) if len(labels) else 0,
        "box": asdict(target.box),
        "db_path": str(db_path.resolve()),
        "scores_csv": str(scores_csv.resolve()),
    }

    if not len(scores):
        return metrics

    bedroc_alpha = getattr(args, "bedroc_alpha", None) or 20.0
    ranking_scores = -affinities
    metrics["bedroc_alpha"] = bedroc_alpha
    if metrics["n_actives"] and metrics["n_decoys"]:
        metrics["roc_auc"] = roc_auc(labels, ranking_scores)
        metrics["ef1"] = enrichment_factor(labels, ranking_scores, 0.01)
        metrics["ef5"] = enrichment_factor(labels, ranking_scores, 0.05)
        metrics["bedroc"] = bedroc_like(labels, ranking_scores, bedroc_alpha)
        metrics["log_auc"] = log_auc(labels, ranking_scores)
    return metrics


def write_single_target_summary(path: Path, metrics: dict[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"targets": [metrics]}, indent=2), encoding="utf-8")
    return path


def write_per_site_outputs(
    *,
    target: TargetSpec,
    target_output: Path,
    args: argparse.Namespace,
    db_path: Path,
    label_map: dict[Path, int],
) -> tuple[list[dict[str, object]], Path | None]:
    per_site_scores = collect_scores_from_db_by_site(db_path, label_map)
    if not per_site_scores:
        return [], None

    sites_root = target_output / "sites"
    site_metrics: list[dict[str, object]] = []
    for site_id, scores in sorted(per_site_scores.items(), key=lambda item: site_sort_key(item[0])):
        site_dir = sites_root / site_id
        scores_csv = site_dir / "scores.csv"
        write_scores_csv(scores_csv, scores)
        metrics = build_metrics_from_scores(
            target=target,
            scores=scores,
            args=args,
            db_path=db_path,
            scores_csv=scores_csv,
            expected_ligands=len(label_map),
            site_aggregation="single_site",
            site_id=site_id,
        )
        summary_json = write_single_target_summary(site_dir / "summary.json", metrics)
        metrics["summary_json"] = str(summary_json.resolve())
        metrics["result_dir"] = str(site_dir.resolve())
        site_metrics.append(metrics)

    manifest_path = target_output / "site_metrics.json"
    manifest_path.write_text(json.dumps(site_metrics, indent=2), encoding="utf-8")
    return site_metrics, manifest_path


def score_target(target: TargetSpec, target_output: Path, args: argparse.Namespace) -> dict[str, object]:
    paths, label_map = stage_target_inputs(
        target=target,
        seed=args.seed,
        center_source=args.center_source,
    )
    setup_log = run_setup_for_workspace(paths=paths, args=args, target_output=target_output)

    receptor_pdbqt = find_prepared_receptor_pdbqt(paths["macro"])
    centers_tsv_src = paths["macro"] / receptor_pdbqt.stem / "centers.tsv"
    site_evaluation: dict[str, object] | None = None
    site_summary_json: Path | None = None
    site_summary_csv: Path | None = None

    if args.center_source == "auto":
        centers_tsv_src = generate_auto_sites(paths)
        try:
            site_evaluation, site_summary_json, site_summary_csv = evaluate_auto_sites(
                target=target,
                centers_tsv_path=centers_tsv_src,
                target_output=target_output,
                hit_threshold_a=float(args.site_hit_threshold),
            )
        except Exception as exc:
            print(f"  [WARN] Site evaluation failed: {exc}")

    db_path = run_docking_for_workspace(target=target, paths=paths, target_output=target_output)

    # Back up auto-detected binding site centers for reproducibility
    centers_tsv_dst = target_output / "centers.tsv"
    if centers_tsv_src.exists():
        shutil.copy2(str(centers_tsv_src), str(centers_tsv_dst))
        print(f"  Backed up centers.tsv --> {centers_tsv_dst}")

    scores = collect_scores_from_db(db_path, label_map)
    per_site_metrics, per_site_manifest = write_per_site_outputs(
        target=target,
        target_output=target_output,
        args=args,
        db_path=db_path,
        label_map=label_map,
    )

    scores_csv = target_output / "scores.csv"
    write_scores_csv(scores_csv, scores)
    metrics = build_metrics_from_scores(
        target=target,
        scores=scores,
        args=args,
        db_path=db_path,
        scores_csv=scores_csv,
        expected_ligands=len(label_map),
        site_aggregation="best_of_sites",
    )
    metrics["setup_log"] = str(setup_log.resolve())
    metrics["dock_log"] = str((target_output / "dock_v02.log").resolve())
    metrics["centers_tsv"] = str(centers_tsv_dst.resolve()) if centers_tsv_dst.exists() else None
    if site_evaluation is not None:
        metrics["site_evaluation"] = site_evaluation
    if site_summary_json is not None:
        metrics["site_summary_json"] = str(site_summary_json.resolve())
    if site_summary_csv is not None:
        metrics["site_summary_csv"] = str(site_summary_csv.resolve())
    if per_site_metrics:
        metrics["per_site_metrics"] = per_site_metrics
    if per_site_manifest is not None:
        metrics["per_site_manifest"] = str(per_site_manifest.resolve())
    return metrics


def aggregate_metrics(per_target: list[dict[str, object]]) -> dict[str, float]:
    summary: dict[str, float] = {}
    for key in ("roc_auc", "ef1", "ef5", "bedroc", "log_auc"):
        values = [float(entry[key]) for entry in per_target if key in entry]
        if values:
            summary[f"{key}_mean"] = float(sum(values) / len(values))
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a DUD-E benchmark through setup.py and dock_v02.py.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dataset-root",
        default=str(DEFAULT_DATASET_ROOT),
        help="Root directory containing one subdirectory per benchmark target.",
    )
    parser.add_argument(
        "--target",
        action="append",
        help=f"Target subdirectory to benchmark. Defaults to {', '.join(DEFAULT_TARGETS)}.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Directory where benchmark workspaces, logs, CSVs, and summary JSON will be written.",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "gpu", "cpu", "cuda", "opencl"],
        default="cpu",
        help="Mode forwarded to docking/setup.py.",
    )
    parser.add_argument(
        "--center-source",
        choices=["auto", "crystal"],
        default="auto",
        help="How UltiDock should obtain centers.tsv for this benchmark arm.",
    )
    parser.add_argument(
        "--receptor-prepare-command",
        help="Template used to convert non-PDBQT receptors before canonicalization.",
    )
    parser.add_argument(
        "--ligand-prepare-command",
        help="Template used to convert non-PDBQT ligands before normalization.",
    )
    parser.add_argument("--seed", type=int, default=0, help="Shared deterministic seed.")
    parser.add_argument("--vina-cpu", type=int, default=1, help="Vina CPU count written via setup.py.")
    parser.add_argument(
        "--vina-exhaustiveness",
        type=int,
        default=8,
        help="Vina exhaustiveness written via setup.py.",
    )
    parser.add_argument(
        "--vina-num-modes",
        type=int,
        default=1,
        help="Vina num_modes written via setup.py.",
    )
    parser.add_argument(
        "--autosites",
        type=int,
        default=6,
        help="Number of auto-detected sites to request when --center-source=auto.",
    )
    parser.add_argument(
        "--box-padding",
        type=float,
        default=4.0,
        help="Padding added on each side when inferring the docking box from a reference ligand.",
    )
    parser.add_argument(
        "--box-min-size",
        type=float,
        default=18.0,
        help="Minimum box size per axis when inferring the docking box from a reference ligand.",
    )
    parser.add_argument(
        "--bedroc-alpha",
        type=float,
        help="If provided, also compute a BEDROC-style early enrichment score.",
    )
    parser.add_argument(
        "--site-hit-threshold",
        type=float,
        default=DEFAULT_HIT_THRESHOLD_A,
        help="Distance threshold in Angstrom for cavity-finder success metrics.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    dataset_root = Path(args.dataset_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    target_names = args.target or list(DEFAULT_TARGETS)
    per_target: list[dict[str, object]] = []
    failed: list[tuple[str, str]] = []

    for i, target_name in enumerate(target_names, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(target_names)}] {target_name.upper()}")
        print(f"{'='*60}")
        try:
            target = load_target_spec(dataset_root, target_name, args)
            target_output = output_dir / target_name
            target_output.mkdir(parents=True, exist_ok=True)
            metrics = score_target(target=target, target_output=target_output, args=args)
            per_target.append(metrics)
            print(f"  [OK] {target_name}: ROC-AUC = {metrics.get('roc_auc', '?')}")
        except Exception as exc:
            msg = str(exc)
            print(f"  [FAIL] {target_name}: {msg}")
            failed.append((target_name, msg))

    summary = {
        "dataset_root": str(dataset_root),
        "output_dir": str(output_dir),
        "targets": per_target,
        "aggregate": aggregate_metrics(per_target),
        "seed": args.seed,
        "mode": args.mode,
        "center_source": args.center_source,
        "n_completed": len(per_target),
        "n_failed": len(failed),
        "failed": [{"target": n, "reason": r} for n, r in failed],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nSummary: {summary_path}")

    if failed:
        print(f"\nFailed targets ({len(failed)}):")
        for name, reason in failed:
            print(f"  [FAIL] {name}: {reason}")

    if per_target:
        agg = summary["aggregate"]
        print(f"\nAggregate (over {len(per_target)} targets):")
        for k, v in agg.items():  # type: ignore[union-attr]
            print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")


if __name__ == "__main__":
    main()
