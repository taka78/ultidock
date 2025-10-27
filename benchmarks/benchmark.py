#!/usr/bin/env python3
"""Benchmark the full Ultidock pipeline.

This utility stages a benchmarking dataset inside ``docking/``, runs the
standard ``run.py`` pipeline under resource monitoring, and aggregates
post-run artefacts (SQLite summaries, grid accuracy checks) into a JSON
report.  It is intentionally thin: the docking workflow itself remains in
``docking/`` – this script only orchestrates it for reproducible studies.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sqlite3
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:  # pragma: no cover - optional dependency for richer metrics
    import psutil  # type: ignore
except Exception:  # pragma: no cover - fall back to coarse timing only
    psutil = None  # type: ignore

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_ROOT = REPO_ROOT / "docking"
DEFAULT_DB_PATH = DOCKING_ROOT / "RESULTS_DIR" / "ultidock_results.db"


@dataclass
class LigandSpec:
    """Ligand entry within a benchmarking dataset."""

    name: str
    source: Path
    reference: Optional[Path] = None


@dataclass
class BenchmarkDataset:
    """Description of a benchmarking dataset."""

    name: str
    receptor: Path
    ligands: List[LigandSpec]
    mode: str = "cpu"
    description: Optional[str] = None
    run_args: List[str] = field(default_factory=list)

    @property
    def slug(self) -> str:
        slug = [c.lower() if c.isalnum() else "-" for c in self.name.strip()]
        cleaned = "".join(slug).strip("-")
        return cleaned or "dataset"


@dataclass
class AccuracyRecord:
    ligand: str
    reference_path: Path
    closest_site: str
    min_distance_A: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "ligand": self.ligand,
            "reference_path": str(self.reference_path),
            "closest_site": self.closest_site,
            "min_distance_A": self.min_distance_A,
        }


@dataclass
class BenchmarkResult:
    dataset: BenchmarkDataset
    run_command: List[str]
    returncode: int
    wall_time_s: float
    peak_rss_mb: Optional[float]
    cpu_user_s: Optional[float]
    cpu_system_s: Optional[float]
    db_metrics: Dict[str, object]
    accuracy_records: List[AccuracyRecord]
    log_path: Path

    def to_dict(self) -> Dict[str, object]:
        payload = {
            "dataset": self.dataset.name,
            "mode": self.dataset.mode,
            "run_command": self.run_command,
            "returncode": self.returncode,
            "wall_time_s": self.wall_time_s,
            "peak_rss_mb": self.peak_rss_mb,
            "cpu_user_s": self.cpu_user_s,
            "cpu_system_s": self.cpu_system_s,
            "db_metrics": self.db_metrics,
            "log_path": str(self.log_path),
            "accuracy_records": [rec.to_dict() for rec in self.accuracy_records],
        }
        return payload


def _resolve_path(base: Path, maybe_relative: str) -> Path:
    candidate = Path(maybe_relative)
    if not candidate.is_absolute():
        candidate = (base / candidate).resolve()
    return candidate


def load_dataset(path: Path) -> BenchmarkDataset:
    data = json.loads(path.read_text())
    base = path.parent

    receptor = _resolve_path(base, data["receptor"])
    ligands = []
    for entry in data.get("ligands", []):
        source_key = entry.get("source") or entry.get("path")
        if not source_key:
            raise ValueError(f"Ligand entry missing 'source' in {path}")
        source_path = _resolve_path(base, source_key)
        reference = entry.get("reference")
        ligands.append(
            LigandSpec(
                name=entry.get("name") or Path(source_key).stem,
                source=source_path,
                reference=_resolve_path(base, reference) if reference else None,
            )
        )

    dataset = BenchmarkDataset(
        name=data["name"],
        receptor=receptor,
        ligands=ligands,
        mode=data.get("mode", "cpu"),
        description=data.get("description"),
        run_args=list(data.get("run_args", [])),
    )
    return dataset


def stage_workspace(dataset: BenchmarkDataset) -> Dict[str, Path]:
    """Clean docking/ and copy the dataset's receptor/ligands into place."""

    subprocess.run(["python3", "clean.py", "-y", "--all"], check=True, cwd=DOCKING_ROOT)

    paths = {
        "macro": DOCKING_ROOT / "MACRO_MOL_DIR",
        "ligands": DOCKING_ROOT / "LIGANDS_DIR",
        "results": DOCKING_ROOT / "RESULTS_DIR",
        "analysis": DOCKING_ROOT / "ANALYSIS_DIR",
        "autodock": DOCKING_ROOT / "AUTODOCK_GPU_DIR",
        "vina": DOCKING_ROOT / "VINA_DIR",
        "docking": DOCKING_ROOT / "DOCKING_DIR",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)

    receptor_target = paths["macro"] / dataset.receptor.name
    shutil.copy2(dataset.receptor, receptor_target)

    for ligand in dataset.ligands:
        shutil.copy2(ligand.source, paths["ligands"] / Path(ligand.source).name)

    # Provide expected bin directories so setup can reuse compiled binaries if present.
    (paths["autodock"] / "bin").mkdir(parents=True, exist_ok=True)
    (paths["autodock"] / "autogrid").mkdir(parents=True, exist_ok=True)
    (paths["vina"] / "bin").mkdir(parents=True, exist_ok=True)

    paths["receptor_target"] = receptor_target
    return paths


def _gather_process_tree_memory(proc: "psutil.Process") -> float:
    rss = 0.0
    try:
        rss += proc.memory_info().rss
    except psutil.Error:  # pragma: no cover - psutil may race with exit
        return 0.0
    for child in proc.children(recursive=True):
        try:
            rss += child.memory_info().rss
        except psutil.Error:
            continue
    return rss


def _collect_cpu_times(proc: "psutil.Process") -> tuple[Optional[float], Optional[float]]:
    user = 0.0
    system = 0.0
    try:
        cpu = proc.cpu_times()
        user += cpu.user
        system += cpu.system
    except psutil.Error:
        return None, None

    for child in proc.children(recursive=True):
        try:
            cpu = child.cpu_times()
        except psutil.Error:
            continue
        user += cpu.user
        system += cpu.system
    return user, system


def run_pipeline(dataset: BenchmarkDataset, workspace: Dict[str, Path], output_dir: Path) -> BenchmarkResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / f"{dataset.slug}-pipeline.log"

    cmd = [
        "python3",
        "run.py",
        "--mode",
        dataset.mode,
        "--LIGANDS_DIR",
        str(workspace["ligands"].resolve()),
        "--DOCKING_DIR",
        str(workspace["docking"].resolve()),
        "--ANALYSIS_DIR",
        str(workspace["analysis"].resolve()),
        "--VINA_DIR",
        str(workspace["vina"].resolve()),
        "--AUTODOCK_GPU_DIR",
        str(workspace["autodock"].resolve()),
        "--MACRO_MOL_DIR",
        str(workspace["macro"].resolve()),
        "--RESULTS_DIR",
        str(workspace["results"].resolve()),
    ] + dataset.run_args

    start = time.perf_counter()
    print(f"[INFO] Capturing pipeline output to {log_path}")

    with open(log_path, "w", encoding="utf-8") as log_file:
        proc = subprocess.Popen(
            cmd,
            cwd=DOCKING_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        peak_rss_mb: Optional[float] = None
        cpu_user_s: Optional[float] = None
        cpu_system_s: Optional[float] = None

        if psutil is not None:
            try:
                ps_proc = psutil.Process(proc.pid)
            except psutil.Error:
                ps_proc = None
        else:
            ps_proc = None

        while True:
            ret = proc.poll()
            if ps_proc is not None:
                rss_bytes = _gather_process_tree_memory(ps_proc)
                if rss_bytes:
                    peak_rss_mb = max(peak_rss_mb or 0.0, rss_bytes / (1024**2))
            if ret is not None:
                proc.wait()
                break
            time.sleep(0.5)

        if ps_proc is not None:
            cpu_user_s, cpu_system_s = _collect_cpu_times(ps_proc)

    wall_time_s = time.perf_counter() - start
    returncode = proc.returncode if proc.returncode is not None else -1

    db_metrics = collect_db_metrics(DEFAULT_DB_PATH)
    accuracy_records = compute_grid_accuracy(dataset, workspace)

    return BenchmarkResult(
        dataset=dataset,
        run_command=cmd,
        returncode=returncode,
        wall_time_s=wall_time_s,
        peak_rss_mb=peak_rss_mb,
        cpu_user_s=cpu_user_s,
        cpu_system_s=cpu_system_s,
        db_metrics=db_metrics,
        accuracy_records=accuracy_records,
        log_path=log_path,
    )


def collect_db_metrics(db_path: Path) -> Dict[str, object]:
    if not db_path.exists():
        return {}

    metrics: Dict[str, object] = {"path": str(db_path)}
    try:
        with sqlite3.connect(db_path) as conn:
            cur = conn.cursor()
            total = cur.execute("SELECT COUNT(*) FROM docking_results").fetchone()[0]
            metrics["records"] = int(total)
            if total:
                row = cur.execute(
                    "SELECT MIN(\"binding_affinity (kcal/mol)\"), MAX(\"binding_affinity (kcal/mol)\"), AVG(\"binding_affinity (kcal/mol)\") FROM docking_results"
                ).fetchone()
                metrics["affinity_min"] = row[0]
                metrics["affinity_max"] = row[1]
                metrics["affinity_avg"] = row[2]
    except sqlite3.Error as exc:  # pragma: no cover - defensive
        metrics["error"] = str(exc)
    return metrics


def parse_centers(centers_path: Path) -> Dict[str, tuple[float, float, float]]:
    centers: Dict[str, tuple[float, float, float]] = {}
    if not centers_path.exists():
        return centers

    header: List[str] | None = None
    with centers_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(",")
            if header is None:
                header = [p.strip().lower() for p in parts]
                continue
            row = {k: v for k, v in zip(header, parts)}
            site_key = row.get("site") or row.get("site_id") or row.get("id")
            if not site_key:
                continue
            try:
                x = float(row.get("center_x") or row.get("x") or row.get("cx"))
                y = float(row.get("center_y") or row.get("y") or row.get("cy"))
                z = float(row.get("center_z") or row.get("z") or row.get("cz"))
            except (TypeError, ValueError):
                continue
            centers[site_key] = (x, y, z)
    return centers


def compute_centroid_from_pdbqt(path: Path) -> Optional[tuple[float, float, float]]:
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
    except FileNotFoundError:
        return None

    if not xs:
        return None
    return (
        sum(xs) / len(xs),
        sum(ys) / len(ys),
        sum(zs) / len(zs),
    )


def compute_grid_accuracy(dataset: BenchmarkDataset, workspace: Dict[str, Path]) -> List[AccuracyRecord]:
    receptor_copy = workspace.get("receptor_target")
    if receptor_copy is None:
        return []

    centers_dir = receptor_copy.parent / receptor_copy.stem
    centers_path = centers_dir / "centers.tsv"
    centers = parse_centers(centers_path)
    if not centers:
        return []

    records: List[AccuracyRecord] = []
    for ligand in dataset.ligands:
        if ligand.reference is None:
            continue
        centroid = compute_centroid_from_pdbqt(ligand.reference)
        if centroid is None:
            continue
        best_site = None
        best_dist = math.inf
        for site_id, center in centers.items():
            dist = math.dist(centroid, center)
            if dist < best_dist:
                best_dist = dist
                best_site = site_id
        if best_site is None:
            continue
        records.append(
            AccuracyRecord(
                ligand=ligand.name,
                reference_path=ligand.reference,
                closest_site=best_site,
                min_distance_A=best_dist,
            )
        )
    return records


def write_summary(result: BenchmarkResult, output_dir: Path) -> Path:
    summary_path = output_dir / f"{result.dataset.slug}-summary.json"
    summary_path.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return summary_path


def run(args: argparse.Namespace) -> List[BenchmarkResult]:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results: List[BenchmarkResult] = []
    for dataset_path in args.dataset:
        dataset = load_dataset(Path(dataset_path).resolve())
        print(f"[INFO] Benchmarking dataset: {dataset.name} (mode={dataset.mode})")
        workspace = stage_workspace(dataset)
        result = run_pipeline(dataset, workspace, output_dir)
        write_summary(result, output_dir)
        results.append(result)
        if result.returncode != 0:
            print(f"[WARN] Pipeline exited with code {result.returncode}")
        else:
            print(f"[OK] Completed in {result.wall_time_s:.1f}s; summary → {output_dir}")
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark the Ultidock pipeline end-to-end",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dataset",
        "-d",
        action="append",
        required=True,
        help="Path to a benchmarking dataset JSON file (can be specified multiple times).",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Directory where benchmark logs and summaries will be stored.",
    )
    return parser


def main(argv: Optional[Iterable[str]] = None) -> List[BenchmarkResult]:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    main()