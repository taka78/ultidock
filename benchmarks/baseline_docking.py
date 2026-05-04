#!/usr/bin/env python3
"""
Standalone Baseline Vina Docking Script.

Usage:
  python3 benchmarks/baseline_docking.py --receptor <path> --ligands <dir>
  python3 benchmarks/baseline_docking.py --receptor <path> --ligands <dir> --target cdk2

Center modes:
  --center-mode centroid   Use geometric centroid of the receptor (default)
  --center-mode crystal    Use centroid of a co-crystal ligand (requires --crystal-ligand)
  --center-mode manual     Use explicit --center-x/y/z values

Output structure (when --target is provided):
  benchmarks/results/dude/<target>/baseline_<center_mode>/
    ├── baseline_results.db
    ├── poses/
    ├── scores.csv
    └── summary.json

- Runs Vina for all .pdbqt ligands in the ligands directory in parallel.
- Parses binding affinity from the output PDBQT file (REMARK VINA RESULT).
- Logs results into an SQLite database matching the main UltiDock schema.
- Writes a summary.json compatible with statistical_analysis.py.
"""

import argparse
import concurrent.futures
import csv
import gzip
import json
import math
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_ROOT = REPO_ROOT / "benchmarks" / "results" / "dude"


def get_receptor_centroid(pdbqt_path: Path) -> tuple[float, float, float]:
    """Calculate the geometric center (centroid) of ATOM/HETATM records in a PDBQT file."""
    x_coords, y_coords, z_coords = [], [], []
    with pdbqt_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
                except ValueError:
                    continue
    if not x_coords:
        raise ValueError(f"No valid ATOM/HETATM coordinates found in {pdbqt_path}")

    cx = sum(x_coords) / len(x_coords)
    cy = sum(y_coords) / len(y_coords)
    cz = sum(z_coords) / len(z_coords)
    return cx, cy, cz


def get_crystal_ligand_centroid(mol2_path: Path) -> tuple[float, float, float]:
    """Calculate centroid from ATOM records in a .mol2 or .mol2.gz file."""
    x_coords, y_coords, z_coords = [], [], []

    opener = gzip.open if str(mol2_path).endswith(".gz") else open
    in_atoms = False
    with opener(mol2_path, "rt", errors="ignore") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            elif stripped.startswith("@<TRIPOS>"):
                in_atoms = False
                continue
            if in_atoms and stripped:
                parts = stripped.split()
                if len(parts) >= 5:
                    try:
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                        x_coords.append(x)
                        y_coords.append(y)
                        z_coords.append(z)
                    except ValueError:
                        continue

    if not x_coords:
        raise ValueError(f"No valid atom coordinates found in {mol2_path}")

    cx = sum(x_coords) / len(x_coords)
    cy = sum(y_coords) / len(y_coords)
    cz = sum(z_coords) / len(z_coords)
    return cx, cy, cz


def init_db(db_path: Path) -> None:
    """Initialize the SQLite database with the same schema as docking_results."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        conn.execute('''
            CREATE TABLE IF NOT EXISTS docking_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ligand_name TEXT,
                "binding_affinity (kcal/mol)" REAL,
                "rmsd_lb (Å)" REAL,
                "rmsd_ub (Å)" REAL,
                "docking_file" TEXT,
                binding_site TEXT,
                docking_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                cuda_setup_s REAL,
                rest_setup_s REAL,
                total_setup_s REAL,
                job_wait_s REAL,
                job_duration_s REAL,
                run_time_s REAL,
                processing_time_s REAL,
                run_plus_processing_s REAL,
                shutdown_s REAL,
                gpu_memory_mb REAL
            )
        ''')


def save_to_db(db_path: Path, result: dict) -> None:
    with sqlite3.connect(db_path, timeout=30.0) as conn:
        conn.execute('''
            INSERT INTO docking_results (
                ligand_name,
                "binding_affinity (kcal/mol)",
                "rmsd_lb (Å)",
                "rmsd_ub (Å)",
                docking_file,
                binding_site,
                job_duration_s
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (
            result["ligand_name"],
            result["affinity"],
            result["rmsd_lb"],
            result["rmsd_ub"],
            str(result["docking_file"]),
            result["binding_site"],
            result["duration_s"]
        ))


def run_vina(
    ligand_path: Path,
    receptor_path: Path,
    cx: float, cy: float, cz: float,
    box_size: float,
    out_path: Path,
    vina_bin: str = "vina",
    exhaustiveness: int = 8,
    cpu: int = 1,
    binding_site_label: str = "baseline",
) -> dict:
    """Run Vina and extract the top binding affinity from the output PDBQT."""
    start_t = time.time()

    cmd = [
        vina_bin,
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--center_x", f"{cx:.3f}",
        "--center_y", f"{cy:.3f}",
        "--center_z", f"{cz:.3f}",
        "--size_x", f"{box_size:.1f}",
        "--size_y", f"{box_size:.1f}",
        "--size_z", f"{box_size:.1f}",
        "--out", str(out_path),
        "--exhaustiveness", str(exhaustiveness),
        "--cpu", str(cpu)
    ]

    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        duration = time.time() - start_t

        # Parse affinity and RMSD from the output PDBQT file directly (Model 1)
        affinity = None
        rmsd_lb = 0.0
        rmsd_ub = 0.0

        if out_path.exists():
            with out_path.open("r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                affinity = float(parts[3])
                                rmsd_lb = float(parts[4])
                                rmsd_ub = float(parts[5])
                            except ValueError:
                                pass
                        break  # Only need the first (best) model

        return {
            "success": True,
            "ligand_name": ligand_path.stem,
            "affinity": affinity,
            "rmsd_lb": rmsd_lb,
            "rmsd_ub": rmsd_ub,
            "docking_file": out_path,
            "binding_site": binding_site_label,
            "duration_s": duration,
            "log": proc.stdout
        }
    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "ligand_name": ligand_path.stem,
            "error": e.stdout,
            "docking_file": out_path
        }
    except Exception as e:
        return {
            "success": False,
            "ligand_name": ligand_path.stem,
            "error": str(e),
            "docking_file": out_path
        }


def main():
    parser = argparse.ArgumentParser(
        description="Baseline parallel Vina docking script.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--receptor", required=True, type=Path, help="Path to receptor .pdbqt")
    parser.add_argument("--ligands", required=True, type=Path, help="Directory containing ligand .pdbqt files")
    parser.add_argument(
        "--target", type=str, default=None,
        help="Target name (e.g. cdk2). When set, output goes to "
             "<results-root>/<target>/baseline_<center-mode>/ automatically.",
    )
    parser.add_argument(
        "--results-root", type=Path, default=DEFAULT_RESULTS_ROOT,
        help="Root directory for structured benchmark results.",
    )
    parser.add_argument("--output_db", type=Path, default=None, help="SQLite DB path (auto-derived when --target is set)")
    parser.add_argument("--output_dir", type=Path, default=None, help="Directory for output poses (auto-derived when --target is set)")
    parser.add_argument("--vina_bin", default="vina", help="Vina executable name/path")
    parser.add_argument("--workers", type=int, default=os.cpu_count(), help="Number of parallel ligands")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina exhaustiveness")
    parser.add_argument("--seed", type=int, default=0, help="Vina random seed for reproducibility")
    parser.add_argument("--box-size", type=float, default=25.0, help="Cubic box side length in Angstroms")

    # Center mode
    parser.add_argument(
        "--center-mode",
        choices=["centroid", "crystal", "manual"],
        default="centroid",
        help="How to determine box center: 'centroid' (receptor geometric center), "
             "'crystal' (co-crystal ligand centroid), 'manual' (explicit coordinates)"
    )
    parser.add_argument("--crystal-ligand", type=Path, help="Path to crystal_ligand.mol2 (required for --center-mode crystal)")
    parser.add_argument("--center-x", type=float, help="Manual center X coordinate")
    parser.add_argument("--center-y", type=float, help="Manual center Y coordinate")
    parser.add_argument("--center-z", type=float, help="Manual center Z coordinate")

    args = parser.parse_args()

    if not args.receptor.exists():
        sys.exit(f"Error: Receptor not found at {args.receptor}")
    if not args.ligands.is_dir():
        sys.exit(f"Error: Ligands directory not found at {args.ligands}")

    # Resolve Vina binary (allow both PATH lookup and absolute paths)
    vina_bin = args.vina_bin
    if not os.path.isabs(vina_bin) and not shutil.which(vina_bin):
        sys.exit(f"Error: Vina executable '{vina_bin}' not found in PATH.")

    # Derive structured output paths when --target is set
    if args.target:
        method_label = f"baseline_{args.center_mode}"
        result_dir = args.results_root / args.target / method_label
        result_dir.mkdir(parents=True, exist_ok=True)
        if args.output_db is None:
            args.output_db = result_dir / "baseline_results.db"
        if args.output_dir is None:
            args.output_dir = result_dir / "poses"
    else:
        # Legacy fallback for unstructured usage
        if args.output_db is None:
            args.output_db = Path("docking/RESULTS_DIR/baseline_results.db")
        if args.output_dir is None:
            args.output_dir = Path("docking/RESULTS_DIR/baseline_poses")
        result_dir = args.output_dir.parent

    args.output_dir.mkdir(parents=True, exist_ok=True)
    init_db(args.output_db)

    # 1. Determine center
    if args.center_mode == "crystal":
        if not args.crystal_ligand:
            sys.exit("Error: --crystal-ligand is required when --center-mode is 'crystal'")
        if not args.crystal_ligand.exists():
            sys.exit(f"Error: Crystal ligand not found at {args.crystal_ligand}")
        cx, cy, cz = get_crystal_ligand_centroid(args.crystal_ligand)
        site_label = "baseline_crystal"
        print(f"[INIT] Crystal ligand centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    elif args.center_mode == "manual":
        if args.center_x is None or args.center_y is None or args.center_z is None:
            sys.exit("Error: --center-x, --center-y, --center-z are required for --center-mode manual")
        cx, cy, cz = args.center_x, args.center_y, args.center_z
        site_label = "baseline_manual"
        print(f"[INIT] Manual center: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    else:  # centroid
        cx, cy, cz = get_receptor_centroid(args.receptor)
        site_label = "baseline_centroid"
        print(f"[INIT] Receptor centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    box_size = args.box_size
    print(f"[INIT] Box size: {box_size:.1f}×{box_size:.1f}×{box_size:.1f} Å")

    # 2. Find ligands
    ligand_files = list(args.ligands.glob("*.pdbqt"))
    if not ligand_files:
        sys.exit(f"Error: No .pdbqt files found in {args.ligands}")
    print(f"[INIT] Found {len(ligand_files)} ligands to dock using {args.workers} workers.")

    # 3. Parallel Docking
    completed = 0
    failed = 0
    start_time = time.time()

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for lig in ligand_files:
            out_file = args.output_dir / f"{lig.stem}_out.pdbqt"
            fut = executor.submit(
                run_vina,
                ligand_path=lig,
                receptor_path=args.receptor,
                cx=cx, cy=cy, cz=cz,
                box_size=box_size,
                out_path=out_file,
                vina_bin=vina_bin,
                exhaustiveness=args.exhaustiveness,
                cpu=1,  # each worker gets 1 thread to maximize overall throughput
                binding_site_label=site_label,
            )
            futures[fut] = lig

        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            completed += 1
            if res["success"]:
                if res["affinity"] is not None:
                    save_to_db(args.output_db, res)
                    print(f"[{completed}/{len(ligand_files)}] Docked {res['ligand_name']} -> {res['affinity']} kcal/mol")
                else:
                    failed += 1
                    print(f"[{completed}/{len(ligand_files)}] {res['ligand_name']} failed: Could not parse affinity")
            else:
                failed += 1
                error_snippet = str(res.get('error', ''))[:200].replace('\n', ' ')
                print(f"[{completed}/{len(ligand_files)}] {res['ligand_name']} failed: {error_snippet}")

    total_time = time.time() - start_time
    succeeded = completed - failed
    print("-" * 50)
    print(f"BASELINE DOCKING COMPLETE")
    print(f"Center mode: {args.center_mode} ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"Box size:    {box_size:.1f} Å")
    print(f"Total time:  {total_time:.1f} s")
    print(f"Total docked:{completed}")
    print(f"Succeeded:   {succeeded}")
    print(f"Failed:      {failed}")
    print(f"Results DB:  {args.output_db.resolve()}")

    # ── Write scores.csv + summary.json matching the UltiDock benchmark schema ──
    # Collect all docked results from the DB
    with sqlite3.connect(args.output_db) as conn:
        rows = conn.execute(
            'SELECT ligand_name, "binding_affinity (kcal/mol)", '
            '"rmsd_lb (Å)", "rmsd_ub (Å)", docking_file, binding_site '
            'FROM docking_results'
        ).fetchall()

    # Build label map from filename prefixes (active__ / decoy__)
    score_rows: list[dict] = []
    for ligand_name, affinity, rmsd_lb, rmsd_ub, docking_file, binding_site in rows:
        if affinity is None:
            continue
        # Determine label from filename prefix
        if ligand_name and ligand_name.startswith("active__"):
            label = 1
        elif ligand_name and ligand_name.startswith("decoy__"):
            label = 0
        else:
            label = -1  # unknown — excluded from metrics
        score_rows.append({
            "ligand_id": ligand_name or "",
            "label": label,
            "best_affinity_kcal_mol": float(affinity),
            "best_binding_site": binding_site or site_label,
            "dock_rows": 1,
            "prepared_path": docking_file or "",
        })

    # Keep only best (most negative) affinity per ligand
    best_by_ligand: dict[str, dict] = {}
    for row in score_rows:
        lid = row["ligand_id"]
        if lid not in best_by_ligand or row["best_affinity_kcal_mol"] < best_by_ligand[lid]["best_affinity_kcal_mol"]:
            best_by_ligand[lid] = row
    score_rows = list(best_by_ligand.values())

    # Write scores.csv
    scores_csv_path = result_dir / "scores.csv"
    with scores_csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ligand_id", "label", "best_affinity_kcal_mol", "best_binding_site", "dock_rows", "prepared_path"])
        for row in score_rows:
            writer.writerow([
                row["ligand_id"], row["label"],
                f"{row['best_affinity_kcal_mol']:.4f}",
                row["best_binding_site"], row["dock_rows"],
                row["prepared_path"],
            ])
    print(f"Scores CSV: {scores_csv_path}")

    # Compute metrics (only for labelled ligands)
    labelled = [r for r in score_rows if r["label"] in (0, 1)]
    metrics: dict = {
        "target": args.target or args.receptor.stem,
        "mode": f"baseline_{args.center_mode}",
        "n_actives": sum(1 for r in labelled if r["label"] == 1),
        "n_decoys": sum(1 for r in labelled if r["label"] == 0),
        "n_total": len(labelled),
        "n_docked": succeeded,
        "n_failed": failed,
        "total_time_s": round(total_time, 1),
        "box": {
            "center_x": cx, "center_y": cy, "center_z": cz,
            "size_x": box_size, "size_y": box_size, "size_z": box_size,
        },
        "db_path": str(args.output_db.resolve()),
        "scores_csv": str(scores_csv_path.resolve()),
    }

    if labelled and metrics["n_actives"] > 0 and metrics["n_decoys"] > 0:
        labels_arr = np.array([r["label"] for r in labelled], dtype=int)
        affinities_arr = np.array([r["best_affinity_kcal_mol"] for r in labelled], dtype=float)
        ranking = -affinities_arr  # more negative = better → negate for ranking

        # Import scoring functions from the benchmark module
        try:
            from benchmarks.dude_docking_benchmark import roc_auc, enrichment_factor, bedroc_like, log_auc
        except ImportError:
            from dude_docking_benchmark import roc_auc, enrichment_factor, bedroc_like, log_auc

        metrics["roc_auc"] = roc_auc(labels_arr, ranking)
        metrics["ef1"] = enrichment_factor(labels_arr, ranking, 0.01)
        metrics["ef5"] = enrichment_factor(labels_arr, ranking, 0.05)
        metrics["bedroc_alpha"] = 20.0
        metrics["bedroc"] = bedroc_like(labels_arr, ranking, 20.0)
        metrics["log_auc"] = log_auc(labels_arr, ranking)

        print(f"  ROC-AUC = {metrics['roc_auc']:.4f}")
        print(f"  EF1%    = {metrics['ef1']:.2f}")
        print(f"  EF5%    = {metrics['ef5']:.2f}")
        print(f"  BEDROC  = {metrics['bedroc']:.4f}")
        print(f"  LogAUC  = {metrics['log_auc']:.4f}")
    else:
        print("  [WARN] Cannot compute enrichment metrics — missing active/decoy labels.")

    # Write summary.json (same schema as dude_docking_benchmark.py)
    summary = {"targets": [metrics]}
    summary_path = result_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Summary:    {summary_path}")


if __name__ == "__main__":
    main()
