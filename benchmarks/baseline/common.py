"""Shared helpers for receptor-only baseline cavity-finder benchmarks."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BENCHMARK_ROOT.parent
DEFAULT_DATASET_ROOT = BENCHMARK_ROOT / "datasets"
DEFAULT_THRESHOLDS = (4.0, 5.0, 10.0)

for import_root in (str(REPO_ROOT), str(BENCHMARK_ROOT)):
    if import_root not in sys.path:
        sys.path.insert(0, import_root)

import cavity_recovery_benchmark as cavity  # noqa: E402
from molguard.io.receptor_prep import prepare_receptor_pdbqt, receptor_pdbqt_output_name  # noqa: E402


@dataclass(frozen=True)
class PredictedSite:
    """A method-proposed binding-site center."""

    site_id: str
    center: tuple[float, float, float]
    score: float | None = None
    source: str = ""


@dataclass(frozen=True)
class PredictionResult:
    """Raw method output plus the receptor file actually given to the method."""

    sites: list[PredictedSite]
    raw_output_dir: Path
    method_receptor: Path


Predictor = Callable[[Path, Path, argparse.Namespace], PredictionResult]


def default_jobs() -> int:
    for env_name in ("ULTIDOCK_BENCHMARK_WORKERS", "ULTIDOCK_WORKERS"):
        raw_value = os.environ.get(env_name)
        if raw_value:
            try:
                return max(1, int(raw_value))
            except ValueError:
                pass
    cpu_count = os.cpu_count() or 1
    return max(1, min(4, cpu_count // 2 or 1))


def add_common_args(
    parser: argparse.ArgumentParser,
    *,
    method_name: str,
    default_output_dir: Path,
) -> None:
    parser.add_argument(
        "--dataset-root",
        default=str(DEFAULT_DATASET_ROOT),
        help="Root directory containing DUD-E target subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(default_output_dir),
        help=f"Directory for {method_name} baseline outputs.",
    )
    parser.add_argument(
        "--targets",
        help="Comma-separated target names, or 'all'. Default: all discoverable targets.",
    )
    parser.add_argument("--max-targets", type=int, help="Run only the first N selected targets.")
    parser.add_argument(
        "--jobs",
        type=int,
        default=default_jobs(),
        help="Number of targets to process concurrently.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Deterministic receptor-prep seed.")
    parser.add_argument("--autosites", type=int, default=6, help="Number of sites to evaluate.")
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
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate method outputs even when centers.tsv already exists.",
    )
    parser.add_argument(
        "--receptor-prepare-command",
        help="Optional receptor conversion command template with {input} and {output}.",
    )
    parser.add_argument(
        "--baseline-receptor-model",
        choices=["ultidock-prepared", "original"],
        default="ultidock-prepared",
        help=(
            "Receptor model passed to the baseline finder. "
            "'ultidock-prepared' uses the same molguard/Meeko receptor prep as Ultidock, "
            "then converts the prepared PDBQT back to PDB records for tools that need PDB."
        ),
    )
    parser.add_argument("--dry-run", action="store_true", help="Print selected targets and exit.")


def run_command(command: list[str], *, cwd: Path) -> None:
    cwd.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            command,
            cwd=cwd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        stdout = (exc.stdout or "").strip()
        stderr = (exc.stderr or "").strip()
        chunks = [
            f"command failed: {' '.join(command)}",
            f"exit code: {exc.returncode}",
        ]
        if stdout:
            chunks.append(f"stdout:\n{stdout}")
        if stderr:
            chunks.append(f"stderr:\n{stderr}")
        raise RuntimeError("\n\n".join(chunks)) from exc


def split_tool_command(command: str) -> list[str]:
    """Split a tool command and absolutize path-like executables."""
    parts = shlex.split(command)
    if not parts:
        raise ValueError("tool command cannot be empty")
    executable = Path(parts[0]).expanduser()
    if "/" in parts[0] or "\\" in parts[0]:
        parts[0] = str(executable.resolve())
    return parts


def materialize_pdb_receptor(receptor_input: Path, work_dir: Path) -> Path:
    """Copy or unpack a receptor PDB into a method-local work directory."""
    work_dir.mkdir(parents=True, exist_ok=True)
    suffix = cavity._effective_suffix(receptor_input)  # type: ignore[attr-defined]
    if suffix != ".pdb":
        raise ValueError(
            f"baseline methods require a PDB receptor; got {receptor_input.name}"
        )

    staged = work_dir / f"{receptor_input.stem.removesuffix('.pdb')}.pdb"
    if receptor_input.suffix.lower() == ".gz":
        import gzip

        with gzip.open(receptor_input, "rb") as src, staged.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(receptor_input, staged)
    return staged


def _infer_element_from_atom_name(atom_name_field: str) -> str:
    field = atom_name_field[:4].ljust(4)
    token = "".join(ch for ch in field.strip() if ch.isalpha())
    if not token:
        return ""
    two_letter = {
        "BR",
        "CA",
        "CL",
        "CO",
        "CU",
        "FE",
        "MG",
        "MN",
        "NA",
        "NI",
        "ZN",
    }
    if field[0].isspace():
        return token[0].upper()
    token = token.upper()
    if len(token) >= 2 and token[:2] in two_letter:
        return token[:2].capitalize()
    return token[0].upper()


def pdb_from_pdbqt(pdbqt_path: Path, pdb_path: Path) -> Path:
    """Write a PDB-like coordinate file from a prepared receptor PDBQT."""
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    out_lines: list[str] = []
    for raw in pdbqt_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        record = raw[:66].ljust(76)
        element = raw[76:78].strip() if len(raw) >= 78 else ""
        if not element:
            element = _infer_element_from_atom_name(raw[12:16])
        out_lines.append(f"{record}{element:>2s}")
    if not out_lines:
        raise ValueError(f"prepared receptor contains no atom records: {pdbqt_path}")
    out_lines.append("END")
    pdb_path.write_text("\n".join(out_lines) + "\n", encoding="ascii")
    return pdb_path


def prepare_method_receptor(
    *,
    receptor_input: Path,
    target_output: Path,
    method_work_dir: Path,
    args: argparse.Namespace,
) -> Path:
    """
    Stage the receptor that will be passed to a baseline method.

    By default the baseline receives the same receptor model as Ultidock:
    molguard/Meeko-prepared PDBQT converted back into PDB coordinate records.
    """
    method_work_dir.mkdir(parents=True, exist_ok=True)
    if args.baseline_receptor_model == "original":
        return materialize_pdb_receptor(receptor_input, method_work_dir)

    prepared_dir = target_output / "prepared"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    prepared_pdbqt = prepared_dir / receptor_pdbqt_output_name(receptor_input)
    prepared_pdb = prepared_dir / f"{prepared_pdbqt.stem}.pdb"
    if args.force or not prepared_pdbqt.exists() or not prepared_pdb.exists():
        prepare_receptor_pdbqt(
            input_path=receptor_input,
            output_path=prepared_pdbqt,
            prepare_command=args.receptor_prepare_command,
            seed=int(args.seed),
            timestamp="BASELINE_BENCHMARK",
        )
        pdb_from_pdbqt(prepared_pdbqt, prepared_pdb)

    staged = method_work_dir / prepared_pdb.name
    if args.force or not staged.exists():
        shutil.copy2(prepared_pdb, staged)
    return staged


def normalize_sites(sites: list[PredictedSite], limit: int) -> list[PredictedSite]:
    normalized: list[PredictedSite] = []
    for index, site in enumerate(sites[:limit], start=1):
        normalized.append(
            PredictedSite(
                site_id=f"S{index}",
                center=site.center,
                score=site.score,
                source=site.source,
            )
        )
    return normalized


def read_existing_centers(path: Path, limit: int) -> list[PredictedSite]:
    parsed_rows = cavity.parse_scored_centers_tsv(path)
    sites = [
        PredictedSite(
            site_id=str(row["site_id"]),
            center=row["center"],
            score=row.get("fitness_score"),
            source=str(path),
        )
        for row in sorted(parsed_rows, key=lambda row: cavity._site_sort_key(str(row["site_id"])))
    ]
    return normalize_sites(sites, limit)


def write_centers_tsv(path: Path, target_name: str, sites: list[PredictedSite]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["target", "site_id", "center_x", "center_y", "center_z", "method_score", "source"]
        )
        for site in sites:
            writer.writerow(
                [
                    target_name,
                    site.site_id,
                    f"{site.center[0]:.4f}",
                    f"{site.center[1]:.4f}",
                    f"{site.center[2]:.4f}",
                    "" if site.score is None else f"{site.score:.6g}",
                    site.source,
                ]
            )


def _csv_value(value: Any) -> Any:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.4f}"
        return ""
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fieldnames})


def write_predictions_tsv(
    path: Path,
    site_rows: list[dict[str, Any]],
    *,
    method_name: str,
) -> None:
    """Write normalized prediction centers for cross-dataset evaluators."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "target_id",
        "method",
        "rank",
        "site_id",
        "center_x",
        "center_y",
        "center_z",
        "score",
        "source",
    ]
    ordered = sorted(
        site_rows,
        key=lambda row: (str(row["target"]), int(row["rank_by_method"]), str(row["site_id"])),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in ordered:
            writer.writerow(
                {
                    "target_id": row["target"],
                    "method": method_name,
                    "rank": int(row["rank_by_method"]),
                    "site_id": row["site_id"],
                    "center_x": _csv_value(row["center_x"]),
                    "center_y": _csv_value(row["center_y"]),
                    "center_z": _csv_value(row["center_z"]),
                    "score": _csv_value(row.get("method_score")),
                    "source": row.get("source", ""),
                }
            )


def evaluate_sites(
    *,
    method_name: str,
    target_name: str,
    reference_ligand: Path,
    receptor_input: Path,
    centers_tsv: Path,
    raw_output_dir: Path,
    method_receptor: Path | None,
    sites: list[PredictedSite],
    thresholds: list[float],
    box_size: float,
    autosites: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    reference_atoms = cavity.calculate_reference_atoms(reference_ligand)
    crystal_center = cavity.calculate_reference_center(reference_ligand)
    site_rows: list[dict[str, Any]] = []
    box_label = cavity.threshold_label(box_size)

    for site_order, site in enumerate(sites, start=1):
        distance = cavity.euclidean_distance(crystal_center, site.center)
        dcc = cavity.distance_to_closest_reference_atom(site.center, reference_atoms)
        delta = cavity.vector_delta(crystal_center, site.center)
        row: dict[str, Any] = {
            "target": target_name,
            "method": method_name,
            "site_id": site.site_id,
            "site_order": site_order,
            "rank_by_method": site_order,
            "method_score": site.score,
            "distance_a": distance,
            "dcc_a": dcc,
            "center_x": site.center[0],
            "center_y": site.center[1],
            "center_z": site.center[2],
            "delta_x": delta[0],
            "delta_y": delta[1],
            "delta_z": delta[2],
            "required_cube_side_a": cavity.required_cube_side(delta),
            f"within_{box_label}_cube": cavity.within_cube(delta, box_size),
            "source": site.source,
        }
        for threshold in thresholds:
            row[f"hit_at_{cavity.threshold_label(threshold)}"] = distance <= threshold
            row[f"dcc_hit_at_{cavity.threshold_label(threshold)}"] = dcc <= threshold
        site_rows.append(row)

    if not site_rows:
        raise ValueError(f"{target_name}: no predicted sites")

    best = min(site_rows, key=lambda row: float(row["distance_a"]))
    best_dcc = min(site_rows, key=lambda row: float(row["dcc_a"]))
    rank1 = min(site_rows, key=lambda row: int(row["rank_by_method"]))
    rank3_pool = [row for row in site_rows if int(row["rank_by_method"]) <= 3]
    rank3_best = min(rank3_pool, key=lambda row: float(row["dcc_a"]))
    summary: dict[str, Any] = {
        "target": target_name,
        "status": "ok",
        "method": method_name,
        "autosites": autosites,
        "n_sites": len(site_rows),
        "avg_predicted_sites": len(site_rows),
        "crystal_center_x": crystal_center[0],
        "crystal_center_y": crystal_center[1],
        "crystal_center_z": crystal_center[2],
        "reference_ligand_atom_count": len(reference_atoms),
        "best_site_id": best["site_id"],
        "best_site_order": best["site_order"],
        "best_distance_a": best["distance_a"],
        "best_required_cube_side_a": best["required_cube_side_a"],
        f"best_within_{box_label}_cube": best[f"within_{box_label}_cube"],
        "best_dcc_site_id": best_dcc["site_id"],
        "best_dcc_rank_by_method": best_dcc["rank_by_method"],
        "best_dcc_a": best_dcc["dcc_a"],
        "rank1_site_id": rank1["site_id"],
        "rank1_method_score": rank1["method_score"],
        "rank1_distance_a": rank1["distance_a"],
        "rank1_dcc_a": rank1["dcc_a"],
        "rank3_best_site_id": rank3_best["site_id"],
        "rank3_best_rank_by_method": rank3_best["rank_by_method"],
        "rank3_best_dcc_a": rank3_best["dcc_a"],
        "reference_ligand": str(reference_ligand.resolve()),
        "receptor_input": str(receptor_input.resolve()),
        "centers_tsv": str(centers_tsv.resolve()),
        "raw_output_dir": str(raw_output_dir.resolve()),
        "method_receptor": str(method_receptor.resolve()) if method_receptor else "",
        "error": "",
    }
    for threshold in thresholds:
        label = cavity.threshold_label(threshold)
        summary[f"any_site_success_at_{label}"] = bool(best[f"hit_at_{label}"])
        summary[f"rank1_dcc_success_at_{label}"] = bool(rank1[f"dcc_hit_at_{label}"])
        summary[f"rank3_dcc_success_at_{label}"] = any(
            bool(row[f"dcc_hit_at_{label}"]) for row in rank3_pool
        )
        summary[f"any_site_dcc_success_at_{label}"] = bool(best_dcc[f"dcc_hit_at_{label}"])

    evaluation = {
        "target": target_name,
        "method": method_name,
        "reference_ligand": str(reference_ligand.resolve()),
        "receptor_input": str(receptor_input.resolve()),
        "centers_tsv": str(centers_tsv.resolve()),
        "raw_output_dir": str(raw_output_dir.resolve()),
        "method_receptor": str(method_receptor.resolve()) if method_receptor else "",
        "crystal_center": {
            "x": crystal_center[0],
            "y": crystal_center[1],
            "z": crystal_center[2],
        },
        "reference_ligand_atom_count": len(reference_atoms),
        "ranking": "rank_by_method follows the method score/order; site labels are method-local names",
        "thresholds_a": thresholds,
        "box_size_a": box_size,
        "summary": summary,
        "sites": site_rows,
    }
    return summary, site_rows, evaluation


def failure_row(
    *,
    method_name: str,
    target_name: str,
    autosites: int,
    error: str,
    reference_ligand: Path | None = None,
    receptor_input: Path | None = None,
    centers_tsv: Path | None = None,
    raw_output_dir: Path | None = None,
    method_receptor: Path | None = None,
) -> dict[str, Any]:
    return {
        "target": target_name,
        "status": "failed",
        "method": method_name,
        "autosites": autosites,
        "reference_ligand": str(reference_ligand.resolve()) if reference_ligand else "",
        "receptor_input": str(receptor_input.resolve()) if receptor_input else "",
        "centers_tsv": str(centers_tsv.resolve()) if centers_tsv else "",
        "raw_output_dir": str(raw_output_dir.resolve()) if raw_output_dir else "",
        "method_receptor": str(method_receptor.resolve()) if method_receptor else "",
        "error": error,
    }


def write_run_outputs(
    *,
    output_dir: Path,
    summary_rows: list[dict[str, Any]],
    site_rows: list[dict[str, Any]],
    thresholds: list[float],
    box_size: float,
    metadata: dict[str, Any],
) -> None:
    box_label = cavity.threshold_label(box_size)
    threshold_labels = [cavity.threshold_label(value) for value in thresholds]
    summary_fields = [
        "target",
        "status",
        "method",
        "autosites",
        "n_sites",
        "avg_predicted_sites",
        "crystal_center_x",
        "crystal_center_y",
        "crystal_center_z",
        "reference_ligand_atom_count",
        "best_site_id",
        "best_site_order",
        "best_distance_a",
        "best_required_cube_side_a",
        f"best_within_{box_label}_cube",
        "best_dcc_site_id",
        "best_dcc_rank_by_method",
        "best_dcc_a",
        "rank1_site_id",
        "rank1_method_score",
        "rank1_distance_a",
        "rank1_dcc_a",
        "rank3_best_site_id",
        "rank3_best_rank_by_method",
        "rank3_best_dcc_a",
    ]
    for label in threshold_labels:
        summary_fields.append(f"any_site_success_at_{label}")
        summary_fields.append(f"rank1_dcc_success_at_{label}")
        summary_fields.append(f"rank3_dcc_success_at_{label}")
        summary_fields.append(f"any_site_dcc_success_at_{label}")
    summary_fields.extend(
        [
            "reference_ligand",
            "receptor_input",
            "method_receptor",
            "centers_tsv",
            "raw_output_dir",
            "error",
        ]
    )

    site_fields = [
        "target",
        "method",
        "site_id",
        "site_order",
        "rank_by_method",
        "method_score",
        "distance_a",
        "dcc_a",
        "center_x",
        "center_y",
        "center_z",
        "delta_x",
        "delta_y",
        "delta_z",
        "required_cube_side_a",
        f"within_{box_label}_cube",
    ]
    for label in threshold_labels:
        site_fields.append(f"hit_at_{label}")
        site_fields.append(f"dcc_hit_at_{label}")
    site_fields.append("source")

    _write_csv(output_dir / "summary.csv", summary_rows, summary_fields)
    _write_csv(output_dir / "sites.csv", site_rows, site_fields)
    write_predictions_tsv(output_dir / "predictions.tsv", site_rows, method_name=str(metadata["method"]))
    summary_json = {
        "metadata": metadata,
        "n_targets": len(summary_rows),
        "n_ok": sum(1 for row in summary_rows if row.get("status") == "ok"),
        "n_failed": sum(1 for row in summary_rows if row.get("status") != "ok"),
        "targets": summary_rows,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary_json, indent=2), encoding="utf-8")


def run_one_target(
    *,
    method_name: str,
    target_name: str,
    index: int,
    total: int,
    args: argparse.Namespace,
    dataset_root: Path,
    output_dir: Path,
    thresholds: list[float],
    predictor: Predictor,
) -> dict[str, Any]:
    print(f"\n[{index}/{total}] {target_name}")
    target_output = output_dir / target_name
    target_output.mkdir(parents=True, exist_ok=True)
    centers_tsv = target_output / "centers.tsv"
    receptor_input: Path | None = None
    reference_ligand: Path | None = None
    raw_output_dir: Path | None = None
    method_receptor: Path | None = None

    try:
        receptor_input, reference_ligand = cavity.load_target_files(dataset_root, target_name)
        if centers_tsv.exists() and not args.force:
            sites = read_existing_centers(centers_tsv, int(args.autosites))
            raw_output_dir = target_output / "raw" / method_name
            print(f"[centers] reusing {centers_tsv} with {len(sites)} row(s)")
        else:
            prediction = predictor(receptor_input, target_output, args)
            sites = prediction.sites
            raw_output_dir = prediction.raw_output_dir
            method_receptor = prediction.method_receptor
            sites = normalize_sites(sites, int(args.autosites))
            if not sites:
                raise ValueError(f"{target_name}: {method_name} produced no sites")
            write_centers_tsv(centers_tsv, target_name, sites)
            print(f"[centers] wrote {len(sites)} --> {centers_tsv}")

        summary, site_rows, evaluation = evaluate_sites(
            method_name=method_name,
            target_name=target_name,
            reference_ligand=reference_ligand,
            receptor_input=receptor_input,
            centers_tsv=centers_tsv,
            raw_output_dir=raw_output_dir,
            method_receptor=method_receptor,
            sites=sites,
            thresholds=thresholds,
            box_size=float(args.box_size),
            autosites=int(args.autosites),
        )
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
                    "method": method_name,
                    "error": error,
                    "traceback": traceback.format_exc(),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return {
            "summary": failure_row(
                method_name=method_name,
                target_name=target_name,
                autosites=int(args.autosites),
                error=error,
                reference_ligand=reference_ligand,
                receptor_input=receptor_input,
                centers_tsv=centers_tsv,
                raw_output_dir=raw_output_dir,
                method_receptor=method_receptor,
            ),
            "site_rows": [],
        }


def ordered_rows(
    results_by_index: dict[int, dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    site_rows: list[dict[str, Any]] = []
    for index in sorted(results_by_index):
        result = results_by_index[index]
        summary_rows.append(result["summary"])
        site_rows.extend(result.get("site_rows", []))
    return summary_rows, site_rows


def run_baseline(
    *,
    method_name: str,
    args: argparse.Namespace,
    predictor: Predictor,
) -> None:
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if args.autosites < 1:
        raise SystemExit("--autosites must be >= 1")

    dataset_root = Path(args.dataset_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    thresholds = cavity.parse_thresholds(args.thresholds)
    target_names = cavity.parse_targets(args.targets, dataset_root)
    if args.max_targets is not None:
        target_names = target_names[: args.max_targets]
    if not target_names:
        raise SystemExit(f"No targets found in {dataset_root}")
    args.jobs = min(args.jobs, len(target_names))

    print(f"Method:       {method_name}")
    print(f"Dataset root: {dataset_root}")
    print(f"Output dir:   {output_dir}")
    print(f"Targets:      {len(target_names)}")
    print(f"Jobs:         {args.jobs}")
    for name in target_names:
        print(f"  - {name}")
    if args.dry_run:
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "method": method_name,
        "dataset_root": str(dataset_root),
        "output_dir": str(output_dir),
        "targets": target_names,
        "autosites": args.autosites,
        "jobs": args.jobs,
        "thresholds_a": thresholds,
        "box_size_a": args.box_size,
        "git_commit": cavity.git_value("rev-parse", "HEAD"),
        "git_branch": cavity.git_value("branch", "--show-current"),
        "git_status_short": cavity.git_status_short(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    results_by_index: dict[int, dict[str, Any]] = {}

    def persist_completed_outputs() -> None:
        summary_rows, all_site_rows = ordered_rows(results_by_index)
        write_run_outputs(
            output_dir=output_dir,
            summary_rows=summary_rows,
            site_rows=all_site_rows,
            thresholds=thresholds,
            box_size=float(args.box_size),
            metadata=metadata,
        )

    if args.jobs == 1:
        for index, target_name in enumerate(target_names, start=1):
            results_by_index[index] = run_one_target(
                method_name=method_name,
                target_name=target_name,
                index=index,
                total=len(target_names),
                args=args,
                dataset_root=dataset_root,
                output_dir=output_dir,
                thresholds=thresholds,
                predictor=predictor,
            )
            persist_completed_outputs()
    else:
        print(f"\n[parallel] using {args.jobs} target worker(s)")
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs,
            thread_name_prefix=f"{method_name}-worker",
        ) as pool:
            futures = {
                pool.submit(
                    run_one_target,
                    method_name=method_name,
                    target_name=target_name,
                    index=index,
                    total=len(target_names),
                    args=args,
                    dataset_root=dataset_root,
                    output_dir=output_dir,
                    thresholds=thresholds,
                    predictor=predictor,
                ): index
                for index, target_name in enumerate(target_names, start=1)
            }
            for future in concurrent.futures.as_completed(futures):
                index = futures[future]
                results_by_index[index] = future.result()
                persist_completed_outputs()
                print(f"[progress] completed {len(results_by_index)}/{len(target_names)} target(s)")

    print(f"\nSummary CSV: {output_dir / 'summary.csv'}")
    print(f"Sites CSV:   {output_dir / 'sites.csv'}")
