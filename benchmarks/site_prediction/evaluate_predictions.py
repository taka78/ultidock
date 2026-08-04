#!/usr/bin/env python3
"""Evaluate normalized site-prediction outputs against normalized labels."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_NORMALIZED_ROOT = REPO_ROOT / "benchmarks" / "site_prediction" / "datasets" / "normalized"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmarks" / "results" / "site_prediction" / "evaluation"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.site_prediction.evaluation.metrics import (  # noqa: E402
    evaluate_target_standard,
    success_rate,
)
from benchmarks.site_prediction.schema import (  # noqa: E402
    PredictedSite,
    prediction_rows_from_tsv,
    target_from_normalized_dir,
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


def group_predictions(
    predictions: list[PredictedSite],
) -> dict[tuple[str, str], list[PredictedSite]]:
    grouped: dict[tuple[str, str], list[PredictedSite]] = defaultdict(list)
    for prediction in predictions:
        grouped[(prediction.method, prediction.target_id)].append(prediction)
    return dict(grouped)


def evaluate_predictions(
    *,
    dataset: str,
    normalized_root: Path,
    predictions_tsv: Path,
    threshold_a: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    predictions = prediction_rows_from_tsv(predictions_tsv)
    grouped = group_predictions(predictions)
    per_site_rows: list[dict[str, Any]] = []
    per_target_rows: list[dict[str, Any]] = []

    for (method, target_id), target_predictions in sorted(grouped.items()):
        target = target_from_normalized_dir(dataset, normalized_root / dataset / target_id)
        evaluated = evaluate_target_standard(
            target_id=target_id,
            labels=target.labels,
            predictions=target_predictions,
            threshold_a=threshold_a,
        )
        top_n_hits = evaluated["top_n"]
        top_n_plus_2_hits = evaluated["top_n_plus_2"]
        all_sites_hits = evaluated["all_sites"]
        per_target_rows.append(
            {
                "dataset": dataset,
                "target_id": target_id,
                "method": method,
                "n_reference_sites": len(target.labels),
                "n_predictions": len(target_predictions),
                "top_n_success_rate": success_rate(top_n_hits),
                "top_n_plus_2_success_rate": success_rate(top_n_plus_2_hits),
                "all_sites_success_rate": success_rate(all_sites_hits),
            }
        )
        for protocol, hits in evaluated.items():
            for hit in hits:
                per_site_rows.append(
                    {
                        "dataset": dataset,
                        "target_id": target_id,
                        "method": method,
                        "protocol": protocol,
                        "true_site_id": hit.true_site_id,
                        "top_k": hit.top_k,
                        "hit": hit.hit,
                        "dca_a": hit.dca_a,
                        "matched_site_id": hit.matched_site_id,
                        "matched_rank": hit.matched_rank,
                    }
                )

    by_method: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in per_site_rows:
        by_method[str(row["method"])].append(row)
    target_rows_by_method: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in per_target_rows:
        target_rows_by_method[str(row["method"])].append(row)

    summary_rows: list[dict[str, Any]] = []
    for method, target_rows in sorted(target_rows_by_method.items()):
        rows = by_method.get(method, [])
        # Keep zero-reference structures visible for provenance, but exclude
        # them from site-level rates because they have no eligible ground truth.
        n_evaluable_targets = sum(int(row["n_reference_sites"]) > 0 for row in target_rows)
        for protocol in ("top_n", "top_n_plus_2", "all_sites"):
            protocol_rows = [row for row in rows if row["protocol"] == protocol]
            if not protocol_rows:
                continue
            summary_rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "protocol": protocol,
                    "threshold_a": threshold_a,
                    "n_targets": len(target_rows),
                    "n_evaluable_targets": n_evaluable_targets,
                    "n_reference_sites": len(protocol_rows),
                    "success_rate": sum(1 for row in protocol_rows if row["hit"]) / len(protocol_rows),
                }
            )
    return summary_rows, per_target_rows, per_site_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate normalized site-prediction predictions.tsv output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. coach420 or holo4k.")
    parser.add_argument(
        "--normalized-root",
        default=str(DEFAULT_NORMALIZED_ROOT),
        help="Root containing normalized/<dataset>/<target_id> directories.",
    )
    parser.add_argument("--predictions-tsv", required=True, help="Normalized predictions TSV.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Evaluation output dir.")
    parser.add_argument("--threshold", type=float, default=4.0, help="DCA hit threshold in Angstrom.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir).resolve()
    normalized_root = Path(args.normalized_root).resolve()
    predictions_tsv = Path(args.predictions_tsv).resolve()
    summary_rows, per_target_rows, per_site_rows = evaluate_predictions(
        dataset=args.dataset,
        normalized_root=normalized_root,
        predictions_tsv=predictions_tsv,
        threshold_a=float(args.threshold),
    )

    label_protocols: set[str] = set()
    for row in per_target_rows:
        metadata_path = normalized_root / args.dataset / str(row["target_id"]) / "metadata.json"
        if not metadata_path.is_file():
            continue
        target_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        protocol = str(target_metadata.get("label_protocol") or "").strip()
        if protocol:
            label_protocols.add(protocol)
    reference_site_counts = {
        str(row["target_id"]): int(row["n_reference_sites"])
        for row in per_target_rows
    }

    metadata = {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "dataset": args.dataset,
        "normalized_root": str(normalized_root),
        "predictions_tsv": str(predictions_tsv),
        "threshold_a": float(args.threshold),
        "distance_metric": "DCA: predicted center to nearest reference-ligand atom",
        "label_protocols": sorted(label_protocols),
        "n_targets": len(reference_site_counts),
        "n_evaluable_targets": sum(count > 0 for count in reference_site_counts.values()),
        "n_reference_sites": sum(reference_site_counts.values()),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_csv(
        output_dir / "summary_by_method.csv",
        summary_rows,
        [
            "dataset",
            "method",
            "protocol",
            "threshold_a",
            "n_targets",
            "n_evaluable_targets",
            "n_reference_sites",
            "success_rate",
        ],
    )
    _write_csv(
        output_dir / "per_target.csv",
        per_target_rows,
        [
            "dataset",
            "target_id",
            "method",
            "n_reference_sites",
            "n_predictions",
            "top_n_success_rate",
            "top_n_plus_2_success_rate",
            "all_sites_success_rate",
        ],
    )
    _write_csv(
        output_dir / "per_site.csv",
        per_site_rows,
        [
            "dataset",
            "target_id",
            "method",
            "protocol",
            "true_site_id",
            "top_k",
            "hit",
            "dca_a",
            "matched_site_id",
            "matched_rank",
        ],
    )
    print(f"Summary CSV: {output_dir / 'summary_by_method.csv'}")
    print(f"Per-target CSV: {output_dir / 'per_target.csv'}")
    print(f"Per-site CSV: {output_dir / 'per_site.csv'}")


if __name__ == "__main__":
    main()
