"""Helpers for resumable benchmark-arm summaries."""

from __future__ import annotations

import json
import statistics
from pathlib import Path
from typing import Any, Iterable


SCREENING_METRIC_KEYS = ("roc_auc", "ef1", "ef5", "bedroc", "log_auc")


def method_output_paths(output_dir: Path, method_label: str) -> dict[str, Path]:
    """Return aggregate output paths for a benchmark arm."""
    return {
        "summary": output_dir / f"{method_label}__summary.json",
        "status": output_dir / f"{method_label}__status.json",
        "runtime_summary": output_dir / f"{method_label}__runtime_summary.json",
        "failure_summary": output_dir / f"{method_label}__failure_summary.json",
    }


def _json_load(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return data


def _numeric_value(raw: Any) -> float | None:
    if isinstance(raw, (int, float)):
        return float(raw)
    return None


def load_metrics_from_method_dir(method_dir: Path) -> dict[str, Any] | None:
    """Load a per-target metrics record, supporting legacy summary.json files."""
    metrics_path = method_dir / "metrics.json"
    if metrics_path.exists():
        return _json_load(metrics_path)

    summary_path = method_dir / "summary.json"
    if not summary_path.exists():
        return None

    summary = _json_load(summary_path)
    targets = summary.get("targets")
    if isinstance(targets, list) and targets:
        first = targets[0]
        if isinstance(first, dict):
            return dict(first)
    return None


def load_failure_from_method_dir(method_dir: Path) -> dict[str, Any] | None:
    failure_path = method_dir / "failure.json"
    if not failure_path.exists():
        return None
    return _json_load(failure_path)


def collect_method_state(
    *,
    output_dir: Path,
    method_label: str,
    target_names: Iterable[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    """Collect completed, failed, and pending targets for a benchmark arm."""
    completed: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    pending: list[str] = []

    seen: set[str] = set()
    for target_name in target_names:
        name = str(target_name).strip().lower()
        if not name or name in seen:
            continue
        seen.add(name)

        method_dir = output_dir / name / method_label
        metrics = load_metrics_from_method_dir(method_dir)
        if metrics is not None:
            metrics = dict(metrics)
            metrics.setdefault("target", name)
            completed.append(metrics)
            continue

        failure = load_failure_from_method_dir(method_dir)
        if failure is not None:
            failure = dict(failure)
            failure.setdefault("target", name)
            failed.append(failure)
            continue

        pending.append(name)

    completed.sort(key=lambda entry: str(entry.get("target", "")))
    failed.sort(key=lambda entry: str(entry.get("target", "")))
    pending.sort()
    return completed, failed, pending


def aggregate_screening_metrics(per_target: list[dict[str, Any]]) -> dict[str, float]:
    aggregate: dict[str, float] = {}
    for key in SCREENING_METRIC_KEYS:
        values = [
            numeric
            for entry in per_target
            if (numeric := _numeric_value(entry.get(key))) is not None
        ]
        if values:
            aggregate[f"{key}_mean"] = float(sum(values) / len(values))
    return aggregate


def aggregate_site_validation(per_target: list[dict[str, Any]]) -> dict[str, Any]:
    evaluations = [
        entry["site_evaluation"]
        for entry in per_target
        if isinstance(entry.get("site_evaluation"), dict)
    ]
    if not evaluations:
        return {}

    best_distances = [
        numeric
        for evaluation in evaluations
        if (numeric := _numeric_value(evaluation.get("best_distance_a"))) is not None
    ]
    best_site_orders = [
        numeric
        for evaluation in evaluations
        if (numeric := _numeric_value(evaluation.get("best_site_order"))) is not None
    ]
    any_site_hits = [
        bool(evaluation.get("any_site_success"))
        for evaluation in evaluations
        if isinstance(evaluation.get("any_site_success"), bool)
    ]

    summary: dict[str, Any] = {
        "n_targets_with_site_evaluation": len(evaluations),
    }
    if best_distances:
        summary["mean_best_distance_a"] = float(sum(best_distances) / len(best_distances))
        summary["median_best_distance_a"] = float(statistics.median(best_distances))
        summary["min_best_distance_a"] = float(min(best_distances))
        summary["max_best_distance_a"] = float(max(best_distances))
    if best_site_orders:
        summary["mean_best_site_order"] = float(sum(best_site_orders) / len(best_site_orders))
        summary["median_best_site_order"] = float(statistics.median(best_site_orders))
    if any_site_hits:
        summary["any_site_success_count"] = int(sum(any_site_hits))
        summary["any_site_success_rate"] = float(sum(any_site_hits) / len(any_site_hits))
    return summary


def build_runtime_summary(
    per_target: list[dict[str, Any]],
    failed_targets: list[dict[str, Any]],
) -> dict[str, Any]:
    success_times = [
        {"target": str(entry.get("target", "")), "status": "completed", "elapsed_s": numeric}
        for entry in per_target
        if (numeric := _numeric_value(entry.get("runner_elapsed_s"))) is not None
    ]
    failure_times = [
        {"target": str(entry.get("target", "")), "status": str(entry.get("status", "failed")), "elapsed_s": numeric}
        for entry in failed_targets
        if (numeric := _numeric_value(entry.get("elapsed_s"))) is not None
    ]
    observed = success_times + failure_times

    summary: dict[str, Any] = {
        "n_completed_with_runtime": len(success_times),
        "n_failed_with_runtime": len(failure_times),
        "n_targets_with_runtime": len(observed),
        "per_target_runtime_s": sorted(observed, key=lambda item: item["target"]),
    }
    if observed:
        all_values = [item["elapsed_s"] for item in observed]
        summary["total_observed_runtime_s"] = float(sum(all_values))
        summary["mean_observed_runtime_s"] = float(sum(all_values) / len(all_values))
        summary["median_observed_runtime_s"] = float(statistics.median(all_values))
    if success_times:
        values = [item["elapsed_s"] for item in success_times]
        summary["total_completed_runtime_s"] = float(sum(values))
        summary["mean_completed_runtime_s"] = float(sum(values) / len(values))
    if failure_times:
        values = [item["elapsed_s"] for item in failure_times]
        summary["total_failed_runtime_s"] = float(sum(values))
        summary["mean_failed_runtime_s"] = float(sum(values) / len(values))
    return summary


def build_failure_summary(failed_targets: list[dict[str, Any]]) -> dict[str, Any]:
    by_status: dict[str, int] = {}
    for entry in failed_targets:
        status = str(entry.get("status", "failed"))
        by_status[status] = by_status.get(status, 0) + 1

    failures = [
        {
            "target": str(entry.get("target", "")),
            "status": str(entry.get("status", "failed")),
            "reason": str(entry.get("reason", "unknown")),
            "elapsed_s": _numeric_value(entry.get("elapsed_s")),
            "log_path": entry.get("log_path"),
        }
        for entry in failed_targets
    ]
    failures.sort(key=lambda item: item["target"])

    return {
        "n_failed": len(failed_targets),
        "by_status": by_status,
        "failures": failures,
    }


def write_method_outputs(
    *,
    dataset_root: Path,
    output_dir: Path,
    method_label: str,
    target_names: Iterable[str],
    mode: str,
    center_source: str,
    seed: int,
    invocation_wall_time_s: float | None = None,
) -> dict[str, Any]:
    """Rebuild arm-level summaries from per-target artifacts already on disk."""
    completed, failed, pending = collect_method_state(
        output_dir=output_dir,
        method_label=method_label,
        target_names=target_names,
    )
    aggregate = aggregate_screening_metrics(completed)
    site_validation = aggregate_site_validation(completed)
    runtime_summary = build_runtime_summary(completed, failed)
    failure_summary = build_failure_summary(failed)

    tracked_targets = sorted({str(name).strip().lower() for name in target_names if str(name).strip()})
    summary: dict[str, Any] = {
        "dataset_root": str(dataset_root),
        "output_dir": str(output_dir),
        "method_label": method_label,
        "mode": mode,
        "center_source": center_source,
        "seed": seed,
        "n_targets": len(tracked_targets),
        "n_completed": len(completed),
        "n_failed": len(failed),
        "n_pending": len(pending),
        "failed": failure_summary["failures"],
        "targets": completed,
        "aggregate": aggregate,
        "site_validation": site_validation,
    }
    if invocation_wall_time_s is not None:
        summary["last_invocation_wall_time_s"] = float(invocation_wall_time_s)

    status = {
        "method_label": method_label,
        "mode": mode,
        "center_source": center_source,
        "n_targets": len(tracked_targets),
        "n_completed": len(completed),
        "n_failed": len(failed),
        "n_pending": len(pending),
        "completed_targets": [str(entry.get("target", "")) for entry in completed],
        "failed_targets": [str(entry.get("target", "")) for entry in failed],
        "pending_targets": pending,
    }

    paths = method_output_paths(output_dir, method_label)
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    paths["status"].write_text(json.dumps(status, indent=2), encoding="utf-8")
    paths["runtime_summary"].write_text(json.dumps(runtime_summary, indent=2), encoding="utf-8")
    paths["failure_summary"].write_text(json.dumps(failure_summary, indent=2), encoding="utf-8")

    return {
        "summary": summary,
        "status": status,
        "runtime_summary": runtime_summary,
        "failure_summary": failure_summary,
        "paths": paths,
    }
