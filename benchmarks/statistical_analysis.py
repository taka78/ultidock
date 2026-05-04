#!/usr/bin/env python3
"""Statistical comparison of DUD-E benchmark arms.

Usage:
  python3 benchmarks/statistical_analysis.py \\
    --summaries arm_a.json arm_b.json [arm_c.json ...] \\
    --labels ultidock crystal blind \\
    --output-dir benchmarks/results/dude/stats

Reads per-target metrics from each arm's summary.json, aligns targets,
and runs pairwise Wilcoxon signed-rank tests across all metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

try:
    from scipy.stats import wilcoxon
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

METRICS = ("roc_auc", "ef1", "ef5", "bedroc", "log_auc")


def load_arm(summary_path: Path) -> dict[str, dict[str, float]]:
    """Return {target_name: {metric: value}} for one arm."""
    data = json.loads(summary_path.read_text("utf-8"))
    result = {}
    for entry in data.get("targets", []):
        target = entry["target"]
        result[target] = {m: float(entry[m]) for m in METRICS if m in entry}
    return result


def pairwise_test(
    arm_a: dict[str, dict[str, float]],
    arm_b: dict[str, dict[str, float]],
    label_a: str,
    label_b: str,
) -> list[dict]:
    """Run Wilcoxon signed-rank on shared targets for each metric."""
    shared = sorted(set(arm_a) & set(arm_b))
    if len(shared) < 5:
        print(f"  [WARN] Only {len(shared)} shared targets between {label_a} and {label_b}; "
              "Wilcoxon requires ≥5 for meaningful results.")

    results = []
    for metric in METRICS:
        vals_a = [arm_a[t].get(metric) for t in shared]
        vals_b = [arm_b[t].get(metric) for t in shared]

        # Filter out Nones
        paired = [(a, b) for a, b in zip(vals_a, vals_b) if a is not None and b is not None]
        if len(paired) < 5:
            results.append({
                "comparison": f"{label_a} vs {label_b}",
                "metric": metric,
                "n_targets": len(paired),
                "mean_a": None,
                "mean_b": None,
                "delta": None,
                "statistic": None,
                "p_value": None,
                "significant_005": None,
                "note": "Insufficient paired samples (<5)"
            })
            continue

        a_vals = np.array([p[0] for p in paired])
        b_vals = np.array([p[1] for p in paired])

        # Wilcoxon requires at least one non-zero difference
        diffs = a_vals - b_vals
        if np.all(diffs == 0):
            stat, p_val = 0.0, 1.0
        elif HAS_SCIPY:
            stat, p_val = wilcoxon(a_vals, b_vals)
        else:
            stat, p_val = float("nan"), float("nan")

        results.append({
            "comparison": f"{label_a} vs {label_b}",
            "metric": metric,
            "n_targets": len(paired),
            "mean_a": float(np.mean(a_vals)),
            "mean_b": float(np.mean(b_vals)),
            "delta": float(np.mean(diffs)),
            "statistic": float(stat),
            "p_value": float(p_val),
            "significant_005": bool(p_val < 0.05),
        })
    return results


def print_results_table(all_results: list[dict]) -> None:
    """Pretty-print the comparison results."""
    header = f"{'Comparison':<30} {'Metric':<10} {'N':>3} {'Mean A':>8} {'Mean B':>8} {'Δ':>8} {'p-value':>10} {'Sig?':>5}"
    print(header)
    print("-" * len(header))
    for r in all_results:
        if r.get("note"):
            print(f"{r['comparison']:<30} {r['metric']:<10} {r['n_targets']:>3}  {'— insufficient data —'}")
            continue
        sig_star = "  *" if r["significant_005"] else ""
        print(
            f"{r['comparison']:<30} {r['metric']:<10} {r['n_targets']:>3} "
            f"{r['mean_a']:>8.4f} {r['mean_b']:>8.4f} {r['delta']:>+8.4f} "
            f"{r['p_value']:>10.4f}{sig_star}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Pairwise statistical comparison of DUD-E benchmark arms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--summaries",
        nargs="+",
        required=True,
        type=Path,
        help="Paths to summary.json files for each arm.",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        required=True,
        help="Human-readable labels for each arm (same order as --summaries).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/results/dude/stats"),
        help="Directory to write statistical comparison results.",
    )
    args = parser.parse_args()

    if len(args.summaries) != len(args.labels):
        sys.exit("Error: --summaries and --labels must have the same length.")
    if len(args.summaries) < 2:
        sys.exit("Error: Need at least 2 arms to compare.")
    if not HAS_SCIPY:
        print("WARNING: scipy not installed. p-values will be NaN.", file=sys.stderr)

    # Load all arms
    arms = {}
    for path, label in zip(args.summaries, args.labels):
        if not path.exists():
            sys.exit(f"Error: Summary not found at {path}")
        arms[label] = load_arm(path)
        print(f"Loaded {label}: {len(arms[label])} targets")

    # Pairwise comparisons
    all_results = []
    for (lab_a, arm_a), (lab_b, arm_b) in combinations(arms.items(), 2):
        print(f"\n{'='*60}")
        print(f"  {lab_a} vs {lab_b}")
        print(f"{'='*60}")
        results = pairwise_test(arm_a, arm_b, lab_a, lab_b)
        all_results.extend(results)
        print_results_table(results)

    # Write CSV
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "statistical_comparison.csv"
    fieldnames = ["comparison", "metric", "n_targets", "mean_a", "mean_b", "delta",
                  "statistic", "p_value", "significant_005", "note"]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nResults written to: {csv_path}")

    # Write JSON
    json_path = args.output_dir / "statistical_comparison.json"
    json_path.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
    print(f"JSON written to: {json_path}")


if __name__ == "__main__":
    main()
