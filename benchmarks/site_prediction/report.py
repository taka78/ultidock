#!/usr/bin/env python3
"""Generate Markdown/HTML summary pages for site-prediction benchmark outputs."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmarks" / "results" / "site_prediction" / "report"


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _table(rows: list[dict[str, Any]], columns: list[str], *, limit: int = 50) -> str:
    if not rows or not columns:
        return "_No rows found._\n"
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines) + "\n"


def _markdown_to_html(markdown: str) -> str:
    body = []
    for line in markdown.splitlines():
        if line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("| "):
            body.append(f"<p><code>{html.escape(line)}</code></p>")
        elif line.strip():
            body.append(f"<p>{html.escape(line)}</p>")
        else:
            body.append("")
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Site Prediction Benchmark Report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:1040px;margin:2rem auto;"
        "line-height:1.45}code{background:#f6f6f6;padding:.15rem .25rem}</style>"
        "</head><body>"
        + "\n".join(body)
        + "</body></html>\n"
    )


def _ablation_summary_rows(
    summary: list[dict[str, str]],
    per_target: list[dict[str, str]],
) -> list[dict[str, str]]:
    methods = sorted(
        {
            str(row.get("method") or "")
            for row in summary
            if str(row.get("method") or "").startswith("cav-emps-")
        }
    )
    rows = []
    for method in methods:
        rates = {
            str(row.get("protocol") or ""): str(row.get("success_rate") or "")
            for row in summary
            if row.get("method") == method
        }
        prediction_counts = [
            int(row.get("n_predictions") or 0)
            for row in per_target
            if row.get("method") == method and int(row.get("n_reference_sites") or 0) > 0
        ]
        rows.append(
            {
                "configuration": method.removeprefix("cav-emps-"),
                "method": method,
                "top_n": rates.get("top_n", ""),
                "top_n_plus_2": rates.get("top_n_plus_2", ""),
                "all_emitted_max_6": rates.get("all_sites", ""),
                "mean_predictions_per_target": (
                    f"{sum(prediction_counts) / len(prediction_counts):.4f}"
                    if prediction_counts
                    else ""
                ),
            }
        )
    return rows


def build_report(evaluation_dir: Path, output_dir: Path) -> dict[str, Path]:
    summary = _read_csv(evaluation_dir / "summary_by_method.csv")
    per_target = _read_csv(evaluation_dir / "per_target.csv")
    per_site = _read_csv(evaluation_dir / "per_site.csv")
    metadata_path = evaluation_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}

    failures = [
        row
        for row in per_site
        if str(row.get("hit", "")).strip().lower() in {"0", "false", "no"}
    ]
    ablation_summary = _ablation_summary_rows(summary, per_target)
    ablation_section = []
    if ablation_summary:
        ablation_section = [
            "## CaV-EMPS Ablation Summary",
            "",
            "All configurations for one target reuse the same AutoGrid C/e/d arrays.",
            "`all_emitted_max_6` is ranking-independent recall over every emitted center; "
            "the no-rescue configuration may emit fewer than six.",
            "",
            _table(
                ablation_summary,
                [
                    "configuration",
                    "top_n",
                    "top_n_plus_2",
                    "all_emitted_max_6",
                    "mean_predictions_per_target",
                ],
            ),
            "",
        ]
    markdown = [
        "# Site Prediction Benchmark Report",
        "",
        "## Framing",
        "",
        "This report compares binding-site proposal methods, not docking enrichment.",
        "`cav-emps` is CaV-EMPS: Cavity detection via Electrostatic Map Pocket Scoring.",
        "DCA is measured from predicted pocket centers to the nearest reference-ligand atom.",
        "`all_sites` is ranking-independent portfolio recall: it asks whether any emitted site recovers the reference site.",
        "",
        "## Metadata",
        "",
        _table([metadata], list(metadata.keys()) if metadata else []),
        "",
        "## Method Summary",
        "",
        _table(
            summary,
            [
                "dataset",
                "method",
                "protocol",
                "threshold_a",
                "n_targets",
                "n_evaluable_targets",
                "n_targets_with_predictions",
                "n_empty_targets",
                "n_reference_sites",
                "success_rate",
            ],
        ),
        "",
        *ablation_section,
        "## Per Target",
        "",
        _table(
            per_target,
            [
                "dataset",
                "target_id",
                "method",
                "prediction_status",
                "n_reference_sites",
                "n_predictions",
                "top_n_success_rate",
                "top_n_plus_2_success_rate",
                "all_sites_success_rate",
            ],
            limit=200,
        ),
        "",
        "## Failed Site Recoveries",
        "",
        _table(
            failures,
            ["dataset", "target_id", "method", "protocol", "true_site_id", "top_k", "dca_a", "matched_site_id", "matched_rank"],
            limit=200,
        ),
        "",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    if ablation_summary:
        with (output_dir / "ablation_summary.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(ablation_summary[0]))
            writer.writeheader()
            writer.writerows(ablation_summary)
    report_md = output_dir / "site_benchmark_report.md"
    report_html = output_dir / "site_benchmark_report.html"
    report_md.write_text("\n".join(markdown), encoding="utf-8")
    report_html.write_text(_markdown_to_html(report_md.read_text(encoding="utf-8")), encoding="utf-8")
    return {"report_md": report_md, "report_html": report_html}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a site-prediction benchmark report from evaluator CSV outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--evaluation-dir", required=True, type=Path)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outputs = build_report(args.evaluation_dir.resolve(), args.output_dir.resolve())
    print(f"Markdown: {outputs['report_md']}")
    print(f"HTML:     {outputs['report_html']}")


if __name__ == "__main__":
    main()
