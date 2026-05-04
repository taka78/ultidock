#!/usr/bin/env python3
"""Scientific plots for DUD-E benchmark results.

Reads the summary.json produced by dude_docking_benchmark.py and generates:
  - Per-target ROC curves (from scores.csv)
  - Aggregate ROC-AUC bar chart
  - EF1% / EF5% bar chart
  - Summary table (printed + CSV)

Requires: matplotlib, numpy  (both are already used by the benchmark)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ── data loaders ──────────────────────────────────────────────────────────

def load_summary(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text("utf-8"))


def load_scores_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (labels, ranking_scores) from a scores.csv."""
    _LABEL_MAP = {"active": 1, "decoy": 0, "1": 1, "0": 0}
    labels: list[int] = []
    affinities: list[float] = []
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raw = row["label"].strip().lower()
            labels.append(_LABEL_MAP.get(raw, int(raw)))
            affinities.append(float(row["best_affinity_kcal_mol"]))
    lab = np.asarray(labels, dtype=int)
    aff = np.asarray(affinities, dtype=float)
    return lab, -aff  # more negative = better binding --> negate for ranking


# ── ROC helpers ───────────────────────────────────────────────────────────

def compute_roc_curve(
    labels: np.ndarray, scores: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return (fpr, tpr) arrays for plotting."""
    order = np.argsort(-scores, kind="mergesort")
    sorted_labels = labels[order]

    positives = labels.sum()
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])

    tps = np.cumsum(sorted_labels)
    fps = np.cumsum(1 - sorted_labels)

    tpr = np.concatenate([[0.0], tps / positives])
    fpr = np.concatenate([[0.0], fps / negatives])
    return fpr, tpr


# ── plotting ──────────────────────────────────────────────────────────────

def _apply_style(ax: "plt.Axes") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)


def plot_roc_per_target(
    summary: dict[str, Any], output_dir: Path, *, max_per_figure: int = 12
) -> list[Path]:
    """Plot ROC curves grouped into multi-panel figures."""
    targets = summary.get("targets", [])
    if not targets:
        return []

    paths: list[Path] = []
    for chunk_start in range(0, len(targets), max_per_figure):
        chunk = targets[chunk_start : chunk_start + max_per_figure]
        n = len(chunk)
        cols = min(n, 4)
        rows = (n + cols - 1) // cols

        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.5 * rows), squeeze=False)
        fig.suptitle("DUD-E ROC Curves", fontsize=14, fontweight="bold", y=1.01)

        for idx, entry in enumerate(chunk):
            r, c = divmod(idx, cols)
            ax = axes[r][c]
            name = entry["target"]
            auc_val = entry.get("roc_auc", float("nan"))

            csv_path = Path(str(entry.get("scores_csv", "")))
            if csv_path.exists():
                labels, scores = load_scores_csv(csv_path)
                fpr, tpr = compute_roc_curve(labels, scores)
                ax.plot(fpr, tpr, color="#2563eb", linewidth=1.5)
            ax.plot([0, 1], [0, 1], "--", color="#94a3b8", linewidth=0.8)
            ax.set_title(f"{name.upper()} (AUC={auc_val:.3f})", fontsize=10)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
            ax.set_xlabel("FPR", fontsize=8)
            ax.set_ylabel("TPR", fontsize=8)
            _apply_style(ax)

        # hide unused subplots
        for idx in range(n, rows * cols):
            r, c = divmod(idx, cols)
            axes[r][c].set_visible(False)

        fig.tight_layout()
        page = chunk_start // max_per_figure + 1
        out_path = output_dir / f"roc_curves_page{page}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(out_path)

    return paths


def plot_auc_bar_chart(summary: dict[str, Any], output_dir: Path) -> Path | None:
    """Horizontal bar chart of ROC-AUC per target, sorted."""
    targets = summary.get("targets", [])
    if not targets:
        return None

    names = [t["target"].upper() for t in targets]
    aucs = [float(t.get("roc_auc", 0)) for t in targets]

    # sort descending
    order = np.argsort(aucs)[::-1]
    names = [names[i] for i in order]
    aucs = [aucs[i] for i in order]

    h = max(4, len(names) * 0.35)
    fig, ax = plt.subplots(figsize=(8, h))
    y = np.arange(len(names))

    colors = ["#16a34a" if a >= 0.7 else "#eab308" if a >= 0.5 else "#dc2626" for a in aucs]
    ax.barh(y, aucs, color=colors, edgecolor="white", height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("ROC-AUC", fontsize=11)
    ax.set_title("DUD-E: ROC-AUC per Target", fontsize=13, fontweight="bold")
    ax.axvline(0.5, color="#94a3b8", linestyle="--", linewidth=0.8, label="random")
    ax.set_xlim(0, 1.05)
    ax.invert_yaxis()
    _apply_style(ax)
    ax.legend(fontsize=8)

    fig.tight_layout()
    out_path = output_dir / "roc_auc_bar.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_ef_bar_chart(summary: dict[str, Any], output_dir: Path) -> Path | None:
    """Grouped bar chart of EF1% and EF5% per target."""
    targets = summary.get("targets", [])
    if not targets:
        return None

    names = [t["target"].upper() for t in targets]
    ef1 = [float(t.get("ef1", 0)) for t in targets]
    ef5 = [float(t.get("ef5", 0)) for t in targets]

    # sort by EF1% descending
    order = np.argsort(ef1)[::-1]
    names = [names[i] for i in order]
    ef1 = [ef1[i] for i in order]
    ef5 = [ef5[i] for i in order]

    h = max(4, len(names) * 0.4)
    fig, ax = plt.subplots(figsize=(9, h))
    y = np.arange(len(names))
    bar_h = 0.35

    ax.barh(y - bar_h / 2, ef1, height=bar_h, color="#2563eb", label="EF 1%")
    ax.barh(y + bar_h / 2, ef5, height=bar_h, color="#7c3aed", label="EF 5%")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Enrichment Factor", fontsize=11)
    ax.set_title("DUD-E: Enrichment Factors per Target", fontsize=13, fontweight="bold")
    ax.axvline(1.0, color="#94a3b8", linestyle="--", linewidth=0.8, label="random (EF=1)")
    ax.invert_yaxis()
    ax.legend(fontsize=9)
    _apply_style(ax)

    fig.tight_layout()
    out_path = output_dir / "enrichment_factors.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_summary_table(summary: dict[str, Any], output_dir: Path) -> Path | None:
    """Render a pretty table figure of per-target metrics."""
    targets = summary.get("targets", [])
    if not targets:
        return None

    cols = ["Target", "Actives", "Decoys", "ROC-AUC", "EF 1%", "EF 5%"]
    has_bedroc = any("bedroc" in t for t in targets)
    if has_bedroc:
        cols.append("BEDROC")

    rows: list[list[str]] = []
    for t in targets:
        row = [
            t["target"].upper(),
            str(t.get("n_actives", "")),
            str(t.get("n_decoys", "")),
            f"{t.get('roc_auc', 0):.4f}",
            f"{t.get('ef1', 0):.2f}",
            f"{t.get('ef5', 0):.2f}",
        ]
        if has_bedroc:
            row.append(f"{t.get('bedroc', 0):.4f}")
        rows.append(row)

    # add aggregate row
    agg = summary.get("aggregate", {})
    agg_row = [
        "MEAN",
        "",
        "",
        f"{agg.get('roc_auc_mean', 0):.4f}",
        f"{agg.get('ef1_mean', 0):.2f}",
        f"{agg.get('ef5_mean', 0):.2f}",
    ]
    if has_bedroc:
        agg_row.append(f"{agg.get('bedroc_mean', 0):.4f}")
    rows.append(agg_row)

    # matplotlib table figure
    n_rows = len(rows) + 1  # +1 for header
    fig_h = max(2.5, 0.35 * n_rows)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.axis("off")

    table = ax.table(
        cellText=rows,
        colLabels=cols,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.3)

    # style header
    for j in range(len(cols)):
        cell = table[0, j]
        cell.set_facecolor("#1e293b")
        cell.set_text_props(color="white", fontweight="bold")

    # highlight mean row
    for j in range(len(cols)):
        cell = table[len(rows), j]
        cell.set_facecolor("#f1f5f9")
        cell.set_text_props(fontweight="bold")

    ax.set_title("DUD-E Benchmark Summary", fontsize=13, fontweight="bold", pad=15)
    fig.tight_layout()
    out_path = output_dir / "summary_table.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # also write CSV
    csv_path = output_dir / "summary_table.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        writer.writerows(rows)

    return out_path


# ── main ──────────────────────────────────────────────────────────────────

def generate_all_plots(summary_path: Path, output_dir: Path | None = None) -> list[Path]:
    """Generate all benchmark plots.  Returns list of created files."""
    if not HAS_MPL:
        print("[WARN] matplotlib not installed – skipping plot generation.")
        print("       Install with: pip install matplotlib")
        return []

    summary = load_summary(summary_path)
    if output_dir is None:
        output_dir = summary_path.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []

    rocs = plot_roc_per_target(summary, output_dir)
    created.extend(rocs)

    auc_bar = plot_auc_bar_chart(summary, output_dir)
    if auc_bar:
        created.append(auc_bar)

    ef_bar = plot_ef_bar_chart(summary, output_dir)
    if ef_bar:
        created.append(ef_bar)

    tbl = plot_summary_table(summary, output_dir)
    if tbl:
        created.append(tbl)

    for p in created:
        print(f"  --> {p}")
    return created


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate scientific plots from DUD-E benchmark results.",
    )
    parser.add_argument(
        "summary_json",
        help="Path to summary.json produced by dude_docking_benchmark.py",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for plot files (default: <summary_dir>/plots)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary_path = Path(args.summary_json).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else None

    if not summary_path.exists():
        print(f"ERROR: {summary_path} does not exist", file=sys.stderr)
        sys.exit(1)

    plots = generate_all_plots(summary_path, output_dir)
    print(f"\nGenerated {len(plots)} plot(s).")


if __name__ == "__main__":
    main()
