#!/usr/bin/env python3
"""Run a DUD-E binding-site recovery baseline with P2Rank."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path
from typing import Any

BASELINE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BASELINE_ROOT.parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmarks" / "results" / "baseline" / "p2rank"

if str(BASELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(BASELINE_ROOT))

from common import (  # noqa: E402
    PredictedSite,
    PredictionResult,
    add_common_args,
    prepare_method_receptor,
    run_baseline,
    run_command,
    split_tool_command,
)


def _float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def _row_value(row: dict[str, str], *names: str) -> str | None:
    lowered = {key.strip().lower(): value for key, value in row.items()}
    for name in names:
        value = lowered.get(name)
        if value not in (None, ""):
            return value
    return None


def parse_p2rank_predictions(csv_path: Path) -> list[PredictedSite]:
    sites: list[tuple[int, int, PredictedSite]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row_index, row in enumerate(reader, start=1):
            x = _float_or_none(_row_value(row, "center_x", "x", "cx"))
            y = _float_or_none(_row_value(row, "center_y", "y", "cy"))
            z = _float_or_none(_row_value(row, "center_z", "z", "cz"))
            if x is None or y is None or z is None:
                continue

            rank_value = _float_or_none(_row_value(row, "rank"))
            rank = int(rank_value) if rank_value is not None else row_index
            score = _float_or_none(_row_value(row, "score", "probability"))
            source_name = _row_value(row, "name", "pocket", "id") or csv_path.name
            sites.append(
                (
                    rank,
                    row_index,
                    PredictedSite(
                        site_id=f"S{row_index}",
                        center=(x, y, z),
                        score=score,
                        source=f"{csv_path.name}:{source_name}",
                    ),
                )
            )
    return [site for _, _, site in sorted(sites, key=lambda item: (item[0], item[1]))]


def find_predictions_csv(raw_dir: Path) -> Path:
    preferred = sorted(raw_dir.rglob("*_predictions.csv"))
    if preferred:
        return preferred[0]
    for path in sorted(raw_dir.rglob("*.csv")):
        try:
            if parse_p2rank_predictions(path):
                return path
        except Exception:
            continue
    raise FileNotFoundError(f"P2Rank predictions CSV not found under {raw_dir}")


def predict_p2rank(
    receptor_input: Path,
    target_output: Path,
    args: argparse.Namespace,
) -> PredictionResult:
    raw_dir = target_output / "raw" / "p2rank"
    if args.force and raw_dir.exists():
        shutil.rmtree(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    staged_receptor = prepare_method_receptor(
        receptor_input=receptor_input,
        target_output=target_output,
        method_work_dir=raw_dir / "input",
        args=args,
    )
    command = [
        *split_tool_command(args.p2rank_cmd),
        "-f",
        str(staged_receptor),
        "-o",
        str(raw_dir),
        *args.p2rank_extra_arg,
    ]
    run_command(command, cwd=raw_dir)
    predictions_csv = find_predictions_csv(raw_dir)
    sites = parse_p2rank_predictions(predictions_csv)
    return PredictionResult(sites=sites, raw_output_dir=raw_dir, method_receptor=staged_receptor)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate P2Rank site centers against DUD-E crystal-ligand centers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_args(parser, method_name="p2rank", default_output_dir=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--p2rank-cmd",
        default="prank predict",
        help="P2Rank command prefix. Use an absolute path if prank is not on PATH.",
    )
    parser.add_argument(
        "--p2rank-extra-arg",
        action="append",
        default=[],
        help="Extra argument appended to the P2Rank command. Repeat for multiple args.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_baseline(method_name="p2rank", args=args, predictor=predict_p2rank)


if __name__ == "__main__":
    main()
