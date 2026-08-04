#!/usr/bin/env python3
"""Run a DUD-E binding-site recovery baseline with fpocket."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

BASELINE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BASELINE_ROOT.parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmarks" / "results" / "baseline" / "fpocket"

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

POCKET_RE = re.compile(r"pocket(\d+)", re.IGNORECASE)
SCORE_RE = re.compile(r"^\s*(?:score|drug\s*score|druggability\s*score)\s*[:=]\s*(-?\d+(?:\.\d+)?)", re.I)


def _pocket_number(path: Path) -> int:
    match = POCKET_RE.search(path.name)
    if match:
        return int(match.group(1))
    return 10**9


def _parse_coord_line(line: str) -> tuple[float, float, float] | None:
    if line.startswith(("ATOM", "HETATM")):
        try:
            return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            pass
    parts = line.split()
    if len(parts) >= 8:
        try:
            return (float(parts[5]), float(parts[6]), float(parts[7]))
        except ValueError:
            return None
    return None


def center_from_coord_file(path: Path) -> tuple[float, float, float] | None:
    coords: list[tuple[float, float, float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        coord = _parse_coord_line(line)
        if coord is not None:
            coords.append(coord)
    if not coords:
        return None
    count = float(len(coords))
    return (
        sum(coord[0] for coord in coords) / count,
        sum(coord[1] for coord in coords) / count,
        sum(coord[2] for coord in coords) / count,
    )


def parse_score(info_path: Path) -> float | None:
    if not info_path.exists():
        return None
    for line in info_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = SCORE_RE.match(line)
        if not match:
            continue
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None


def pocket_coord_candidates(pockets_dir: Path, pocket_number: int) -> list[Path]:
    candidates = [
        pockets_dir / f"pocket{pocket_number}_atm.pdb",
        pockets_dir / f"pocket{pocket_number}_vert.pqr",
    ]
    candidates.extend(sorted(pockets_dir.glob(f"pocket{pocket_number}_*.pdb")))
    candidates.extend(sorted(pockets_dir.glob(f"pocket{pocket_number}_*.pqr")))
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in candidates:
        if path in seen or not path.exists():
            continue
        seen.add(path)
        unique.append(path)
    return unique


def parse_fpocket_output(fpocket_output: Path) -> list[PredictedSite]:
    pockets_dir = fpocket_output / "pockets"
    if not pockets_dir.is_dir():
        raise FileNotFoundError(f"fpocket pockets directory not found: {pockets_dir}")

    pocket_numbers = sorted(
        {
            _pocket_number(path)
            for path in pockets_dir.glob("pocket*_*.*")
            if _pocket_number(path) < 10**9
        }
    )
    sites: list[tuple[float, int, PredictedSite]] = []
    for pocket_number in pocket_numbers:
        center = None
        source = ""
        for coord_path in pocket_coord_candidates(pockets_dir, pocket_number):
            center = center_from_coord_file(coord_path)
            if center is not None:
                source = str(coord_path.relative_to(fpocket_output))
                break
        if center is None:
            continue

        score = parse_score(pockets_dir / f"pocket{pocket_number}_info.txt")
        sort_score = score if score is not None else -float(pocket_number)
        sites.append(
            (
                sort_score,
                pocket_number,
                PredictedSite(
                    site_id=f"S{pocket_number}",
                    center=center,
                    score=score,
                    source=source,
                ),
            )
        )
    return [site for _, _, site in sorted(sites, key=lambda item: (-item[0], item[1]))]


def find_fpocket_output(raw_dir: Path, staged_receptor: Path) -> Path:
    expected = raw_dir / f"{staged_receptor.stem}_out"
    if expected.is_dir():
        return expected
    outputs = sorted(path for path in raw_dir.glob("*_out") if path.is_dir())
    if outputs:
        return outputs[0]
    raise FileNotFoundError(f"fpocket output directory not found under {raw_dir}")


def predict_fpocket(
    receptor_input: Path,
    target_output: Path,
    args: argparse.Namespace,
) -> PredictionResult:
    raw_dir = target_output / "raw" / "fpocket"
    if args.force and raw_dir.exists():
        shutil.rmtree(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    staged_receptor = prepare_method_receptor(
        receptor_input=receptor_input,
        target_output=target_output,
        method_work_dir=raw_dir,
        args=args,
    )
    command = [
        *split_tool_command(args.fpocket_cmd),
        "-f",
        staged_receptor.name,
        *args.fpocket_extra_arg,
    ]
    run_command(command, cwd=raw_dir)
    fpocket_output = find_fpocket_output(raw_dir, staged_receptor)
    sites = parse_fpocket_output(fpocket_output)
    return PredictionResult(
        sites=sites,
        raw_output_dir=fpocket_output,
        method_receptor=staged_receptor,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate fpocket site centers against DUD-E crystal-ligand centers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_args(parser, method_name="fpocket", default_output_dir=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--fpocket-cmd",
        default="fpocket",
        help="fpocket command. Use an absolute path if fpocket is not on PATH.",
    )
    parser.add_argument(
        "--fpocket-extra-arg",
        action="append",
        default=[],
        help="Extra argument appended to the fpocket command. Repeat for multiple args.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_baseline(method_name="fpocket", args=args, predictor=predict_fpocket)


if __name__ == "__main__":
    main()
