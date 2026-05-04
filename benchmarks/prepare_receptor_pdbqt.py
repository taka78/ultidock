#!/usr/bin/env python3
"""Prepare one receptor as a deterministic PDBQT."""

from __future__ import annotations

import argparse
from pathlib import Path

from molguard.io.receptor_prep import prepare_receptor_pdbqt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a receptor into a canonical PDBQT for benchmarking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", help="Input receptor file (.pdbqt/.pdb/.mol2).")
    parser.add_argument("output", help="Output receptor PDBQT path.")
    parser.add_argument(
        "--prepare-command",
        help=(
            "External conversion command template. Use placeholders {input}, {output}, "
            "and optionally {seed}. Example: "
            "\"mk_prepare_receptor.py -i {input} -o {output}\""
        ),
    )
    parser.add_argument("--seed", type=int, default=0, help="Seed forwarded to the template.")
    parser.add_argument(
        "--timestamp",
        default="BENCHMARK",
        help="Stable canonicalization stamp embedded in the output REMARK.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    digest = prepare_receptor_pdbqt(
        input_path=Path(args.input),
        output_path=Path(args.output),
        prepare_command=args.prepare_command,
        seed=args.seed,
        timestamp=args.timestamp,
    )
    print(digest)


if __name__ == "__main__":
    main()
