#!/usr/bin/env python3
"""Prepare a batch of ligands as deterministic PDBQT files."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    from benchmarks.dude_prep import prepare_ligand_inputs
except ModuleNotFoundError:
    from dude_prep import prepare_ligand_inputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert ligands into normalized PDBQT files for benchmarking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Input ligand paths. Accepts directories plus .pdbqt or .mol2 files. "
            "Multi-record MOL2 files are split into one prepared ligand per record."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where prepared ligand PDBQT files will be written.",
    )
    parser.add_argument(
        "--prepare-command",
        help=(
            "External conversion command template. Use placeholders {input}, {output}, "
            "and optionally {seed}. Example: "
            "\"mk_prepare_ligand.py -i {input} -o {output}\""
        ),
    )
    parser.add_argument("--seed", type=int, default=0, help="Seed forwarded to the template.")
    parser.add_argument(
        "--manifest",
        help="Optional JSON output listing source -> prepared file mappings.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outputs = prepare_ligand_inputs(
        inputs=[Path(path) for path in args.inputs],
        output_dir=Path(args.output_dir),
        prepare_command=args.prepare_command,
        seed=args.seed,
        manifest_path=Path(args.manifest) if args.manifest else None,
    )
    print(f"prepared={len(outputs)}")


if __name__ == "__main__":
    main()
