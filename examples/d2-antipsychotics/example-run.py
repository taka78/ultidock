#!/usr/bin/env python3
"""Dock the bundled three-ligand comparison against the D2 receptor."""

import argparse
import json
from pathlib import Path
import sys

EXAMPLE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXAMPLE_DIR.parent))

from common import run_pipeline, stage_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("site_method", nargs="?", default="cav-emps",
                        choices=("cav-emps", "fpocket", "p2rank"))
    parser.add_argument("--mode", default="cpu", choices=("cpu", "gpu", "auto", "cuda", "opencl"))
    parser.add_argument("--dry-run", action="store_true", help="List inputs and settings without docking.")
    args = parser.parse_args()
    dataset = json.loads((EXAMPLE_DIR / "dataset.json").read_text())
    receptor = EXAMPLE_DIR / dataset["receptor"]
    ligands = [EXAMPLE_DIR / item["source"] for item in dataset["ligands"]]
    for path in [receptor, *ligands]:
        if not path.is_file():
            parser.error(f"Missing example input: {path.name}")
    print(f"D2 receptor: {receptor.name}")
    print("Ligands: " + ", ".join(path.name for path in ligands))
    print(f"Site method: {args.site_method}; engine mode: {args.mode}")
    if args.dry_run:
        print("Would stage only these inputs in a new example workspace, then run docking.")
        return
    paths = stage_inputs(EXAMPLE_DIR, receptor, ligands)
    run_pipeline(paths, mode=args.mode, site_method=args.site_method)
    print(f"Results: {paths['results']}")
    print(f"Best-pose PDBQT files: {paths['docking']}")


if __name__ == "__main__":
    main()
