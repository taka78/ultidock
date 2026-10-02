#!/usr/bin/env python3
"""Run the SERT escitalopram example using the shared helpers."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXAMPLE_DIR.parent))

from common import run_pipeline, stage_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the SERT escitalopram docking example.")
    parser.add_argument(
        "site_method",
        nargs="?",
        choices=("cav-emps", "p2rank", "fpocket"),
        default="cav-emps",
        help="Method used to propose docking sites (default: cav-emps).",
    )
    parser.add_argument(
        "--mode",
        choices=("auto", "gpu", "cpu", "cuda", "opencl"),
        default="gpu",
        help="Docking backend (default: gpu; use cpu when no GPU runtime is available).",
    )
    args = parser.parse_args()

    receptor = EXAMPLE_DIR / "5i6x_edited.pdbqt"
    if not receptor.exists():
        fallback = EXAMPLE_DIR / "5i6x.pdbqt"
        if fallback.exists():
            receptor = fallback
        else:
            raise FileNotFoundError("Expected 5i6x_edited.pdbqt or 5i6x.pdbqt in the example folder")

    ligand = EXAMPLE_DIR / "escitalopram-e.pdbqt"
    if not ligand.exists():
        raise FileNotFoundError("Missing ligand file escitalopram-e.pdbqt")

    workspace_paths = stage_inputs(EXAMPLE_DIR, receptor, [ligand])
    run_pipeline(workspace_paths, mode=args.mode, site_method=args.site_method)


if __name__ == "__main__":
    main()
