#!/usr/bin/env python3
"""Run the SERT escitalopram example using the shared helpers."""

from __future__ import annotations

import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXAMPLE_DIR.parent))

from common import run_pipeline, stage_inputs


def main() -> None:
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
    run_pipeline(workspace_paths, mode="cpu")


if __name__ == "__main__":
    main()