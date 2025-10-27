#!/usr/bin/env python3
"""Run the GABAA benzodiazepine example with the shared setup pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXAMPLE_DIR.parent))

from common import run_pipeline, stage_inputs


def main() -> None:
    receptor = EXAMPLE_DIR / "4COF_edited.pdbqt"
    if not receptor.exists():
        fallback = EXAMPLE_DIR / "4COF.pdbqt"
        if fallback.exists():
            receptor = fallback
        else:
            raise FileNotFoundError("Expected 4COF_edited.pdbqt or 4COF.pdbqt in the example folder")

    ligands = [
        EXAMPLE_DIR / "aspirine-e.pdbqt",
        EXAMPLE_DIR / "ibuprofen-e.pdbqt",
        EXAMPLE_DIR / "morphine-e.pdbqt",
    ]

    for ligand in ligands:
        if not ligand.exists():
            raise FileNotFoundError(f"Missing ligand file: {ligand.name}")

    workspace_paths = stage_inputs(EXAMPLE_DIR, receptor, ligands)
    run_pipeline(workspace_paths, mode="gpu")


if __name__ == "__main__":
    main()