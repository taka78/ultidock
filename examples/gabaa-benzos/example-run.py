#!/usr/bin/env python3
import os
import sys
import shutil
import importlib
from pathlib import Path
import subprocess


# ---- paths ----
SCRIPT_DIR = Path(__file__).resolve().parent
print(f"script dir: {SCRIPT_DIR}")
ROOT_DIR = SCRIPT_DIR.parent.parent                                  # repo root
print(f"root dir: {ROOT_DIR}")
DOCKING_DIR = ROOT_DIR / "docking"
print(f"docking dir: {DOCKING_DIR}")
LIGANDS_DIR = DOCKING_DIR / "LIGANDS_DIR"
MACRO_MOL_DIR = DOCKING_DIR / "MACRO_MOL_DIR"
sys.path.append(DOCKING_DIR)

# example inputs (pick existing filenames)
RECEPTOR_SRC = SCRIPT_DIR /"4COF_edited.pdbqt"
print(f"looking for receptor PDBQT at {RECEPTOR_SRC}")
if not RECEPTOR_SRC.exists():
    alt = SCRIPT_DIR / "4COF.pdbqt"
    if alt.exists():
        RECEPTOR_SRC = alt
    else:
        raise FileNotFoundError("Receptor PDBQT not found (expected 4COF_edited.pdbqt or 4COF.pdbqt)")

LIGAND_SRC = SCRIPT_DIR /"aspirine-e.pdbqt"
if not LIGAND_SRC.exists():
    raise FileNotFoundError("Ligand PDBQT not found (expected aspirine-e.pdbqt)") #placeholder, needs to be replaced

# ---- make dirs ----
LIGANDS_DIR.mkdir(parents=True, exist_ok=True)
MACRO_MOL_DIR.mkdir(parents=True, exist_ok=True)

# ---- copy files into working dirs ----
shutil.copy2(str(RECEPTOR_SRC), str(MACRO_MOL_DIR / "4COF-edited.pdbqt"))
shutil.copy2(str(LIGAND_SRC),   str(LIGANDS_DIR   / "aspirine-e.pdbqt"))
shutil.copy2(str(LIGAND_SRC),   str(LIGANDS_DIR   / "ibuprofen-e.pdbqt"))
shutil.copy2(str(LIGAND_SRC),   str(LIGANDS_DIR   / "morphine-e.pdbqt"))
print(ROOT_DIR)
print(f"[ok] receptor → {MACRO_MOL_DIR/'4COF-edited.pdbqt'}")
print(f"[ok] ligand   → {LIGANDS_DIR/'aspirine-e.pdbqt'}")
print(f"[ok] ligand   → {LIGANDS_DIR/'ibuprofen-e.pdbqt'}")
print(f"[ok] ligand   → {LIGANDS_DIR/'morphine-e.pdbqt'}")

subprocess.run(["python3", "setup.py", "--example"], check=True, cwd=DOCKING_DIR)

subprocess.run(["python3", "dock_v02.py"], check=True, cwd=DOCKING_DIR)
subprocess.run(["python3", "analyse_docking_results.py"], check=True, cwd=DOCKING_DIR)