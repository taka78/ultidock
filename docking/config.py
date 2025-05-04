# config.py
# Auto-generated config.py
import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

LIGANDS_DIR = os.path.join(BASE_DIR, "LIGANDS_DIR")
DOCKING_DIR = os.path.join(BASE_DIR, "DOCKING_DIR")
ANALYSIS_DIR = os.path.join(BASE_DIR, "ANALYSIS_DIR")
VINA_DIR = os.path.join(BASE_DIR, "VINA_DIR")
VINA_GPU_DIR = os.path.join(BASE_DIR, "VINA_GPU_DIR")
MACRO_MOL_DIR = os.path.join(BASE_DIR, "MACRO_MOL_DIR")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
GPU_TYPE = "AMD"
DB_PATH = os.path.join(RESULTS_DIR, "ultidock_results.db")
