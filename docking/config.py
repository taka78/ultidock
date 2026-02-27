# config.py
# Auto-generated config.py
import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

LIGANDS_DIR = '/home/taka/development/ultidock/docking/LIGANDS_DIR'
DOCKING_DIR = '/home/taka/development/ultidock/docking/DOCKING_DIR'
ANALYSIS_DIR = '/home/taka/development/ultidock/docking/ANALYSIS_DIR'
VINA_DIR = '/home/taka/development/ultidock/docking/VINA_DIR'
AUTODOCK_GPU_DIR = '/home/taka/development/ultidock/docking/AUTODOCK_GPU_DIR'
MACRO_MOL_DIR = '/home/taka/development/ultidock/docking/MACRO_MOL_DIR'
RESULTS_DIR = '/home/taka/development/ultidock/docking/RESULTS_DIR'
GPU_TYPE = 'CPU'
DB_PATH = '/home/taka/development/ultidock/docking/RESULTS_DIR/ultidock_results.db'
NUMWI = '128'
GRID_MODE = "centers"      # ligand | residues | centers | blind
GRID_SPACING = 0.375
GRID_MARGIN = 5.0         # Å
GRID_CAP = 150.0           # Å cap per axis for blind mode
AUTO_GRID_BIN = os.path.join(AUTODOCK_GPU_DIR, "autogrid", "autogrid4")
CENTERS_TSV  = os.path.join(MACRO_MOL_DIR, "centers.tsv")  # path or None
REF_LIGAND_PDB = None    # path to co-crystal/ref ligand if GRID_MODE="ligand"
HOTSPOT_NMS_MINSEP_A = 2.0
R_MIN_CAVITY_A = 3.0   # minimum inscribed sphere radius for cavity acceptance
SURFACE_SHELL__MIN_A = 2.0  # min/max distance from protein surface for surface pockets
SURFACE_SHELL__MAX_A = 20.0
SURFACE_NMS_MINSEP_A = 5        # voxels for non-max suppression of surface pockets
MAX_CENTER_DIST_A = 10.0         
CONTACT_SHELL_A = 4.0           # voxels ≤5 Å from surface count as “contact”
HOTSPOT_BOX_ANGLE = 35       # minimum box side length (Å)
MIN_SURFACE_FRAC = 0.01        # ~0.2% of box must be near-surface
AUTOSITES = 6
