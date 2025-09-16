#dock_v02.py
import os, stat
import threading
import numpy as np
from threading import Semaphore
import queue
import glob
import subprocess
import time
from pathlib import Path
import psutil
import uuid
import gc
import xml.etree.ElementTree as ET
import shutil  
from threading import Barrier
from make_grids import (
    HotspotGPFGenerator,
    autogenerate_centers_tsv,
    ensure_grids_multi_centers,
    ensure_whole_protein_maps,
)
from config import (
    LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, AUTODOCK_GPU_DIR, MACRO_MOL_DIR,
    DB_PATH, GPU_TYPE, RESULTS_DIR, NUMWI, GRID_MODE, GRID_MARGIN, GRID_CAP, CENTERS_TSV, REF_LIGAND_PDB,
    GRID_SPACING, AUTO_GRID_BIN, AUTOSITES, R_MIN_CAVITY_A
)
from db_manager import DockingDatabaseManager

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
BINARY_PATH = os.path.join(AUTODOCK_GPU_DIR, "bin", "autodock_gpu_" + NUMWI + "wi")
COMPILER_SCRIPT = os.path.join(os.path.dirname(__file__), "autodock-gpu-compiler.sh")

BACKEND = (GPU_TYPE or "CPU").upper()
IS_CUDA = BACKEND == "CUDA"
IS_OPENCL = BACKEND == "OPENCL"
IS_GPU = IS_CUDA or IS_OPENCL

def list_nvidia_gpus():
    """Return a list of GPU indices [0,1,...]. Falls back to [0] if unknown."""
    if shutil.which("nvidia-smi") is None:
        return [0]
    try:
        out = subprocess.check_output(["nvidia-smi", "-L"], text=True)
        idxs = []
        for line in out.strip().splitlines():
            # lines look like: "GPU 0: NVIDIA GeForce RTX 5070 (UUID: ...)"
            if line.startswith("GPU "):
                num = line.split()[1].rstrip(":")
                idxs.append(int(num))
        return idxs or [0]
    except Exception:
        return [0]

GPU_IDS = list_nvidia_gpus()
# BoundedSemaphore protects against accidental double-release
GPU_SLOTS_PER_DEV = int(os.environ.get("GPU_SLOTS_PER_DEV", "1"))
GPU_SEMAPHORES = {i: threading.BoundedSemaphore(GPU_SLOTS_PER_DEV) for i in GPU_IDS}
_GPU_RR_LOCK = threading.Lock()
_GPU_RR = 0

def next_gpu_id():
    """Round-robin GPU selection; safe even if GPU_IDS == [0] fallback."""
    global _GPU_RR
    with _GPU_RR_LOCK:
        if not GPU_IDS:
            # extremely defensive; shouldn't happen because list_nvidia_gpus returns [0] on failure
            return 0
        gid = GPU_IDS[_GPU_RR % len(GPU_IDS)]
        _GPU_RR = (_GPU_RR + 1) % (10_000_000)  # avoid unbounded growth
        return gid

def prepare_sites_for_docking(receptor_pdbqt: str, macro_dir: str):
    """
    Returns: list[{
      'site_id': 'S1', 'center': (cx,cy,cz), 'npts': (nx,ny,nz),
      'spacing': float, 'fld_path': str, 'out_dir': str
    }]
    """
    macro_dir = Path(macro_dir); macro_dir.mkdir(parents=True, exist_ok=True)
    centers_tsv_path = Path(CENTERS_TSV)

    # If user provided centers.tsv (GRID_MODE == "centers"), just build per-site maps.
    if (GRID_MODE or "").lower() == "centers" and centers_tsv_path.exists():
        sites = ensure_grids_multi_centers(
            receptor_pdbqt=receptor_pdbqt,
            out_root=str(macro_dir),
            centers_tsv=str(centers_tsv_path),
            autogrid4_bin=str(Path(AUTODOCK_GPU_DIR) / "autogrid" / "autogrid4"),
        )
        return sites

    # Otherwise, auto-detect with your old algorithm and write centers.tsv, then per-site maps.
    # Hybrid = internal cavities first, then C/E/D hotspots as fallback.
    generator = HotspotGPFGenerator(
        autogrid4_bin=str(Path(AUTODOCK_GPU_DIR) / "autogrid" / "autogrid4")
    )
    sites = generator.prepare_centers_and_grids(
        receptor_pdbqt=receptor_pdbqt,
        out_root=str(macro_dir),
        centers_tsv_path=str(centers_tsv_path),
        mode="hybrid",               # "internal" | "maps" | "hybrid"
        n_sites=6,                   # or tune
        whole_spacing=float(GRID_SPACING),
        whole_cap_ang=float(GRID_CAP),      # 30–100 Å per axis cap for whole maps
        hotspot_box_ang=18.0,        # per-site half box for maps mode
        tau_rel=0.58,                # hotspot relative threshold
        min_sep_A=7.0,               # NMS min separation
    )
    return sites

class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        print(f"[INIT] ligands discovered: {len(self.FILES)} in {LIGANDS_DIR}")

        # === receptor selection using your variables ===
        macro_candidates = sorted(glob.glob(os.path.join(MACRO_MOL_DIR, "*.pdbqt")))
        if not macro_candidates:
            raise FileNotFoundError(f"No receptor .pdbqt found in {MACRO_MOL_DIR}")

        # keep your naming style
        self.MACRO_MOL_DIR = MACRO_MOL_DIR
        self.MACRO_MOL = macro_candidates[0]                  # <— path to receptor file
        self.RECEPTOR_STEM = Path(self.MACRO_MOL).stem

        # grids land under MACRO_MOL_DIR/<receptor>/grids
        self.GRID_ROOT = MACRO_MOL_DIR
        os.makedirs(self.GRID_ROOT, exist_ok=True)

        # optional: only if you ever use residues mode
        def _RES_PRED(line: str) -> bool:
            chain = line[21].strip()
            resi  = line[22:26].strip()
            resn  = line[17:20].strip()
            try: idx = int(resi)
            except: return False
            # example: orthosteric pocket on chain A; include key residues you trust
            return (chain == "A" and 120 <= idx <= 145) or (resn in {"HIS","ASP","SER"} and chain == "A")


        # Build/ensure whole-protein maps and per-site maps in-place under MACRO_MOL_DIR
        ref_lig = REF_LIGAND_PDB if (REF_LIGAND_PDB and os.path.exists(REF_LIGAND_PDB) and REF_LIGAND_PDB.lower().endswith(".pdbqt")) else None
        self.SITES = prepare_sites_for_docking(
            receptor_pdbqt=self.MACRO_MOL,
            macro_dir=self.MACRO_MOL_DIR
        )
        for s in self.SITES:
            fld = s["fld_path"]
            if not os.path.exists(fld):
                raise FileNotFoundError(f"[preflight] Missing FLD: {fld}\nLog: {os.path.join(s['out_dir'], 'grid.glg')}")

                    # normalize to the same structure used by centers-mode

        # NOTE: barrier/event were not used; keeping them None avoids accidental waits later
        self.barrier = None
        self.event = None

        # Initialize database manager and queue for async writes
        cpu_count = os.cpu_count() or 1
        cpus_per_vina = 2
        max_parallel = max(1, cpu_count // cpus_per_vina)
        print(f"[INIT] allowing up to {max_parallel} concurrent Vina runs")
        self.vina_sem = Semaphore(max_parallel)

        self.db_manager = DockingDatabaseManager(DB_PATH)
        self.db_queue = queue.Queue()
        self.db_thread = threading.Thread(target=self._db_worker, daemon=True)
        self.db_thread.start()


    def _db_worker(self):
        print("DB worker started")
        buffer = []
        seen_none = False

        while True:
            try:
                print("DB worker waiting for item...")
                item = self.db_queue.get(timeout=10)
                print(f"DB worker received: {item}")

                if item == "INIT":
                    print("Received warm-up INIT. Touching DB...")
                    self.db_manager.insert_bulk([])  # Safe no-op
                    self.db_queue.task_done()
                    continue

                elif item == "SHUTDOWN":
                    print("Shutdown signal received. Flushing and exiting DB thread.")
                    if buffer:
                        self.db_manager.insert_bulk(buffer)
                    self.db_queue.task_done()
                    break
                # Item is valid
                added = 0
                if isinstance(item, list):
                    for rec in item:
                        if isinstance(rec, tuple) and len(rec) == 5:
                            buffer.append(rec); added += 1
                        else:
                            print(f"[DB] skipping bad sub-record: {rec}")
                elif isinstance(item, tuple) and len(item) == 5:
                    buffer.append(item); added = 1
                else:
                    print(f"[DB] skipping malformed record: {item}")

                self.db_queue.task_done()  # Very important!

                if len(buffer) >= 500:
                    try:
                        print(f"[DB] flushing {len(buffer)} records")
                        self.db_manager.insert_bulk(buffer)
                        buffer.clear()
                    except Exception as e:
                        # don't die; log and continue
                        print(f"[DB] insert_bulk failed on batch of {len(buffer)}: {e}")
                        buffer.clear()

            except queue.Empty:
                if buffer:
                    try:
                        print(f"[DB] timeout flush: writing {len(buffer)} records")
                        self.db_manager.insert_bulk(buffer)
                    except Exception as e:
                        print(f"[DB] timeout flush failed on {len(buffer)} records: {e}")
                    buffer.clear()


    def memory_monitor(self, threshold):
        while True:
            available_memory = psutil.virtual_memory().available
            if available_memory > threshold:
                break
            time.sleep(1)

    def process_files(self):
        threads = []
        completed_threads = 0
        lock = threading.Lock()

        def thread_callback():
            nonlocal completed_threads
            with lock:
                completed_threads += 1

        for ligand_file in self.FILES:
            thread = ProcessFileThread(
                [ligand_file],
                self.barrier,
                self.event,
                thread_callback,
                self.db_queue,
                self.vina_sem,
                gpu_id= next_gpu_id() 
            )
            thread._parent = self   
            threads.append(thread)

        # 1 KB threshold is effectively always true; keep but log once
        threshold = 1024
        if psutil.virtual_memory().available <= threshold:
            print("[WARN] extremely low available memory reported; waiting...")
            while psutil.virtual_memory().available <= threshold:
                time.sleep(1)

        for thread in threads:
            thread.start()
            print(f"[RUN] launched {len(threads)} worker threads")


        while completed_threads < len(threads):
            time.sleep(1)

        time.sleep(2)  # Give threads time to finish putting into the queue

        self.db_queue.put("INIT")
        self.db_queue.put("SHUTDOWN")
        self.db_queue.join()           # Wait for all DB queue tasks to complete
        self.db_thread.join()          # Wait for the DB thread to fully exit
        self.db_manager.close()        # Only now it's safe to close DB

        gc.collect()
        self.check_memory()

    def check_memory(self):
        bunch = []
        for i in range(0, len(self.FILES), 6):
            current_bunch = self.FILES[i:i+6]
            bunch.extend(current_bunch)
            memory_info = psutil.Process().memory_info()
            print(f"Memory usage: {memory_info.rss / (1024 * 1024)} MB")
            break  # ← prevent barrier from ever being reached


    def run(self):
        self.process_files()


'''
###this is a placeholder for the GPFGenerator class, which is not used in the current version of the script.
    def extract_atom_types_from_pdbqt(self, ligand_file):
        atom_types = set()
        try:
            with open(ligand_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_type = line[77:79].strip()  # AutoDock-style atom type
                        if atom_type:
                            atom_types.add(atom_type)
        except Exception as e:
            print(f"[ERROR] Failed to extract atom types: {e}")
        return " ".join(sorted(atom_types))
'''


class ProcessFileThread(threading.Thread):
    
    def __init__(self, bunch, barrier, event, callback, db_queue, vina_sem, gpu_id=0):
        super().__init__()
        self.bunch = bunch
        self.barrier = barrier
        self.event = event
        self.callback = callback
        self.db_queue = db_queue
        self.vina_sem = vina_sem
        self.gpu_id = gpu_id
        self.MACRO_MOL_DIR = MACRO_MOL_DIR

    '''LEGACY CODE
    # Function to parse ligand's .pdbqt file and extract atomic coordinates
    def calculate_grid_center_and_size(self, ligand_file):
        x_coords, y_coords, z_coords = [], [], []
        
        with open(ligand_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])  # Extract x, y, z coordinates
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
        
        # Calculate geometric center (mean of x, y, z coordinates)
        center_x = np.mean(x_coords)
        center_y = np.mean(y_coords)
        center_z = np.mean(z_coords)
        
        # Calculate grid size (based on the bounding box)
        size_x = np.max(x_coords) - np.min(x_coords) + 10  # Add padding to grid size
        size_y = np.max(y_coords) - np.min(y_coords) + 10
        size_z = np.max(z_coords) - np.min(z_coords) + 10
        
        return (center_x, center_y, center_z), (size_x, size_y, size_z)
    '''
            
    def parse_vina_output_file(self, filepath, receptor_name, ligand_file):
        # Wait briefly for filesystem
        for _ in range(3):
            if os.path.exists(filepath):
                break
            time.sleep(0.05)
        else:
            print(f"parse_vina_output_file: File not found after retries: {filepath}")
            return []

        results = []
        current_model = None
        current_affinity = None
        current_rmsd_lb = None
        current_rmsd_ub = None
        current_ligand_id = None
        inside_model_block = False

        def finalize():
            nonlocal current_model, current_affinity, current_rmsd_lb, current_rmsd_ub, current_ligand_id, inside_model_block
            if inside_model_block and (current_model is not None) and (current_affinity is not None):
                lid = current_ligand_id or Path(ligand_file).stem or "ligand"
                tag = f"{lid}-Model{current_model}-{receptor_name}"
                results.append((tag, current_affinity, current_rmsd_lb, current_rmsd_ub, ligand_file))
            # reset
            current_model = None
            current_affinity = None
            current_rmsd_lb = None
            current_rmsd_ub = None
            current_ligand_id = None
            inside_model_block = False

        with open(filepath, "r", encoding="utf-8", errors="ignore") as fp:
            for line in fp:
                line = line.strip()

                if line.startswith("MODEL"):
                    # start new model (do not finalize here; we finalize on ENDMDL)
                    inside_model_block = True
                    current_affinity = None
                    current_rmsd_lb = None
                    current_rmsd_ub = None
                    current_ligand_id = None
                    try:
                        current_model = int(line.split()[1])
                    except (IndexError, ValueError):
                        current_model = None

                elif line.startswith("REMARK VINA RESULT:") and inside_model_block:
                    parts = line.split()
                    if len(parts) >= 6 and parts[0:3] == ['REMARK', 'VINA', 'RESULT:']:
                        try:
                            current_affinity = float(parts[3])
                            current_rmsd_lb = float(parts[4])
                            current_rmsd_ub = float(parts[5])
                        except ValueError:
                            pass

                elif line.startswith("REMARK  Name") and inside_model_block:
                    # Some writers use: "REMARK  Name = LIGAND"
                    parts = line.split()
                    if len(parts) >= 4:
                        current_ligand_id = parts[-1]  # safer than [3]

                elif line.startswith("ENDMDL") and inside_model_block:
                    finalize()

        # EOF without ENDMDL: finalize any open block
        if inside_model_block:
            finalize()

        return results


    def parse_adgpu_xml(self, xml_path, receptor_name, ligand_file_fallback):
        """
        Parse AutoDock-GPU XML and return a list of tuples:
        (tag, affinity, rmsd_lb, rmsd_ub, ligand_file)
        tag format: <ligand_id>-Model<run_id>-<receptor_name>
        """
        results = []
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()

            # ligand id from <ligand> element (fallback to provided path)
            ligand_path = root.findtext("ligand") or ligand_file_fallback
            ligand_id = Path(ligand_path).stem

            runs = root.find("runs")
            if runs is None:
                return results

            for run in runs.findall("run"):
                # run id becomes the "Model" number in your schema
                run_id_text = run.get("id")
                try:
                    run_id = int(run_id_text) if run_id_text is not None else 1
                except ValueError:
                    run_id = 1

                # binding energy
                e_txt = run.findtext("free_NRG_binding")
                if e_txt is None:
                    # skip malformed entries
                    continue
                try:
                    affinity = float(e_txt)
                except ValueError:
                    continue

                tag = f"{ligand_id}-Model{run_id}-{receptor_name}"
                # AD-GPU XML doesn't include RMSD values
                results.append((tag, affinity, None, None, ligand_file_fallback))

        except Exception as e:
            print(f"[XML-parse] Failed on {xml_path}: {e}")

        return results


    def run(self):
        # Calculate grid center and size dynamically for each ligand
        parent = getattr(self, "_parent", None)
        if parent is None:
            print("[WARN] Missing parent reference; cannot access precomputed grids.")
            return
        # use sites prepared in DockingProcessor.__init__
        sites = parent.SITES
        receptor_name = Path(parent.MACRO_MOL).stem
        try:
            buffer = []
            BATCH_SIZE = 500
            for ligand_file in self.bunch:
                # ensure output dir exists
                Path(DOCKING_DIR).mkdir(parents=True, exist_ok=True)
                lig_stem = Path(ligand_file).stem
                for receptor_path in sorted(Path(MACRO_MOL_DIR).glob("*.pdbqt")):
                    receptor_pdbqt = str(receptor_path)
                    sites = prepare_sites_for_docking(receptor_pdbqt, macro_dir=MACRO_MOL_DIR)
                for site in sites:
                    cx, cy, cz = site["center"]
                    nx, ny, nz = site["npts"]
                    sp = float(site["spacing"])
                    fld_base = site["fld_path"]                 # <-- pass this to AD-GPU --ffile

                    size_x = (nx - 1) * sp                      # <-- Vina box from npts/spacing
                    size_y = (ny - 1) * sp
                    size_z = (nz - 1) * sp
                    site_id = site["site_id"]
                    center  = tuple(float(x) for x in site["center"])
                    spacing = float(site["spacing"])
                    npts    = np.array(site["npts"], dtype=float)
                    box_xyz = tuple((npts * spacing).tolist())
                    fld_file= site["fld_path"]

                    out_stem  = Path(DOCKING_DIR) / f"{site_id}__{lig_stem}_{uuid.uuid4().hex[:8]}"
                    gpu_pdbqt = f"{out_stem}_out.pdbqt"
                    xml_out   = f"{out_stem}.xml"
                    vina_out  = f"{out_stem}.pdbqt"

                    # Per-GPU semaphore (one job per GPU at a time)
                    sem = GPU_SEMAPHORES.get(self.gpu_id, next(iter(GPU_SEMAPHORES.values())))
                    env = os.environ.copy()
                    env["OMP_NUM_THREADS"] = "1"
                    env["MKL_NUM_THREADS"] = "1"
                    env["OPENBLAS_NUM_THREADS"] = "1"
                    site_dir = os.path.dirname(fld_file)              # <-- folder with FLD + .map files
                    fld_base = os.path.basename(fld_file)
                    print("Using grid:", fld_file)
                    # If you need hard isolation per thread, uncomment:
                    # env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)
                    acquired_gpu = False
                    acquired_cpu = False
                    print(ligand_file)

                    try:
                        if GPU_TYPE == "NVIDIA" or GPU_TYPE == "CUDA":
                            sem.acquire()
                            acquired_gpu = True
                            result = subprocess.run(
                                [
                                    f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_cuda_{NUMWI}wi",
                                    "--lfile",   str(Path(ligand_file).resolve()),
                                    "--ffile",   fld_base,
                                    "--nrun",    "512",
                                    "--nev",   "5000000",            # no early stopping
                                    "--gbest",   "5",            # write <resnam>_out.pdbqt (optional but nice)
                                    "--xmloutput","1",           # ensure XML is produced
                                    "--resnam",  str(out_stem),  # basename; outputs land in DOCKING_DIR
                                    "--devnum",  str(int(self.gpu_id) + 1),  # AD-GPU is 1-indexed
                                    "--lsrat", "80.0",         # increase local search rate
                                    "--lsmet", "sw",       # use Solis-Wets local search
                                    "--initswgens", "300",  # More initial poses
                                    "--stopstd", "0.02",     # Tighter convergence
                                    "--psize", "300",      # Population size
                                    "--dang", "30",         # Torsional angle change (degrees)
                                    "--dmov","1.0",         # Torsional step size (angstroms)
                                    "--mrat", "5",         # Mutation rate
                                    "--autostop", "0",        # Disable autostop
                                    "--lsit", "1000"        # Increase max iterations for local search
                                ],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                cwd=site_dir,
                                env=env,
                                timeout=10000
                            )
                            print(f"[AutoDock-GPU stdout]\n{result.stdout}")
                            print(f"[AutoDock-GPU stderr]\n{result.stderr}")

                            if result.returncode != 0:
                                print(f"Docking failed for: {ligand_file}")
                                print(f"STDOUT:\n{result.stdout.strip()}")
                                print(f"STDERR:\n{result.stderr.strip()}")
                                continue
                            # Ensure output files are created
                        
                        elif GPU_TYPE == "AMD" or GPU_TYPE == "OPENCL":
                            sem.acquire()
                            acquired_gpu = True
                            # --- AutoDock-GPU (XML-first) ---
                            result = subprocess.run(
                                [
                                    f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_ocl_{NUMWI}wi",
                                    "--lfile",   str(Path(ligand_file).resolve()),
                                    "--ffile",   fld_base,
                                    "--nrun",    "512",
                                    "--gbest",   "5",            # write <resnam>_out.pdbqt (optional but nice)
                                    "--xmloutput","1",           # ensure XML is produced
                                    "--resnam",  str(out_stem),  # basename; outputs land in DOCKING_DIR
                                    "--devnum",  str(int(self.gpu_id) + 1),  # AD-GPU is 1-indexed
                                    "--lsrat", "50.0",         # increase local search rate
                                    "--lsmet", "sw"       # use Solis-Wets local search
                                ],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                cwd=site_dir,
                                env=env,
                                timeout=10000
                            )
                            print(f"[AutoDock-GPU stdout]\n{result.stdout}")
                            print(f"[AutoDock-GPU stderr]\n{result.stderr}")

                            if result.returncode != 0:
                                print(f"Docking failed for: {ligand_file}")
                                print(f"STDOUT:\n{result.stdout.strip()}")
                                print(f"STDERR:\n{result.stderr.strip()}")
                                continue
                            # Ensure output files are created

                        elif GPU_TYPE == "CPU" or GPU_TYPE == "VINA": # CPU/Vina fallback, i mean, why not?
                            # --- Vina CPU fallback (PDBQT parsing) ---
                            # CPU/Vina path uses ONLY the CPU limiter
                            self.vina_sem.acquire()
                            acquired_cpu = True
                            print(f"[VINA/{site_id}] center:", tuple(f"{v:.3f}" for v in center),
                                " size:", tuple(f"{v:.3f}" for v in box_xyz))
                            result = subprocess.run(
                                [
                                    f"{VINA_DIR}/bin/vina",
                                    "--receptor", str(parent.MACRO_MOL),
                                    "--ligand",   str(ligand_file),
                                    "--center_x", f"{cx:.3f}",
                                    "--center_y", f"{cy:.3f}",
                                    "--center_z", f"{cz:.3f}",
                                    "--size_x",   f"{size_x:.3f}",
                                    "--size_y",   f"{size_y:.3f}",
                                    "--size_z",   f"{size_z:.3f}",
                                    "--cpu", "2",
                                    "--out",      str(vina_out)
                                ],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                timeout=12000
                            )
                            if result.returncode != 0:
                                print(f"Vina failed")   
                                print(f"STDOUT:\n{result.stdout.strip()}")
                                print(f"STDERR:\n{result.stderr.strip()}")
                                continue

                        # Ensure output files are created
                    except subprocess.TimeoutExpired:
                        print(f"Docking timed out for: {ligand_file}")
                        sem.release()
                        self.vina_sem.release()
                        continue
                    finally:
                        if acquired_gpu:
                            sem.release()
                        if acquired_cpu:
                            self.vina_sem.release()

                    print(f"[{threading.current_thread().name}] Docking {ligand_file}")
                    if GPU_TYPE in ("NVIDIA", "AMD", "OPENCL", "CUDA"):
                        parsed_results = self.parse_adgpu_xml(xml_out, receptor_name, ligand_file)
                    else:
                        parsed_results = self.parse_vina_output_file(vina_out, receptor_name, ligand_file)

                    # parsed_results now unified shape: (tag, affinity, rmsd_lb, rmsd_ub, ligand_file)
                    if not parsed_results:
                        print(f"No valid docking data for: {ligand_file}")
                        continue
                    buffer.extend(parsed_results)


                    if len(buffer) >= BATCH_SIZE:
                        self.db_queue.put(buffer.copy())
                        buffer.clear()

                    # memory log
                    memory_info = psutil.Process().memory_info()
                    print(f"Memory usage: {memory_info.rss / (1024 * 1024):.2f} MB")


                    # Set event if memory usage is above the limit 
                    if memory_info.rss > 2252800000000: #Set your own limit
                        print(f"Memory usage exceeded the limit.")
                        #self.event.set()
                        #self.barrier.wait()
            if buffer:
                self.db_queue.put(buffer)

        except Exception as e:
            print(e)

        gc.collect()
        # Delete the files
        del self.bunch
        # Notify the main thread that this thread has completed its work
        self.callback()


if __name__ == "__main__":
    start_time = time.time()
    processor = DockingProcessor()
    processor.run()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))
###i can be stupid 