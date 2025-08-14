#dock_v02.py
import os, stat
import threading
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
import numpy as np
import shutil  
from threading import Barrier
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, AUTODOCK_GPU_DIR, MACRO_MOL_DIR, DB_PATH, GPU_TYPE, RESULTS_DIR
from db_manager import DockingDatabaseManager

def _ensure_vina_exec():
    for exe in ("vina", "vina_split"):
        path = os.path.join(VINA_DIR, "bin", exe)
        try:
            st = os.stat(path)
            # if owner-x isn’t already set, add it (resulting in at least 0o755)
            if not (st.st_mode & stat.S_IXUSR):
                os.chmod(path, st.st_mode | stat.S_IXUSR)
        except FileNotFoundError:
            print(f"Warning: {exe} not found at {path}")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
BINARY_PATH = os.path.join(AUTODOCK_GPU_DIR, "bin", "autodock_gpu_128wi")
COMPILER_SCRIPT = os.path.join(os.path.dirname(__file__), "autodock-gpu-compiler.sh")

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

class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        print(f"[INIT] ligands discovered: {len(self.FILES)} in {LIGANDS_DIR}")

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

class GPFGenerator:
    def calculate_grid_center_and_size(self, receptor_file, padding=10.0):
        ##in the past i was mixing macro molecules and ligands, which is a bad idea
        # This function calculates the geometric center and size of the grid based on the receptor file.
        x_coords, y_coords, z_coords = [], [], []
        with open(receptor_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        x_coords.append(x)
                        y_coords.append(y)
                        z_coords.append(z)
                    except ValueError:
                        continue

        center_x = np.mean(x_coords)
        center_y = np.mean(y_coords)
        center_z = np.mean(z_coords)

        size_x = np.max(x_coords) - np.min(x_coords) + padding
        size_y = np.max(y_coords) - np.min(y_coords) + padding
        size_z = np.max(z_coords) - np.min(z_coords) + padding

        return (center_x, center_y, center_z), (size_x, size_y, size_z)

    def write_gpf_file(self, receptor_file, center, size, output_gpf, spacing=0.375):
        ###yup, i realised that generating all maps file in one go is more computationally efficient
        # than parsing all ligands for needed atom types
        # and then generating maps for each ligand separately.
        # This is because AutoDock can generate all maps in one go, and then use them
        # for all ligands separately, as needed. which is much faster than generating maps for each ligand separately.
        DEFAULT_AUTODOCK_ATOM_TYPES = [
            "A",  # aliphatic carbon
            "C",  # aromatic carbon
            "HD", # hydrogen donor
            "N",  # nitrogen
            "NA", # nonpolar nitrogen
            "OA", # oxygen acceptor
            "S",  # sulfur
            "SA", # sulfur acceptor (optional, rarely used)
            "Cl", "Br", "F", "I", # halogens
            "Zn", "Mg", "Ca", "Fe", "Mn"  # metals (optional based on use-case)
        ]
        try:
            # Ensure grid point counts are odd
            npts = [int(dim / spacing) | 1 for dim in size]  # force odd with bitwise OR

            base_name = os.path.splitext(os.path.basename(receptor_file))[0]
            fld_filename = f"{base_name}.maps.fld"
            output_gpf_path = os.path.join(MACRO_MOL_DIR, os.path.basename(output_gpf))

            with open(output_gpf_path, 'w') as f:
                f.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
                f.write(f"gridfld {fld_filename}\n")
                f.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
                f.write(f"spacing {spacing}\n")
                f.write(f"receptor {os.path.basename(receptor_file)}\n")
                f.write(f"ligand_types {' '.join(DEFAULT_AUTODOCK_ATOM_TYPES)}\n")

                for atom in DEFAULT_AUTODOCK_ATOM_TYPES:
                    f.write(f"map {base_name}.{atom}.map\n")

                f.write(f"elecmap {base_name}.e.map\n")
                f.write(f"dsolvmap {base_name}.d.map\n")
                f.write("dielectric -0.1465\n")


            print(f"[INFO] GPF file written: {output_gpf_path}")
        except Exception as e:
            print(f"[ERROR] Failed to write GPF file: {e}")
            raise

    def run_autogrid(self, gpf_path):
        log_path = gpf_path.replace(".gpf", ".glg")
        try:
            subprocess.run(
                [f"{AUTODOCK_GPU_DIR}/autogrid/autogrid4", "-p", os.path.basename(gpf_path), "-l", os.path.basename(log_path)],
                check=True,
                stderr=subprocess.PIPE,
                text=True,
                cwd=MACRO_MOL_DIR,
                timeout=1200
            )
            print(f"[INFO] AutoGrid finished. Log written to: {log_path}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] AutoGrid execution failed:\n{e.stderr}")
            raise

    def create_gpf(self, receptor_file, output_gpf="grid_params.gpf", spacing=0.375):
        center, size = self.calculate_grid_center_and_size(receptor_file)
        self.write_gpf_file(receptor_file, center, size, output_gpf, spacing)
        self.run_autogrid(os.path.join(MACRO_MOL_DIR, os.path.basename(output_gpf)))


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
    
    def parse_vina_output_file(self, filepath, receptor_name,ligand_file):
    
        #  Retry a few times to let file system catch up
        for _ in range(3):
            if os.path.exists(filepath):
                break
            time.sleep(0.05)  # Small delay (50ms) yay
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

        with open(filepath, 'r') as fp:
            for line in fp:
                line = line.strip()

                if line.startswith("MODEL"):
                    # Start new block
                    if current_model is not None and current_ligand_id and current_affinity is not None:
                        tag = f"{current_ligand_id}-{receptor_name}-M{current_model}"
                        results.append((tag, current_affinity, current_rmsd_lb, current_rmsd_ub))

                    # Reset state
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
                            continue

                elif line.startswith("REMARK  Name") and inside_model_block:
                    parts = line.split()
                    if len(parts) >= 4:
                        current_ligand_id = parts[3]

                elif line.startswith("ENDMDL") and inside_model_block:
                    # Finalize this model block
                    if current_model is not None and current_ligand_id and current_affinity is not None:
                        tag = f"{current_ligand_id}-Model{current_model}-{receptor_name}"
                        results.append((tag, current_affinity, current_rmsd_lb, current_rmsd_ub, ligand_file))

                    # Reset all block data
                    current_model = None
                    current_affinity = None
                    current_rmsd_lb = None
                    current_rmsd_ub = None
                    current_ligand_id = None
                    inside_model_block = False

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
        macro_mols = glob.glob(f"{MACRO_MOL_DIR}/*.pdbqt")
        if not macro_mols:
            print("No macro molecule found in", MACRO_MOL_DIR)
            return
        macro_mol = macro_mols[0]
        receptor_name = macro_mol.split("/")[-1]
        fld_files = glob.glob(os.path.join(MACRO_MOL_DIR, "*.maps.fld"))
        if not fld_files:
            print("No .maps.fld found in", MACRO_MOL_DIR)
            return
        fld_file = fld_files[0]
        try:
            buffer = []
            BATCH_SIZE = 500
            grid_center, grid_size = self.calculate_grid_center_and_size(f"{macro_mol}")
            for ligand_file in self.bunch:
                # ensure output dir exists
                Path(DOCKING_DIR).mkdir(parents=True, exist_ok=True)

                lig_stem  = Path(ligand_file).stem
                out_stem  = Path(DOCKING_DIR) / f"{lig_stem}_{uuid.uuid4().hex[:8]}"
                gpu_pdbqt = f"{out_stem}_out.pdbqt"   # AD-GPU best pose (if --gbest 1)
                xml_out   = f"{out_stem}.xml"     # AD-GPU XML
                vina_out  = f"{out_stem}.pdbqt"       # Vina CPU output target

                # Per-GPU semaphore (one job per GPU at a time)
                sem = GPU_SEMAPHORES.get(self.gpu_id, next(iter(GPU_SEMAPHORES.values())))

                env = os.environ.copy()
                env["OMP_NUM_THREADS"] = "1"
                env["MKL_NUM_THREADS"] = "1"
                env["OPENBLAS_NUM_THREADS"] = "1"
                # If you need hard isolation per thread, uncomment:
                # env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)
                acquired_gpu = False
                acquired_cpu = False
                print(ligand_file)

                try:
                    if GPU_TYPE == "NVIDIA":
                        sem.acquire()
                        acquired_gpu = True
                        # --- AutoDock-GPU (XML-first) ---
                        result = subprocess.run(
                            [
                                f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_128wi",
                                "--lfile",   str(Path(ligand_file).resolve()),
                                "--ffile",   str(Path(fld_file).resolve()),
                                "--nrun",    "50",
                                "--gbest",   "5",            # write <resnam>_out.pdbqt (optional but nice)
                                "--xmloutput","1",           # ensure XML is produced
                                "--resnam",  str(out_stem),  # basename; outputs land in DOCKING_DIR
                                "--devnum",  str(int(self.gpu_id) + 1)  # AD-GPU is 1-indexed
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            cwd=MACRO_MOL_DIR,
                            env=env,
                            timeout=100
                        )
                        print(f"[AutoDock-GPU stdout]\n{result.stdout}")
                        print(f"[AutoDock-GPU stderr]\n{result.stderr}")

                        if result.returncode != 0:
                            print(f"Docking failed for: {ligand_file}")
                            print(f"STDOUT:\n{result.stdout.strip()}")
                            print(f"STDERR:\n{result.stderr.strip()}")
                            continue
                        # Ensure output files are created

                    else:
                        # --- Vina CPU fallback (PDBQT parsing) ---
                        # CPU/Vina path uses ONLY the CPU limiter
                        self.vina_sem.acquire()
                        acquired_cpu = True
                        result = subprocess.run(
                            [
                                f"{VINA_DIR}/bin/vina",
                                "--receptor", str(macro_mol),
                                "--ligand",   str(ligand_file),
                                "--center_x", str(grid_center[0]),
                                "--center_y", str(grid_center[1]),
                                "--center_z", str(grid_center[2]),
                                "--size_x",   str(grid_size[0]),
                                "--size_y",   str(grid_size[1]),
                                "--size_z",   str(grid_size[2]),
                                "--cpu",      "2",
                                "--out",      str(vina_out)
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            timeout=1200
                        )
                        if result.returncode != 0:
                            print(f"Vina failed for: {ligand_file}")
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
                if GPU_TYPE == "NVIDIA":
                    parsed_results = self.parse_adgpu_xml(xml_out, receptor_name, ligand_file)
                elif GPU_TYPE == "AMD":
                    parsed_results = self.parse_adgpu_xml(gpu_pdbqt, receptor_name, ligand_file) #definetly gonna change this
                elif GPU_TYPE == "CPU":
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
    gpf_gen = GPFGenerator()
    fld_path = gpf_gen.create_gpf(f"{MACRO_MOL_DIR}/4h10_edited-autodock-with-remark.pdbqt")  # returns .fld path
    processor = DockingProcessor()
    processor.run()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))
###i can be stupid 