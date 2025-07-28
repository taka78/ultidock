#dock_v02.py
import os, stat
import threading
from threading import Semaphore
import queue
import glob
import subprocess
import time
import psutil
import uuid
import gc
import numpy as np
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


class DockingProcessor:

    def __init__(self):
        _ensure_vina_exec()
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        self.barrier = Barrier(6)
        self.event = threading.Event()
        # Initialize database manager and queue for async writes
        cpu_count = os.cpu_count() or 1
        cpus_per_vina = 2
        max_parallel = max(1, cpu_count // cpus_per_vina)
        print(f"Allowing up to {max_parallel} concurrent Vina runs")
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
                if isinstance(item, list):
                    for rec in item:
                        if isinstance(rec, tuple) and len(rec) == 5:
                            buffer.append(rec)
                        else:
                            print(f"Skipping bad sub-record in list: {rec}")
                elif isinstance(item, tuple) and len(item) == 4:
                    buffer.append(item)
                else:
                    print(f"Skipping malformed record: {item}")

                self.db_queue.task_done()  #Very important!

                if len(buffer) >= 500:
                    print(f"Flushing 500 records to DB")
                    self.db_manager.insert_bulk(buffer)
                    buffer.clear()

            except queue.Empty:
                if buffer:
                    print(f"Timeout flush: Writing {len(buffer)} records to DB")
                    self.db_manager.insert_bulk(buffer)
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
                self.vina_sem      # ← new argument
            )
            threads.append(thread)

        threshold = 1024  # available memory threshold in bytes
        while True:
            if psutil.virtual_memory().available > threshold:
                break
            time.sleep(1)

        for thread in threads:
            thread.start()

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


class GPFGenerator:

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

    def calculate_grid_center_and_size(self, ligand_file, padding=10.0):
        x_coords, y_coords, z_coords = [], [], []
        with open(ligand_file, 'r') as f:
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

    def write_gpf_file(self, receptor_file, atom_types, center, size, output_gpf, spacing=0.375):
        try:
            # Ensure grid point counts are ODD (AutoDock likes that)
            npts_x = int(size[0] / spacing)
            npts_y = int(size[1] / spacing)
            npts_z = int(size[2] / spacing)

            npts_x = npts_x + 1 if npts_x % 2 == 0 else npts_x
            npts_y = npts_y + 1 if npts_y % 2 == 0 else npts_y
            npts_z = npts_z + 1 if npts_z % 2 == 0 else npts_z

            base_name = os.path.splitext(os.path.basename(receptor_file))[0]
            fld_filename = f"{base_name}.maps.fld"

            with open(output_gpf, 'w') as f:
                f.write(f"npts {npts_x} {npts_y} {npts_z}\n")  # MUST come before gridcenter
                f.write(f"gridfld {fld_filename}\n")         # MUST come before gridcenter
                f.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
                f.write(f"spacing {spacing}\n")
                f.write(f"receptor {MACRO_MOL_DIR}/{os.path.basename(receptor_file)}\n")
                f.write(f"ligand_types {atom_types}\n")

                for atom in atom_types.split():
                    f.write(f"map {base_name}.{atom}.map\n")

                f.write(f"elecmap {base_name}.e.map\n")
                f.write(f"dsolvmap {base_name}.d.map\n")
                f.write("dielectric -0.1465\n")
                f.write("torsdof 0\n")  # TODO : Update if torsion degrees can be parsed DEFINETLY WILL DO THAT

            print(f"[INFO] GPF file created successfully: {output_gpf}")

        except Exception as e:
            print(f"[ERROR] Failed to write GPF: {e}")
            raise


    def create_gpf(self, receptor_file, ligand_template, output_gpf="grid_params.gpf", spacing=0.375):
        atom_types = self.extract_atom_types_from_pdbqt(ligand_template)
        center, size = self.calculate_grid_center_and_size(ligand_template)
        self.write_gpf_file(receptor_file, atom_types, center, size, output_gpf, spacing)

class ProcessFileThread(threading.Thread):
    
    def __init__(self, bunch, barrier, event, callback, db_queue, vina_sem):
        super().__init__()
        self.bunch = bunch
        self.barrier = barrier
        self.event = event
        self.callback = callback
        self.db_queue = db_queue
        self.vina_sem = vina_sem



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


    def run(self):
        # Calculate grid center and size dynamically for each ligand
        macro_mols = glob.glob(f"{MACRO_MOL_DIR}/*.pdbqt")
        if not macro_mols:
            print("No macro molecule found in", MACRO_MOL_DIR)
            return
        macro_mol = macro_mols[0]
        receptor_name = macro_mol.split("/")[-1]
        fld_files = glob.glob(os.path.join(MACRO_MOL_DIR, "*.maps.fld"))
        fld_file = fld_files[0] if fld_files else None
        try:
            BATCH_SIZE = 500
            grid_center, grid_size = self.calculate_grid_center_and_size(f"{macro_mol}")
            for ligand_file in self.bunch:
                print(ligand_file)
                output_file = f"{DOCKING_DIR}/{ligand_file[-26:]}_{uuid.uuid4()}.pdbqt" # will that last though?
                self.vina_sem.acquire()  # Acquire the semaphore before running Vina
                try:
                    if GPU_TYPE == "NVIDIA":
                        # Run AutoDock-GPU
                        result = subprocess.run([
                            f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_128wi",
                            "--lfile", f"{ligand_file}",
                            "--ffile", f"{fld_file}",
                            "--nrun", "20"
                            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=MACRO_MOL_DIR, timeout=1200)
                        print(f"[AutoDock-GPU stdout]\n{result.stdout}") #debugging time
                        print(f"[AutoDock-GPU stderr]\n{result.stderr}") #debugging time
                    else:
                        result = subprocess.run([
                            f"{VINA_DIR}/bin/vina",
                            "--receptor", f"{macro_mol}",
                            "--ligand", f"{ligand_file}",
                            "--center_x", str(grid_center[0]),
                            "--center_y", str(grid_center[1]),
                            "--center_z", str(grid_center[2]),
                            "--size_x", str(grid_size[0]),
                            "--size_y", str(grid_size[1]),
                            "--size_z", str(grid_size[2]),
                            "--cpu", "2",
                            "--out", output_file
                        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=1200)
                finally:
                    self.vina_sem.release()
                print(f"[{threading.current_thread().name}] Docking {ligand_file}")


                # Check for subprocess failure
                if result.returncode != 0:
                    print(f"Vina failed for: {ligand_file}")
                    print(f"STDOUT:\n{result.stdout.strip()}")
                    print(f"STDERR:\n{result.stderr.strip()}")
                    continue  # skip to the next ligand

                # Wait until the output file really exists
                timeout = 5
                while not os.path.exists(output_file) and timeout > 0:
                    time.sleep(0.1)
                    timeout -= 0.1

                if not os.path.exists(output_file):
                    print(f"Output file missing for: {ligand_file}")
                    continue

                parsed_results = self.parse_vina_output_file(output_file, receptor_name, ligand_file)

                if not parsed_results:
                    print(f"No valid docking data in: {output_file}")
                    continue

                buffer.extend(parsed_results)

                if len(buffer) >= BATCH_SIZE:
                    self.db_queue.put(buffer.copy())
                    buffer.clear()
 
                # Check memory usage
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