#dock_v02.py
import os
import threading
import queue
import glob
import subprocess
import time
import psutil
import uuid
import gc
import numpy as np
from threading import Barrier
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, MACRO_MOL_DIR, DB_PATH, GPU_TYPE, RESULTS_DIR
from db_manager import DockingDatabaseManager



class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        self.barrier = Barrier(6)
        self.event = threading.Event()
        # Initialize database manager and queue for async writes
        self.db_manager = DockingDatabaseManager(DB_PATH)
        self.db_queue = queue.Queue()
        self.db_thread = threading.Thread(target=self._db_worker, daemon=True)
        self.db_thread.start()


    def _db_worker(self):
        print("🧵 DB worker started")
        buffer = []
        seen_none = False

        while True:
            try:
                print("🕒 DB worker waiting for item...")
                item = self.db_queue.get(timeout=10)
                print(f"📥 DB worker received: {item}")

                if item == "INIT":
                    print("🧊 Received warm-up INIT. Touching DB...")
                    self.db_manager.insert_bulk([])  # Safe no-op
                    self.db_queue.task_done()
                    continue

                elif item == "SHUTDOWN":
                    print("🛑 Shutdown signal received. Flushing and exiting DB thread.")
                    if buffer:
                        self.db_manager.insert_bulk(buffer)
                    self.db_queue.task_done()
                    break

                # Item is valid
                if isinstance(item, list):
                    for rec in item:
                        if isinstance(rec, tuple) and len(rec) == 4:
                            buffer.append(rec)
                        else:
                            print(f"❌ Skipping bad sub-record in list: {rec}")
                elif isinstance(item, tuple) and len(item) == 4:
                    buffer.append(item)
                else:
                    print(f"❌ Skipping malformed record: {item}")

                self.db_queue.task_done()  # ✅ Very important!

                if len(buffer) >= 500:
                    print(f"📤 Flushing 500 records to DB")
                    self.db_manager.insert_bulk(buffer)
                    buffer.clear()

            except queue.Empty:
                if buffer:
                    print(f"⏳ Timeout flush: Writing {len(buffer)} records to DB")
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
            thread = ProcessFileThread([ligand_file], self.barrier, self.event, thread_callback, self.db_queue)
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

    def run_autogrid(receptor_file, grid_center, grid_size):
        try:
            subprocess.run([
                "autogrid4",
                "-p", "grid_parameters.gpf",  # .gpf file with grid parameters
                "-l", f"{receptor_file}.glg"  # .glg log file
            ])
            print(f"AutoGrid completed for receptor {receptor_file}")
        except Exception as e:
            print(f"Error running AutoGrid: {e}")

    # Function to generate .gpf file for AutoGrid
    def create_gpf(receptor_file, grid_center, grid_size):
        gpf_content = f"""\nnpts {grid_size[0]} {grid_size[1]} {grid_size[2]}  # Grid size\ngridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]}  # Grid center\nreceptor {receptor_file}  # Receptor PDBQT file\nspacing 0.375  # Grid spacing\n"""
        with open("grid_parameters.gpf", "w") as gpf_file:
            gpf_file.write(gpf_content)

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

class ProcessFileThread(threading.Thread):
    
    def __init__(self, bunch, barrier, event, callback, db_queue):
        super().__init__()
        self.bunch = bunch
        self.barrier = barrier
        self.event = event
        self.callback = callback
        self.db_queue = db_queue



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
    
    def parse_vina_output_file(self, filepath, receptor_name):
    
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
                        results.append((tag, current_affinity, current_rmsd_lb, current_rmsd_ub))

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
        try:
            buffer = []
            BATCH_SIZE = 500
            grid_center, grid_size = self.calculate_grid_center_and_size(f"{macro_mol}")
            for ligand_file in self.bunch:
                print(ligand_file)
                output_file = f"{DOCKING_DIR}/{ligand_file[-26:]}_{uuid.uuid4()}.pdbqt" # will that last though?
                subprocess.run([
                    f"{VINA_DIR}/bin/vina",
                    "--receptor", f"{macro_mol}",
                    "--ligand", f"{ligand_file}",
                    "--center_x", str(grid_center[0]),
                    "--center_y", str(grid_center[1]),
                    "--center_z", str(grid_center[2]),
                    "--size_x", str(grid_size[0]),
                    "--size_y", str(grid_size[1]),
                    "--size_z", str(grid_size[2]),
                    "--cpu", "6", # Number of CPU cores to use- CHANGE IT YOU LITTLE M SERIES MAC USERS
                    "--out", output_file #i said it's not gonna last
                ])

                parsed_results = self.parse_vina_output_file(output_file, receptor_name)

                if not parsed_results:
                    print(f"⚠️ No valid docking data in: {output_file}")
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
                    self.event.set()
                    self.barrier.wait()
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
    processor = DockingProcessor()
    processor.run()
###im stupid
#     