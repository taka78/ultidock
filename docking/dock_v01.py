import threading
import glob
import subprocess
import time
import psutil
import resource
import uuid
import gc
import numpy as np
from threading import Barrier
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, MACRO_MOL_DIR



class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        self.barrier = Barrier(6)
        self.event = threading.Event()

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

        for i in range(0, len(self.FILES), 100):
            bunch = self.FILES[i:i+100]
            thread = ProcessFileThread(bunch, self.barrier, self.event, thread_callback)
            threads.append(thread)

        threshold = 1024
        while True:
            available_memory = psutil.virtual_memory().available
            if available_memory > threshold:
                break
            time.sleep(1)

        for thread in threads:
            thread.start()

        while completed_threads < len(threads):
            time.sleep(1)

        gc.collect()

        # After threads have completed, check memory
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
        gpf_content = f"""\
        npts {grid_size[0]} {grid_size[1]} {grid_size[2]}  # Grid size
        gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]}  # Grid center
        receptor {receptor_file}  # Receptor PDBQT file
        spacing 0.375  # Grid spacing
        """
        with open("grid_parameters.gpf", "w") as gpf_file:
            gpf_file.write(gpf_content)


    # Function to parse ligand's .pdbqt file and extract atomic coordinates
    def calculate_grid_center_and_size(ligand_file):
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



    def check_memory(self):
        bunch = []
        for i in range(0, len(self.FILES), 6):
            current_bunch = self.FILES[i:i+6]
            bunch.extend(current_bunch)
            memory_info = psutil.Process().memory_info()
            print(f"Memory usage: {memory_info.rss / (1024 * 1024)} MB")

            if memory_info.rss > 22528000:
                print("Memory usage exceeded the limit.")
                self.event.set()
                self.barrier.wait()
                break

        print(f"Bunch loaded.")

    def run(self):
        self.process_files()

class ProcessFileThread(threading.Thread):
    def __init__(self, bunch, barrier, event, callback):
        super().__init__()
        self.bunch = bunch
        self.barrier = barrier
        self.event = event
        self.callback = callback

    def run(self):
        # Calculate grid center and size dynamically for each ligand
        grid_center, grid_size = calculate_grid_center_and_size(f)
        
        try:
            for f in self.bunch:
                subprocess.run([
                    f"{VINA_DIR}",
                    "--receptor", f"{MACRO_MOL_DIR}",
                    "--ligand", f"{f}",
                    "--center_x", str(grid_center[0]),
                    "--center_y", str(grid_center[1]),
                    "--center_z", str(grid_center[2]),
                    "--size_x", str(grid_size[0]),
                    "--size_y", str(grid_size[1]),
                    "--size_z", str(grid_size[2]),
                    "--cpu", "6",
                    "--out", f"{f}_{uuid.uuid4()}.pdbqt"
                ])

                # Check memory usage
                memory_info = psutil.Process().memory_info()
                print(f"Memory usage: {memory_info.rss / (1024 * 1024)} GB")

                # Set event if memory usage is above the limit
                if memory_info.rss > 2252800000000:
                    print("Memory usage exceeded the limit.")
                    self.event.set()
                    self.barrier.wait()

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
