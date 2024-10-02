import threading
import glob
import subprocess
import time
import psutil
import resource
import uuid
import os
import pandas as pd
import gc
from threading import Barrier
import multiprocessing
import numpy as np
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, MACRO_MOL_DIR



all_files = glob.glob("f{LIGANDS_DIR}/*.pdbqt")

def chunk_list(data, chunk_size):
    new_data = []
    for filename in all_files:
        basename = os.path.basename(filename)
        new_data.append(basename)
    chunks = [new_data[i:i + chunk_size] for i in range(0, len(new_data), chunk_size)]
    return chunks

#create dataframe
def create_dataframe(chunks):
    dataframes = []
    for chunk in chunks:
        dataframe = pd.DataFrame(chunk).T
        dataframes.append(dataframe)
    final_dataframe = pd.concat(dataframes, ignore_index=True)
    return final_dataframe

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

# Modify the dock function to use dynamic grid size and center
def dock(f):
    print(f"Docking {f}")
    print(f"Memory usage: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2} MB")
    
    ligand_file = f"{LIGANDS_DIR}/{f}"
    
    # Calculate grid center and size dynamically for each ligand
    grid_center, grid_size = calculate_grid_center_and_size(ligand_file)
    
    try:
        subprocess.run([
            f"{VINA_DIR}",
            "--receptor", f"{MACRO_MOL_DIR}/*.pdbqt",
            "--ligand", ligand_file,
            "--center_x", str(grid_center[0]),
            "--center_y", str(grid_center[1]),
            "--center_z", str(grid_center[2]),
            "--size_x", str(grid_size[0]),
            "--size_y", str(grid_size[1]),
            "--size_z", str(grid_size[2]),
            "--cpu", "6",
            "--out", f"{f}_{uuid.uuid4()}.pdbqt"
        ])
        print(f"Docking {f} completed")
    except Exception as e:
        print(f"Error docking {f}: {e}")

def main():
    chunks = chunk_list(all_files, 1000)
    print(chunks)
    all_file_bunches = create_dataframe(chunks)
    first_bunch = all_file_bunches.iloc[0]
    
    # Create a pool of workers to process ligands in parallel
    pool = multiprocessing.Pool()
    pool.map(dock, first_bunch)

if __name__ == "__main__":
    main()
##there is 944951 ligands in the ligand-prep folder
###select the first row of the dataframe




