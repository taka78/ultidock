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



all_files = glob.glob("/mnt/g/ligands-prep/*.pdbqt")

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


def dock(f):
    print(f"Docking {f}")
    print(f"Memory usage: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2} MB")
    try:
        subprocess.run([
            "/mnt/c/Users/taha_/OneDrive/Belgeler/docking/autodock_vina_1_1_2_linux_x86/bin/vina",
            "--receptor", "/mnt/g/4h10_edited-autodock-with-remark.pdbqt",
            "--ligand", f"/mnt/g/ligands-prep/{f}",
            "--center_x", "20",
            "--center_y", "-15",
            "--center_z", "0",
            "--size_x", "30",
            "--size_y", "30",
            "--size_z", "30",
            "--cpu", "6",
            "--out", f"{f}_{uuid.uuid4()}.pdbqt"
        ])
        print(f"Docking {f} completed")
    except Exception as e:
        print(f"Error docking {f}: {e}")



def main():
    chunks = chunk_list(all_files, 1000)
    all_file_bunches = create_dataframe(chunks)
    first_bunch = all_file_bunches.iloc[0]
    pool = multiprocessing.Pool()
    pool.map(dock, first_bunch)

if __name__ == "__main__":
    main()

##there is 944951 ligands in the ligand-prep folder
###select the first row of the dataframe




