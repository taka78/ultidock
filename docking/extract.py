import threading
import time
import glob
import shutil
import subprocess
import os
import gzip
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR

class ProcessFileThread(threading.Thread):
    def __init__(self, f, extraction_barrier):
        super().__init__()
        self.f = f  # File path
        self.extraction_barrier = extraction_barrier

    def run(self):
        # First, extract the archive
        self.extract_all()
        # Wait until all files are extracted
        self.extraction_barrier.wait()  # Sync point after extraction
        # Then, process with vina_split
        self.process_file()

    def extract_all(self):
        # Ensure the file has the .gz extension
        if not self.f.endswith(".gz"):
            print(f"{self.f} is not a .gz file.")
            return

        # Get the path without the .gz extension for the output file
        extracted_file_path = self.f[:-3]

        # Open the .gz file and extract its contents
        try:
            with gzip.open(self.f, 'rb') as gz_file:
                with open(extracted_file_path, 'wb') as extracted_file:
                    shutil.copyfileobj(gz_file, extracted_file)
            print(f"Extracted: {extracted_file_path}")
            
            # After successful extraction, delete the .gz file
            os.remove(self.f)
            print(f"Deleted archive: {self.f}")

        except Exception as e:
            print(f"Failed to extract {self.f}: {e}")

    def process_file(self):
        extracted_file_path = self.f[:-3]  # Path without the .gz extension
        
        if not extracted_file_path.endswith(".pdbqt"):
            print(f"{extracted_file_path} is not a ligand file.")
            return
        try:
            """Processes the file using vina_split."""
            print(f"Processing {extracted_file_path}")
            subprocess.run([f"{VINA_DIR}/bin/vina_split", "--input", f"{extracted_file_path}"])
            os.remove(extracted_file_path)
        except Exception as e:
            print(e)
            pass

def main():
    print(LIGANDS_DIR)

    # Phase 1: Find all .gz files
    FILES = glob.glob(f"{LIGANDS_DIR}/*.gz")
    
    # Create a barrier that will block the process_file phase until all threads have finished extracting
    extraction_barrier = threading.Barrier(len(FILES) + 1)  # +1 includes the main thread

    # Create a thread for each file (both for extraction and split process)
    threads = []
    for f in FILES:
        t = ProcessFileThread(f, extraction_barrier)
        threads.append(t)
        t.start()

    # Main thread waits on the barrier as well
    extraction_barrier.wait()  # This ensures the main thread participates in the barrier

    # Wait for all threads to complete
    for t in threads:
        t.join()

if __name__ == "__main__":
    main()