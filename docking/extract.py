import threading
import time
import glob
import shutil
import subprocess
import gzip
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR

def delete_the_archives(self):
    """Processes the file `f` and then deletes it."""


class ProcessFileThread(threading.Thread):
    def __init__(self, f, barrier, event):
        super().__init__()
        self.f = f
        self.barrier = barrier
        self.event = event

    def run(self):
        self.extract_all()
        self.barrier.wait()
        self.event.set()
        self.process_file()
        self.barrier.wait()
        self.event.set()

    def extract_all(self):
        print(self.f)
        try:
            with gzip.open(f'{self.f}',"rb") as f_in:
                with open(self.f) as f_out:
                    f_out.write(f"{self.f}.a")
        except Exception as e:
            print(f"{self.f}", "occured error", e)
            pass

    def process_file(self):
        """Processes the file `f` and then deletes it."""
        print(f"Processing {self.f}")
        subprocess.run([f"{VINA_DIR}/bin/vina_split", "--input", f"{self.f}"])
        
        # Delete the file

    
def main():
    print(LIGANDS_DIR)
    FILES = glob.glob(f"{LIGANDS_DIR}/ligands_raw/*")
    #print(FILES)

    # Create a barrier object
    barrier = threading.Barrier(len(FILES)+1)

    # Create an event object to signal to the main thread that all of the threads have finished processing their files
    event = threading.Event()

    # Create a thread for each file
    threads = []
    for f in FILES:
        t = ProcessFileThread(f, barrier, event)
        threads.append(t)
        t.start()

    # Wait for all of the threads to finish processing their files
    for t in threads:
        t.join()

    # Wait for the event to be set before moving on
    event.wait()

    for f in FILES:
        #print(f"Deleting archive {f}")
        subprocess.run(["sudo", "rm", "-f", f"{f}"])
    

if __name__ == "__main__":
    main()