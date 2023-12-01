import threading
import time
import glob
import subprocess

class ProcessFileThread(threading.Thread):
    def __init__(self, f, barrier, event):
        super().__init__()
        self.f = f
        self.barrier = barrier
        self.event = event

    def run(self):
        self.process_file()
        self.barrier.wait()
        self.event.set()

    def process_file(self):
        """Processes the file `f` and then deletes it."""
        print(f"Processing {self.f}")
        subprocess.run(["/mnt/c/Users/taha_/OneDrive/Belgeler/docking/autodock_vina_1_1_2_linux_x86/bin/vina_split", "--input", f"{self.f}"])
        # Delete the file

def main():
    FILES = glob.glob("/mnt/g/ligands-prep/*")

    # Create a barrier object
    barrier = threading.Barrier(len(FILES))

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

    # Remove all processed files
    for f in FILES:
        subprocess.run(["sudo rm", f"{f}"])

if __name__ == "__main__":
    main()