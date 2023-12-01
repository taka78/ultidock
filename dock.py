import threading
import glob
import subprocess
import time
import psutil
import resource
import uuid
import gc
from threading import Barrier


class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob("*.pdbqt")
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
        try:
            for f in self.bunch:
                subprocess.run([
                    "/mnt/c/Users/taha_/OneDrive/Belgeler/docking/autodock_vina_1_1_2_linux_x86/bin/vina",
                    "--receptor", "/mnt/g/4h10_edited-autodock-with-remark.pdbqt",
                    "--ligand", f"{f}",
                    "--center_x", "20",
                    "--center_y", "-15",
                    "--center_z", "0",
                    "--size_x", "30",
                    "--size_y", "30",
                    "--size_z", "30",
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
