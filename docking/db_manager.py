# db_manager.py
import os
import sqlite3
import threading

# Pull the results directory from config
try:
    from config import RESULTS_DIR
except ImportError:
    RESULTS_DIR = "results"  # fallback if config is missing (for testability)

class DockingDatabaseManager:
    def __init__(self, db_filename='ultidock_results.db'):
        # Ensure the directory exists
        os.makedirs(RESULTS_DIR, exist_ok=True)

        self.db_path = os.path.join(RESULTS_DIR, db_filename)
        self.lock = threading.Lock()

        # Connect with multi-thread support
        self.connection = sqlite3.connect(self.db_path, check_same_thread=False)
        self.cursor = self.connection.cursor()

        self._initialize_db()

    def _initialize_db(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS docking_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ligand_name TEXT,
                binding_affinity REAL,
                rmsd_lb REAL,
                rmsd_ub REAL,
                docking_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        self.connection.commit()

    def insert_docking_result(self, ligand_name, binding_affinity, rmsd_lb, rmsd_ub):
        with self.lock:
            try:
                self.cursor.execute('''
                    INSERT INTO docking_results (ligand_name, binding_affinity, rmsd_lb, rmsd_ub)
                    VALUES (?, ?, ?, ?)
                ''', (ligand_name, binding_affinity, rmsd_lb, rmsd_ub))
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while inserting docking result: {e}")
                self.connection.rollback()

    def close(self):
        if self.connection:
            self.connection.close()
