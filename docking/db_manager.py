# db_manager.py
import os
import sqlite3
import threading

# Pull the results directory from config
try:
    from config import RESULTS_DIR
except ImportError:
    RESULTS_DIR = "RESULTS"  # fallback if config is missing (for testability)

class DockingDatabaseManager:
    def __init__(self, db_filename='ultidock_results.db'):
        # Ensure the directory exists
        os.makedirs(RESULTS_DIR, exist_ok=True)

        self.db_path = os.path.join(RESULTS_DIR, db_filename)
        self.lock = threading.Lock()

        # Connect with multi-thread support
        self.connection = sqlite3.connect(self.db_path, check_same_thread=False)
        self.cursor = self.connection.cursor()

        # Performance tuning
        self.connection.execute("PRAGMA journal_mode=WAL;")
        self.connection.execute("PRAGMA synchronous = NORMAL;")

        self._initialize_db()

    def _initialize_db(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS docking_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ligand_name TEXT,
                "binding_affinity (kcal/mol)" REAL,
                "rmsd_lb (\u00c5)" REAL,
                "rmsd_ub (\u00c5)" REAL,
                "docking_file" TEXT,
                docking_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')

        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS simulation_metadata (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ligand_name TEXT,
                ligand_file TEXT,
                complex_gro TEXT,
                topol_top TEXT,
                index_file TEXT,
                em_log TEXT,
                nvt_xtc TEXT,
                rmsd_xvg TEXT,
                com_xvg TEXT,
                hbonds_xvg TEXT,
                rmsd_equilibrated REAL,
                com_initial REAL,
                com_equilibrated REAL,
                com_shift REAL,
                hbonds INTEGER,
                em_success BOOLEAN,
                topology_pass BOOLEAN,
                umbrella_ready BOOLEAN DEFAULT 0
            )
        ''')##yup im cooking eheheh

        self.connection.commit()

    def insert_docking_result(self, ligand_name, binding_affinity, rmsd_lb, rmsd_ub, docking_file):
        with self.lock:
            try:
                self.cursor.execute('''
                    INSERT INTO docking_results (
                        ligand_name,
                        "binding_affinity (kcal/mol)",
                        "rmsd_lb (\u00c5)",
                        "rmsd_ub (\u00c5)",
                        "docking_file"
                    )
                    VALUES (?, ?, ?, ?, ?)
                ''', (ligand_name, binding_affinity, rmsd_lb, rmsd_ub, docking_file))
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while inserting docking result: {e}")
                self.connection.rollback()

    def insert_bulk(self, records):
        """Batch-insert a list of (ligand_name, affinity, rmsd_lb, rmsd_ub, docking_file)."""
        if not records:
            print("\u2139\ufe0f insert_bulk: received empty list, skipping.")
            return

        with self.lock:
            try:
                self.cursor.executemany(
                    '''
                    INSERT INTO docking_results (
                        ligand_name,
                        "binding_affinity (kcal/mol)",
                        "rmsd_lb (\u00c5)",
                        "rmsd_ub (\u00c5)",
                        "docking_file"
                    )
                    VALUES (?, ?, ?, ?, ?)
                    ''',
                    records
                )
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred in bulk insert: {e}")
                self.connection.rollback()

    def insert_simulation_metadata(self, metadata):
        """Insert or update simulation metadata by ligand_name"""
        with self.lock:
            try:
                self.cursor.execute('''
                    INSERT INTO simulation_metadata (
                        ligand_name, ligand_file, complex_gro, topol_top, index_file,
                        em_log, nvt_xtc, rmsd_xvg, com_xvg, hbonds_xvg,
                        rmsd_equilibrated, com_initial, com_equilibrated, com_shift,
                        hbonds, em_success, topology_pass, umbrella_ready
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', metadata)
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while inserting simulation metadata: {e}")
                self.connection.rollback()

    def close(self):
        if self.connection:
            self.connection.close()
