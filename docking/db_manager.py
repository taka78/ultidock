# db_manager.py
import os
import sqlite3
import threading

class DockingDatabaseManager:
    def __init__(self, db_path='docking_results.db'):
        self.db_path = db_path
        # Create a lock for thread-safe database operations
        self.lock = threading.Lock()
        # Allow the connection to be shared across threads
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
        """Insert a new docking result into the database safely."""
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