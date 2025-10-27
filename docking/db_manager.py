# db_manager.py
import os
import sqlite3
import threading

# Pull the results directory from config
try:
    from config import RESULTS_DIR
except ImportError:
    RESULTS_DIR = "RESULTS_DIR"  # fallback if config is missing (for testability)

class DockingDatabaseManager:
    def __init__(self, db_filename='ultidock_results.db'):
        # Ensure the directory exists
        os.makedirs(RESULTS_DIR, exist_ok=True)

        self.db_path = os.path.join(RESULTS_DIR, db_filename)
        self.lock = threading.Lock()

        # Connect with multi-thread support
        self.connection = sqlite3.connect(self.db_path, check_same_thread=False)
        self.connection.row_factory = sqlite3.Row
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
                "rmsd_lb (Å)" REAL,
                "rmsd_ub (Å)" REAL,
                "docking_file" TEXT,
                binding_site TEXT,
                docking_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        self.connection.commit()
        self._ensure_binding_site_column()
        self._ensure_metrics_columns()

        self._metric_columns = [
            "cuda_setup_s",
            "rest_setup_s",
            "total_setup_s",
            "job_wait_s",
            "job_duration_s",
            "run_time_s",
            "processing_time_s",
            "run_plus_processing_s",
            "shutdown_s",
            "gpu_memory_mb",
        ]

    def _ensure_binding_site_column(self):
        """Add the binding_site column if an older database is missing it."""
        self.cursor.execute("PRAGMA table_info(docking_results)")
        columns = {row[1] for row in self.cursor.fetchall()}
        if "binding_site" not in columns:
            try:
                self.cursor.execute('ALTER TABLE docking_results ADD COLUMN binding_site TEXT')
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while altering docking_results schema: {e}")
                self.connection.rollback()

    def _ensure_metrics_columns(self):
        """Add timing/memory metric columns to docking_results when missing."""
        metric_columns = {
            "cuda_setup_s": "REAL",
            "rest_setup_s": "REAL",
            "total_setup_s": "REAL",
            "job_wait_s": "REAL",
            "job_duration_s": "REAL",
            "run_time_s": "REAL",
            "processing_time_s": "REAL",
            "run_plus_processing_s": "REAL",
            "shutdown_s": "REAL",
            "gpu_memory_mb": "REAL",
        }

        try:
            self.cursor.execute("PRAGMA table_info(docking_results)")
            existing = {row[1] for row in self.cursor.fetchall()}
        except sqlite3.Error as e:
            print(f"An error occurred while introspecting docking_results: {e}")
            return

        for column, column_type in metric_columns.items():
            if column in existing:
                continue
            try:
                self.cursor.execute(
                    f"ALTER TABLE docking_results ADD COLUMN {column} {column_type}"
                )
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while adding column {column}: {e}")
                self.connection.rollback()

    def insert_docking_result(self, ligand_name, binding_affinity, rmsd_lb, rmsd_ub, docking_file, binding_site):
        with self.lock:
            try:
                self.cursor.execute('''
                    INSERT INTO docking_results (
                        ligand_name,
                        "binding_affinity (kcal/mol)",
                        "rmsd_lb (Å)",
                        "rmsd_ub (Å)",
                        "docking_file",
                        binding_site
                    )
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (ligand_name, binding_affinity, rmsd_lb, rmsd_ub, docking_file, binding_site))
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred while inserting docking result: {e}")
                self.connection.rollback()

    def insert_bulk(self, records):
        """Batch-insert a list of (ligand_name, affinity, rmsd_lb, rmsd_ub, docking_file, binding_site)."""
        if not records:
            print("ℹinsert_bulk: received empty list, skipping.")
            return

        with self.lock:
            try:
                self.cursor.executemany(
                    '''
                    INSERT INTO docking_results (
                        ligand_name,
                        "binding_affinity (kcal/mol)",
                        "rmsd_lb (Å)",
                        "rmsd_ub (Å)",
                        "docking_file",
                        binding_site
                    )
                    VALUES (?, ?, ?, ?, ?, ?)
                    ''',
                    records
                )
                self.connection.commit()
            except sqlite3.Error as e:
                print(f"An error occurred in bulk insert: {e}")
                self.connection.rollback()

    def update_metrics_by_ligand_prefix(self, ligand_prefix, binding_site, metrics):
        """Update metrics for all rows whose ligand_name starts with the prefix."""
        if not metrics or not ligand_prefix:
            return 0

        assignments = ", ".join(f"{key} = ?" for key in metrics)
        values = list(metrics.values())

        pattern = f"{ligand_prefix}%"
        query = f"UPDATE docking_results SET {assignments} WHERE ligand_name LIKE ?"
        values.append(pattern)

        if binding_site:
            query += " AND binding_site = ?"
            values.append(binding_site)

        with self.lock:
            try:
                self.cursor.execute(query, values)
                self.connection.commit()
                return self.cursor.rowcount
            except sqlite3.Error as e:
                print(
                    f"An error occurred while updating metrics for prefix {ligand_prefix}: {e}"
                )
                self.connection.rollback()
                return 0

    def count_matches_by_prefix(self, ligand_prefix, binding_site=None):
        """Return how many rows share the ligand_name prefix (and optional site)."""
        if not ligand_prefix:
            return 0

        pattern = f"{ligand_prefix}%"

        with self.lock:
            try:
                if binding_site:
                    self.cursor.execute(
                        "SELECT COUNT(*) FROM docking_results "
                        "WHERE ligand_name LIKE ? AND binding_site = ?",
                        (pattern, binding_site),
                    )
                else:
                    self.cursor.execute(
                        "SELECT COUNT(*) FROM docking_results "
                        "WHERE ligand_name LIKE ?",
                        (pattern,),
                    )
                row = self.cursor.fetchone()
                return row[0] if row else 0
            except sqlite3.Error as e:
                print(
                    f"An error occurred while counting matches for prefix {ligand_prefix}: {e}"
                )
                return 0

    def fetch_rows_by_ligand_prefix(self, ligand_prefix, binding_site=None):
        """Return metric-bearing rows for ligand names starting with the prefix."""
        if not ligand_prefix:
            return []

        pattern = f"{ligand_prefix}%"
        column_sql = ", ".join(self._metric_columns)
        query = (
            f"SELECT id, ligand_name, binding_site, {column_sql} "
            "FROM docking_results WHERE ligand_name LIKE ?"
        )
        params = [pattern]

        if binding_site:
            query += " AND binding_site = ?"
            params.append(binding_site)

        with self.lock:
            try:
                self.cursor.execute(query, params)
                rows = self.cursor.fetchall() or []
            except sqlite3.Error as e:
                print(
                    f"An error occurred while fetching rows for prefix {ligand_prefix}: {e}"
                )
                return []

        return [dict(row) for row in rows]

    # Backwards-compatible wrappers -------------------------------------------------

    def update_metrics(self, docking_file, binding_site, metrics):
        """Deprecated shim that treats docking_file as a ligand prefix."""
        return self.update_metrics_by_ligand_prefix(docking_file, binding_site, metrics)

    def count_matches(self, docking_file, binding_site=None):
        """Deprecated shim that treats docking_file as a ligand prefix."""
        return self.count_matches_by_prefix(docking_file, binding_site)

    def close(self):
        if self.connection:
            self.connection.close()