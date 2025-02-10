import sqlite3
import mmap
import glob
import re
from pathlib import Path




conn = sqlite3.connect("docking_results.db")

db_path = "docking_results.db"
path = Path("docking_results.db")
folder_path = path.parent.absolute()

cur = conn.cursor()

cur.execute('''
CREATE TABLE IF NOT EXISTS docking_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    zinc_id TEXT,
    model INTEGER,
    affinity REAL,
    rmsd_lb REAL,
    rmsd_ub REAL
);
''')
conn.commit()

files = glob.glob(f"{folder_path}/docking-outputs/*.pdbqt")

def process_pdbqt_file(file_path):
    with open(file_path, "r") as f:
        mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        content = mmapped_file.read().decode()

    # Extract relevant information using regex (avoids line-by-line reads)
    models = re.findall(r"MODEL (\d+)", content)
    affinities = re.findall(r"REMARK VINA RESULT:\s+([-.\d]+)", content)
    rmsd_lb = re.findall(r"REMARK VINA RESULT:\s+[-.\d]+\s+([-.\d]+)", content)
    rmsd_ub = re.findall(r"REMARK VINA RESULT:\s+[-.\d]+\s+[-.\d]+\s+([-.\d]+)", content)
    zinc_id_match = re.search(r"REMARK\s+Name\s+=\s+(\S+)", content)

    if not zinc_id_match:
        return []  # Skip if no ZINC ID found

    zinc_id = zinc_id_match.group(1)

    # Convert parsed data to a structured list of tuples for SQLite insertion
    data = [
        (zinc_id, int(models[i]), float(affinities[i]), float(rmsd_lb[i]), float(rmsd_ub[i]))
        for i in range(len(models))
    ]

    return data

# Process all files and bulk insert into SQLite
all_data = []
for file in files:
    print(file)
    all_data.extend(process_pdbqt_file(file))

# Bulk insert into SQLite
cur.executemany(
    "INSERT INTO docking_results (zinc_id, model, affinity, rmsd_lb, rmsd_ub) VALUES (?, ?, ?, ?, ?)",
    all_data
)
conn.commit()
conn.close()
