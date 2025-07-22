#!/usr/bin/env python3
"""
analyze_docking_results.py
========================
Query the SQLite docking_results DB and filter:
  • affinity < threshold
  • rmsd_lb < threshold
  • rmsd_ub < threshold
  • ligand_id appears >1 times
  • model > min_model

Write results by default to a CSV file named "<date>-docking-results.csv" (placed next to the database) and also print all rows to the console.

Usage:
    python analyze_docking_results.py \
        [--db <db_path>] \
        [--affinity <value>] \
        [--rmsd-lb <value>] \
        [--rmsd-ub <value>] \
        [--min-model <value>] \
        [--out <output_file>]

If --out is provided, it overrides the default filename and format is chosen by extension.
"""
import time
import argparse
import sqlite3
import os
from datetime import date
from config import DB_PATH
import pandas as pd

# Determine default output directory based on DB_PATH
_db_dir = os.path.dirname(DB_PATH) or '.'

# Default thresholds
try:
    from config import (
        DEFAULT_AFFINITY,
        DEFAULT_RMSD_LB,
        DEFAULT_RMSD_UB,
        DEFAULT_MIN_MODEL,
    )
except ImportError:
    DEFAULT_AFFINITY = -6.0        # kcal/mol // you should change this according to how much chemically active your macromolecule.
    DEFAULT_RMSD_LB = 3.0         # Å // you should change this according to how big your macromolecule's docking site is.
    DEFAULT_RMSD_UB = 8.0         # Å // you should change this according to how big your macromolecule's docking site is.
    DEFAULT_MIN_MODEL = 2          # integer // you should change this according to how picky you are.

# Default output filename: YYYY-MM-DD-docking-results.csv
DEFAULT_OUT = os.path.join(_db_dir, f"{time.strftime('%Y-%m-%d-%H-%M-%S')}-docking-results.csv")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze docking results stored in SQLite"
    )
    parser.add_argument(
        "--db", default=DB_PATH,
        help="Path to SQLite database file (default from config)"
    )
    parser.add_argument(
        "--affinity", type=float, default=DEFAULT_AFFINITY,
        help=f"Maximum binding affinity (kcal/mol) [default: {DEFAULT_AFFINITY}]"
    )
    parser.add_argument(
        "--rmsd-lb", type=float, default=DEFAULT_RMSD_LB,
        help=f"Maximum RMSD lower bound (Å) [default: {DEFAULT_RMSD_LB}]"
    )
    parser.add_argument(
        "--rmsd-ub", type=float, default=DEFAULT_RMSD_UB,
        help=f"Maximum RMSD upper bound (Å) [default: {DEFAULT_RMSD_UB}]"
    )
    parser.add_argument(
        "--min-model", type=int, default=DEFAULT_MIN_MODEL,
        help=f"Minimum model index (exclusive) [default: {DEFAULT_MIN_MODEL}]"
    )
    parser.add_argument(
        "--out", default=DEFAULT_OUT,
        help=f"Output file path (default: {DEFAULT_OUT})"
    )
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    query = """
    WITH filtered AS (
      SELECT
        ligand_name,
        "binding_affinity (kcal/mol)" AS affinity,
        "rmsd_lb (Å)" AS rmsd_lb,
        "rmsd_ub (Å)" AS rmsd_ub,
        "docking_file",
        CAST(
          substr(
            ligand_name,
            instr(ligand_name, 'Model') + 5,
            instr(substr(ligand_name, instr(ligand_name, 'Model')+5), '-') - 1
          ) AS INTEGER
        ) AS model,
        substr(ligand_name, 1, instr(ligand_name, '-') - 1) AS zinc_id
      FROM docking_results
      WHERE
        "binding_affinity (kcal/mol)" < ?
        AND "rmsd_lb (Å)" < ?
        AND "rmsd_ub (Å)" < ?
    ), dupes AS (
      SELECT zinc_id
      FROM filtered
      GROUP BY zinc_id
      HAVING COUNT(*) > 1
    )
    SELECT *
    FROM filtered
    WHERE zinc_id IN (SELECT zinc_id FROM dupes)
      AND model > ?
    ORDER BY affinity;
    """

    cur.execute(query, (
        args.affinity,
        args.rmsd_lb,
        args.rmsd_ub,
        args.min_model,
    ))
    rows = [dict(r) for r in cur.fetchall()]
    df = pd.DataFrame(rows)

    # Write to output file
    out_path = args.out
    ext = os.path.splitext(out_path)[1].lower()
    if ext in ('.xlsx', '.xls'):
        df.to_excel(out_path, index=False)
        print(f"🔹 Results written to Excel: {out_path}")
    else:
        df.to_csv(out_path, index=False)
        print(f"🔹 Results written to CSV: {out_path}")

    # Always print all results to screen
    if not df.empty:
        print(df.to_string(index=False))
    else:
        print("No rows matched the filters.")

if __name__ == '__main__':
    main()