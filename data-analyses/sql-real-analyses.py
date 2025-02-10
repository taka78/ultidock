import sqlite3
import pandas as pd
from pathlib import Path

# Connect to SQLite database
db_path = "docking_results.db"
path = Path("docking_results.db")
folder_path = path.parent.absolute()
conn = sqlite3.connect(db_path)
cur = conn.cursor()

cur.execute("SELECT COUNT(*) FROM docking_results;")
total_records = cur.fetchone()[0]
print(f"Total records in docking_results: {total_records}")



cur.execute("SELECT MIN(affinity), MAX(affinity) FROM docking_results;")
min_affinity, max_affinity = cur.fetchone()
print(f"Affinity values range from {min_affinity} to {max_affinity}")



#choosing the affinity values less than -6 kcal/mol and Dist from RMSD l.b. less than 2Å and Best mode RMSD u.b. less than 3Å and creating a new dataframe with them
query = """
WITH great_ones AS (
    SELECT * 
    FROM docking_results 
    WHERE affinity < -7
      AND rmsd_lb < 2
      AND rmsd_ub < 3
),
fabulous_ones AS (
    SELECT * 
    FROM great_ones 
    WHERE zinc_id IN (
        SELECT zinc_id 
        FROM great_ones 
        GROUP BY zinc_id 
        HAVING COUNT(*) > 1
    )
)
SELECT *
FROM fabulous_ones
WHERE model > 2;
"""

# Execute the query and load the results into a Pandas DataFrame
we_need_this = pd.read_sql_query(query, conn)

# Close the database connection
conn.close()

# Print the results
print(we_need_this)

# Optionally, save the filtered results
we_need_this.to_csv("filtered_results.csv", index=False)
we_need_this.to_json("filtered_results.json", orient="records")
