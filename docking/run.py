#integration coming
import os
from config import BASE_DIR

# List your scripts in one place
scripts = [
    "setup.py",
    "extract.py",
    "dock_v02.py",
    "analyse_docking_results.py"
]

for script in scripts:
    path = os.path.join(BASE_DIR, script)
    with open(path) as f:
        exec(f.read(), globals())
