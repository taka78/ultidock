import os
import importlib.util

# First run setup.py (which creates config.py), yup im stupid shut up.
with open("setup.py") as f:
    exec(f.read(), globals())

# Now that setup.py has run, we can safely import config
import sys
import importlib
config = importlib.import_module("config")

# Scripts that require config
scripts = [
    "extract.py",
    "dock_v02.py",
    "analyse_docking_results.py"
]

for script in scripts:
    path = os.path.join(config.BASE_DIR, script)
    with open(path) as f:
        exec(f.read(), globals())
