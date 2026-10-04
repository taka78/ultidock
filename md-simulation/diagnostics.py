"""MD dependency checks using only the Python standard library."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
from pathlib import Path
import shutil


def check_dependencies(gmx="gmx", acpype="acpype", obabel="obabel", *, workspace=None) -> bool:
    ready = True
    for name, executable in (("GROMACS", gmx), ("ACPYPE", acpype), ("Open Babel", obabel),
                             ("AmberTools antechamber", "antechamber"),
                             ("AmberTools parmchk2", "parmchk2"), ("AmberTools tleap", "tleap"),
                             ("AmberTools sqm", "sqm")):
        location = shutil.which(executable)
        status = "OK" if location else "WARN"
        detail = location or f"{executable} not found; required for MD"
        print(f"  [{status}] {name:30s} {detail}")
        ready &= location is not None
    for name in ("rdkit", "numpy", "scipy"):
        spec = importlib.util.find_spec(name)
        if spec is None:
            print(f"  [WARN] Python {name:23s} not installed; required for MD")
            ready = False
            continue
        try:
            version = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            version = "version unknown"
        print(f"  [OK] Python {name:23s} {version} ({spec.origin})")
    workspace = workspace or Path(__file__).resolve().parent / "workspace"
    print(f"  MD workspace: {workspace}")
    print(f"  MD dependencies: {'ready' if ready else 'incomplete; see md-simulation/README.md'}")
    return ready


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--acpype", default="acpype")
    parser.add_argument("--obabel", default="obabel")
    args = parser.parse_args(argv)
    return 0 if check_dependencies(args.gmx, args.acpype, args.obabel) else 1


if __name__ == "__main__":
    raise SystemExit(main())
