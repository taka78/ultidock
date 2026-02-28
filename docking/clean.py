#!/usr/bin/env python3
"""
Ultidock cleaner: nuke build artifacts and run outputs in one go.

I'm too lazy to write this to relative paths, so it assumes you run it from the docking/ directory. Maybe later I will fix this.
Usage: python3 clean.py [options]

Default behavior = DRY RUN (prints what it would delete).
Use -y / --yes to actually delete.

Categories:
  build    -> autodock-gpu builds (bin/, *.o/*.a etc) + autogrid binary
  results  -> docking outputs (DOCKING_DIR, ANALYSIS_DIR)
  maps     -> generated maps & grid files in MACRO_MOL_DIR (*.map, *.fld, *.gpf, *.glg)
  pycache  -> __pycache__ and *.pyc everywhere
  db       -> results DB file (off by default)
  ligands  -> LIGANDS_DIR contents (off by default)

Shortcuts:
  --all    -> build + results + maps + pycache + db + ligands
"""

from __future__ import annotations
import argparse
import os
import sys
import shutil
from pathlib import Path
import fnmatch

SCRIPT_DIR = Path(__file__).resolve().parent
DOCKING_DIR_DEFAULT = SCRIPT_DIR  
CONFIG_FILE = DOCKING_DIR_DEFAULT / "config.py"
SENTINEL = ".ultidock_sentinel"
PROTECT_NAMES = {
    "Makefile", "MAKEFILE",          # dispatcher or root Makefile
    "Makefile.Cuda", "Makefile.OpenCL",
    "CMakeLists.txt",                # also generally important
}

def _has_sentinel(d: Path) -> bool:
    try:
        return (d / SENTINEL).is_file()
    except Exception:
        return False


def _load_config():
    cfg = {}
    # try to import docking/config.py
    sys.path.insert(0, str(DOCKING_DIR_DEFAULT))
    try:
        import config as C  # type: ignore
        cfg.update({
            "LIGANDS_DIR": Path(getattr(C, "LIGANDS_DIR")),
            "DOCKING_DIR": Path(getattr(C, "DOCKING_DIR")),
            "ANALYSIS_DIR": Path(getattr(C, "ANALYSIS_DIR")),
            "VINA_DIR": Path(getattr(C, "VINA_DIR")),
            "AUTODOCK_GPU_DIR": Path(getattr(C, "AUTODOCK_GPU_DIR")),
            "MACRO_MOL_DIR": Path(getattr(C, "MACRO_MOL_DIR")),
            "RESULTS_DIR": Path(getattr(C, "RESULTS_DIR")),
            "DB_PATH": Path(getattr(C, "DB_PATH")),
            "NUMWI": int(getattr(C, "NUMWI", 128)),
        })
    except Exception:
        # fallbacks relative to docking/
        b = DOCKING_DIR_DEFAULT
        cfg.update({
            "LIGANDS_DIR": b / "LIGANDS_DIR",
            "DOCKING_DIR": b / "DOCKING_DIR",
            "ANALYSIS_DIR": b / "ANALYSIS_DIR",
            "VINA_DIR": b / "VINA_DIR",
            "AUTODOCK_GPU_DIR": b / "AUTODOCK_GPU_DIR",
            "MACRO_MOL_DIR": b / "MACRO_MOL_DIR",
            "RESULTS_DIR": b / "RESULTS_DIR",
            "DB_PATH": (b / "RESULTS_DIR" / "ultidock_results.db"),
            "NUMWI": 128,
        })
    return cfg

CFG = _load_config()


def clean_config(*, dry: bool):
    # optional: require sentinel in docking/ root for extra safety
    if not _has_sentinel(DOCKING_DIR_DEFAULT):
        print("[SKIP] config: missing sentinel in docking/ — nothing deleted.")
        return

    if not _ensure_inside_repo(CONFIG_FILE):
        print("[SKIP] config: path outside repo guard.")
        return

    if CONFIG_FILE.exists() or CONFIG_FILE.is_symlink():
        _safe_delete(CONFIG_FILE, dry=dry)
    else:
        print("[INFO] config.py not found; nothing to do.")

    # also remove compiled bytecode to avoid stale imports
    _glob_delete(DOCKING_DIR_DEFAULT / "__pycache__", ["config.*.pyc"], dry=dry, only_files=True)


def _safe_delete(path: Path, *, dry: bool):
    # Absolute guard: never remove protected names
    if path.name in PROTECT_NAMES:
        print(f"[SKIP] protected: {path}")
        return
    try:
        if path.is_symlink() or path.is_file():
            print(f"DEL file: {path}")
            if not dry:
                path.unlink(missing_ok=True)
        elif path.is_dir():
            print(f"DEL dir : {path}")
            if not dry:
                shutil.rmtree(path, ignore_errors=True)
    except Exception as e:
        print(f"[WARN] failed to delete {path}: {e}")


def _glob_delete(root: Path, patterns: list[str], *, dry: bool, only_files=False, exclude: set[str] | None=None):
    if not root.exists():
        return
    exclude = exclude or set()
    for pat in patterns:
        for p in root.rglob(pat):
            if p.name in PROTECT_NAMES or p.name in exclude:
                print(f"[SKIP] protected/excluded: {p}")
                continue
            if only_files and p.is_dir():
                continue
            _safe_delete(p, dry=dry)


def _ensure_inside_repo(p: Path) -> bool:
    # basic guardrail: don’t allow deleting outside the repo tree
    try:
        return SCRIPT_DIR in p.resolve().parents or p.resolve() == SCRIPT_DIR
    except Exception:
        return False

def clean_build(*, dry: bool):
    ad = CFG["AUTODOCK_GPU_DIR"]
    # 1) binaries & links
    _glob_delete(ad / "bin", ["autodock_gpu*", "adgpu_analysis*"], dry=dry)
    # 2) object/libs and generated headers
    pats = ["*.o", "*.obj", "*.a", "*.so", "*.dylib", "*.dll", "*.exp", "*.lib",
            "host/inc/stringify.h"]
    _glob_delete(ad, pats, dry=dry)
    # 3) make/cmake caches
    _glob_delete(ad, ["CMakeCache.txt", "cmake-build-*", "build", "config.log", "config.status"], dry=dry)

    # AutoGrid (if present)
    ag = ad / "autogrid"
    _glob_delete(ag, ["autogrid4", "*.o", "*.a", "*.so", "Makefile", "config.log", "config.status"], dry=dry)

def clean_results(*, dry: bool):
    _glob_delete(CFG["DOCKING_DIR"], ["*"], dry=dry)
    _glob_delete(CFG["ANALYSIS_DIR"], ["*"], dry=dry)

def clean_config(*, dry: bool):
    # remove this block:
    # if not _has_sentinel(DOCKING_DIR_DEFAULT):
    #     print("[SKIP] config: missing sentinel in docking/ — nothing deleted.")
    #     return
    print(CONFIG_FILE)
    if not _ensure_inside_repo(CONFIG_FILE):
        print("[SKIP] config: path outside repo guard.")
        return

    if CONFIG_FILE.exists() or CONFIG_FILE.is_symlink():
        _safe_delete(CONFIG_FILE, dry=dry)

    else:
        print("[INFO] config.py not found; nothing to do.")

    _glob_delete(DOCKING_DIR_DEFAULT / "__pycache__", ["config.*.pyc"], dry=dry, only_files=True)


def clean_maps(*, dry: bool):
    mm = CFG["MACRO_MOL_DIR"]
    if not mm.exists():
        return
    # Delete everything in MACRO_MOL_DIR that is NOT a receptor .pdbqt.
    # Generated outputs include:
    #   - per-receptor subdirs  (e.g. 5i6x_edited/ with S1/ … S6/ inside)
    #   - clash distance grids  (clash_dist.npz.npy, clash_dist.meta.txt)
    #   - hotspot centres       (centers.tsv)
    #   - autogrid products     (*.map, *.fld, *.gpf, *.glg)
    # All of these sit inside mm as children; we walk one level and skip .pdbqt.
    for child in sorted(mm.iterdir()):
        if child.suffix.lower() == ".pdbqt":
            print(f"[KEEP] receptor: {child.name}")
            continue
        _safe_delete(child, dry=dry)


def clean_pycache(*, dry: bool):
    _glob_delete(SCRIPT_DIR, ["__pycache__"], dry=dry)
    _glob_delete(SCRIPT_DIR, ["*.pyc", "*.pyo"], dry=dry, only_files=True)

def clean_db(*, dry: bool):
    _safe_delete(CFG["DB_PATH"], dry=dry)

def clean_ligands(*, dry: bool):
    # nukes the ligands folder contents (keep the folder)
    lig = CFG["LIGANDS_DIR"]
    if lig.exists():
        for p in lig.iterdir():
            _safe_delete(p, dry=dry)

def parse_args():
    ap = argparse.ArgumentParser(description="Ultidock cleaner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--build", action="store_true", help="Clean AutoDock-GPU & AutoGrid build artifacts")
    ap.add_argument("--results", action="store_true", help="Clean docking & analysis outputs")
    ap.add_argument("--maps", action="store_true", help="Clean generated map/grid files in MACRO_MOL_DIR")
    ap.add_argument("--pycache", action="store_true", help="Clean __pycache__ and *.pyc")
    ap.add_argument("--db", action="store_true", help="Delete results database file")
    ap.add_argument("--ligands", action="store_true", help="Delete all files in LIGANDS_DIR")
    ap.add_argument("--all", action="store_true", help="Do everything (build + results + maps + pycache + db + ligands)")
    ap.add_argument("-y", "--yes", action="store_true", help="Actually delete (otherwise perform a dry run)")
    ap.add_argument("--config", action="store_true", help="Delete docking/config.py (and its __pycache__ entry)")
    return ap.parse_args()

def main():
    args = parse_args()
    dry = not args.yes

    if args.all:
        args.build = args.results = args.maps = args.pycache = args.db = args.ligands = args.config = True

    # If no flags, default to a safe sweep that won’t destroy DB or ligands.
    if not any([args.build, args.results, args.maps, args.pycache, args.db, args.ligands]):
        args.build = args.results = args.maps = args.pycache = True

    # Guardrails: ensure targets are inside repo
    for k, v in CFG.items():
        if isinstance(v, Path) and k not in ("NUMWI",):
            if not _ensure_inside_repo(v) and v.exists():
                print(f"[ABORT] Refusing to touch path outside repo: {k} -> {v}")
                return 2

    print(f"== Ultidock clean  (dry-run={dry}) ==")
    print(f"Config:\n  AUTODOCK_GPU_DIR={CFG['AUTODOCK_GPU_DIR']}\n  DOCKING_DIR={CFG['DOCKING_DIR']}\n  ANALYSIS_DIR={CFG['ANALYSIS_DIR']}\n  MACRO_MOL_DIR={CFG['MACRO_MOL_DIR']}\n  RESULTS_DIR={CFG['RESULTS_DIR']}\n  LIGANDS_DIR={CFG['LIGANDS_DIR']}\n  DB_PATH={CFG['DB_PATH']}\n")

    if args.build:   clean_build(dry=dry)
    if args.results: clean_results(dry=dry)
    if args.maps:    clean_maps(dry=dry)
    if args.pycache: clean_pycache(dry=dry)
    if args.db:      clean_db(dry=dry)
    if args.ligands: clean_ligands(dry=dry)
    if args.config:  clean_config(dry=dry)

    clean_config(dry=dry)  # always clean config.py

    if dry:
        print("\n(DRY RUN) Nothing deleted. Re-run with -y to apply.")
    else:
        print("\nDone.")

if __name__ == "__main__":
    raise SystemExit(main())
