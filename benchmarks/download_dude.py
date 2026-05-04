#!/usr/bin/env python3
"""Download DUD-E benchmark datasets from dude.docking.org.

Each target directory will contain:
  receptor.pdb           – the protein structure
  crystal_ligand.mol2    – the co-crystal reference ligand
  actives_final.mol2.gz  – active ligands (compressed multi-record MOL2)
  decoys_final.mol2.gz   – decoy ligands (compressed multi-record MOL2)
  benchmark.json         – metadata consumed by dude_docking_benchmark.py
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Complete list of all 102 DUD-E targets (gene-name slug used on the server)
# ---------------------------------------------------------------------------
ALL_DUDE_TARGETS: tuple[str, ...] = (
    "aa2ar", "abl1", "ace", "aces", "ada", "ada17", "adrb1", "adrb2",
    "akt1", "akt2", "aldr", "ampc", "andr", "aofb", "bace1", "braf",
    "cah2", "casp3", "cdk2", "comt", "cp2c9", "cp3a4", "csf1r", "cxcr4",
    "def", "dhi1", "dpp4", "drd3", "dyr", "egfr", "esr1", "esr2",
    "fa10", "fa7", "fabp4", "fak1", "fgfr1", "fkb1a", "fnta", "fpps",
    "gcr", "glcm", "gria2", "grik1", "hdac2", "hdac8", "hivint", "hivpr",
    "hivrt", "hmdh", "hs90a", "hxk4", "igf1r", "inha", "ital", "jak2",
    "kif11", "kit", "kith", "kpcb", "lck", "lkha4", "mapk2", "mcr",
    "met", "mk01", "mk10", "mk14", "mmp13", "mp2k1", "nos1", "nram",
    "pa2ga", "parp1", "pde5a", "pgh1", "pgh2", "plk1", "pnph", "ppara",
    "ppard", "pparg", "prgr", "ptn1", "pur2", "pygm", "pyrd", "reni",
    "rock1", "rxra", "sahh", "src", "tgfr1", "thb", "thrb", "try1",
    "tryb1", "tysy", "urok", "vgfr2", "wee1", "xiap",
)

BASE_URL = "http://dude.docking.org/targets"

# Files to download per target (server filename --> local filename)
TARGET_FILES = {
    "receptor.pdb": "receptor.pdb",
    "crystal_ligand.mol2": "crystal_ligand.mol2",
    "actives_final.mol2.gz": "actives_final.mol2.gz",
    "decoys_final.mol2.gz": "decoys_final.mol2.gz",
}


def download_file(url: str, dest: Path, *, force: bool = False) -> bool:
    """Download *url* to *dest*.  Returns True on success."""
    if dest.exists() and not force:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        print(f"  ↓ {url}")
        urllib.request.urlretrieve(url, str(dest))
        return True
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        print(f"  [FAIL] FAILED {url}: {exc}")
        if dest.exists():
            dest.unlink()
        return False


def write_benchmark_json(target_dir: Path) -> None:
    """Write a benchmark.json that dude_docking_benchmark.py can consume."""
    meta = {
        "receptor": "receptor.pdb",
        "reference_ligand": "crystal_ligand.mol2",
        "actives": "actives_final.mol2.gz",
        "decoys": "decoys_final.mol2.gz",
    }
    path = target_dir / "benchmark.json"
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


def download_target(name: str, dataset_root: Path, *, force: bool = False) -> bool:
    """Download one DUD-E target.  Returns True on full success."""
    target_dir = dataset_root / name
    target_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    for remote_name, local_name in TARGET_FILES.items():
        url = f"{BASE_URL}/{name}/{remote_name}"
        if not download_file(url, target_dir / local_name, force=force):
            ok = False
    if ok:
        write_benchmark_json(target_dir)
    return ok


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download DUD-E benchmark datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--targets",
        help=(
            "Comma-separated list of target names to download, "
            "or 'all' for all 102 targets. "
            "Default: all."
        ),
        default="all",
    )
    parser.add_argument(
        "--dataset-root",
        default=str(Path(__file__).resolve().parent / "datasets"),
        help="Root directory where target subdirectories will be created.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_targets",
        help="Just print the list of available targets and exit.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.list_targets:
        for name in ALL_DUDE_TARGETS:
            print(name)
        print(f"\n{len(ALL_DUDE_TARGETS)} targets total")
        return

    if args.targets.strip().lower() == "all":
        targets = list(ALL_DUDE_TARGETS)
    else:
        targets = [t.strip().lower() for t in args.targets.split(",") if t.strip()]
        unknown = [t for t in targets if t not in ALL_DUDE_TARGETS]
        if unknown:
            print(f"ERROR: unknown target(s): {', '.join(unknown)}", file=sys.stderr)
            print(f"  Available: {', '.join(ALL_DUDE_TARGETS)}", file=sys.stderr)
            sys.exit(1)

    dataset_root = Path(args.dataset_root).resolve()
    print(f"Downloading {len(targets)} DUD-E target(s) --> {dataset_root}\n")

    succeeded = 0
    failed_targets: list[str] = []
    for i, name in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] {name}")
        if download_target(name, dataset_root, force=args.force):
            succeeded += 1
        else:
            failed_targets.append(name)

    print(f"\n{'='*50}")
    print(f"Downloaded: {succeeded}/{len(targets)}")
    if failed_targets:
        print(f"Failed: {', '.join(failed_targets)}")
        sys.exit(1)
    print("All downloads completed successfully.")


if __name__ == "__main__":
    main()
