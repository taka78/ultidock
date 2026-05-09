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
import shutil
import sys
import time
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
RETRYABLE_HTTP_STATUS = {408, 429, 500, 502, 503, 504, 520, 522, 524}

# Files to download per target (server filename --> local filename)
TARGET_FILES = {
    "receptor.pdb": "receptor.pdb",
    "crystal_ligand.mol2": "crystal_ligand.mol2",
    "actives_final.mol2.gz": "actives_final.mol2.gz",
    "decoys_final.mol2.gz": "decoys_final.mol2.gz",
}


def _is_retryable_download_error(exc: Exception) -> bool:
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code in RETRYABLE_HTTP_STATUS
    return isinstance(exc, (urllib.error.URLError, TimeoutError, OSError))


def _download_once(url: str, dest: Path, *, timeout: float) -> None:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Ultidock-DUD-E-downloader/1.0"},
    )
    tmp_dest = dest.with_name(f"{dest.name}.part")
    if tmp_dest.exists():
        tmp_dest.unlink()
    with urllib.request.urlopen(request, timeout=timeout) as response:
        with tmp_dest.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    tmp_dest.replace(dest)


def download_file(
    url: str,
    dest: Path,
    *,
    force: bool = False,
    retries: int = 5,
    retry_delay: float = 5.0,
    max_retry_delay: float = 60.0,
    timeout: float = 60.0,
) -> bool:
    """Download *url* to *dest*.  Returns True on success."""
    if dest.exists() and not force:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  ↓ {url}")
    attempts = max(1, retries + 1)
    for attempt in range(1, attempts + 1):
        try:
            _download_once(url, dest, timeout=timeout)
            if attempt > 1:
                print(f"    [OK] recovered on attempt {attempt}/{attempts}")
            return True
        except (urllib.error.URLError, urllib.error.HTTPError, OSError, TimeoutError) as exc:
            if dest.exists():
                dest.unlink()
            tmp_dest = dest.with_name(f"{dest.name}.part")
            if tmp_dest.exists():
                tmp_dest.unlink()
            retryable = _is_retryable_download_error(exc)
            if attempt >= attempts or not retryable:
                print(f"  [FAIL] FAILED {url}: {exc}")
                return False
            wait_s = min(max_retry_delay, retry_delay * (2 ** (attempt - 1)))
            print(
                f"    [retry {attempt}/{retries}] {exc}; "
                f"waiting {wait_s:.1f}s"
            )
            time.sleep(wait_s)
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


def download_target(
    name: str,
    dataset_root: Path,
    *,
    force: bool = False,
    retries: int = 5,
    retry_delay: float = 5.0,
    max_retry_delay: float = 60.0,
    timeout: float = 60.0,
) -> bool:
    """Download one DUD-E target.  Returns True on full success."""
    target_dir = dataset_root / name
    target_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    for remote_name, local_name in TARGET_FILES.items():
        url = f"{BASE_URL}/{name}/{remote_name}"
        if not download_file(
            url,
            target_dir / local_name,
            force=force,
            retries=retries,
            retry_delay=retry_delay,
            max_retry_delay=max_retry_delay,
            timeout=timeout,
        ):
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
        "--retries",
        type=int,
        default=5,
        help="Retry count per file for transient HTTP/network failures.",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=5.0,
        help="Initial retry delay in seconds; later attempts use exponential backoff.",
    )
    parser.add_argument(
        "--max-retry-delay",
        type=float,
        default=60.0,
        help="Maximum delay between retries in seconds.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Per-request timeout in seconds.",
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
    if args.retries < 0:
        raise SystemExit("--retries must be >= 0")
    if args.retry_delay < 0:
        raise SystemExit("--retry-delay must be >= 0")
    if args.max_retry_delay < 0:
        raise SystemExit("--max-retry-delay must be >= 0")
    if args.timeout <= 0:
        raise SystemExit("--timeout must be > 0")

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
        if download_target(
            name,
            dataset_root,
            force=args.force,
            retries=args.retries,
            retry_delay=args.retry_delay,
            max_retry_delay=args.max_retry_delay,
            timeout=args.timeout,
        ):
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
