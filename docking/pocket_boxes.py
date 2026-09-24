"""Run local pocket finders and translate their sites into Ultidock boxes."""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.baseline.common import pdb_from_pdbqt  # noqa: E402
from benchmarks.baseline.fpocket.run_fpocket_baseline import (  # noqa: E402
    find_fpocket_output,
    parse_fpocket_output,
)
from benchmarks.baseline.p2rank.run_p2rank_baseline import (  # noqa: E402
    find_predictions_csv,
    parse_p2rank_predictions,
)

LOCAL_BINARIES = {
    "fpocket": REPO_ROOT / "external" / "fpocket" / "bin" / "fpocket",
    "p2rank": REPO_ROOT / "external" / "bin" / "prank",
}


def resolve_binary(method: str, override: Path | None = None) -> Path:
    path = (override or LOCAL_BINARIES[method]).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(
            f"{method} executable not found: {path}. Install it locally under external/ "
            "or pass --tool with an executable path."
        )
    return path


def create_pocket_boxes(
    *,
    method: str,
    receptor_pdbqt: Path,
    output_tsv: Path,
    work_dir: Path,
    box_size: float = 35.0,
    max_sites: int = 6,
    tool: Path | None = None,
) -> int:
    """Write fixed-size boxes centered on the finder's ranked pocket centers."""
    if method not in LOCAL_BINARIES:
        raise ValueError(f"Unknown pocket method: {method}")
    if not math.isfinite(box_size) or box_size <= 0:
        raise ValueError("box_size must be a positive finite number")
    if max_sites < 1:
        raise ValueError("max_sites must be at least 1")
    receptor_pdbqt = receptor_pdbqt.expanduser().resolve()
    if not receptor_pdbqt.is_file() or receptor_pdbqt.suffix.lower() != ".pdbqt":
        raise ValueError(f"Expected a receptor PDBQT file: {receptor_pdbqt}")
    binary = resolve_binary(method, tool)
    work_dir = work_dir.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = pdb_from_pdbqt(receptor_pdbqt, work_dir / f"{receptor_pdbqt.stem}.pdb")

    if method == "fpocket":
        raw_dir = work_dir / "fpocket"
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
        raw_dir.mkdir()
        staged = raw_dir / receptor_pdb.name
        shutil.copy2(receptor_pdb, staged)
        subprocess.run([str(binary), "-f", staged.name], cwd=raw_dir, check=True)
        sites = parse_fpocket_output(find_fpocket_output(raw_dir, staged))
    else:
        raw_dir = work_dir / "p2rank"
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
        raw_dir.mkdir()
        subprocess.run(
            [str(binary), "predict", "-f", str(receptor_pdb), "-o", str(raw_dir)],
            cwd=raw_dir,
            check=True,
        )
        sites = parse_p2rank_predictions(find_predictions_csv(raw_dir))

    if not sites:
        raise ValueError(f"{method} returned no pocket centers for {receptor_pdbqt}")

    spacing = 0.375
    npts = max(1, int(round(box_size / spacing)))
    if npts % 2 == 0:
        npts += 1
    output_tsv = output_tsv.expanduser().resolve()
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8") as handle:
        handle.write(
            "# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tscore\tfamily\n"
        )
        handle.write(f"# meta policy=provided site_count={min(len(sites), max_sites)} ranking_score={method}\n")
        for rank, site in enumerate(sites[:max_sites], start=1):
            x, y, z = site.center
            score = "" if site.score is None else f"{site.score:.4f}"
            handle.write(
                f"{receptor_pdbqt.stem}\tS{rank}\t{x:.3f}\t{y:.3f}\t{z:.3f}\t"
                f"{npts}\t{npts}\t{npts}\t{spacing:.3f}\t{score}\t{method}\n"
            )
    return min(len(sites), max_sites)
