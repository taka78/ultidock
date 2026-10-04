"""Recover scored coordinates, independently of the legacy analysis database."""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
import xml.etree.ElementTree as ET


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@dataclass
class Pose:
    ligand: str
    receptor: str
    binding_site: str
    engine: str
    model: int
    score: float
    artifact: str
    artifact_sha256: str
    coordinates: list[str]
    receptor_sha256: str | None = None
    ligand_sha256: str | None = None

    def record(self) -> dict:
        return asdict(self)


def models(path: Path, *, dlg: bool = False) -> dict[int, list[str]]:
    blocks, block, number = {}, [], None
    for raw in path.read_text().splitlines():
        if dlg:
            if not raw.startswith("DOCKED: "):
                continue
            raw = raw[len("DOCKED: "):]
        if raw.startswith("MODEL"):
            if number is not None:
                raise ValueError(f"Nested MODEL in {path}")
            number, block = int(raw.split()[1]), [raw]
        elif number is not None:
            block.append(raw)
            if dlg and re.match(r"USER\s+Run\s*=", raw):
                number = int(raw.split("=")[1])
            if raw.startswith("ENDMDL"):
                if number in blocks:
                    raise ValueError(f"Duplicate model/run {number} in {path}")
                blocks[number] = block
                number, block = None, []
    if number is not None:
        raise ValueError(f"Truncated MODEL {number} in {path}")
    return blocks


def identity(path: Path, receptor: str, legacy: bool):
    sidecar = path.with_suffix(".md.json")
    if sidecar.exists():
        data = json.loads(sidecar.read_text())
        if data["receptor_id"] != receptor:
            return None
        if data["output_sha256"] != digest(path):
            raise ValueError(f"Docking artifact changed since its provenance was recorded: {path}")
        if path.suffix == ".xml" and "coordinates_sha256" in data:
            coordinates = path.with_suffix(".dlg")
            if not coordinates.is_file() or digest(coordinates) != data["coordinates_sha256"]:
                raise ValueError(f"Docking coordinates changed or missing since provenance was recorded: {coordinates}")
        return (data["ligand_id"], str(data["binding_site"]), data["receptor_sha256"],
                data["source_ligand_sha256"])
    match = re.fullmatch(r"(?:(.+)__)?S?(\d+)__(.+)_[0-9a-f]{8}", path.stem)
    if not match:
        raise ValueError(f"No unambiguous ligand/site provenance for {path.name}")
    owner, site, ligand = match.groups()
    if owner is None and not legacy:
        raise ValueError("Legacy filenames lack a receptor ID; use --legacy-single-receptor "
                         "only after checking this directory contains one receptor")
    if owner is not None and owner != receptor:
        return None
    return ligand, site, None, None


def select_poses(directory: Path, receptor: str, engine: str, top: int = 5,
                 legacy: bool = False, manifest: Path | None = None) -> list[Pose]:
    """Rank distinct ligand IDs within one receptor and one scoring engine."""
    if top < 1:
        raise ValueError("top must be positive")
    suffix = ".pdbqt" if engine == "vina" else ".xml"
    candidates = []
    paths = sorted(directory.glob(f"*{suffix}"))
    if manifest is not None:
        # An integrated run must never select an older score from a shared directory.
        records = json.loads(manifest.read_text())["outputs"]
        paths = []
        for record in records:
            path = Path(record["path"]).resolve()
            if path.parent != directory.resolve():
                raise ValueError(f"Docking manifest artifact is outside {directory}: {path}")
            if not path.is_file() or digest(path) != record["sha256"]:
                raise ValueError(f"Docking manifest artifact changed or missing: {path}")
            if path.suffix == suffix:
                paths.append(path)
        paths = sorted(set(paths))
    for path in paths:
        if path.stem.endswith("-best"):
            continue
        info = identity(path, receptor, legacy)
        if info is None:
            continue
        ligand, site, receptor_hash, ligand_hash = info
        if engine == "vina":
            blocks = models(path)
            scores = {}
            for number, block in blocks.items():
                values = [float(line.split()[3]) for line in block
                          if line.startswith("REMARK VINA RESULT:")]
                if len(values) != 1:
                    raise ValueError(f"Missing/duplicate score for model {number}: {path}")
                scores[number] = values[0]
        else:
            xml = ET.parse(path).getroot()
            scores = {}
            for run in xml.findall("./runs/run"):
                number = int(run.attrib["id"])
                if number in scores:
                    raise ValueError(f"Duplicate run {number}: {path}")
                scores[number] = float(run.findtext("free_NRG_binding"))
            dlg_path = path.with_suffix(".dlg")
            if not dlg_path.exists():
                raise ValueError(f"Exact AutoDock-GPU coordinates require {dlg_path}; "
                                 "a -best.pdbqt file is insufficient")
            blocks = models(dlg_path, dlg=True)
        if not scores or any(not math.isfinite(score) for score in scores.values()):
            raise ValueError(f"Missing or non-finite docking scores: {path}")
        number = min(scores, key=lambda n: (scores[n], n))
        if number not in blocks:
            raise ValueError(f"Scored run {number} has no coordinates: {path}")
        if engine == "adgpu":
            energies = [float(line.split("=")[1].split()[0]) for line in blocks[number]
                        if "Estimated Free Energy of Binding" in line]
            if len(energies) != 1 or abs(energies[0] - scores[number]) > 0.011:
                raise ValueError(f"XML/DLG energy mismatch for run {number}: {path}")
        coordinates = [line for line in blocks[number] if line.startswith(("ATOM  ", "HETATM"))]
        if not coordinates:
            raise ValueError(f"Empty pose {number}: {path}")
        artifact = path if engine == "vina" else path.with_suffix(".dlg")
        candidates.append(Pose(ligand, receptor, site, engine, number, scores[number],
                               str(artifact.resolve()), digest(artifact), coordinates,
                               receptor_hash, ligand_hash))
    hashes = {pose.receptor_sha256 for pose in candidates if pose.receptor_sha256}
    if len(hashes) > 1:
        raise ValueError("Different receptor preparations share this receptor ID; split the runs")
    for ligand in {pose.ligand for pose in candidates}:
        hashes = {pose.ligand_sha256 for pose in candidates
                  if pose.ligand == ligand and pose.ligand_sha256}
        if len(hashes) > 1:
            raise ValueError(f"Different input molecules share ligand ID {ligand}; split the runs")
    ordered = sorted(candidates, key=lambda p: (p.score, p.ligand, p.binding_site,
                                               p.model, p.artifact))
    unique = {}
    for pose in ordered:
        unique.setdefault(pose.ligand, pose)
    result = list(unique.values())[:top]
    if not result:
        raise ValueError(f"No scored {engine} poses for receptor {receptor} in {directory}")
    return result
