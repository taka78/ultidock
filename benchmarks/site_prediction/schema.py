"""Shared data contracts for site-prediction benchmarks."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

Coordinate = tuple[float, float, float]
SITE_PREDICTION_METHODS = ("cav-emps", "fpocket", "p2rank")
METHOD_ALIASES = {
    "ultidock": "cav-emps",
    "cavemps": "cav-emps",
    "cav_emps": "cav-emps",
    "cav-emps": "cav-emps",
    "fpocket": "fpocket",
    "p2rank": "p2rank",
}


def normalize_method_name(method: str) -> str:
    """Return the canonical site-benchmark method ID."""

    key = method.strip().lower()
    return METHOD_ALIASES.get(key, key)


@dataclass(frozen=True)
class BindingSiteLabel:
    """One reference ligand-binding site used for post-hoc evaluation."""

    site_id: str
    ligand_id: str
    center: Coordinate
    ligand_atoms: tuple[Coordinate, ...]
    source: str = ""


@dataclass(frozen=True)
class TargetRecord:
    """One normalized benchmark target."""

    dataset: str
    target_id: str
    receptor_input: Path
    labels: tuple[BindingSiteLabel, ...]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class PredictedSite:
    """One method-predicted site center."""

    method: str
    target_id: str
    rank: int
    site_id: str
    center: Coordinate
    score: float | None = None
    source: str = ""


def _coordinate(value: Any, *, field: str) -> Coordinate:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{field} must be a 3-number coordinate")
    try:
        return (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a 3-number coordinate") from exc


def labels_from_json(path: Path) -> tuple[BindingSiteLabel, ...]:
    """Load normalized binding-site labels from ``labels.json``."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    raw_sites = payload.get("binding_sites")
    if not isinstance(raw_sites, list):
        raise ValueError(f"{path}: labels.json must contain a binding_sites list")

    labels: list[BindingSiteLabel] = []
    for index, raw_site in enumerate(raw_sites, start=1):
        if not isinstance(raw_site, dict):
            raise ValueError(f"{path}: binding_sites[{index}] must be an object")
        raw_atoms = raw_site.get("ligand_atoms")
        if not isinstance(raw_atoms, list) or not raw_atoms:
            raise ValueError(f"{path}: binding_sites[{index}].ligand_atoms must be non-empty")
        label = BindingSiteLabel(
            site_id=str(raw_site.get("site_id") or f"L{index}"),
            ligand_id=str(raw_site.get("ligand_id") or f"ligand_{index}"),
            center=_coordinate(raw_site.get("center"), field=f"binding_sites[{index}].center"),
            ligand_atoms=tuple(
                _coordinate(atom, field=f"binding_sites[{index}].ligand_atoms[]")
                for atom in raw_atoms
            ),
            source=str(raw_site.get("source") or ""),
        )
        labels.append(label)
    return tuple(labels)


def target_from_normalized_dir(dataset: str, target_dir: Path) -> TargetRecord:
    """Load one normalized target directory."""

    receptor = target_dir / "receptor_input.pdb"
    labels_path = target_dir / "labels.json"
    metadata_path = target_dir / "metadata.json"
    if not receptor.is_file():
        raise FileNotFoundError(f"missing normalized receptor: {receptor}")
    if not labels_path.is_file():
        raise FileNotFoundError(f"missing normalized labels: {labels_path}")
    metadata = (
        json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata_path.is_file()
        else {}
    )
    return TargetRecord(
        dataset=dataset,
        target_id=target_dir.name,
        receptor_input=receptor,
        labels=labels_from_json(labels_path),
        metadata=metadata,
    )


def prediction_rows_from_tsv(path: Path) -> list[PredictedSite]:
    """Read the normalized method prediction table."""

    import csv

    rows: list[PredictedSite] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(reader, start=1):
            rows.append(
                PredictedSite(
                    method=normalize_method_name(str(row.get("method") or "")),
                    target_id=str(row.get("target_id") or ""),
                    rank=int(row.get("rank") or index),
                    site_id=str(row.get("site_id") or f"S{index}"),
                    center=(
                        float(row["center_x"]),
                        float(row["center_y"]),
                        float(row["center_z"]),
                    ),
                    score=float(row["score"]) if row.get("score") not in (None, "") else None,
                    source=str(row.get("source") or ""),
                )
            )
    return rows
