#!/usr/bin/env python3
"""Download and normalize COACH420/HOLO4K-style datasets for site prediction.

Reference labels are deliberately not defined as "every non-water HETATM
residue." That shortcut fragments multi-residue ligands, admits small or
distant cofactors, and changes both the Top-n cutoff and its denominator. The
normalizer therefore reproduces the ligand eligibility parameters distributed
with P2Rank's COACH420/HOLO4K benchmark data.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW_ROOT = REPO_ROOT / "benchmarks" / "site_prediction" / "datasets" / "raw"
DEFAULT_NORMALIZED_ROOT = REPO_ROOT / "benchmarks" / "site_prediction" / "datasets" / "normalized"
DEFAULT_MANIFEST_ROOT = REPO_ROOT / "benchmarks" / "site_prediction" / "manifests"
DEFAULT_SOURCE_ROOT = DEFAULT_RAW_ROOT / "p2rank-datasets"
P2RANK_DATASET_ARCHIVES = (
    "https://github.com/rdk/p2rank-datasets/archive/refs/heads/master.zip",
    "https://github.com/rdk/p2rank-datasets/archive/refs/heads/main.zip",
)

# This is benchmark ground truth, not a method-tuning knob. Changing these
# values creates a different label corpus and invalidates direct comparisons
# with runs generated under the P2Rank-compatible protocol.
P2RANK_IGNORED_HET_GROUPS = {
    "HOH",
    "DOD",
    "WAT",
    "NAG",
    "MAN",
    "UNK",
    "GLC",
    "ABA",
    "MPD",
    "GOL",
    "SO4",
    "PO4",
}
P2RANK_LIGAND_CLUSTERING_DISTANCE_A = 1.7
P2RANK_MIN_LIGAND_HEAVY_ATOMS = 5
P2RANK_LIGAND_PROTEIN_CONTACT_DISTANCE_A = 4.0
P2RANK_LIGAND_CENTER_PROTEIN_DISTANCE_A = 5.5
LABEL_PROTOCOL = "p2rank_relevant_ligands_v1"


def _looks_like_p2rank_dataset_root(path: Path) -> bool:
    return path.is_dir() and (
        any(path.rglob("coach420.ds"))
        or any(path.rglob("holo4k.ds"))
        or any(path.rglob("*coach420*.ds"))
        or any(path.rglob("*holo4k*.ds"))
    )


def ensure_dataset_source(
    *,
    source_root: Path | None,
    raw_root: Path = DEFAULT_RAW_ROOT,
    download: bool = True,
    dry_run: bool = False,
) -> Path:
    """Return a local p2rank-datasets checkout, downloading it if needed."""

    raw_root = raw_root.resolve()
    target_root = (source_root.resolve() if source_root else raw_root / "p2rank-datasets")
    if _looks_like_p2rank_dataset_root(target_root):
        return target_root
    if source_root is not None and source_root.exists():
        raise FileNotFoundError(f"{source_root} exists but does not look like rdk/p2rank-datasets")
    if not download:
        raise FileNotFoundError(
            f"p2rank-datasets not found at {target_root}. Re-run with --download "
            "or pass --source-root to an existing checkout."
        )
    if dry_run:
        return target_root

    raw_root = raw_root.resolve()
    raw_root.mkdir(parents=True, exist_ok=True)
    archive_path = raw_root / "p2rank-datasets.zip"
    last_error: Exception | None = None
    for url in P2RANK_DATASET_ARCHIVES:
        try:
            print(f"[download] {url}")
            urllib.request.urlretrieve(url, archive_path)
            break
        except Exception as exc:
            last_error = exc
    else:
        raise RuntimeError(
            "Could not download rdk/p2rank-datasets. Check network access or "
            "clone it manually and pass --source-root."
        ) from last_error

    extract_root = raw_root / "_p2rank_datasets_extract"
    if extract_root.exists():
        shutil.rmtree(extract_root)
    extract_root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as archive:
        archive.extractall(extract_root)

    candidates = [path for path in extract_root.iterdir() if path.is_dir()]
    candidates.extend(path for path in extract_root.rglob("*") if path.is_dir())
    source_candidate = next((path for path in candidates if _looks_like_p2rank_dataset_root(path)), None)
    if source_candidate is None:
        raise RuntimeError(f"Downloaded archive did not contain COACH420/HOLO4K dataset files: {archive_path}")
    if target_root.exists():
        shutil.rmtree(target_root)
    target_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source_candidate), str(target_root))
    shutil.rmtree(extract_root, ignore_errors=True)
    return target_root


def _pdb_candidates_from_ds(source_root: Path, dataset: str) -> list[Path]:
    ds_files = sorted(source_root.rglob(f"{dataset}.ds"))
    if not ds_files:
        ds_files = sorted(source_root.rglob(f"*{dataset}*.ds"))
    candidates: list[Path] = []
    for ds_file in ds_files:
        for raw in ds_file.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            for token in line.replace(",", " ").split():
                if token.lower().endswith((".pdb", ".pdb.gz")):
                    path = Path(token)
                    candidates.append(path if path.is_absolute() else (ds_file.parent / path))
    return candidates


def _scan_pdbs(source_root: Path, dataset: str) -> list[Path]:
    candidates = _pdb_candidates_from_ds(source_root, dataset)
    if not candidates:
        dataset_dirs = [path for path in source_root.rglob("*") if path.is_dir() and dataset in path.name.lower()]
        roots = dataset_dirs or [source_root]
        for root in roots:
            candidates.extend(sorted(root.rglob("*.pdb")))
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved.exists() and resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return sorted(unique)


def _coord(line: str) -> tuple[float, float, float] | None:
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        return None


def _ligand_key(line: str) -> tuple[str, str, str, str]:
    return (
        line[17:20].strip(),
        line[21].strip() or "_",
        line[22:26].strip(),
        line[26].strip() or "",
    )


def _element(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    atom_name = line[12:16].strip().upper()
    while atom_name and atom_name[0].isdigit():
        atom_name = atom_name[1:]
    return atom_name[:1]


def _is_heavy_atom(line: str) -> bool:
    return _element(line) not in {"H", "D"}


def _atom_identity(line: str) -> tuple[str, str, str, str, str, str]:
    return (
        line[:6].strip(),
        line[12:16].strip(),
        line[17:20].strip(),
        line[21].strip() or "_",
        line[22:26].strip(),
        line[26].strip() or "",
    )


def _occupancy(line: str) -> float:
    try:
        return float(line[54:60])
    except ValueError:
        return 0.0


def _select_altloc_records(lines: list[str]) -> list[str]:
    """Select one deterministic conformer for each PDB atom identity."""

    selected: dict[tuple[str, str, str, str, str, str], tuple[tuple[int, float, int], int, str]] = {}
    for order, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        altloc = line[16].strip().upper()
        priority = (
            1 if not altloc else 0,
            _occupancy(line),
            1 if altloc == "A" else 0,
        )
        key = _atom_identity(line)
        current = selected.get(key)
        if current is None or priority > current[0]:
            selected[key] = (priority, order, line)
    return [entry[2] for entry in sorted(selected.values(), key=lambda entry: entry[1])]


def _logical_pdb_lines(lines: list[str]) -> tuple[list[str], int]:
    """Split malformed ``TERHETATM`` records present in the source corpus."""

    logical: list[str] = []
    repaired = 0
    for line in lines:
        if line.startswith("TERHETATM"):
            logical.extend(("TER", line[3:]))
            repaired += 1
        else:
            logical.append(line)
    return logical, repaired


def _squared_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((a[axis] - b[axis]) ** 2 for axis in range(3))


def _within_distance(
    query: list[tuple[float, float, float]],
    reference: list[tuple[float, float, float]],
    threshold_a: float,
) -> bool:
    threshold_sq = float(threshold_a) ** 2
    return any(_squared_distance(a, b) <= threshold_sq for a in query for b in reference)


def _cluster_ligand_groups(
    groups: list[tuple[tuple[str, str, str, str], list[tuple[float, float, float]]]],
    threshold_a: float,
) -> list[list[int]]:
    """Join covalently connected HET groups using P2Rank's distance cutoff."""

    if not groups:
        return []
    parent = list(range(len(groups)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    cell_size = float(threshold_a)
    threshold_sq = cell_size**2
    buckets: dict[tuple[int, int, int], list[tuple[int, tuple[float, float, float]]]] = defaultdict(list)
    for group_index, (_, atoms) in enumerate(groups):
        for atom in atoms:
            cell = tuple(math.floor(value / cell_size) for value in atom)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for other_index, other_atom in buckets.get(
                            (cell[0] + dx, cell[1] + dy, cell[2] + dz), []
                        ):
                            if other_index != group_index and _squared_distance(atom, other_atom) <= threshold_sq:
                                union(group_index, other_index)
            buckets[cell].append((group_index, atom))

    components: dict[int, list[int]] = defaultdict(list)
    for group_index in range(len(groups)):
        components[find(group_index)].append(group_index)
    return sorted(components.values(), key=lambda component: min(component))


def _ligand_id(keys: list[tuple[str, str, str, str]]) -> str:
    return "+".join(
        f"{resname}:{chain}:{resseq}{icode}"
        for resname, chain, resseq, icode in keys
    )


def _normalized_pdb_data(pdb_path: Path, *, dataset: str) -> tuple[list[str], list[dict[str, object]], dict[str, object]]:
    receptor_lines: list[str] = []
    source_lines = pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    raw_lines, repaired_concatenated_records = _logical_pdb_lines(source_lines)
    selected_atoms = _select_altloc_records(raw_lines)
    protein_atoms: list[tuple[float, float, float]] = []
    ligand_groups: dict[tuple[str, str, str, str], list[tuple[float, float, float]]] = defaultdict(list)
    ignored_groups: set[tuple[str, str, str, str]] = set()

    for raw in raw_lines:
        if raw.startswith("ATOM"):
            receptor_lines.append(raw)
        elif raw.startswith(("TER", "END")):
            receptor_lines.append(raw)

    for raw in selected_atoms:
        coord = _coord(raw)
        if coord is None or not _is_heavy_atom(raw):
            continue
        if raw.startswith("ATOM"):
            protein_atoms.append(coord)
        elif raw.startswith("HETATM"):
            resname = raw[17:20].strip().upper()
            # Exclude ignored groups before clustering so a glycan/cofactor
            # cannot bridge otherwise independent ligand components.
            if resname in P2RANK_IGNORED_HET_GROUPS:
                ignored_groups.add(_ligand_key(raw))
                continue
            ligand_groups[_ligand_key(raw)].append(coord)

    if not receptor_lines:
        raise ValueError(f"{pdb_path}: no ATOM records")
    if not protein_atoms:
        raise ValueError(f"{pdb_path}: no protein heavy atoms")

    ordered_groups = sorted(ligand_groups.items())
    components = _cluster_ligand_groups(
        ordered_groups,
        P2RANK_LIGAND_CLUSTERING_DISTANCE_A,
    )

    labels: list[dict[str, object]] = []
    excluded_small = 0
    excluded_no_contact = 0
    excluded_distant_center = 0
    for component in components:
        # One connected component is one reference ligand. Treating each
        # residue as a site inflates the label count and changes Top-n itself.
        keys = [ordered_groups[group_index][0] for group_index in component]
        atoms = [
            atom
            for group_index in component
            for atom in ordered_groups[group_index][1]
        ]
        if len(atoms) < P2RANK_MIN_LIGAND_HEAVY_ATOMS:
            excluded_small += 1
            continue
        center = [sum(atom[axis] for atom in atoms) / len(atoms) for axis in range(3)]
        if not _within_distance(atoms, protein_atoms, P2RANK_LIGAND_PROTEIN_CONTACT_DISTANCE_A):
            excluded_no_contact += 1
            continue
        if not _within_distance([tuple(center)], protein_atoms, P2RANK_LIGAND_CENTER_PROTEIN_DISTANCE_A):
            excluded_distant_center += 1
            continue
        index = len(labels) + 1
        labels.append(
            {
                "site_id": f"L{index}",
                "ligand_id": _ligand_id(keys),
                "center": center,
                "ligand_atoms": atoms,
                "source": f"{pdb_path.name} P2Rank-compatible relevant HETATM ligand",
                "het_groups": [
                    {
                        "resname": resname,
                        "chain": chain,
                        "resseq": resseq,
                        "icode": icode,
                    }
                    for resname, chain, resseq, icode in keys
                ],
                "heavy_atom_count": len(atoms),
            }
        )

    metadata = {
        "dataset": dataset,
        "target_id": pdb_path.stem.lower(),
        "source_pdb": str(pdb_path.resolve()),
        "normalization": (
            "ATOM records as receptor_input; P2Rank-compatible relevant HETATM "
            "components as ligand labels"
        ),
        "label_protocol": LABEL_PROTOCOL,
        "label_parameters": {
            "ignore_het_groups": sorted(P2RANK_IGNORED_HET_GROUPS),
            "ligand_clustering_distance_a": P2RANK_LIGAND_CLUSTERING_DISTANCE_A,
            "min_ligand_heavy_atoms": P2RANK_MIN_LIGAND_HEAVY_ATOMS,
            "ligand_protein_contact_distance_a": P2RANK_LIGAND_PROTEIN_CONTACT_DISTANCE_A,
            "ligand_center_protein_distance_a": P2RANK_LIGAND_CENTER_PROTEIN_DISTANCE_A,
            "altloc_selection": "blank_then_highest_occupancy_then_A",
        },
        "label_filter_diagnostics": {
            "repaired_concatenated_ter_hetatm_records": repaired_concatenated_records,
            "ignored_het_groups": len(ignored_groups),
            "candidate_het_groups": len(ordered_groups),
            "clustered_ligand_components": len(components),
            "excluded_small_components": excluded_small,
            "excluded_no_contact_components": excluded_no_contact,
            "excluded_distant_center_components": excluded_distant_center,
        },
        "n_reference_sites": len(labels),
    }
    return receptor_lines, labels, metadata


def normalize_one_pdb(
    pdb_path: Path,
    output_dir: Path,
    *,
    dataset: str,
    labels_only: bool = False,
) -> dict[str, object]:
    receptor_lines, labels, metadata = _normalized_pdb_data(pdb_path, dataset=dataset)

    target_id = pdb_path.stem.lower()
    target_dir = output_dir / dataset / target_id
    target_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = target_dir / "receptor_input.pdb"
    if labels_only:
        # Label repair must never perturb the receptor seen by prediction
        # methods; that would turn reevaluation into a new prediction run.
        if not receptor_path.is_file():
            raise FileNotFoundError(f"labels-only normalization requires existing receptor: {receptor_path}")
        metadata["receptor_input_preserved"] = True
    else:
        receptor_path.write_text("\n".join(receptor_lines) + "\n", encoding="utf-8")
    (target_dir / "labels.json").write_text(
        json.dumps(
            {
                "label_protocol": LABEL_PROTOCOL,
                "binding_sites": labels,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (target_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def prepare_dataset(
    *,
    dataset: str,
    source_root: Path,
    normalized_root: Path,
    manifest_root: Path,
    max_targets: int | None = None,
    dry_run: bool = False,
    labels_only: bool = False,
) -> list[dict[str, object]]:
    """Normalize one dataset and return manifest rows."""

    pdbs = _scan_pdbs(source_root, dataset)
    if max_targets is not None:
        pdbs = pdbs[:max_targets]
    if dry_run:
        print(f"Dataset: {dataset}")
        print(f"Source root: {source_root}")
        print(f"PDB candidates: {len(pdbs)}")
        for path in pdbs[:20]:
            print(f"  - {path}")
        return []
    if not pdbs:
        raise FileNotFoundError(f"No PDB files found for {dataset} under {source_root}")

    if labels_only:
        # Validate the complete corpus before the first in-place write. This
        # prevents a failed repair from mixing old and new label protocols.
        preflight_failures: list[str] = []
        for pdb_path in pdbs:
            try:
                _normalized_pdb_data(pdb_path, dataset=dataset)
                receptor_path = normalized_root / dataset / pdb_path.stem.lower() / "receptor_input.pdb"
                if not receptor_path.is_file():
                    raise FileNotFoundError(f"missing existing receptor: {receptor_path}")
            except Exception as exc:
                preflight_failures.append(f"{pdb_path}: {exc}")
        if preflight_failures:
            details = "\n".join(f"  - {failure}" for failure in preflight_failures)
            raise RuntimeError(
                "Labels-only preflight failed; no labels were regenerated:\n"
                f"{details}"
            )

    manifest_rows = []
    for pdb_path in pdbs:
        try:
            row = normalize_one_pdb(
                pdb_path,
                normalized_root,
                dataset=dataset,
                labels_only=labels_only,
            )
        except Exception as exc:
            print(f"[FAIL] {pdb_path}: {exc}", file=sys.stderr)
            continue
        if row is not None:
            manifest_rows.append(row)
            print(f"[OK] {row['target_id']} ({row['n_reference_sites']} site(s))")

    manifest_root.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_root / f"{dataset}.jsonl"
    manifest_path.write_text(
        "\n".join(json.dumps(row, sort_keys=True) for row in manifest_rows) + "\n",
        encoding="utf-8",
    )
    print(f"Normalized targets: {len(manifest_rows)}")
    print(f"Manifest: {manifest_path}")
    return manifest_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download and normalize COACH420/HOLO4K from rdk/p2rank-datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--dataset", choices=["coach420", "holo4k"], required=True)
    parser.add_argument("--source-root", type=Path, help="Existing p2rank-datasets checkout.")
    parser.add_argument("--raw-root", default=str(DEFAULT_RAW_ROOT), type=Path)
    parser.add_argument("--normalized-root", default=str(DEFAULT_NORMALIZED_ROOT), type=Path)
    parser.add_argument("--manifest-root", default=str(DEFAULT_MANIFEST_ROOT), type=Path)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--download", dest="download", action="store_true", default=True)
    parser.add_argument("--no-download", dest="download", action="store_false")
    parser.add_argument(
        "--labels-only",
        action="store_true",
        help="Regenerate labels/metadata while preserving existing receptor_input.pdb files.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    source_root = ensure_dataset_source(
        source_root=args.source_root,
        raw_root=args.raw_root,
        download=bool(args.download),
        dry_run=bool(args.dry_run),
    )
    if args.dry_run and not source_root.exists():
        print(f"Dataset source would be downloaded to: {source_root}")
        return
    prepare_dataset(
        dataset=args.dataset,
        source_root=source_root.resolve(),
        normalized_root=args.normalized_root.resolve(),
        manifest_root=args.manifest_root.resolve(),
        max_targets=args.max_targets,
        dry_run=args.dry_run,
        labels_only=bool(args.labels_only),
    )


if __name__ == "__main__":
    main()
