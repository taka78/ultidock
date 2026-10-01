#!/usr/bin/env python3
"""Run and report the preregistered PDB 8DD2 GABAA site-recovery case study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


EXAMPLE_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXAMPLE_DIR.parents[1]
PREREGISTRATION_PATH = EXAMPLE_DIR / "preregistration.json"
PREDICT_SCRIPT = REPO_ROOT / "benchmarks" / "site_prediction" / "predict_sites.py"
DATASET_ID = "gabaa_8dd2"
TARGET_ID = "8dd2"
RCSB_PDB_URL = "https://files.rcsb.org/download/8DD2.pdb"
DEFAULT_AUTOGRID4 = (
    REPO_ROOT / "docking" / "AUTODOCK_GPU_DIR" / "autogrid" / "autogrid4"
)
DEFAULT_FPOCKET = REPO_ROOT / "external" / "bin" / "fpocket"
DEFAULT_P2RANK = REPO_ROOT / "external" / "bin" / "prank"
BUFFER_VALUES = (0.0, 2.0, 4.0)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _load_preregistration() -> dict[str, Any]:
    return json.loads(PREREGISTRATION_PATH.read_text(encoding="utf-8"))


def _download_pdb(destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            print(f"[download] 8DD2 from RCSB (attempt {attempt}/3)")
            request = urllib.request.Request(
                RCSB_PDB_URL,
                headers={"User-Agent": "Ultidock-GABAA-case-study/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                with temporary.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
            temporary.replace(destination)
            return
        except Exception as exc:  # pragma: no cover - network behavior is external
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(2**attempt)
    raise RuntimeError(f"Could not download {RCSB_PDB_URL}: {last_error}")


def _stage_source_pdb(
    output_dir: Path,
    supplied_pdb: Path | None,
    *,
    force: bool,
    allow_source_mismatch: bool,
    preregistration: dict[str, Any],
) -> Path:
    source_pdb = output_dir / "source" / "8DD2.pdb"
    source_pdb.parent.mkdir(parents=True, exist_ok=True)
    if supplied_pdb is not None:
        supplied_pdb = supplied_pdb.expanduser().resolve()
        if not supplied_pdb.is_file():
            raise FileNotFoundError(f"8DD2 PDB not found: {supplied_pdb}")
        if supplied_pdb != source_pdb and (
            force or not source_pdb.exists() or _sha256(supplied_pdb) != _sha256(source_pdb)
        ):
            shutil.copy2(supplied_pdb, source_pdb)
    elif force or not source_pdb.exists():
        _download_pdb(source_pdb)

    expected = str(preregistration["structure"]["rcsb_revision_sha256"])
    observed = _sha256(source_pdb)
    if observed != expected and not allow_source_mismatch:
        raise RuntimeError(
            "8DD2 source hash does not match the preregistered RCSB revision: "
            f"expected {expected}, found {observed}. Pass --allow-source-mismatch "
            "only after documenting the PDB revision difference."
        )
    return source_pdb


def _pdb_float(line: str, start: int, end: int, default: float = 0.0) -> float:
    try:
        return float(line[start:end])
    except ValueError:
        return default


def _atom_selection_key(line: str) -> tuple[str, str, str, str, str]:
    return (
        line[21:22],
        line[22:26],
        line[26:27],
        line[17:20],
        line[12:16],
    )


def _altloc_priority(line: str) -> tuple[float, int, str]:
    altloc = line[16:17]
    preferred = 2 if altloc == " " else 1 if altloc == "A" else 0
    return (_pdb_float(line, 54, 60), preferred, chr(255 - ord(altloc or " ")))


def _selected_receptor_atoms(lines: Iterable[str], chains: set[str]) -> list[str]:
    """Select one deterministic conformer while retaining only receptor ATOM records."""
    selected: dict[tuple[str, str, str, str, str], tuple[tuple[float, int, str], int, str]] = {}
    for order, line in enumerate(lines):
        if not line.startswith("ATOM  ") or line[21:22] not in chains:
            continue
        key = _atom_selection_key(line)
        candidate = (_altloc_priority(line), -order, line)
        if key not in selected or candidate[:2] > selected[key][:2]:
            selected[key] = candidate
    kept = [entry[2] for entry in sorted(selected.values(), key=lambda item: -item[1])]
    normalized = []
    for line in kept:
        if line[16:17] != " ":
            line = line[:16] + " " + line[17:]
        normalized.append(line)
    return normalized


def _atom_coordinate(line: str) -> tuple[float, float, float]:
    return (
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54]),
    )


def _ligand_lines_for_site(lines: Iterable[str], definition: dict[str, Any]) -> list[str]:
    result = []
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != definition["pdb_resname"]:
            continue
        if line[21:22] != definition["pdb_chain"]:
            continue
        if int(line[22:26]) != int(definition["pdb_resseq"]):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        result.append(line)
    return result


def _coordinate_center(atoms: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    count = float(len(atoms))
    return tuple(sum(atom[axis] for atom in atoms) / count for axis in range(3))  # type: ignore[return-value]


def prepare_case_study(
    output_dir: Path,
    source_pdb: Path,
    preregistration: dict[str, Any],
) -> tuple[Path, Path]:
    """Create ligand-free predictor input and a physically separate truth directory."""
    raw_lines = source_pdb.read_text(encoding="utf-8", errors="ignore").splitlines()
    receptor_chains = set(preregistration["structure"]["predictor_chains"])
    receptor_atoms = _selected_receptor_atoms(raw_lines, receptor_chains)
    found_chains = {line[21:22] for line in receptor_atoms}
    if found_chains != receptor_chains:
        raise RuntimeError(
            f"Expected receptor chains {sorted(receptor_chains)}, found {sorted(found_chains)}"
        )

    target_dir = output_dir / "normalized" / DATASET_ID / TARGET_ID
    target_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = target_dir / "receptor_input.pdb"
    receptor_output = [
        "REMARK 900 ULTIDOCK GABAA 8DD2 BLIND SITE-RECOVERY INPUT",
        "REMARK 900 ONLY HUMAN RECEPTOR CHAINS A B C D E ARE PRESENT",
        "REMARK 900 ALL HETATM, GLYCAN, FAB AND GROUND-TRUTH LIGAND RECORDS REMOVED",
    ]
    previous_chain = None
    for line in receptor_atoms:
        chain = line[21:22]
        if previous_chain is not None and chain != previous_chain:
            receptor_output.append("TER")
        receptor_output.append(line)
        previous_chain = chain
    receptor_output.extend(["TER", "END"])
    receptor_path.write_text("\n".join(receptor_output) + "\n", encoding="ascii")

    ground_truth_dir = output_dir / "ground_truth"
    ligand_dir = ground_truth_dir / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    labels = []
    for definition in preregistration["ground_truth_sites"]:
        ligand_lines = _ligand_lines_for_site(raw_lines, definition)
        if not ligand_lines:
            raise RuntimeError(
                "Could not extract preregistered ligand "
                f"{definition['pdb_resname']}:{definition['pdb_chain']}:{definition['pdb_resseq']}"
            )
        atoms = [_atom_coordinate(line) for line in ligand_lines]
        expected_count = 7 if definition["pdb_resname"] == "ABU" else 23
        if len(atoms) != expected_count:
            raise RuntimeError(
                f"{definition['site_id']}: expected {expected_count} heavy atoms, found {len(atoms)}"
            )
        ligand_path = ligand_dir / f"{definition['site_id']}.pdb"
        ligand_path.write_text("\n".join([*ligand_lines, "END"]) + "\n", encoding="ascii")
        labels.append(
            {
                **definition,
                "ligand_id": (
                    f"{definition['pdb_resname']}:{definition['pdb_chain']}:"
                    f"{definition['pdb_resseq']}"
                ),
                "center": [round(value, 6) for value in _coordinate_center(atoms)],
                "ligand_atoms": [[round(value, 6) for value in atom] for atom in atoms],
                "atom_count": len(atoms),
                "source": "PDB 8DD2 deposited ligand coordinates; withheld from predictors",
                "extracted_ligand": str(ligand_path.resolve()),
            }
        )

    labels_path = ground_truth_dir / "labels.json"
    _write_json(
        labels_path,
        {
            "case_study_id": preregistration["case_study_id"],
            "source_pdb": str(source_pdb.resolve()),
            "source_pdb_sha256": _sha256(source_pdb),
            "binding_sites": labels,
        },
    )

    chain_atoms: dict[str, int] = {chain: 0 for chain in sorted(receptor_chains)}
    chain_residues: dict[str, set[tuple[str, str, str]]] = {
        chain: set() for chain in sorted(receptor_chains)
    }
    for line in receptor_atoms:
        chain = line[21:22]
        chain_atoms[chain] += 1
        chain_residues[chain].add((line[17:20], line[22:26], line[26:27]))
    metadata = {
        "dataset": DATASET_ID,
        "target_id": TARGET_ID,
        "pdb_id": "8DD2",
        "source_pdb_sha256": _sha256(source_pdb),
        "receptor_input_sha256": _sha256(receptor_path),
        "predictor_has_ground_truth_ligands": False,
        "predictor_record_types": ["ATOM"],
        "receptor_chains": sorted(receptor_chains),
        "chain_subunits": preregistration["structure"]["subunits"],
        "chain_atom_counts": chain_atoms,
        "chain_residue_counts": {
            chain: len(residues) for chain, residues in chain_residues.items()
        },
        "excluded": preregistration["structure"]["excluded_from_predictor"],
        "ground_truth_path_withheld_from_predictors": str(labels_path.resolve()),
        "prepared_at_utc": _utc_now(),
    }
    _write_json(target_dir / "metadata.json", metadata)
    return receptor_path, labels_path


def _method_list(value: str) -> list[str]:
    methods = [part.strip().lower() for part in value.split(",") if part.strip()]
    invalid = [method for method in methods if method not in {"cav-emps", "fpocket", "p2rank"}]
    if invalid:
        raise ValueError(f"Unknown method(s): {', '.join(invalid)}")
    return list(dict.fromkeys(methods))


def _command_text(command: list[str]) -> str:
    return shlex.join(command)


def _run_recorded(command: list[str], *, cwd: Path, log_path: Path) -> None:
    started = _utc_now()
    print(f"\n$ {_command_text(command)}", flush=True)
    completed = subprocess.run(command, cwd=cwd, check=False)
    record = {
        "started_at_utc": started,
        "finished_at_utc": _utc_now(),
        "cwd": str(cwd.resolve()),
        "command": command,
        "command_text": _command_text(command),
        "returncode": completed.returncode,
    }
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)


def run_predictions(
    output_dir: Path,
    args: argparse.Namespace,
    preregistration: dict[str, Any],
) -> dict[str, str]:
    normalized_root = output_dir / "normalized"
    prediction_root = output_dir / "predictions"
    command_log = output_dir / "commands.jsonl"
    common = [
        sys.executable,
        str(PREDICT_SCRIPT),
        "--dataset",
        DATASET_ID,
        "--normalized-root",
        str(normalized_root),
        "--output-dir",
        str(prediction_root),
        "--targets",
        TARGET_ID,
        "--autosites",
        str(preregistration["prediction_protocol"]["portfolio_size"]),
        "--jobs",
        str(max(1, int(args.jobs))),
        "--seed",
        str(args.seed),
    ]
    if args.force:
        common.append("--force")

    failures = {}
    for method in args.methods:
        method_output = prediction_root / method
        if args.force and method_output.exists():
            shutil.rmtree(method_output)
        command = [*common, "--method", method]
        if method == "cav-emps":
            command.extend(
                [
                    "--autogrid4-bin",
                    str(Path(args.autogrid4_bin).expanduser()),
                    "--keep-artifacts",
                    "--cav-emps-ablation-profiles",
                    "all",
                ]
            )
            if args.receptor_prepare_command:
                command.extend(["--receptor-prepare-command", args.receptor_prepare_command])
        elif method == "fpocket":
            command.extend(["--fpocket-cmd", args.fpocket_cmd])
        elif method == "p2rank":
            command.extend(["--p2rank-cmd", args.p2rank_cmd])
        try:
            _run_recorded(command, cwd=REPO_ROOT, log_path=command_log)
        except subprocess.CalledProcessError as exc:
            failures[method] = (
                f"prediction command exited with status {exc.returncode}: "
                f"{_command_text(command)}"
            )
            print(f"[FAIL] {method}: {failures[method]}", file=sys.stderr)
    failure_path = output_dir / "prediction_failures.json"
    if failures:
        _write_json(failure_path, failures)
    else:
        failure_path.unlink(missing_ok=True)
    return failures


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _expected_configurations(
    methods: list[str], preregistration: dict[str, Any]
) -> list[str]:
    configurations = []
    if "cav-emps" in methods:
        configurations.extend(
            f"cav-emps-{profile}"
            for profile in preregistration["prediction_protocol"]["cav_emps_profiles"]
        )
    configurations.extend(method for method in ("fpocket", "p2rank") if method in methods)
    return configurations


def _collect_prediction_rows(
    output_dir: Path,
    methods: list[str],
    preregistration: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    rows = []
    missing: dict[str, str] = {}
    prediction_root = output_dir / "predictions"
    failure_path = output_dir / "prediction_failures.json"
    if failure_path.is_file():
        missing.update(json.loads(failure_path.read_text(encoding="utf-8")))
    for method in methods:
        path = prediction_root / method / "predictions.tsv"
        method_rows = _read_tsv(path)
        if not path.is_file():
            missing.setdefault(method, f"missing predictions file: {path}")
        elif not method_rows:
            missing.setdefault(method, f"empty predictions file: {path}")
        for row in method_rows:
            if row.get("target_id") != TARGET_ID:
                continue
            normalized = dict(row)
            normalized["rank"] = int(row["rank"])
            for key in ("center_x", "center_y", "center_z"):
                normalized[key] = float(row[key])
            normalized["score"] = (
                None if row.get("score") in (None, "") else float(row["score"])
            )
            normalized["prediction_file"] = str(path.resolve())
            rows.append(normalized)

    rows.sort(key=lambda row: (str(row["method"]), int(row["rank"]), str(row["site_id"])))
    expected = _expected_configurations(methods, preregistration)
    emitted = {str(row["method"]) for row in rows}
    for configuration in expected:
        if configuration not in emitted:
            parent = "cav-emps" if configuration.startswith("cav-emps-") else configuration
            missing.setdefault(configuration, missing.get(parent, "no prediction rows emitted"))
    return rows, missing


def _parse_centers_tsv(path: Path) -> tuple[list[dict[str, Any]], dict[str, str]]:
    header: list[str] | None = None
    metadata: dict[str, str] = {}
    rows = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if raw.startswith("# meta "):
            for token in raw[len("# meta ") :].split():
                if "=" in token:
                    key, value = token.split("=", 1)
                    metadata[key] = value
            continue
        if raw.startswith("# receptor"):
            header = [part.strip().lstrip("# ") for part in raw.split("\t")]
            continue
        if not raw or raw.startswith("#"):
            continue
        if header is None:
            raise ValueError(f"{path}: data row appears before centers header")
        values = raw.split("\t")
        row = {key: value for key, value in zip(header, values)}
        for key in (
            "cx",
            "cy",
            "cz",
            "spacing",
            "r_peak",
            "F",
            "raw_F",
            "legacy_F",
            "common_physics_score",
            "common_physics_raw",
            "common_edt_A",
        ):
            if row.get(key) not in (None, ""):
                row[key] = float(row[key])
        for key in ("nx", "ny", "nz"):
            if row.get(key) not in (None, ""):
                row[key] = int(row[key])
        rows.append(row)
    return rows, metadata


def _cav_centers_by_configuration(
    prediction_rows: list[dict[str, Any]],
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, dict[str, str]]]:
    details: dict[tuple[str, str], dict[str, Any]] = {}
    metadata: dict[str, dict[str, str]] = {}
    paths = {
        (str(row["method"]), Path(str(row["source"])))
        for row in prediction_rows
        if str(row["method"]).startswith("cav-emps-") and row.get("source")
    }
    for method, path in sorted(paths, key=lambda item: (item[0], str(item[1]))):
        if not path.is_file():
            continue
        center_rows, center_meta = _parse_centers_tsv(path)
        metadata[method] = center_meta
        for row in center_rows:
            details[(method, str(row["site_id"]))] = row
    return details, metadata


def _parse_fld(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="ignore")

    def comment_values(name: str) -> list[float]:
        match = re.search(rf"^{re.escape(name)}\s+(.+)$", text, flags=re.MULTILINE | re.I)
        if not match:
            raise ValueError(f"{path}: missing {name}")
        return [float(value) for value in match.group(1).split()]

    def dimension(name: str) -> int:
        match = re.search(rf"^{name}\s*=\s*(\d+)", text, flags=re.MULTILINE | re.I)
        if not match:
            raise ValueError(f"{path}: missing {name}")
        return int(match.group(1))

    spacing = comment_values("#SPACING")[0]
    center = tuple(comment_values("#CENTER")[:3])
    nelements = tuple(int(value) for value in comment_values("#NELEMENTS")[:3])
    shape = tuple(dimension(f"dim{axis}") for axis in (1, 2, 3))
    expected_shape = tuple(value + 1 for value in nelements)
    if shape != expected_shape:
        raise ValueError(f"{path}: dimensions {shape} != NELEMENTS+1 {expected_shape}")
    origin = tuple(
        center[axis] - 0.5 * nelements[axis] * spacing for axis in range(3)
    )
    return {
        "spacing": spacing,
        "center": center,
        "nelements": nelements,
        "shape": shape,
        "origin": origin,
    }


def _nearest_grid_index(
    center: tuple[float, float, float],
    origin: tuple[float, float, float],
    spacing: float,
    shape: tuple[int, int, int],
) -> tuple[int, int, int]:
    return tuple(
        max(0, min(shape[axis] - 1, int(round((center[axis] - origin[axis]) / spacing))))
        for axis in range(3)
    )  # type: ignore[return-value]


def _flat_fortran_index(index: tuple[int, int, int], shape: tuple[int, int, int]) -> int:
    i, j, k = index
    return i + shape[0] * (j + shape[1] * k)


def _sample_ascii_map(
    path: Path,
    shape: tuple[int, int, int],
    wanted_indices: set[int],
) -> dict[int, float]:
    """Stream an AutoGrid map once and retain only requested Fortran-order voxels."""
    samples: dict[int, float] = {}
    value_index = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped[0] not in "0123456789-+.":
                continue
            for token in stripped.split():
                try:
                    value = float(token)
                except ValueError:
                    continue
                if value_index in wanted_indices:
                    samples[value_index] = value
                value_index += 1
    expected = math.prod(shape)
    if value_index != expected:
        raise ValueError(f"{path}: expected {expected} map values, found {value_index}")
    missing = wanted_indices - samples.keys()
    if missing:
        raise ValueError(f"{path}: failed to sample map indices {sorted(missing)}")
    return samples


def _combined_map_samples(
    output_dir: Path,
    prediction_rows: list[dict[str, Any]],
    cav_metadata: dict[str, dict[str, str]],
) -> tuple[dict[int, dict[str, Any]], str]:
    method = "cav-emps-combined-final"
    combined = sorted(
        (row for row in prediction_rows if row["method"] == method),
        key=lambda row: int(row["rank"]),
    )[:6]
    if not combined:
        return {}, "combined-final predictions are missing"
    map_root = output_dir / "predictions" / "cav-emps" / TARGET_ID / "maps"
    fld_candidates = sorted(map_root.rglob("*.fld"))
    if len(fld_candidates) != 1:
        return {}, f"expected one FLD under {map_root}, found {len(fld_candidates)}"
    fld = fld_candidates[0]
    fld_meta = _parse_fld(fld)
    center_meta = cav_metadata.get(method, {})
    filenames = {
        "C": center_meta.get("map_C_filename", ""),
        "e": center_meta.get("map_E_filename", ""),
        "d": center_meta.get("map_D_filename", ""),
    }
    if not all(filenames.values()):
        return {}, "combined-final centers metadata does not identify C/e/d map files"
    map_paths = {}
    for channel, filename in filenames.items():
        candidates = sorted(map_root.rglob(filename))
        if len(candidates) != 1:
            return {}, f"expected one {channel} map named {filename}, found {len(candidates)}"
        map_paths[channel] = candidates[0]

    rank_to_flat = {}
    rank_to_index = {}
    for row in combined:
        center = (float(row["center_x"]), float(row["center_y"]), float(row["center_z"]))
        index = _nearest_grid_index(
            center,
            fld_meta["origin"],
            fld_meta["spacing"],
            fld_meta["shape"],
        )
        rank = int(row["rank"])
        rank_to_index[rank] = index
        rank_to_flat[rank] = _flat_fortran_index(index, fld_meta["shape"])
    wanted = set(rank_to_flat.values())
    channel_samples = {
        channel: _sample_ascii_map(path, fld_meta["shape"], wanted)
        for channel, path in map_paths.items()
    }
    result = {}
    for rank, flat_index in rank_to_flat.items():
        result[rank] = {
            "map_grid_i": rank_to_index[rank][0],
            "map_grid_j": rank_to_index[rank][1],
            "map_grid_k": rank_to_index[rank][2],
            "map_C_raw": channel_samples["C"][flat_index],
            "map_e_raw": channel_samples["e"][flat_index],
            "map_d_raw": channel_samples["d"][flat_index],
            "map_C_file": str(map_paths["C"].resolve()),
            "map_e_file": str(map_paths["e"].resolve()),
            "map_d_file": str(map_paths["d"].resolve()),
        }
    return result, ""


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[axis] - b[axis]) ** 2 for axis in range(3)))


def _box_ligand_metrics(box: dict[str, Any], label: dict[str, Any]) -> dict[str, Any]:
    center = tuple(float(value) for value in box["center"])
    side = tuple(float(value) for value in box["side"])
    ligand_center = tuple(float(value) for value in label["center"])
    atoms = [tuple(float(value) for value in atom) for atom in label["ligand_atoms"]]
    axis_extent = [
        max(abs(atom[axis] - center[axis]) for atom in atoms) for axis in range(3)
    ]
    face_clearance = min(side[axis] / 2.0 - axis_extent[axis] for axis in range(3))
    centroid_inside = all(
        abs(ligand_center[axis] - center[axis]) <= side[axis] / 2.0
        for axis in range(3)
    )
    result = {
        "centroid_inside": int(centroid_inside),
        "face_clearance_A": face_clearance,
        "center_to_ligand_centroid_A": _distance(center, ligand_center),
        "dca_A": min(_distance(center, atom) for atom in atoms),
        "required_side_x_A": 2.0 * axis_extent[0],
        "required_side_y_A": 2.0 * axis_extent[1],
        "required_side_z_A": 2.0 * axis_extent[2],
    }
    for buffer_A in BUFFER_VALUES:
        key = f"full_containment_{int(buffer_A)}A"
        result[key] = int(face_clearance + 1e-9 >= buffer_A)
    return result


def _build_boxes(
    prediction_rows: list[dict[str, Any]],
    cav_details: dict[tuple[str, str], dict[str, Any]],
    standardized_side_A: float,
) -> list[dict[str, Any]]:
    boxes = []
    for row in prediction_rows:
        if int(row["rank"]) > 6:
            continue
        method = str(row["method"])
        center = (float(row["center_x"]), float(row["center_y"]), float(row["center_z"]))
        common = {
            "method": method,
            "rank": int(row["rank"]),
            "site_id": str(row["site_id"]),
            "center": center,
            "score": row.get("score"),
            "source": row.get("source", ""),
        }
        boxes.append(
            {
                **common,
                "box_variant": "standardized_35.625A",
                "side": (standardized_side_A,) * 3,
            }
        )
        if method.startswith("cav-emps-"):
            detail = cav_details.get((method, str(row["site_id"])))
            if detail:
                spacing = float(detail["spacing"])
                boxes.append(
                    {
                        **common,
                        "box_variant": "native",
                        "side": (
                            float(detail["nx"]) * spacing,
                            float(detail["ny"]) * spacing,
                            float(detail["nz"]) * spacing,
                        ),
                    }
                )
    return boxes


def _summary_for_boxes(
    method: str,
    variant: str,
    boxes: list[dict[str, Any]],
    labels: list[dict[str, Any]],
    missing_reason: str = "",
) -> dict[str, Any]:
    boxes = sorted(boxes, key=lambda box: int(box["rank"]))
    centroid_covered = set()
    contained = {buffer_A: set() for buffer_A in BUFFER_VALUES}
    min_dca = {label["site_id"]: math.inf for label in labels}
    redundant = 0
    off_target = 0
    already_covered = set()
    for box in boxes:
        box_sites = set()
        for label in labels:
            metrics = _box_ligand_metrics(box, label)
            site_id = label["site_id"]
            if metrics["centroid_inside"]:
                centroid_covered.add(site_id)
            if metrics["full_containment_0A"]:
                box_sites.add(site_id)
            for buffer_A in BUFFER_VALUES:
                if metrics[f"full_containment_{int(buffer_A)}A"]:
                    contained[buffer_A].add(site_id)
            min_dca[site_id] = min(min_dca[site_id], float(metrics["dca_A"]))
        if not box_sites:
            off_target += 1
        elif not (box_sites - already_covered):
            redundant += 1
        already_covered.update(box_sites)
    finite_dca = [value for value in min_dca.values() if math.isfinite(value)]
    missed = [
        label["pocket_class"]
        for label in labels
        if label["site_id"] not in contained[0.0]
    ]
    return {
        "method": method,
        "box_variant": variant,
        "status": "ok" if boxes else "missing",
        "n_predictions": len(boxes),
        "n_truth_sites": len(labels),
        "centroid_sites_covered": len(centroid_covered),
        "full_ligands_contained_0A": len(contained[0.0]),
        "full_ligands_contained_2A": len(contained[2.0]),
        "full_ligands_contained_4A": len(contained[4.0]),
        "redundant_boxes": redundant,
        "off_target_boxes": off_target,
        "mean_min_dca_A": (
            sum(finite_dca) / len(finite_dca) if finite_dca else ""
        ),
        "max_min_dca_A": max(finite_dca) if finite_dca else "",
        "missed_pocket_classes": "; ".join(sorted(set(missed))),
        "missing_information": missing_reason if not boxes else "",
    }


def _first_containing_rank(
    containment_rows: list[dict[str, Any]],
    method: str,
    variant: str,
    truth_site_id: str,
    buffer_A: float = 0.0,
) -> int | None:
    key = f"full_containment_{int(buffer_A)}A"
    ranks = [
        int(row["rank"])
        for row in containment_rows
        if row["method"] == method
        and row["box_variant"] == variant
        and row["truth_site_id"] == truth_site_id
        and int(row[key]) == 1
    ]
    return min(ranks) if ranks else None


def _closest_row(
    containment_rows: list[dict[str, Any]],
    method: str,
    variant: str,
    truth_site_id: str,
) -> dict[str, Any] | None:
    matches = [
        row
        for row in containment_rows
        if row["method"] == method
        and row["box_variant"] == variant
        and row["truth_site_id"] == truth_site_id
    ]
    return min(matches, key=lambda row: (float(row["dca_A"]), int(row["rank"]))) if matches else None


def _format_float(value: Any, digits: int = 3) -> Any:
    if value in (None, ""):
        return ""
    return f"{float(value):.{digits}f}"


def evaluate_case_study(
    output_dir: Path,
    methods: list[str],
    preregistration: dict[str, Any],
) -> dict[str, Path]:
    labels_payload = json.loads((output_dir / "ground_truth" / "labels.json").read_text())
    labels = labels_payload["binding_sites"]
    prediction_rows, missing = _collect_prediction_rows(output_dir, methods, preregistration)
    cav_details, cav_metadata = _cav_centers_by_configuration(prediction_rows)
    standardized_side_A = float(
        preregistration["prediction_protocol"]["standardized_box_side_A"]
    )
    boxes = _build_boxes(prediction_rows, cav_details, standardized_side_A)
    evaluation_dir = output_dir / "evaluation"
    evaluation_dir.mkdir(parents=True, exist_ok=True)

    prediction_fields = [
        "target_id",
        "method",
        "rank",
        "site_id",
        "center_x",
        "center_y",
        "center_z",
        "score",
        "source",
        "prediction_file",
    ]
    _write_tsv(evaluation_dir / "predictions_all.tsv", prediction_fields, prediction_rows)

    truth_rows = []
    for label in labels:
        truth_rows.append(
            {
                "site_id": label["site_id"],
                "display_name": label["display_name"],
                "ligand_id": label["ligand_id"],
                "ligand": label["ligand"],
                "pocket_class": label["pocket_class"],
                "interface": label["interface"],
                "center_x": label["center"][0],
                "center_y": label["center"][1],
                "center_z": label["center"][2],
                "atom_count": label["atom_count"],
                "preregistered_map_expectation": label["preregistered_map_expectation"],
                "source": label["source"],
            }
        )
    truth_fields = list(truth_rows[0])
    _write_csv(evaluation_dir / "ground_truth_sites.csv", truth_fields, truth_rows)

    containment_rows = []
    for box in boxes:
        for label in labels:
            metrics = _box_ligand_metrics(box, label)
            containment_rows.append(
                {
                    "target_id": TARGET_ID,
                    "method": box["method"],
                    "box_variant": box["box_variant"],
                    "rank": box["rank"],
                    "site_id": box["site_id"],
                    "score": box["score"],
                    "center_x": box["center"][0],
                    "center_y": box["center"][1],
                    "center_z": box["center"][2],
                    "box_side_x_A": box["side"][0],
                    "box_side_y_A": box["side"][1],
                    "box_side_z_A": box["side"][2],
                    "truth_site_id": label["site_id"],
                    "truth_display_name": label["display_name"],
                    "ligand_id": label["ligand_id"],
                    "pocket_class": label["pocket_class"],
                    "interface": label["interface"],
                    **metrics,
                    "prediction_source": box["source"],
                }
            )
    containment_fields = list(containment_rows[0]) if containment_rows else [
        "target_id",
        "method",
        "box_variant",
        "rank",
        "site_id",
        "truth_site_id",
        "missing_information",
    ]
    _write_csv(
        evaluation_dir / "per_box_containment.csv",
        containment_fields,
        containment_rows,
    )

    expected_configs = _expected_configurations(methods, preregistration)
    summaries = []
    for method in expected_configs:
        variants = ["standardized_35.625A"]
        if method.startswith("cav-emps-"):
            variants.insert(0, "native")
        for variant in variants:
            selected_boxes = [
                box for box in boxes if box["method"] == method and box["box_variant"] == variant
            ]
            summaries.append(
                _summary_for_boxes(
                    method,
                    variant,
                    selected_boxes,
                    labels,
                    missing.get(method, ""),
                )
            )
    summary_fields = list(summaries[0])
    _write_csv(evaluation_dir / "portfolio_summary.csv", summary_fields, summaries)
    ablation_rows = [
        row
        for row in summaries
        if row["method"].startswith("cav-emps-")
        and row["box_variant"] == "standardized_35.625A"
    ]
    _write_csv(evaluation_dir / "ablation_summary.csv", summary_fields, ablation_rows)

    map_samples, map_sample_error = _combined_map_samples(
        output_dir,
        prediction_rows,
        cav_metadata,
    )
    combined_method = "cav-emps-combined-final"
    per_site_rows = []
    for label in labels:
        truth_site_id = label["site_id"]
        cav_rank = _first_containing_rank(
            containment_rows,
            combined_method,
            "native",
            truth_site_id,
        )
        cav_match = None
        if cav_rank is not None:
            cav_match = next(
                row
                for row in containment_rows
                if row["method"] == combined_method
                and row["box_variant"] == "native"
                and row["truth_site_id"] == truth_site_id
                and int(row["rank"]) == cav_rank
            )
        else:
            cav_match = _closest_row(
                containment_rows,
                combined_method,
                "native",
                truth_site_id,
            )
        diagnostic_rank = int(cav_match["rank"]) if cav_match else None
        cav_site_id = str(cav_match["site_id"]) if cav_match else ""
        detail = cav_details.get((combined_method, cav_site_id), {})
        sample = map_samples.get(diagnostic_rank or -1, {})
        per_site_rows.append(
            {
                "truth_site_id": truth_site_id,
                "site": label["display_name"],
                "ligand_id": label["ligand_id"],
                "pocket_class": label["pocket_class"],
                "interface": label["interface"],
                "preregistered_map_expectation": label["preregistered_map_expectation"],
                "cav_rank_first_native_full_containment": cav_rank or "",
                "cav_diagnostic_rank": diagnostic_rank or "",
                "cav_site_id": cav_site_id,
                "cav_native_containment_0A": cav_match["full_containment_0A"] if cav_match else "",
                "cav_native_containment_2A": cav_match["full_containment_2A"] if cav_match else "",
                "cav_native_containment_4A": cav_match["full_containment_4A"] if cav_match else "",
                "cav_standardized_rank_0A": _first_containing_rank(
                    containment_rows,
                    combined_method,
                    "standardized_35.625A",
                    truth_site_id,
                )
                or "",
                "cav_dca_A": _format_float(cav_match["dca_A"] if cav_match else ""),
                "map_C_raw_at_cav_center": _format_float(sample.get("map_C_raw", ""), 6),
                "map_e_raw_at_cav_center": _format_float(sample.get("map_e_raw", ""), 6),
                "map_d_raw_at_cav_center": _format_float(sample.get("map_d_raw", ""), 6),
                "EDT_A_at_cav_center": _format_float(detail.get("common_edt_A", ""), 4),
                "cav_family": detail.get("family", ""),
                "cav_portfolio_role": detail.get("portfolio_role", ""),
                "fpocket_rank_standardized_0A": _first_containing_rank(
                    containment_rows,
                    "fpocket",
                    "standardized_35.625A",
                    truth_site_id,
                )
                or "",
                "p2rank_rank_standardized_0A": _first_containing_rank(
                    containment_rows,
                    "p2rank",
                    "standardized_35.625A",
                    truth_site_id,
                )
                or "",
                "missing_information": map_sample_error,
            }
        )
    per_site_fields = list(per_site_rows[0])
    _write_csv(
        evaluation_dir / "per_site_interpretation.csv",
        per_site_fields,
        per_site_rows,
    )

    outputs = {
        "ground_truth": evaluation_dir / "ground_truth_sites.csv",
        "predictions": evaluation_dir / "predictions_all.tsv",
        "containment": evaluation_dir / "per_box_containment.csv",
        "portfolio": evaluation_dir / "portfolio_summary.csv",
        "ablation": evaluation_dir / "ablation_summary.csv",
        "per_site": evaluation_dir / "per_site_interpretation.csv",
    }
    report_outputs = generate_report(
        output_dir,
        labels,
        summaries,
        ablation_rows,
        per_site_rows,
        containment_rows,
        preregistration,
        missing,
    )
    outputs.update(report_outputs)
    _write_visualization(output_dir, prediction_rows, cav_details, labels)
    return outputs


def _markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(title for _, title in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        values = []
        for key, _ in columns:
            value = row.get(key, "")
            if isinstance(value, float):
                value = f"{value:.3f}"
            values.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _html_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    headers = "".join(f"<th>{html.escape(title)}</th>" for _, title in columns)
    body = []
    for row in rows:
        cells = []
        for key, _ in columns:
            value = row.get(key, "")
            if isinstance(value, float):
                value = f"{value:.3f}"
            cells.append(f"<td>{html.escape(str(value))}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def _main_summary_rows(summaries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    wanted = [
        ("cav-emps-combined-final", "native", "CaV-EMPS native boxes"),
        (
            "cav-emps-combined-final",
            "standardized_35.625A",
            "CaV-EMPS standardized",
        ),
        ("fpocket", "standardized_35.625A", "fpocket standardized"),
        ("p2rank", "standardized_35.625A", "P2Rank standardized"),
    ]
    rows = []
    for method, variant, display in wanted:
        match = next(
            (
                row
                for row in summaries
                if row["method"] == method and row["box_variant"] == variant
            ),
            None,
        )
        if match is None:
            match = {
                "method": method,
                "box_variant": variant,
                "status": "not selected",
                "n_predictions": 0,
                "n_truth_sites": 5,
                "centroid_sites_covered": 0,
                "full_ligands_contained_0A": 0,
                "full_ligands_contained_2A": 0,
                "full_ligands_contained_4A": 0,
                "redundant_boxes": 0,
                "off_target_boxes": 0,
                "missed_pocket_classes": "",
                "missing_information": "method was not selected",
            }
        rows.append({**match, "display_method": display})
    return rows


def _coverage_level(
    containment_rows: list[dict[str, Any]],
    method: str,
    variant: str,
    truth_site_id: str,
) -> tuple[str, str]:
    rows = [
        row
        for row in containment_rows
        if row["method"] == method
        and row["box_variant"] == variant
        and row["truth_site_id"] == truth_site_id
    ]
    if any(int(row["full_containment_4A"]) for row in rows):
        return "#1b7f5a", "4 A buffer"
    if any(int(row["full_containment_2A"]) for row in rows):
        return "#55a868", "2 A buffer"
    if any(int(row["full_containment_0A"]) for row in rows):
        return "#e5b642", "0 A containment"
    if any(int(row["centroid_inside"]) for row in rows):
        return "#e99138", "centroid only"
    return "#ba4a4a", "miss"


def _write_containment_svg(
    path: Path,
    labels: list[dict[str, Any]],
    containment_rows: list[dict[str, Any]],
) -> None:
    methods = [
        ("cav-emps-combined-final", "native", "CaV native"),
        ("cav-emps-combined-final", "standardized_35.625A", "CaV 35.625 A"),
        ("fpocket", "standardized_35.625A", "fpocket 35.625 A"),
        ("p2rank", "standardized_35.625A", "P2Rank 35.625 A"),
    ]
    left = 180
    top = 66
    cell_w = 130
    cell_h = 46
    width = left + cell_w * len(labels) + 20
    height = top + cell_h * len(methods) + 70
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#171717}.h{font-size:13px;font-weight:700}.r{font-size:13px}.c{font-size:12px;fill:#fff;font-weight:700}</style>',
        '<text x="16" y="24" class="h">Top-six full-ligand containment</text>',
    ]
    for column, label in enumerate(labels):
        x = left + column * cell_w + cell_w / 2
        short = label["display_name"]
        parts.append(
            f'<text x="{x}" y="52" text-anchor="middle" class="h">{html.escape(short)}</text>'
        )
    for row_index, (method, variant, display) in enumerate(methods):
        y = top + row_index * cell_h
        parts.append(f'<text x="16" y="{y + 29}" class="r">{html.escape(display)}</text>')
        for column, label in enumerate(labels):
            x = left + column * cell_w
            color, level = _coverage_level(
                containment_rows,
                method,
                variant,
                label["site_id"],
            )
            parts.append(
                f'<rect x="{x + 3}" y="{y + 3}" width="{cell_w - 6}" height="{cell_h - 6}" rx="3" fill="{color}"/>'
            )
            parts.append(
                f'<text x="{x + cell_w / 2}" y="{y + 29}" text-anchor="middle" class="c">{html.escape(level)}</text>'
            )
    parts.append(
        f'<text x="16" y="{height - 20}" class="r">Green: buffered containment; yellow: complete at box face; orange: centroid only; red: miss.</text>'
    )
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def generate_report(
    output_dir: Path,
    labels: list[dict[str, Any]],
    summaries: list[dict[str, Any]],
    ablation_rows: list[dict[str, Any]],
    per_site_rows: list[dict[str, Any]],
    containment_rows: list[dict[str, Any]],
    preregistration: dict[str, Any],
    missing: dict[str, str],
) -> dict[str, Path]:
    report_dir = output_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    svg_path = report_dir / "containment_matrix.svg"
    _write_containment_svg(svg_path, labels, containment_rows)
    main_rows = _main_summary_rows(summaries)
    main_columns = [
        ("display_method", "Method"),
        ("status", "Status"),
        ("centroid_sites_covered", "Occupied sites covered /5"),
        ("full_ligands_contained_0A", "Full ligands /5"),
        ("full_ligands_contained_2A", "2 A buffered /5"),
        ("full_ligands_contained_4A", "4 A buffered /5"),
        ("redundant_boxes", "Redundant boxes"),
        ("off_target_boxes", "Off-target boxes"),
        ("missed_pocket_classes", "Missed pocket class"),
    ]
    site_columns = [
        ("site", "Site"),
        ("cav_rank_first_native_full_containment", "CaV rank"),
        ("cav_native_containment_0A", "Native 0 A"),
        ("cav_native_containment_2A", "Native 2 A"),
        ("cav_native_containment_4A", "Native 4 A"),
        ("map_C_raw_at_cav_center", "C raw"),
        ("map_e_raw_at_cav_center", "e raw"),
        ("map_d_raw_at_cav_center", "d raw"),
        ("EDT_A_at_cav_center", "EDT A"),
        ("fpocket_rank_standardized_0A", "fpocket rank"),
        ("p2rank_rank_standardized_0A", "P2Rank rank"),
    ]
    ablation_columns = [
        ("method", "Configuration"),
        ("status", "Status"),
        ("centroid_sites_covered", "Centroids /5"),
        ("full_ligands_contained_0A", "Full 0 A /5"),
        ("full_ligands_contained_2A", "Full 2 A /5"),
        ("full_ligands_contained_4A", "Full 4 A /5"),
        ("redundant_boxes", "Redundant"),
        ("off_target_boxes", "Off-target"),
        ("mean_min_dca_A", "Mean min DCA A"),
    ]
    truth_columns = [
        ("display_name", "Site"),
        ("ligand_id", "Deposited ligand"),
        ("pocket_class", "Class"),
        ("interface", "Interface"),
        ("preregistered_map_expectation", "Preregistered CaV-map hypothesis"),
    ]
    missing_text = (
        "\n".join(f"- `{key}`: {value}" for key, value in sorted(missing.items()))
        if missing
        else "- None."
    )
    markdown = f"""# GABAA 8DD2 Blind Site-Recovery Case Study

Generated: `{_utc_now()}`

## Frozen Protocol

The predictor input contains only human receptor chains A-E from PDB 8DD2.
All deposited ligands, glycans, Fabs and other HETATM records were removed.
The three zolpidem and two GABA molecules were kept outside the normalized
predictor directory and revealed only during this evaluation.

The frozen portfolio contains six predictions per method. The primary endpoint
is complete containment of every deposited ligand atom in an axis-aligned box,
with 0, 2 and 4 Angstrom clearance from every box face. Center-to-ligand-atom
distance is secondary. Standardized comparisons use identical 35.625 Angstrom
boxes. No weights or ranks are tuned from this result.

## Ground Truth and Preregistered Hypotheses

{_markdown_table(labels, truth_columns)}

These C/e/d expectations are preregistered physical hypotheses, not claims
established by the cited structural papers. The controlled ablations are the
test of channel contribution.

## Docking-Box Portfolio Performance

{_markdown_table(main_rows, main_columns)}

![Containment matrix](containment_matrix.svg)

`Occupied sites covered` means ligand-centroid containment. `Full ligands`
means every ligand atom is inside the box. A redundant box contains only sites
already contained by an earlier-ranked box; an off-target box fully contains
none of the five labels.

## Per-Site Interpretation

{_markdown_table(per_site_rows, site_columns)}

The C/e/d columns are raw AutoGrid values at the diagnostic combined-final CaV
center. They are descriptive and are not equivalent to the regional,
transformed field used internally by CaV-EMPS. If no native box fully contains
a ligand, the diagnostic CaV center is the top-six center with the smallest
distance to any ligand atom.

## Controlled CaV-EMPS Ablations

{_markdown_table(ablation_rows, ablation_columns)}

Every profile used the same per-target C/e/d maps. The checked-in suite is
combined-final plus six one-component controls; geometry still acts as a
steric validity mask in the physics-only profile.

## Missing Information

{missing_text}

Missing or empty method output remains explicit and is not silently excluded.

## Interpretation Limits

- CaV-EMPS receives no explicit lipid bilayer. Partly lipid-facing TMD sites
  are therefore a difficult stress test and an acknowledged model limitation.
- The additional gamma2/beta2 diazepam TMD site reported in earlier work may
  help interpret an extra sixth prediction, but it is not one of this
  zolpidem/GABA truth set's five labels.
- This focused mechanistic case study complements, but does not replace, the
  Holo4K containment audit.

## Audit Files

- `../evaluation/ground_truth_sites.csv`
- `../evaluation/predictions_all.tsv`
- `../evaluation/per_box_containment.csv`
- `../evaluation/per_site_interpretation.csv`
- `../evaluation/portfolio_summary.csv`
- `../evaluation/ablation_summary.csv`
- `../provenance.json`

## References

- [RCSB PDB 8DD2](https://www.rcsb.org/structure/8DD2), DOI
  `10.2210/pdb8DD2/pdb`.
- Masiulis et al. (2019), DOI `10.1038/s41586-018-0832-5`.
- Kim et al. (2020), DOI `10.1038/s41586-020-2654-5`.
- Zhu et al. (2022), DOI `10.1038/s41467-022-32212-4`.
- Legesse et al. (2023), DOI `10.1038/s41467-023-40800-1`.
"""
    markdown_path = report_dir / "case_study_report.md"
    markdown_path.write_text(markdown, encoding="utf-8")

    html_path = report_dir / "case_study_report.html"
    html_document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GABAA 8DD2 Blind Site-Recovery Case Study</title>
<style>
body{{font:15px/1.5 system-ui,sans-serif;color:#1c2328;max-width:1280px;margin:32px auto;padding:0 24px}}
h1,h2{{line-height:1.2}} table{{border-collapse:collapse;width:100%;margin:12px 0 28px}}
th,td{{border:1px solid #ccd3d8;padding:7px 9px;text-align:left;vertical-align:top}}
th{{background:#eef2f3}} code{{background:#f2f4f5;padding:2px 4px}} .note{{border-left:4px solid #d09c2d;padding:10px 14px;background:#fff8e8}}
img{{max-width:100%;height:auto}}
</style></head><body>
<h1>GABAA 8DD2 Blind Site-Recovery Case Study</h1>
<p>Generated <code>{html.escape(_utc_now())}</code></p>
<h2>Frozen protocol</h2>
<p>The predictor saw only receptor chains A-E. Three zolpidem and two GABA molecules, all other HETATM records, glycans, and Fab chains were withheld. The primary endpoint is full-ligand containment in the top-six boxes at 0, 2, and 4 Angstrom face buffers. Standardized comparisons use identical 35.625 Angstrom boxes.</p>
<h2>Ground truth and preregistered hypotheses</h2>{_html_table(labels, truth_columns)}
<p class="note">The C/e/d expectations are preregistered physical hypotheses. Literature establishes the occupied sites; controlled ablations test map contribution.</p>
<h2>Docking-box portfolio performance</h2>{_html_table(main_rows, main_columns)}
<img src="containment_matrix.svg" alt="Full-ligand containment matrix">
<h2>Per-site interpretation</h2>{_html_table(per_site_rows, site_columns)}
<p>Raw C/e/d values are sampled at the diagnostic CaV center and are not the regional transformed score itself.</p>
<h2>Controlled CaV-EMPS ablations</h2>{_html_table(ablation_rows, ablation_columns)}
<h2>Missing information</h2><pre>{html.escape(missing_text)}</pre>
<h2>Limitations</h2><p>There is no explicit lipid bilayer. Partly lipid-facing TMD sites are a deliberate stress test. The additional diazepam gamma2/beta2 TMD site may contextualize a sixth prediction but is not a ground-truth zolpidem site. This case study does not replace the Holo4K containment audit.</p>
<h2>References</h2><ul><li><a href="https://www.rcsb.org/structure/8DD2">RCSB PDB 8DD2</a></li><li>Masiulis et al. 2019, DOI 10.1038/s41586-018-0832-5</li><li>Kim et al. 2020, DOI 10.1038/s41586-020-2654-5</li><li>Zhu et al. 2022, DOI 10.1038/s41467-022-32212-4</li><li>Legesse et al. 2023, DOI 10.1038/s41467-023-40800-1</li></ul>
</body></html>"""
    html_path.write_text(html_document, encoding="utf-8")
    return {"report_markdown": markdown_path, "report_html": html_path, "matrix_svg": svg_path}


def _pml_quote(path: Path) -> str:
    return str(path.resolve()).replace("\\", "/").replace('"', '\\"')


def _write_visualization(
    output_dir: Path,
    prediction_rows: list[dict[str, Any]],
    cav_details: dict[tuple[str, str], dict[str, Any]],
    labels: list[dict[str, Any]],
) -> Path:
    visualization_dir = output_dir / "visualization"
    visualization_dir.mkdir(parents=True, exist_ok=True)
    path = visualization_dir / "gabaa_8dd2_case_study.pml"
    receptor = output_dir / "normalized" / DATASET_ID / TARGET_ID / "receptor_input.pdb"
    lines = [
        f'load "{_pml_quote(receptor)}", receptor_8dd2',
        "hide everything",
        "show cartoon, receptor_8dd2",
        "color gray80, receptor_8dd2",
    ]
    for label in labels:
        ligand_path = Path(label["extracted_ligand"])
        object_name = f"truth_{label['site_id']}"
        lines.extend(
            [
                f'load "{_pml_quote(ligand_path)}", {object_name}',
                f"show sticks, {object_name}",
                f"color yellow, {object_name}",
            ]
        )
    method_colors = {
        "cav-emps-combined-final": "cyan",
        "fpocket": "green",
        "p2rank": "magenta",
    }
    for row in prediction_rows:
        method = str(row["method"])
        if method not in method_colors or int(row["rank"]) > 6:
            continue
        safe_method = method.replace("-", "_")
        name = f"{safe_method}_S{int(row['rank'])}"
        lines.append(
            f"pseudoatom {name}, pos=[{row['center_x']},{row['center_y']},{row['center_z']}], label=\"{method} #{row['rank']}\""
        )
        lines.append(f"color {method_colors[method]}, {name}")
        lines.append(f"show spheres, {name}")
    lines.extend(
        [
            "set sphere_scale, 0.65",
            "set label_size, 14",
            "group ground_truth, truth_*",
            "group cav_emps_centers, cav_emps_combined_final_*",
            "group fpocket_centers, fpocket_*",
            "group p2rank_centers, p2rank_*",
            "disable fpocket_centers",
            "disable p2rank_centers",
            "orient receptor_8dd2",
            "zoom receptor_8dd2, 8",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _git_value(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return completed.stdout.strip()


def _write_provenance(
    output_dir: Path,
    source_pdb: Path,
    receptor_path: Path,
    labels_path: Path,
    args: argparse.Namespace,
    preregistration: dict[str, Any],
) -> Path:
    path = output_dir / "provenance.json"
    files = {
        "source_pdb": source_pdb,
        "receptor_input": receptor_path,
        "ground_truth_labels": labels_path,
        "preregistration": PREREGISTRATION_PATH,
        "predict_sites": PREDICT_SCRIPT,
        "make_grids": REPO_ROOT / "docking" / "make_grids.py",
    }
    payload = {
        "case_study_id": preregistration["case_study_id"],
        "generated_at_utc": _utc_now(),
        "python": sys.version,
        "repo_root": str(REPO_ROOT),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "git_status_short": _git_value("status", "--short").splitlines(),
        "methods_requested": args.methods,
        "portfolio_size": preregistration["prediction_protocol"]["portfolio_size"],
        "standardized_box_side_A": preregistration["prediction_protocol"][
            "standardized_box_side_A"
        ],
        "files": {
            name: {"path": str(file.resolve()), "sha256": _sha256(file)}
            for name, file in files.items()
            if file.is_file()
        },
        "commands_log": str((output_dir / "commands.jsonl").resolve()),
        "no_post_hoc_tuning": True,
    }
    _write_json(path, payload)
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the preregistered PDB 8DD2 GABAA blind site-recovery case study.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=EXAMPLE_DIR / "workspace" / "gabaa_8dd2_case_study",
    )
    parser.add_argument("--pdb", type=Path, help="Local archived 8DD2 PDB; otherwise download RCSB.")
    parser.add_argument(
        "--allow-source-mismatch",
        action="store_true",
        help="Allow a PDB whose SHA-256 differs from the preregistered RCSB revision.",
    )
    parser.add_argument("--methods", default="cav-emps,fpocket,p2rank")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--autogrid4-bin", default=str(DEFAULT_AUTOGRID4))
    parser.add_argument(
        "--fpocket-cmd",
        default=str(DEFAULT_FPOCKET) if DEFAULT_FPOCKET.exists() else "fpocket",
    )
    parser.add_argument(
        "--p2rank-cmd",
        default=(f"{DEFAULT_P2RANK} predict" if DEFAULT_P2RANK.exists() else "prank predict"),
    )
    parser.add_argument("--receptor-prepare-command")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.prepare_only and args.evaluate_only:
        raise SystemExit("--prepare-only and --evaluate-only are mutually exclusive")
    try:
        args.methods = _method_list(args.methods)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    output_dir = args.output_dir.expanduser().resolve()
    preregistration = _load_preregistration()

    if args.dry_run:
        print(f"Output: {output_dir}")
        print("Structure: PDB 8DD2, receptor chains A-E only")
        print("Hidden truth: ABU A403, ABU C403, R5R D601, R5R A402, R5R C402")
        print(f"Methods: {', '.join(args.methods)}")
        print("CaV-EMPS profiles: " + ", ".join(preregistration["prediction_protocol"]["cav_emps_profiles"]))
        print("Evaluation: native CaV boxes plus standardized 35.625 A boxes; full containment at 0/2/4 A")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    if args.evaluate_only:
        source_pdb = output_dir / "source" / "8DD2.pdb"
        receptor_path = output_dir / "normalized" / DATASET_ID / TARGET_ID / "receptor_input.pdb"
        labels_path = output_dir / "ground_truth" / "labels.json"
        for required in (source_pdb, receptor_path, labels_path):
            if not required.is_file():
                raise FileNotFoundError(f"--evaluate-only requires {required}")
    else:
        source_pdb = _stage_source_pdb(
            output_dir,
            args.pdb,
            force=bool(args.force),
            allow_source_mismatch=bool(args.allow_source_mismatch),
            preregistration=preregistration,
        )
        receptor_path, labels_path = prepare_case_study(
            output_dir,
            source_pdb,
            preregistration,
        )

    _write_provenance(
        output_dir,
        source_pdb,
        receptor_path,
        labels_path,
        args,
        preregistration,
    )
    if args.prepare_only:
        print(f"Prepared receptor: {receptor_path}")
        print(f"Hidden labels: {labels_path}")
        return

    if not args.evaluate_only:
        if args.force:
            (output_dir / "commands.jsonl").unlink(missing_ok=True)
        prediction_failures = run_predictions(output_dir, args, preregistration)
    else:
        prediction_failures = {}
    outputs = evaluate_case_study(output_dir, args.methods, preregistration)
    _write_provenance(
        output_dir,
        source_pdb,
        receptor_path,
        labels_path,
        args,
        preregistration,
    )
    print(f"\nCase-study report: {outputs['report_markdown']}")
    print(f"HTML report: {outputs['report_html']}")
    print(f"Portfolio summary: {outputs['portfolio']}")
    if prediction_failures:
        print(
            "Partial benchmark: failed method(s): "
            + ", ".join(sorted(prediction_failures)),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
