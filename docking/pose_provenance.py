"""Record the identities of a docking output and its original inputs."""

import hashlib
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def record_docking_run(directory, outputs, failures=None):
    """Persist successful outputs and every skipped screening case, even if all fail."""
    failures = list(failures or [])
    Path(directory).mkdir(parents=True, exist_ok=True)
    target = Path(directory) / ("docking-run-" +
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f") + ".json")
    records = [{"path": str(Path(path).resolve()),
                "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
               for path in sorted(set(outputs))]
    for record in records:
        sidecar = Path(record["path"]).with_suffix(".md.json")
        if sidecar.exists():
            data = json.loads(sidecar.read_text())
            record.update({key: data[key] for key in
                           ("receptor_id", "ligand_id", "binding_site", "engine")})
    status = "failed" if not records else "partial" if failures else "complete"
    failure_report = target.with_suffix(".failures.csv")
    fields = ["stage", "receptor_id", "ligand_id", "binding_site", "error"]
    with failure_report.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(failures)
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps({"schema": 2, "status": status, "outputs": records,
                                    "failures": failures, "failure_report": str(failure_report.resolve()),
                                    "successful_outputs": len(records),
                                    "failed_cases": len(failures)}, indent=2) + "\n")
    temporary.replace(target)
    print(f"[screening] {len(records)} successful docking outputs; {len(failures)} failures. "
          f"Status: {status}. Failure report: {failure_report}", flush=True)
    return target


def record_pose_provenance(output, ligand, receptor, receptor_id, site, engine):
    output, ligand, receptor = map(lambda value: Path(value).resolve(), (output, ligand, receptor))
    def checksum(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    target = output.with_suffix(".md.json")
    data = {
        "schema": 1,
        "engine": engine,
        "receptor_id": receptor_id,
        "binding_site": str(site).removeprefix("S"),
        "ligand_id": ligand.stem,
        "output": str(output),
        "output_sha256": checksum(output),
        "source_ligand": str(ligand),
        "source_ligand_sha256": checksum(ligand),
        "receptor": str(receptor),
        "receptor_sha256": checksum(receptor),
    }
    if engine == "adgpu":
        # The XML holds scores; the matching DLG holds the coordinates used by MD.
        coordinates = output.with_suffix(".dlg")
        if not coordinates.is_file():
            raise ValueError(f"Exact AutoDock-GPU coordinates require {coordinates}")
        data["coordinates_sha256"] = checksum(coordinates)
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(target)
