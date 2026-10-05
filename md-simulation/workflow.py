#!/usr/bin/env python3
"""Ultidock's CLI-only, checkpointed docking-to-GROMACS workflow."""

from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess
import sys

import numpy as np

from chemistry import convert_mol2, posed_molecule, transform_pdb, validate_receptor
from diagnostics import check_dependencies
from poses import digest, select_poses
from protocol import ligand_input, load_protocol, stage_mdp
from systems import (add_topology, embed_bilayer, read_gro, remove_membrane_core_water,
                     restore_ligand, set_water_count, split_itp, topology_charge, write_gro)

ROOT = Path(__file__).resolve().parent


def save(path: Path, data) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def work_path(value=None) -> Path:
    path = (Path(value).resolve() if value else ROOT / "workspace" /
            datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f"))
    if not path.is_relative_to(ROOT) or path == ROOT:
        raise ValueError(f"MD working directories must be under {ROOT}")
    return path


def snapshot(source: Path, target: Path) -> None:
    if source.is_dir():
        if any(p.is_symlink() for p in source.rglob("*")):
            raise ValueError(f"Parameter directory must not contain symlinks: {source}")
        shutil.copytree(source, target)
    else:
        shutil.copy2(source, target)


def preflight(config_path: Path, *, receptor_dir=None, ligands_dir=None, engine=None,
              tools=True, gmx="gmx", acpype="acpype", obabel="obabel") -> dict:
    """Validate the protocol before spending time on docking or creating an MD job."""
    config = load_protocol(config_path)
    validate_receptor(Path(config["receptor_pdb"]), Path(config["docking_receptor_pdbqt"]))
    if engine is not None and config["engine"] != engine:
        raise ValueError(f"MD protocol engine={config['engine']} but docking uses {engine}; "
                         "use mode cpu for vina or a GPU backend for adgpu")
    if receptor_dir is not None:
        receptor = Path(receptor_dir) / f"{config['receptor_id']}.pdbqt"
        if not receptor.is_file():
            raise ValueError(f"MD receptor_id does not match a docking receptor: {receptor}")
        validate_receptor(Path(config["receptor_pdb"]), receptor)
        config["docking_receptor_pdbqt"] = str(receptor.resolve())
    if ligands_dir is not None:
        ligands = sorted(Path(ligands_dir).glob("*.pdbqt"))
        if not ligands:
            raise ValueError(f"No prepared docking ligands in {ligands_dir}")
        for ligand in ligands:
            lines = [line for line in ligand.read_text().splitlines()
                     if line.startswith(("ATOM  ", "HETATM"))]
            # Check chemical identity and atom mapping here; MD placement still comes
            # exclusively from the scored output selected after docking finishes.
            try:
                entry = ligand_input(config, ligand.stem)
                posed_molecule(Path(entry["sdf"]), entry["atom_map"], lines,
                               entry["net_charge"], config["transform"])
            except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
                print(f"[WARN] MD chemical input for {ligand.stem}: {exc}; "
                      "this candidate will be skipped if selected", flush=True)
    for stage in ("ions", "em", "nvt", "npt", "production"):
        stage_mdp(config, stage)
    if tools and not check_dependencies(gmx, acpype, obabel, workspace=ROOT / "workspace"):
        raise ValueError("Install the missing MD dependencies; run ultidock doctor")
    return config


def prepare(config_path: Path, docking_dir: Path, output=None, legacy=False, manifest=None,
            *, config=None) -> Path:
    config = load_protocol(config_path) if config is None else deepcopy(config)
    receptor, docking_receptor = Path(config["receptor_pdb"]), Path(config["docking_receptor_pdbqt"])
    validate_receptor(receptor, docking_receptor)
    selected = select_poses(docking_dir, config["receptor_id"], config["engine"],
                            config["top"], legacy, manifest)
    valid = []
    failures = []
    for pose in selected:
        if pose.receptor_sha256 and pose.receptor_sha256 != digest(docking_receptor):
            raise ValueError("Selected pose was docked to a different receptor preparation")
        try:
            entry = ligand_input(config, pose.ligand)
            mol = posed_molecule(Path(entry["sdf"]), entry["atom_map"], pose.coordinates,
                                 entry["net_charge"], config["transform"])
            valid.append((pose, mol))
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            failures.append({"stage": "md-preparation", "ligand": pose.ligand,
                             "pose": pose.record(), "error": str(exc)})
            print(f"[WARN] Skipping MD for {pose.ligand}: {exc}; continuing", flush=True)
    # Validate every template before creating the run directory.
    templates = {stage: stage_mdp(config, stage)
                 for stage in ("ions", "em", "nvt", "npt", "production")}
    output = work_path(output)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "preparation_failures.json", failures)
    with (output / "preparation_failures.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "ligand", "receptor", "binding_site", "error"])
        writer.writeheader()
        for failure in failures:
            writer.writerow({key: failure[key] for key in ("stage", "ligand", "error")} |
                            {key: failure["pose"].get(key, "") for key in ("receptor", "binding_site")})
    if not valid:
        raise ValueError(f"No selected ligands passed MD preparation; see {output / 'preparation_failures.json'}")
    inputs = output / "inputs"
    inputs.mkdir()
    for key in ("receptor_pdb", "docking_receptor_pdbqt", "water_gro", "force_field_dir"):
        if key in config:
            source = Path(config[key])
            target = inputs / (source.name if key == "force_field_dir" else f"{key}{source.suffix}")
            snapshot(source, target)
            config[key] = target.relative_to(output).as_posix()
    transform_pdb(output / config["receptor_pdb"], inputs / "oriented_receptor.pdb", config["transform"])
    if config["system_type"] == "membrane":
        membrane = config["membrane"]
        for key in ("gro", "topology_dir"):
            source = Path(membrane[key])
            target = inputs / ("bilayer.gro" if key == "gro" else "lipid_parameters")
            snapshot(source, target)
            membrane[key] = target.relative_to(output).as_posix()
        for name in membrane["include_files"]:
            file = (output / membrane["topology_dir"] / name).resolve()
            if not file.is_relative_to(inputs.resolve()) or not file.is_file():
                raise ValueError(f"Invalid lipid include: {name}")
            if "#include" in file.read_text():
                raise ValueError("Supply flattened lipid include files with all parameter dependencies")
            split_itp(file.read_text())
    entries = []
    for rank, (pose, mol) in enumerate(valid, 1):
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", pose.ligand)[:80]
        system = output / f"{rank:02d}_{safe}"
        system.mkdir()
        (system / "selected_pose.pdbqt").write_text("MODEL 1\n" +
            "\n".join(pose.coordinates) + "\nENDMDL\n")
        entry = dict(config["ligands"][pose.ligand])
        source = Path(entry["sdf"])
        target = system / "chemical_source.sdf"
        snapshot(source, target)
        entry["sdf"] = target.relative_to(output).as_posix()
        from rdkit import Chem
        mol = posed_molecule(target, entry["atom_map"], pose.coordinates,
                             entry["net_charge"], config["transform"])
        Chem.MolToMolFile(mol, str(system / "posed_ligand.mol"))
        for stage, template in templates.items():
            (system / f"{stage}.mdp").write_text(template)
        entries.append({"directory": system.name, "chemical_input": entry, "pose": pose.record()})
    # Keep only selected chemical inputs; the original protocol is preserved separately.
    config.pop("ligands")
    snapshot(config_path, output / "submitted_protocol.json")
    if manifest is not None:
        snapshot(manifest, output / "docking_manifest.json")
    immutable = {p.relative_to(output).as_posix(): digest(p) for p in output.rglob("*") if p.is_file()}
    job = {"schema": 1, "config": config, "systems": entries, "immutable_inputs": immutable,
           "preparation_failures": failures,
           "created_utc": datetime.now(timezone.utc).isoformat(),
           "workflow_sha256": {p.name: digest(p) for p in ROOT.glob("*.py")}}
    save(output / "job.json", job)
    save(output / "job.sha256.json", {"job.json": digest(output / "job.json")})
    print(f"Prepared {len(entries)} distinct ligands (requested {config['top']}): {output}")
    return output


def start(config_path: Path, docking_dir: Path, *, through="npt", output=None,
          legacy=False, manifest=None, result_file=None, receptor_dir=None,
          gmx="gmx", acpype="acpype", obabel="obabel"):
    """Select docking results, prepare a job and continue directly into GROMACS."""
    status = {"schema": 1, "docking_dir": str(docking_dir.resolve()),
              "protocol": str(config_path.resolve()), "through": through, "status": "checking"}
    def record():
        if result_file is not None:
            result_file.parent.mkdir(parents=True, exist_ok=True)
            save(result_file, status)
    record()
    try:
        config = preflight(config_path, receptor_dir=receptor_dir, tools=through != "prepare",
                           gmx=gmx, acpype=acpype, obabel=obabel)
        status["status"] = "preparing"
        record()
        jobdir = prepare(config_path, docking_dir, output, legacy, manifest, config=config)
        status["preparation_failures"] = json.loads((jobdir / "preparation_failures.json").read_text())
        status["partial"] = bool(status["preparation_failures"])
        status.update(job_dir=str(jobdir), status="prepared" if through == "prepare" else "running")
        record()
        print(f"MD job: {jobdir}", flush=True)
        if through != "prepare":
            run(jobdir, through, False, gmx=gmx, acpype=acpype, obabel=obabel)
            status["status"] = "complete"
        record()
        return jobdir
    except BaseException as exc:
        status.update(status="failed", error=str(exc))
        record()
        raise
    finally:
        if "job_dir" in status:
            command = ["ultidock", "md", "run", status["job_dir"], "--through",
                       "npt" if through == "prepare" else through]
            for name, value in (("gmx", gmx), ("acpype", acpype), ("obabel", obabel)):
                if value != name:
                    command.extend([f"--{name}", value])
            print(f"Resume this job: {shlex.join(command)}", flush=True)
            if status["status"] == "complete" and through == "npt":
                command[command.index("--through") + 1] = "production"
                print("After reviewing equilibration: " + shlex.join(
                    [*command, "--equilibration-reviewed"]), flush=True)


class Runner:
    def __init__(self, system: Path, *, gmx="gmx", acpype="acpype", obabel="obabel"):
        self.system = system
        self.state_path = system / "state.json"
        self.state = json.loads(self.state_path.read_text()) if self.state_path.exists() else {}
        self.tools = {name: shutil.which(value) or value
                      for name, value in (("gmx", gmx), ("acpype", acpype), ("obabel", obabel))}
        self.env = dict(os.environ, GMX_MAXBACKUP="-1")
        scratch = system / "scratch"
        scratch.mkdir(exist_ok=True)
        self.env["TMPDIR"] = str(scratch)

    def command(self, label, command, stdin=""):
        log = self.system / f"{label}.command.log"
        with (self.system / "commands.jsonl").open("a") as record:
            record.write(json.dumps({"time": datetime.now(timezone.utc).isoformat(),
                                     "label": label, "argv": command, "stdin": stdin}) + "\n")
        print(f"[{self.system.name}] {label}", flush=True)
        with log.open("a") as stream:
            result = subprocess.run(command, cwd=self.system, env=self.env,
                                    input=stdin, text=True, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{label} failed ({result.returncode}); see {log}")

    def gmx(self, label, *args, stdin=""):
        self.command(label, [self.tools["gmx"], *map(str, args)], stdin)

    def stage(self, name, action, outputs):
        if self.state.get(name, {}).get("status") == "complete":
            expected = self.state[name]["outputs"]
            if any(not (self.system / p).is_file() or digest(self.system / p) != value
                   for p, value in expected.items()):
                raise ValueError(f"Completed {name} output is missing/changed; use a new MD job")
            return
        self.state[name] = {**self.state.get(name, {}), "status": "running"}
        save(self.state_path, self.state)
        try:
            action()
            if callable(outputs):
                outputs = outputs()
            for p in outputs:
                if not (self.system / p).is_file():
                    raise RuntimeError(f"{name} did not produce {p}")
            self.state[name] = {"status": "complete", "outputs":
                                {p: digest(self.system / p) for p in outputs}}
            save(self.state_path, self.state)
        except BaseException as exc:
            self.state[name].update(status="failed", error=str(exc))
            save(self.state_path, self.state)
            raise


def build_system(runner: Runner, job: dict, entry: dict):
    system, jobdir, config = runner.system, runner.system.parent, job["config"]
    if runner.state.get("build"):
        # Partial system building edits topologies in place; do not guess at recovery.
        if runner.state["build"]["status"] != "complete":
            raise ValueError("Interrupted/failed system build: correct the inputs and prepare a fresh job")

    def build():
        runner.command("gromacs_version", [runner.tools["gmx"], "--version"])
        runner.command("acpype_version", [runner.tools["acpype"], "-v"])
        runner.command("openbabel_version", [runner.tools["obabel"], "-V"])
        if "force_field_dir" in config:
            snapshot(jobdir / config["force_field_dir"], system / f"{config['force_field']}.ff")
        runner.gmx("protein", "pdb2gmx", "-f", str(jobdir / "inputs/oriented_receptor.pdb"),
                   "-o", "protein.gro", "-p", "topol.top", "-i", "posre_protein.itp",
                   "-ff", config["force_field"], "-water", config["water_model"],
                   stdin=config.get("pdb2gmx_answers", ""))
        runner.command("chemical_conversion", [runner.tools["obabel"], "-imol", "posed_ligand.mol",
                                                "-omol2", "-O", "posed_ligand.mol2"])
        chemical = entry["chemical_input"]
        mol = posed_molecule(jobdir / chemical["sdf"], chemical["atom_map"],
                             entry["pose"]["coordinates"], chemical["net_charge"], config["transform"])
        positions = convert_mol2(system / "posed_ligand.mol2", mol)
        runner.command("parameterize", [runner.tools["acpype"], "-i", "posed_ligand.mol2",
                                         "-b", "LIG", "-a", "gaff2", "-c", "bcc",
                                         "-n", str(chemical["net_charge"]), "-o", "gmx"])
        name = restore_ligand(system / "LIG.acpype/LIG_GMX.itp", system / "LIG.acpype/LIG_GMX.gro",
                              positions, chemical["net_charge"], system)
        protein, box = read_gro(system / "protein.gro")
        ligand, _ = read_gro(system / "ligand.gro")
        # pdb2gmx should retain the common coordinate frame. Detect any inadvertent shift.
        input_xyz = [np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) / 10
                     for line in (jobdir / "inputs/oriented_receptor.pdb").read_text().splitlines()
                     if line.startswith("ATOM  ") and
                     not line[12:16].strip().lstrip("0123456789").startswith("H")]
        from scipy.spatial import cKDTree
        distances, _ = cKDTree(np.array([a.xyz for a in protein])).query(input_xyz)
        if distances.max() > 0.003:
            raise ValueError("pdb2gmx changed protein heavy-atom coordinates")
        atoms = protein + ligand
        globals_, includes, molecules = ["ligand_atomtypes.itp"], ["ligand.itp"], [(name, 1)]
        if config["system_type"] == "membrane":
            lipids, box, counts, report = embed_bilayer(atoms, {
                **config, "membrane": {**config["membrane"],
                    "gro": str(jobdir / config["membrane"]["gro"]),
                    "topology_dir": str(jobdir / config["membrane"]["topology_dir"])}})
            atoms += lipids
            molecules += counts
            save(system / "membrane_composition.json", report)
            for index, source in enumerate(config["membrane"]["include_files"]):
                text = (jobdir / config["membrane"]["topology_dir"] / source).read_text()
                global_text, molecule_text = split_itp(text)
                for suffix, content, destination in (("globals", global_text, globals_),
                                                      ("molecules", molecule_text, includes)):
                    if content.strip():
                        file = f"lipid_{index}_{suffix}.itp"
                        (system / file).write_text(content)
                        destination.append(file)
        add_topology(system / "topol.top", globals_, includes, molecules)
        write_gro(system / "complex.gro", atoms, box)
        if config["system_type"] == "soluble":
            runner.gmx("box", "editconf", "-f", "complex.gro", "-o", "box.gro", "-bt",
                       config["box_shape"], "-d", config["padding_nm"])
        else:
            snapshot(system / "complex.gro", system / "box.gro")
        solvent = (str(jobdir / config["water_gro"]) if "water_gro" in config
                   else "tip4p.gro" if config["water_model"] == "tip4pew" else "spc216.gro")
        runner.gmx("solvate", "solvate", "-cp", "box.gro", "-cs", solvent,
                   "-o", "solvated.gro", "-p", "topol.top")
        if config["system_type"] == "membrane":
            count = remove_membrane_core_water(system / "solvated.gro", config["membrane"])
            set_water_count(system / "topol.top", count)
        runner.gmx("ions_input", "grompp", "-f", "ions.mdp", "-c", "solvated.gro",
                   "-p", "topol.top", "-o", "ions.tpr", "-po", "ions.processed.mdp",
                   "-pp", "preion.top")
        charge = topology_charge((system / "preion.top").read_text())
        pairs = None
        salt = ["-conc", config["salt_molar"]]
        if config["salt_basis"] == "water-count":
            atoms, _ = read_gro(system / "solvated.gro")
            sites = 4 if config["water_model"] in {"tip4pew", "opc"} else 3
            water_count = sum(a.resname == "SOL" for a in atoms) // sites
            pairs = round(config["salt_molar"] * water_count / 55.5)
            salt = ["-np", pairs, "-nn", pairs]
        if not charge and (config["salt_molar"] == 0 or pairs == 0):
            snapshot(system / "solvated.gro", system / "system.gro")
        else:
            runner.gmx("ions", "genion", "-s", "ions.tpr", "-o", "system.gro", "-p", "topol.top",
                       "-pname", "NA", "-nname", "CL", *salt,
                       "-neutral", "-seed", config["seed"], stdin="SOL\n")
        save(system / "salt.json", {"requested_molar": config["salt_molar"],
                                   "basis": config["salt_basis"], "salt_pairs": pairs,
                                   "neutralizing_charge": -charge,
                                   "water_count_reference_molar": 55.5 if pairs is not None else None})
        snapshot(system / "system.gro", system / "restraint_reference.gro")
    def build_outputs():
        files = ["system.gro", "topol.top", "restraint_reference.gro"]
        files += [p.relative_to(system).as_posix() for p in system.glob("*.itp")]
        if "force_field_dir" in config:
            files += [p.relative_to(system).as_posix()
                      for p in (system / f"{config['force_field']}.ff").rglob("*") if p.is_file()]
        return files
    runner.stage("build", build, build_outputs)


def simulate(runner: Runner, config: dict, stage: str, previous: str):
    def action():
        checkpoint = runner.system / f"{stage}.cpt"
        finished = runner.state[stage].get("mdrun_outputs")
        if finished:
            if any(not (runner.system / p).exists() or digest(runner.system / p) != value
                   for p, value in finished.items()):
                raise ValueError("Finished dynamics output changed before energy extraction")
            runner.gmx(f"{stage}_energy", "energy", "-f", f"{stage}.edr", "-o", f"{stage}_energy.xvg",
                       stdin="Temperature\nPressure\nDensity\nPotential\n0\n")
            return
        # If mdrun was interrupted, reuse the exact TPR and its checkpoint.
        if not checkpoint.exists():
            if (runner.system / f"{stage}.gro").exists():
                raise ValueError(f"Untracked {stage} output; prepare a fresh job")
            args = ["grompp", "-f", f"{stage}.mdp", "-c", f"{previous}.gro",
                    "-p", "topol.top", "-o", f"{stage}.tpr", "-po", f"{stage}.processed.mdp"]
            if stage in {"nvt", "npt"}:
                args += ["-r", "restraint_reference.gro"]
            if stage in {"npt", "production"}:
                args += ["-t", f"{previous}.cpt"]
            runner.gmx(f"{stage}_input", *args)
            runner.state[stage]["tpr_sha256"] = digest(runner.system / f"{stage}.tpr")
            save(runner.state_path, runner.state)
        elif not (runner.system / f"{stage}.tpr").exists():
            raise ValueError("Checkpoint has no matching TPR")
        elif runner.state[stage].get("tpr_sha256") != digest(runner.system / f"{stage}.tpr"):
            raise ValueError("TPR changed or does not belong to this checkpointed stage")
        # GROMACS rejects conflicting environment and command-line budgets.
        runner.env["OMP_NUM_THREADS"] = str(config["threads"])
        args = ["mdrun", "-deffnm", stage, "-nt", config["threads"],
                "-ntomp", config["threads"]]
        if checkpoint.exists():
            args += ["-cpi", checkpoint.name, "-append"]
        runner.gmx(stage, *args)
        if stage == "em":
            text = (runner.system / "em.log").read_text()
            if not re.search(r"converged to Fmax", text):
                raise RuntimeError("Energy minimization did not converge; inspect em.log")
        else:
            runner.state[stage]["mdrun_outputs"] = {
                f"{stage}.{suffix}": digest(runner.system / f"{stage}.{suffix}")
                for suffix in ("gro", "tpr", "edr", "log", "cpt")}
            save(runner.state_path, runner.state)
            runner.gmx(f"{stage}_energy", "energy", "-f", f"{stage}.edr", "-o", f"{stage}_energy.xvg",
                       stdin="Temperature\nPressure\nDensity\nPotential\n0\n")
    outputs = [f"{stage}.gro", f"{stage}.tpr", f"{stage}.edr", f"{stage}.log"]
    if stage != "em":
        outputs += [f"{stage}.cpt", f"{stage}_energy.xvg"]
    runner.stage(stage, action, outputs)


def run(jobdir: Path, through: str, reviewed: bool, *, gmx="gmx", acpype="acpype", obabel="obabel"):
    jobdir = work_path(jobdir)
    if through == "production" and not reviewed:
        raise ValueError("Inspect the built complex, membrane composition/pore hydration and NPT "
                         "equilibration first; pass --equilibration-reviewed to start production")
    if not check_dependencies(gmx, acpype, obabel, workspace=ROOT / "workspace"):
        raise ValueError("Install the missing MD dependencies; see md-simulation/README.md")
    with (jobdir / ".run.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("This MD job is already running") from None
        checksum = json.loads((jobdir / "job.sha256.json").read_text())["job.json"]
        if checksum != digest(jobdir / "job.json"):
            raise ValueError("Job configuration changed; prepare a fresh job")
        job = json.loads((jobdir / "job.json").read_text())
        for name, expected in job["workflow_sha256"].items():
            if digest(ROOT / name) != expected:
                raise ValueError("MD workflow code changed; prepare a fresh job")
        for name, expected in job["immutable_inputs"].items():
            path = (jobdir / name).resolve()
            if not path.is_relative_to(jobdir) or not path.is_file() or digest(path) != expected:
                raise ValueError(f"Input changed or missing: {name}; prepare a fresh job")
        stages = ["build", "em", "nvt", "npt", "production"]
        summary = {"through": through, "completed": [], "failed": [],
                   "preparation_failures": job.get("preparation_failures", [])}
        for entry in job["systems"]:
            try:
                runner = Runner(jobdir / entry["directory"], gmx=gmx, acpype=acpype, obabel=obabel)
                build_system(runner, job, entry)
                previous = "system"
                for stage in stages[1:stages.index(through) + 1]:
                    if stage == "production":
                        save(runner.system / "equilibration_review.json", {
                            "reviewed_utc": datetime.now(timezone.utc).isoformat(),
                            "npt_checkpoint_sha256": digest(runner.system / "npt.cpt")})
                    simulate(runner, job["config"], stage, previous)
                    previous = stage
                summary["completed"].append(entry["directory"])
            except Exception as exc:
                summary["failed"].append({"system": entry["directory"], "error": str(exc)})
                print(f"[WARN] MD failed for {entry['directory']}: {exc}; continuing with other systems",
                      flush=True)
            save(jobdir / "run_summary.json", summary)
        if summary["failed"]:
            raise RuntimeError(f"MD finished with {len(summary['failed'])} failed systems and "
                               f"{len(summary['completed'])} completed; see {jobdir / 'run_summary.json'}")
    print(f"Completed through {through}: {jobdir}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="Validate MD inputs and tools before docking")
    check.add_argument("--config", type=Path, required=True)
    check.add_argument("--receptor-dir", type=Path)
    check.add_argument("--ligands-dir", type=Path)
    check.add_argument("--engine", choices=["vina", "adgpu"])
    check.add_argument("--work-dir", type=Path)
    check.add_argument("--inputs-only", action="store_true")
    check.add_argument("--gmx", default="gmx")
    check.add_argument("--acpype", default="acpype")
    check.add_argument("--obabel", default="obabel")
    select = commands.add_parser("select", help="List the best pose of each top ligand")
    select.add_argument("--docking-dir", type=Path, required=True)
    select.add_argument("--receptor-id", required=True)
    select.add_argument("--engine", choices=["vina", "adgpu"], required=True)
    select.add_argument("--top", type=int, default=5)
    select.add_argument("--legacy-single-receptor", action="store_true")
    setup = commands.add_parser("prepare", help="Validate and snapshot a job; start no simulation")
    setup.add_argument("--config", type=Path, required=True)
    setup.add_argument("--docking-dir", type=Path, required=True)
    setup.add_argument("--work-dir", type=Path)
    setup.add_argument("--legacy-single-receptor", action="store_true")
    execute = commands.add_parser("run", help="Continue docking results into MD, or resume a prepared job")
    execute.add_argument("work_dir", type=Path, nargs="?")
    execute.add_argument("--config", type=Path, help="Prepare and run directly from docking results")
    execute.add_argument("--docking-dir", type=Path)
    execute.add_argument("--receptor-dir", type=Path, help="Bind the receptor staged by this docking run")
    execute.add_argument("--work-dir", dest="new_work_dir", type=Path)
    execute.add_argument("--pose-manifest", type=Path, help="Restrict selection to this docking run")
    execute.add_argument("--result-file", type=Path, help="Write the MD job path and handoff status as JSON")
    execute.add_argument("--legacy-single-receptor", action="store_true")
    execute.add_argument("--through", choices=["prepare", "build", "em", "nvt", "npt", "production"], default="npt")
    execute.add_argument("--equilibration-reviewed", action="store_true")
    execute.add_argument("--gmx", default="gmx")
    execute.add_argument("--acpype", default="acpype")
    execute.add_argument("--obabel", default="obabel")
    args = parser.parse_args(argv)
    try:
        if args.command == "check":
            if args.work_dir is not None and work_path(args.work_dir).exists():
                raise ValueError("MD work directory already exists; resume it with ultidock md run PATH")
            preflight(args.config, receptor_dir=args.receptor_dir, ligands_dir=args.ligands_dir,
                      engine=args.engine, tools=not args.inputs_only,
                      gmx=args.gmx, acpype=args.acpype, obabel=args.obabel)
            print("MD preflight passed")
        elif args.command == "select":
            poses = select_poses(args.docking_dir, args.receptor_id, args.engine,
                                  args.top, args.legacy_single_receptor)
            print(json.dumps([p.record() for p in poses], indent=2))
        elif args.command == "prepare":
            prepare(args.config, args.docking_dir, args.work_dir, args.legacy_single_receptor)
        elif args.config is not None:
            if args.work_dir is not None or args.docking_dir is None:
                raise ValueError("Use --config with --docking-dir for a new job; "
                                 "use a positional job directory to resume")
            if args.through == "production" or args.equilibration_reviewed:
                raise ValueError("A new MD job runs through NPT first. Review its equilibration, "
                                 "then resume with --through production --equilibration-reviewed")
            start(args.config, args.docking_dir, through=args.through, output=args.new_work_dir,
                  legacy=args.legacy_single_receptor, manifest=args.pose_manifest,
                  result_file=args.result_file, receptor_dir=args.receptor_dir,
                  gmx=args.gmx, acpype=args.acpype, obabel=args.obabel)
        else:
            if (args.work_dir is None or args.docking_dir is not None or args.new_work_dir is not None
                    or args.pose_manifest is not None or args.result_file is not None
                    or args.receptor_dir is not None
                    or args.legacy_single_receptor or args.through == "prepare"):
                raise ValueError("Use ultidock md run JOB_DIR to resume, or "
                                 "ultidock md run --config PROTOCOL --docking-dir DOCKING_DIR")
            run(args.work_dir, args.through, args.equilibration_reviewed,
                gmx=args.gmx, acpype=args.acpype, obabel=args.obabel)
    except (ValueError, KeyError, OSError, RuntimeError, TypeError) as exc:
        print(f"MD ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
