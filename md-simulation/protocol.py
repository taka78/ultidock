"""Validated Amber/GAFF2 profiles and conservative, editable MD stage templates."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

PROFILES = {
    "amber99sb-ildn": {"tip3p", "tip4pew", "spce"},
    "amber14sb": {"tip3p", "tip4pew"},
    "amber19sb": {"opc", "opc3"},
}


def load_protocol(path: Path) -> dict:
    config = json.loads(path.read_text())
    base = path.resolve().parent
    config.setdefault("force_field", "amber99sb-ildn")
    config.setdefault("water_model", "tip3p")
    config.setdefault("box_shape", "dodecahedron")
    defaults = {"padding_nm": 1.2, "salt_molar": 0.15, "temperature_k": 310,
                "production_ns": 100, "nvt_ps": 100, "npt_ps": 1000,
                "dt_ps": 0.002, "seed": 2026, "top": 5, "threads": 1}
    for key, value in defaults.items():
        config.setdefault(key, value)
        if isinstance(config[key], bool) or not isinstance(config[key], (int, float)):
            raise ValueError(f"{key} must be numeric")
        if not math.isfinite(config[key]) or config[key] < (0 if key == "salt_molar" else 1e-12):
            raise ValueError(f"Invalid {key}")
    if config["dt_ps"] > 0.002:
        raise ValueError("This constrained all-atom protocol supports timesteps up to 2 fs")
    for key in ("seed", "top", "threads"):
        if not isinstance(config[key], int) or config[key] > 2147483647:
            raise ValueError(f"{key} must be a positive 32-bit integer")
    if config["padding_nm"] < 1.2:
        raise ValueError("padding_nm must be at least the 1.2 nm nonbonded cutoff")
    if config["system_type"] not in {"soluble", "membrane"}:
        raise ValueError("system_type must explicitly be soluble or membrane")
    config.setdefault("salt_basis", "water-count" if config["system_type"] == "membrane" else "box-volume")
    if config["salt_basis"] not in {"box-volume", "water-count"}:
        raise ValueError("salt_basis must be box-volume or water-count")
    if config.get("receptor_reviewed") is not True:
        raise ValueError("Set receptor_reviewed=true after reviewing completeness, protonation, "
                         "termini, disulfides and absence of required cofactors")
    if config["engine"] not in {"vina", "adgpu"}:
        raise ValueError("Choose one scoring engine: vina or adgpu")
    ff, water = config["force_field"], config["water_model"]
    if ff not in PROFILES:
        raise ValueError("Automatic GAFF2 parameters currently require amber99sb-ildn, "
                         "amber14sb or amber19sb. CHARMM/CGenFF, OPLS and GROMOS "
                         "need separate compatible parameterization providers")
    if water not in PROFILES[ff]:
        raise ValueError(f"Unsupported water/force-field combination: {ff}/{water}")
    if config["box_shape"] not in {"cubic", "dodecahedron", "triclinic"}:
        raise ValueError("box_shape must be cubic, dodecahedron or triclinic")

    def file(key, obj=config, directory=False):
        if key in obj:
            value = (base / obj[key]).resolve()
            if not (value.is_dir() if directory else value.is_file()):
                raise ValueError(f"Missing input {key}: {value}")
            obj[key] = str(value)

    for key in ("receptor_pdb", "docking_receptor_pdbqt", "water_gro"):
        file(key)
    file("force_field_dir", directory=True)
    if ff != "amber99sb-ildn" and "force_field_dir" not in config:
        raise ValueError(f"Supply force_field_dir for an installed {ff}.ff parameter set")
    if "force_field_dir" in config and Path(config["force_field_dir"]).name != f"{ff}.ff":
        raise ValueError("force_field_dir basename must match force_field + '.ff'")
    if water.startswith("opc") and "water_gro" not in config:
        raise ValueError("OPC models require a matching water_gro solvent coordinate template")
    for ligand, entry in config["ligands"].items():
        file("sdf", entry)
        if not isinstance(entry.get("net_charge"), int) or isinstance(entry["net_charge"], bool):
            raise ValueError(f"Provide an integer net_charge for {ligand}")
        if not isinstance(entry.get("atom_map"), dict):
            raise ValueError(f"Provide an explicit docking-serial -> SDF-index atom_map for {ligand}")

    transform = np.eye(4)
    if config["system_type"] == "membrane":
        membrane = config["membrane"]
        if config["box_shape"] == "dodecahedron":
            raise ValueError("Membrane bilayers require a planar periodic box, not dodecahedron")
        if membrane.get("orientation_reviewed") is not True:
            raise ValueError("Review the protein's membrane orientation and bilayer frame")
        if membrane.get("compatible_force_field") != ff:
            raise ValueError("Bilayer lipid parameters must explicitly match the protein force field")
        file("gro", membrane)
        file("topology_dir", membrane, directory=True)
        if not membrane.get("include_files") or not membrane.get("molecules"):
            raise ValueError("Membrane requires ordered include_files and molecule blocks")
        membrane.setdefault("clash_distance_nm", 0.2)
        membrane["water_sites"] = 4 if water in {"tip4pew", "opc"} else 3
        low, high = membrane["hydrophobic_z_nm"]
        if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in (low, high)):
            raise ValueError("Provide finite membrane hydrophobic_z_nm bounds")
        for pore in membrane.get("water_pores", []):
            values = [*pore["xy_nm"], pore["radius_nm"]]
            if len(values) != 3 or not all(math.isfinite(v) for v in values) or values[-1] <= 0:
                raise ValueError("Invalid water pore cylinder")
        required = {"rvdw", "rcoulomb", "vdw-modifier", "DispCorr"}
        nonbonded = membrane.get("mdp_nonbonded", {})
        if not required <= nonbonded.keys() or not nonbonded.keys() <= required | {"rvdw-switch"}:
            raise ValueError("Supply membrane mdp_nonbonded settings recommended by the lipid "
                             "parameter authors: rvdw, rcoulomb, vdw-modifier, DispCorr; "
                             "optionally rvdw-switch")
        if nonbonded["DispCorr"] not in {"no", "EnerPres", "Ener"} or nonbonded["vdw-modifier"] not in {
            "Potential-shift", "Force-switch", "Potential-switch"}:
            raise ValueError("Unsupported lipid nonbonded settings")
        for key in ("rvdw", "rcoulomb", "rvdw-switch"):
            if key in nonbonded and (not isinstance(nonbonded[key], (int, float)) or
                                    not math.isfinite(nonbonded[key]) or nonbonded[key] <= 0):
                raise ValueError("Invalid lipid cutoff")
        if max(nonbonded["rvdw"], nonbonded["rcoulomb"]) > config["padding_nm"]:
            raise ValueError("Periodic padding must cover the lipid nonbonded cutoff")
        if "switch" in nonbonded["vdw-modifier"].lower() and not (
            0 < nonbonded.get("rvdw-switch", 0) < nonbonded["rvdw"]
        ):
            raise ValueError("Switching requires 0 < rvdw-switch < rvdw")
        clash = membrane["clash_distance_nm"]
        if not isinstance(clash, (int, float)) or not math.isfinite(clash) or clash <= 0:
            raise ValueError("Invalid membrane clash_distance_nm")
        transform = np.asarray(membrane["transform"], dtype=float)
    if transform.shape != (4, 4) or not np.isfinite(transform).all():
        raise ValueError("transform must be a finite 4x4 rigid transform (Angstrom coordinates)")
    rotation = transform[:3, :3]
    if not (np.allclose(transform[3], [0, 0, 0, 1]) and
            np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6) and
            np.isclose(np.linalg.det(rotation), 1, atol=1e-6)):
        raise ValueError("Membrane transform must preserve distances and handedness")
    config["transform"] = transform.tolist()
    return config


def stage_mdp(config: dict, stage: str) -> str:
    """One thermostat bath; no independently thermostatted ligand."""
    if stage not in {"ions", "em", "nvt", "npt", "production"}:
        raise ValueError(f"Unknown MD stage {stage}")
    common = {"cutoff-scheme": "Verlet", "nstlist": 20, "pbc": "xyz",
              "coulombtype": "PME", "rcoulomb": 1.2, "vdwtype": "Cut-off",
              "vdw-modifier": "Potential-shift", "rvdw": 1.2,
              "DispCorr": "EnerPres", "constraints": "h-bonds",
              "constraint-algorithm": "lincs", "lincs-order": 4}
    if config["system_type"] == "membrane":
        common.update(config["membrane"]["mdp_nonbonded"])
    if stage in {"ions", "em"}:
        common.update({"integrator": "steep", "nsteps": 50000, "emtol": 1000, "emstep": 0.01})
        if stage == "ions":
            # This TPR is only a genion input: avoid a charged-system PME warning.
            common["coulombtype"] = "Cut-off"
    else:
        duration_ps = (config["production_ns"] * 1000 if stage == "production"
                       else config[f"{stage}_ps"])
        steps = round(duration_ps / config["dt_ps"])
        if steps < 1 or not math.isclose(steps * config["dt_ps"], duration_ps, abs_tol=1e-8):
            raise ValueError("Stage durations must be integer multiples of dt_ps")
        common.update({"integrator": "md", "dt": config["dt_ps"], "nsteps": steps,
                       "continuation": "no" if stage == "nvt" else "yes",
                       "tcoupl": "v-rescale", "tc-grps": "System", "tau-t": 0.5,
                       "ref-t": config["temperature_k"], "ld-seed": config["seed"],
                       "gen-vel": "yes" if stage == "nvt" else "no", "pcoupl": "no",
                       "nstenergy": 500, "nstlog": 500,
                       "nstxout-compressed": 5000, "compressed-x-grps": "System",
                       "nstxout": 0, "nstvout": 0, "nstfout": 0})
        if stage == "nvt":
            common.update({"gen-temp": config["temperature_k"], "gen-seed": config["seed"]})
        if stage != "production":
            common["define"] = "-DPOSRES -DPOSRES_LIG"
        if stage in {"npt", "production"}:
            membrane = config["system_type"] == "membrane"
            common.update({"pcoupl": "C-rescale", "tau-p": 5,
                           "pcoupltype": "semiisotropic" if membrane else "isotropic",
                           "ref-p": "1 1" if membrane else "1",
                           "compressibility": "4.5e-5 4.5e-5" if membrane else "4.5e-5"})
            if stage == "npt":
                common["refcoord-scaling"] = "com"
    return "\n".join(f"{key} = {value}" for key, value in common.items()) + "\n"
