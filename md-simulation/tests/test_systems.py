import numpy as np
import pytest

from protocol import stage_mdp
from systems import (GroAtom, add_topology, embed_bilayer, read_gro,
                     remove_membrane_core_water, restore_ligand, set_water_count, write_gro)
from systems import topology_charge
from workflow import Runner, simulate


def config(system_type="soluble"):
    return {"system_type": system_type, "temperature_k": 310, "dt_ps": 0.002,
            "seed": 7, "production_ns": 100, "nvt_ps": 100, "npt_ps": 1000,
            "membrane": {"mdp_nonbonded": {"rvdw": 1, "rcoulomb": 1,
                         "vdw-modifier": "Potential-shift", "DispCorr": "EnerPres"}}}


def test_ensembles_duration_restraints_and_velocities():
    settings = config()
    nvt, npt, production = (stage_mdp(settings, s) for s in ("nvt", "npt", "production"))
    assert "gen-vel = yes" in nvt and "gen-temp = 310" in nvt
    assert "POSRES_LIG" in nvt and "POSRES_LIG" in npt and "POSRES" not in production
    assert "gen-vel = no" in npt and "gen-vel = no" in production
    assert "pcoupl = C-rescale" in production and "pcoupltype = isotropic" in production
    assert "nsteps = 50000000" in production
    membrane = stage_mdp(config("membrane"), "npt")
    assert "pcoupltype = semiisotropic" in membrane and "ref-p = 1 1" in membrane
    assert "refcoord-scaling = com" in membrane
    assert "DispCorr = EnerPres" in membrane  # Lipid profile is explicitly supplied.


def test_topology_global_parameters_precede_protein_molecules(tmp_path):
    top = tmp_path / "topol.top"
    top.write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n[ moleculetype ]\n'
                   'Protein 3\n[ atoms ]\n[ system ]\ncomplex\n[ molecules ]\nProtein 1\n')
    add_topology(top, ["ligand_atomtypes.itp"], ["ligand.itp"], [("LIG", 1)])
    text = top.read_text()
    assert text.index("ligand_atomtypes.itp") < text.index("[ moleculetype ]")
    assert text.index("ligand.itp") < text.index("[ system ]")
    assert text.endswith("LIG 1\n")


def test_acpype_reordering_keeps_pose_and_charge(tmp_path):
    itp, gro = tmp_path / "lig.itp", tmp_path / "lig.gro"
    itp.write_text('[ atomtypes ]\nc3 12 0 A 0.3 0.4\n[ moleculetype ]\nLIG 3\n'
                   '[ atoms ]\n1 c3 1 LIG C1 1 0 12\n2 hc 1 LIG H2 1 0 1\n')
    write_gro(gro, [GroAtom(1, "LIG", "C1", np.zeros(3)),
                    GroAtom(1, "LIG", "H2", np.zeros(3))], np.eye(3))
    positions = {"C1": np.array([4, 5, 6]), "H2": np.array([4.1, 5, 6])}
    assert restore_ligand(itp, gro, positions, 0, tmp_path) == "LIG"
    atoms, _ = read_gro(tmp_path / "ligand.gro")
    assert np.allclose(atoms[0].xyz, positions["C1"])
    assert "1 1 1000 1000 1000" in (tmp_path / "posre_ligand.itp").read_text()
    with pytest.raises(ValueError, match="charge"):
        restore_ligand(itp, gro, positions, 1, tmp_path)


def test_membrane_removes_whole_overlapping_molecules_and_keeps_pose(tmp_path):
    top = tmp_path / "lipid.itp"
    top.write_text('[ moleculetype ]\nPOPC 3\n[ atoms ]\n'
                   '1 c 1 POPC C1 1 0 12\n2 c 1 POPC C2 1 0 12\n')
    gro = tmp_path / "seed.gro"
    lipids = [GroAtom(1, "POPC", "C1", np.array([5, 5, 5])),
              GroAtom(1, "POPC", "C2", np.array([5, 5, 5.1])),
              GroAtom(2, "POPC", "C1", np.array([8, 8, 4])),
              GroAtom(2, "POPC", "C2", np.array([8, 8, 4.1]))]
    write_gro(gro, lipids, np.eye(3) * 10)
    solute = [GroAtom(1, "LIG", "C1", np.array([5, 5, 5]))]
    settings = {"box_shape": "triclinic", "padding_nm": 1.2, "membrane": {
        "gro": str(gro), "topology_dir": str(tmp_path), "include_files": [top.name],
        "molecules": [{"name": "POPC", "count": 2, "atom_count": 2}],
        "clash_distance_nm": 0.2, "hydrophobic_z_nm": [4, 6],
        "mdp_nonbonded": {"rvdw": 1.2, "rcoulomb": 1.2}}}
    kept, box, counts, report = embed_bilayer(solute, settings)
    assert len(kept) == 2 and counts == [("POPC", 1)]
    assert report[0]["retained"] == 1 and report[0]["lower_leaflet"] == 1
    assert np.array_equal(solute[0].xyz, [5, 5, 5])
    assert np.array_equal(box, np.eye(3) * 10)


def test_membrane_solvent_filter_keeps_aqueous_and_declared_pore_water(tmp_path):
    gro = tmp_path / "solvated.gro"
    atoms = []
    for i, xyz in enumerate(([2, 2, 2], [2, 2, 5], [5, 5, 5]), 1):
        atoms += [GroAtom(i, "SOL", name, np.array(xyz, dtype=float)) for name in ("OW", "HW1", "HW2")]
    write_gro(gro, atoms, np.eye(3) * 10)
    count = remove_membrane_core_water(gro, {"hydrophobic_z_nm": [4, 6], "water_sites": 3,
                                            "water_pores": [{"xy_nm": [5, 5], "radius_nm": 0.5}]})
    assert count == 2 and len(read_gro(gro)[0]) == 6
    top = tmp_path / "topol.top"
    top.write_text('[ molecules ]\nProtein 1\nSOL 3\n')
    set_water_count(top, count)
    assert top.read_text().rstrip().endswith("SOL 2")


def test_command_failure_stops_pipeline_and_completed_outputs_cannot_change(tmp_path):
    runner = Runner(tmp_path)
    with pytest.raises(RuntimeError, match="failed"):
        runner.stage("em", lambda: runner.command("failure", ["/bin/false"]), ["em.gro"])
    assert runner.state["em"]["status"] == "failed"
    file = tmp_path / "em.gro"
    runner.stage("em", lambda: file.write_text("complete"), ["em.gro"])
    file.write_text("changed")
    with pytest.raises(ValueError, match="missing/changed"):
        runner.stage("em", lambda: pytest.fail("must not rerun"), ["em.gro"])


def test_npt_passes_previous_velocities_and_never_bypasses_warnings(tmp_path):
    runner, calls = Runner(tmp_path), []
    def gmx(label, *args, stdin=""):
        calls.append(args)
        if args[0] == "grompp":
            (tmp_path / "npt.tpr").write_text("tpr")
        elif args[0] == "mdrun":
            for suffix in ("gro", "tpr", "edr", "log", "cpt"):
                (tmp_path / f"npt.{suffix}").write_text("output")
        elif args[0] == "energy":
            (tmp_path / "npt_energy.xvg").write_text("energies")
    runner.gmx = gmx
    settings = {**config(), "threads": 1}
    simulate(runner, settings, "npt", "nvt")
    assert "-t" in calls[0] and "nvt.cpt" in calls[0]
    assert all("-maxwarn" not in args for args in calls)


def test_rounded_amber_defaults_and_integral_charge():
    text = ('[ defaults ]\n1 2 yes 0.5 0.8333\n[ moleculetype ]\nLIG 3\n'
            '[ atoms ]\n1 c3 1 LIG C1 1 1 12\n[ molecules ]\nLIG 2\n')
    assert topology_charge(text) == 2
    with pytest.raises(ValueError, match="incompatible"):
        topology_charge(text.replace("0.8333", "0.5"))
    with pytest.raises(ValueError, match="non-integer"):
        topology_charge(text.replace("C1 1 1 12", "C1 1 0.2 12"))


def test_interrupted_dynamics_reuses_same_tpr_and_checkpoint(tmp_path):
    runner, calls = Runner(tmp_path), []
    fail = [True]
    def gmx(label, *args, stdin=""):
        calls.append(args)
        if args[0] == "grompp":
            (tmp_path / "npt.tpr").write_text("original TPR")
        elif args[0] == "mdrun":
            (tmp_path / "npt.cpt").write_text("checkpoint")
            if fail[0]:
                fail[0] = False
                raise RuntimeError("interrupted")
            for suffix in ("gro", "edr", "log"):
                (tmp_path / f"npt.{suffix}").write_text("finished")
        elif args[0] == "energy":
            (tmp_path / "npt_energy.xvg").write_text("energies")
    runner.gmx = gmx
    with pytest.raises(RuntimeError, match="interrupted"):
        simulate(runner, {**config(), "threads": 1}, "npt", "nvt")
    simulate(runner, {**config(), "threads": 1}, "npt", "nvt")
    assert sum(args[0] == "grompp" for args in calls) == 1
    assert "-cpi" in calls[-2] and "-append" in calls[-2]
