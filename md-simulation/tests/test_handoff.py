from pathlib import Path

import numpy as np
import pytest

from chemistry import dock_atoms, posed_molecule, validate_receptor
from poses import select_poses
from docking.pose_provenance import record_pose_provenance

REPO = Path(__file__).resolve().parents[2]


def atom(serial, element="C", xyz=(1, 2, 3), name=None):
    return (f"ATOM  {serial:5d} {name or element:>4s} LIG A   1    " +
            "".join(f"{v:8.3f}" for v in xyz) + f"  1.00  0.00     0.000 {element}")


def vina(path, scores):
    path.write_text("".join(f"MODEL {i}\nREMARK VINA RESULT: {score} 0 0\n"
                            f"{atom(1, xyz=(i, 2, 3))}\nENDMDL\n"
                            for i, score in enumerate(scores, 1)))


def test_top_five_unique_and_exact_best_model(tmp_path):
    for i in range(7):
        vina(tmp_path / f"protein__S1__lig{i}_12345678.pdbqt", [-i, -i - 0.5])
    vina(tmp_path / "protein__S2__lig6_abcdef01.pdbqt", [-20, -19])
    vina(tmp_path / "other__S1__lig0_12345678.pdbqt", [-100])
    poses = select_poses(tmp_path, "protein", "vina")
    assert [p.ligand for p in poses] == ["lig6", "lig5", "lig4", "lig3", "lig2"]
    assert poses[0].binding_site == "2" and poses[0].model == 1 and poses[0].score == -20
    assert poses[1].model == 2
    assert dock_atoms(poses[1].coordinates)[1][1][0] == 2


def test_real_vina_fixture():
    directory = REPO / "examples/sert-escitalopram/expected-results/vina"
    poses = select_poses(directory, "5i6x_edited", "vina", legacy=True)
    assert len(poses) == 1  # Six sites are still one ligand.
    pose = poses[0]
    assert (pose.ligand, pose.binding_site, pose.model, pose.score) == ("escitalopram-e", "6", 1, -9.644)
    assert np.allclose(dock_atoms(pose.coordinates)[1][1], [-31.354, -21.390, 0.998])


def test_real_gpu_fixture_uses_run_501_not_first_best_pdbqt():
    directory = REPO / "examples/sert-escitalopram/expected-results/AD-Gpu"
    pose = select_poses(directory, "5i6x_edited", "adgpu", legacy=True)[0]
    assert (pose.binding_site, pose.model, pose.score) == ("3", 501, -6.91)
    dlg = Path(pose.artifact).read_text()
    exact_block = dlg.split("DOCKED: USER    Run = 501\n")[1].split("DOCKED: ENDMDL")[0]
    assert all(f"DOCKED: {line}" in exact_block for line in pose.coordinates)


def test_gpu_missing_exact_run_fails(tmp_path):
    xml = tmp_path / "protein__S1__lig_12345678.xml"
    xml.write_text('<autodock_gpu><runs><run id="42"><free_NRG_binding>-9</free_NRG_binding>'
                   '</run></runs></autodock_gpu>')
    with pytest.raises(ValueError, match="Exact AutoDock-GPU"):
        select_poses(tmp_path, "protein", "adgpu")
    xml.with_suffix(".dlg").write_text("DOCKED: MODEL 1\nDOCKED: USER    Run = 1\n"
                                      f"DOCKED: {atom(1)}\nDOCKED: ENDMDL\n")
    with pytest.raises(ValueError, match="has no coordinates"):
        select_poses(tmp_path, "protein", "adgpu")


def test_legacy_ambiguity_and_nonfinite_fail(tmp_path):
    output = tmp_path / "S1__lig_12345678.pdbqt"
    vina(output, [-9])
    with pytest.raises(ValueError, match="Legacy"):
        select_poses(tmp_path, "protein", "vina")
    vina(output, [float("nan")])
    with pytest.raises(ValueError, match="non-finite"):
        select_poses(tmp_path, "protein", "vina", legacy=True)


@pytest.mark.parametrize("change", ["coordinates", "missing"])
def test_gpu_provenance_detects_changed_coordinates_with_unchanged_scores(tmp_path, change):
    receptor, ligand = tmp_path / "receptor.pdbqt", tmp_path / "ligand.pdbqt"
    receptor.write_text(atom(1))
    ligand.write_text(atom(1))
    output = tmp_path / "protein__S1__ligand_12345678.xml"
    output.write_text('<autodock_gpu><runs><run id="42"><free_NRG_binding>-9</free_NRG_binding>'
                      '</run></runs></autodock_gpu>')
    dlg = output.with_suffix(".dlg")
    dlg.write_text("DOCKED: MODEL 1\nDOCKED: USER    Run = 42\n"
                   "DOCKED: USER Estimated Free Energy of Binding = -9 kcal/mol\n"
                   f"DOCKED: {atom(1)}\nDOCKED: ENDMDL\n")
    record_pose_provenance(output, ligand, receptor, "protein", "S1", "adgpu")
    assert select_poses(tmp_path, "protein", "adgpu")[0].model == 42
    if change == "coordinates":
        dlg.write_text(dlg.read_text().replace(atom(1), atom(1, xyz=(20, 30, 40))))
    else:
        dlg.unlink()
    with pytest.raises(ValueError, match="coordinates changed or missing"):
        select_poses(tmp_path, "protein", "adgpu")


def test_gpu_output_without_coordinates_cannot_be_marked_successful(tmp_path):
    output, ligand, receptor = (tmp_path / name for name in ("result.xml", "ligand.pdbqt", "receptor.pdbqt"))
    for path in (output, ligand, receptor):
        path.write_text("fixture")
    with pytest.raises(ValueError, match="Exact AutoDock-GPU coordinates require"):
        record_pose_provenance(output, ligand, receptor, "protein", "S1", "adgpu")
    assert not output.with_suffix(".md.json").exists()


def test_provenance_detects_changed_output_and_mixed_receptors(tmp_path):
    receptor, ligand = tmp_path / "receptor.pdbqt", tmp_path / "input.pdbqt"
    receptor.write_text(atom(1))
    ligand.write_text(atom(1))
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    output = outputs / "protein__S1__input_12345678.pdbqt"
    vina(output, [-9])
    record_pose_provenance(output, ligand, receptor, "protein", "S1", "vina")
    assert select_poses(outputs, "protein", "vina")[0].artifact != str(ligand)
    other = outputs / "protein__S2__input_abcdef01.pdbqt"
    vina(other, [-8])
    receptor.write_text(atom(1, xyz=(4, 5, 6)))
    record_pose_provenance(other, ligand, receptor, "protein", "S2", "vina")
    with pytest.raises(ValueError, match="Different receptor"):
        select_poses(outputs, "protein", "vina")
    output.write_text(output.read_text() + "REMARK altered\n")
    with pytest.raises(ValueError, match="artifact changed"):
        select_poses(outputs, "protein", "vina")


def ethanol_source(path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=17)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return mol


def test_chemical_source_does_not_supply_md_heavy_coordinates(tmp_path):
    sdf = tmp_path / "source.sdf"
    source = ethanol_source(sdf)
    original = source.GetConformer().GetPositions()[:3]
    docked = original + np.array([20, -30, 11])
    lines = [atom(i + 1, source.GetAtomWithIdx(i).GetSymbol(), xyz) for i, xyz in enumerate(docked)]
    transform = np.eye(4)
    transform[:3, 3] = [10, 5, -3]
    result = posed_molecule(sdf, {"1": 0, "2": 1, "3": 2}, lines, 0, transform)
    assert np.allclose(result.GetConformer().GetPositions()[:3],
                       np.round(docked, 3) + transform[:3, 3], atol=1e-6)
    assert result.GetNumAtoms() == source.GetNumAtoms()
    with pytest.raises(ValueError, match="net charge"):
        posed_molecule(sdf, {"1": 0, "2": 1, "3": 2}, lines, 1, transform)
    with pytest.raises(ValueError, match="element mismatch"):
        posed_molecule(sdf, {"1": 0, "2": 2, "3": 1}, lines, 0, transform)
    with pytest.raises(ValueError, match="cover each"):
        posed_molecule(sdf, {"1": 0, "2": 1}, lines, 0, transform)


def test_receptor_frame_mismatch_fails(tmp_path):
    pdb, pdbqt = tmp_path / "protein.pdb", tmp_path / "protein.pdbqt"
    pdb.write_text(atom(1))
    pdbqt.write_text(atom(1, xyz=(2, 2, 3)))
    with pytest.raises(ValueError, match="coordinate frame"):
        validate_receptor(pdb, pdbqt)


def test_polar_hydrogen_pose_is_retained(tmp_path):
    sdf = tmp_path / "source.sdf"
    source = ethanol_source(sdf)
    xyz = source.GetConformer().GetPositions() + [20, 10, 5]
    hydroxyl_h = next(a.GetIdx() for a in source.GetAtomWithIdx(2).GetNeighbors()
                      if a.GetAtomicNum() == 1)
    lines = [atom(i + 1, source.GetAtomWithIdx(i).GetSymbol(), pos)
             for i, pos in enumerate(xyz) if i < 3 or i == hydroxyl_h]
    result = posed_molecule(sdf, {"1": 0, "2": 1, "3": 2}, lines, 0, np.eye(4))
    restored_h = next(a.GetIdx() for a in result.GetAtomWithIdx(2).GetNeighbors()
                      if a.GetAtomicNum() == 1)
    assert np.allclose(result.GetConformer().GetPositions()[restored_h], np.round(xyz[hydroxyl_h], 3))
