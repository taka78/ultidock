"""Opt-in native membrane test with pinned, separately downloaded physical inputs."""
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
import workflow
from systems import read_gro


@pytest.mark.skipif(os.environ.get('ULTIDOCK_MD_NATIVE_TEST') != '1' or
                    not os.environ.get('ULTIDOCK_MD_MEMBRANE_FIXTURE'),
                    reason='Native membrane test requires tools and ULTIDOCK_MD_MEMBRANE_FIXTURE')
def test_native_reviewed_popc_seed_through_short_production(tmp_path, monkeypatch):
    fixture = Path(os.environ['ULTIDOCK_MD_MEMBRANE_FIXTURE']).resolve()
    review = json.loads((fixture / 'review.json').read_text())
    for name, checksum in review['files'].items():
        path = (fixture / name).resolve()
        assert path.is_relative_to(fixture)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum
    root = tmp_path / 'md-simulation'
    root.mkdir()
    for source in workflow.ROOT.glob('*.py'):
        shutil.copy2(source, root / source.name)
    monkeypatch.setattr(workflow, 'ROOT', root)
    jobdir = workflow.start(fixture / 'protocol.json', fixture / 'docking', through='npt')
    system = jobdir / '01_benzene'
    job = json.loads((jobdir / 'job.json').read_text())
    composition = json.loads((system / 'membrane_composition.json').read_text())
    row, = composition
    assert 0 < row['retained'] < row['initial'] == 128
    assert row['upper_leaflet'] > 0 and row['lower_leaflet'] > 0
    assert row['upper_leaflet'] + row['lower_leaflet'] == row['retained']
    atoms, box = read_gro(system / 'system.gro')
    assert sum(a.resname == 'POPC' for a in atoms) == 134 * row['retained']
    config = job['config']
    low, high = config['membrane']['hydrophobic_z_nm']
    assert all(not low < (a.xyz[2] % box[2, 2]) < high for a in atoms
               if a.resname == 'SOL' and a.name.startswith('O'))
    solvated, _ = read_gro(system / 'solvated.gro')
    waters = sum(a.resname == 'SOL' for a in solvated) // 3
    salt = json.loads((system / 'salt.json').read_text())
    pairs = round(config['salt_molar'] * waters / 55.5)
    assert salt['basis'] == 'water-count' and salt['salt_pairs'] == pairs
    neutralizing_charge = int(salt['neutralizing_charge'])
    assert sum(a.resname == 'NA' for a in atoms) == pairs + max(neutralizing_charge, 0)
    assert sum(a.resname == 'CL' for a in atoms) == pairs + max(-neutralizing_charge, 0)
    # The source pose is transformed once into the bilayer frame, then ACPYPE
    # atom ordering is restored before assembly. Check the emitted ligand too.
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(system / 'posed_ligand.mol'), removeHs=False)
    source = np.array([[float(line[i:i+8]) for i in (30, 38, 46)]
                       for line in job['systems'][0]['pose']['coordinates']])
    matrix = np.asarray(config['transform'])
    expected = (matrix @ np.c_[source, np.ones(len(source))].T).T[:, :3]
    np.testing.assert_allclose(mol.GetConformer().GetPositions()[:6], expected, atol=1e-4)
    ligand, _ = read_gro(system / 'ligand.gro')
    heavy = np.array([a.xyz for a in ligand if not a.name.startswith('H')])
    # ACPYPE may reorder atoms: compare coordinate sets within GRO precision.
    assert max(min(np.linalg.norm(x - heavy * 10, axis=1)) for x in expected) < .01
    for stage in ('nvt', 'npt'):
        values = np.loadtxt(system / f'{stage}_energy.xvg', comments=['#', '@'])
        assert np.isfinite(values).all()
        assert 'LINCS WARNING' not in (system / f'{stage}.log').read_text()
    assert 'pcoupltype = semiisotropic' in (system / 'npt.mdp').read_text()
    workflow.run(jobdir, 'production', True)
    state = json.loads((system / 'state.json').read_text())
    assert all(state[stage]['status'] == 'complete' for stage in ('build', 'em', 'nvt', 'npt', 'production'))
    assert (system / 'production.cpt').stat().st_size > 0
    assert 'POSRES' not in (system / 'production.mdp').read_text()
    assert 'LINCS WARNING' not in (system / 'production.log').read_text()
    assert not any('-maxwarn' in json.loads(line)['argv'] for line in (system / 'commands.jsonl').read_text().splitlines())
    print(f'Native membrane evidence: {jobdir}', flush=True)
    print(json.dumps({'composition': composition, 'salt': json.loads((system / 'salt.json').read_text())}), flush=True)
