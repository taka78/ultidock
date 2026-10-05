"""Membrane assembly invariants; native coverage lives in test_membrane_native.py."""
import numpy as np
import pytest
from systems import GroAtom, embed_bilayer, read_gro, remove_membrane_core_water, write_gro


@pytest.fixture
def bilayer(tmp_path):
    (tmp_path / 'lipid.itp').write_text('[ moleculetype ]\nLIP 3\n[ atoms ]\n'
                                      '1 c 1 LIP C1 1 0 12\n2 c 1 LIP C2 1 0 12\n')
    atoms = [GroAtom(i, 'LIP', name, np.array([x, 5., z + dz]))
             for i, (x, z) in enumerate([(5., 5.), (8., 4.), (8., 6.)], 1)
             for name, dz in [('C1', 0.), ('C2', .1)]]
    path = tmp_path / 'bilayer.gro'
    write_gro(path, atoms, np.eye(3) * 10)
    config = {'box_shape': 'triclinic', 'padding_nm': 1.2, 'membrane': {
        'gro': str(path), 'topology_dir': str(tmp_path), 'include_files': ['lipid.itp'],
        'molecules': [{'name': 'LIP', 'count': 3, 'atom_count': 2}],
        'hydrophobic_z_nm': [4, 6], 'clash_distance_nm': .2,
        'mdp_nonbonded': {'rvdw': 1.2, 'rcoulomb': 1.2}}}
    solute = [GroAtom(1, 'ALA', 'CA', np.array([5., 5., 5.]))]
    return config, solute, atoms, path


@pytest.mark.parametrize('fault, message', [
    ('order', 'atom order'), ('extra-atoms', 'exactly the specified'),
    ('short-block', 'exceed coordinate'), ('tilted', 'orthorhombic'),
    ('small-box', 'twice the cutoff'), ('padding', 'periodic padding'),
    ('slab', 'hydrophobic_z_nm'), ('one-leaflet', 'Both membrane leaflets'),
])
def test_invalid_membrane_assembly_is_rejected(bilayer, fault, message):
    config, solute, atoms, path = bilayer
    box = np.eye(3) * 10
    if fault == 'order':
        atoms[0].name = 'WRONG'
    elif fault == 'extra-atoms':
        atoms.append(GroAtom(4, 'SOL', 'OW', np.array([1., 1., 1.])))
    elif fault == 'short-block':
        atoms.pop()
    elif fault == 'tilted':
        box[1, 0] = 1
    elif fault == 'small-box':
        box[0, 0] = 2.4
    elif fault == 'padding':
        solute[0].xyz[0] = .1
    elif fault == 'slab':
        config['membrane']['hydrophobic_z_nm'] = [6, 4]
    elif fault == 'one-leaflet':
        for atom in atoms[4:]:
            atom.xyz[2] -= 2
    write_gro(path, atoms, box)
    with pytest.raises(ValueError, match=message):
        embed_bilayer(solute, config)


def test_periodic_lipid_images_are_removed_as_whole_molecules(bilayer):
    config, solute, atoms, path = bilayer
    for atom in atoms[:2]:
        atom.xyz[0] += 10
    write_gro(path, atoms, np.eye(3) * 10)
    kept, _, counts, report = embed_bilayer(solute, config)
    assert counts == [('LIP', 2)]
    assert {a.residue for a in kept} == {2, 3}
    assert report[0]['upper_leaflet'] == report[0]['lower_leaflet'] == 1
    np.testing.assert_array_equal(solute[0].xyz, [5., 5., 5.])


@pytest.mark.parametrize('sites', [3, 4])
def test_water_filter_uses_periodic_pores_and_keeps_whole_waters(tmp_path, sites):
    names = ['OW', 'HW1', 'HW2'] + (['MW'] if sites == 4 else [])
    atoms = [GroAtom(i, 'SOL', name, np.array(xyz, dtype=float))
             for i, xyz in enumerate([[5, 5, 2], [9.9, 5, 5], [5, 5, 15]], 1)
             for name in names]
    path = tmp_path / 'water.gro'
    write_gro(path, atoms, np.eye(3) * 10)
    count = remove_membrane_core_water(path, {'water_sites': sites, 'hydrophobic_z_nm': [4, 6],
        'water_pores': [{'xy_nm': [.1, 5], 'radius_nm': .25}]})
    kept, _ = read_gro(path)
    assert count == 2 and len(kept) == 2 * sites
    assert {a.residue for a in kept} == {1, 2}


@pytest.mark.parametrize('fault', ['residue', 'oxygen', 'duplicate-name'])
def test_water_filter_rejects_broken_molecule_groups(tmp_path, fault):
    atoms = [GroAtom(1, 'SOL', name, np.array([2., 2., 2.])) for name in ['OW', 'HW1', 'HW2']]
    if fault == 'residue':
        atoms[-1].residue = 2
    elif fault == 'oxygen':
        atoms[0].name = 'HW3'
    else:
        atoms[-1].name = 'HW1'
    path = tmp_path / 'water.gro'
    write_gro(path, atoms, np.eye(3) * 10)
    original = path.read_bytes()
    with pytest.raises(ValueError, match='Solvent atom grouping'):
        remove_membrane_core_water(path, {'water_sites': 3, 'hydrophobic_z_nm': [4, 6]})
    assert path.read_bytes() == original
