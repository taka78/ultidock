"""Build a native software-test fixture from separately downloaded, pinned sources.

No production receptor or docking result is inferred by this helper. The benzene
pose is synthetic; the bilayer and oriented membrane protein are published data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from systems import itp_atoms, read_gro, write_gro

ARCHIVE_MD5 = '45909a83bbed9a0659d039a6633d918a'
OPM_SHA256 = '4dabfbfbfe8ff42eecfa487c77900cf440f96ffbf9db25ead47409736dc59c65'


def prepare(archive: Path, opm: Path, output: Path):
    if hashlib.md5(archive.read_bytes()).hexdigest() != ARCHIVE_MD5:
        raise ValueError('Unexpected Slipids archive checksum')
    if hashlib.sha256(opm.read_bytes()).hexdigest() != OPM_SHA256:
        raise ValueError('Unexpected OPM 1AFO checksum; review upstream changes first')
    output.mkdir(parents=True, exist_ok=False)
    params = output / 'lipid_parameters'
    params.mkdir()
    with tarfile.open(archive) as bundle:
        def read(name):
            return bundle.extractfile('SLipids_2016/' + name).read().decode()
        nonbonded = read('SLipids_FF/ffnonbonded.itp')
        bonded = read('SLipids_FF/ffbonded.itp')
        topology = read('itp_files/POPC.itp')
        (output / 'source-bilayer.gro').write_text(read('boxes/POPC_303K.gro'))
        (output / 'source-README.txt').write_text(read('README.txt'))
    # Namespace lipid atom types so no lipid declaration can override a protein,
    # GAFF2 ligand, solvent or ion parameter. Preserve all numerical parameters.
    section = ''
    types = set()
    for line in nonbonded.splitlines():
        clean = line.split(';')[0].strip()
        if clean.startswith('['):
            section = clean.strip('[] ').lower()
        elif clean and section == 'atomtypes':
            types.add(clean.split()[0])
    arity = {'atomtypes': 1, 'bondtypes': 2, 'constrainttypes': 2, 'angletypes': 3,
             'dihedraltypes': 4, 'pairtypes': 2, 'nonbond_params': 2}
    used_types = {row[1] for row in itp_atoms(topology)[1]}
    def rename(text):
        section, lines = '', []
        for raw in text.splitlines():
            clean, _, comment = raw.partition(';')
            if clean.strip().startswith('['):
                section = clean.strip().strip('[] ').lower()
            elif clean.strip() and not clean.lstrip().startswith('#'):
                fields = clean.split()
                if section in arity and any(t not in used_types | {'X'}
                                            for t in fields[:arity[section]]):
                    # The full distribution includes other lipid species and
                    # duplicate parameters for their unused atom types. Keep
                    # every applicable POPC term, including multi-term torsions.
                    continue
                indexes = [1] if section == 'atoms' else range(arity.get(section, 0))
                for index in indexes:
                    if fields[index] in types:
                        fields[index] = 'S16_' + fields[index]
                clean = ' '.join(fields)
            lines.append(clean + (' ;' + comment if comment else ''))
        return '\n'.join(lines) + '\n'
    (params / 'lipid_globals.itp').write_text(rename(nonbonded) + '\n' + rename(bonded))
    (params / 'POPC.itp').write_text(rename(topology))
    name, rows = itp_atoms(topology)
    atoms, box = read_gro(output / 'source-bilayer.gro')
    lipids = [atom for atom in atoms if atom.resname == 'POPC']
    assert name == 'POPC' and len(rows) == 134 and len(lipids) == 128 * 134
    for offset in range(0, len(lipids), 134):
        assert [a.name for a in lipids[offset:offset+134]] == [row[4] for row in rows]
    phosphorus = np.array([a.xyz[2] for a in lipids if a.name == 'P'])
    midpoint = (phosphorus[phosphorus < box[2, 2] / 2].mean() +
                phosphorus[phosphorus > box[2, 2] / 2].mean()) / 2
    # Preserve the equilibrated XY lattice; extend only the aqueous Z space.
    center_z = 5.5
    for atom in lipids:
        atom.xyz[2] += center_z - midpoint
    box[2, 2] = 11.0
    write_gro(output / 'bilayer.gro', lipids, box)
    lines = [line for line in opm.read_text().splitlines() if line.startswith('ATOM  ')
             and line[76:78].strip() != 'H']
    assert {line[21] for line in lines} == {'A', 'B'} and len(lines) == 570
    # The OPM orientation already puts the membrane normal on Z and its center
    # at zero. Only XY centering and a common translation are applied.
    xyz = np.array([[float(line[i:i+8]) for i in (30, 38, 46)] for line in lines])
    translation = np.r_[box.diagonal()[:2] * 5 - (xyz.min(0)[:2] + xyz.max(0)[:2]) / 2,
                        center_z * 10]
    (output / 'protein.pdb').write_text('\n'.join(lines) + '\nEND\n')
    (output / 'protein.pdbqt').write_text('\n'.join(
        line[:66].ljust(70) + ' 0.000 ' + line[76:78].strip() for line in lines) + '\n')
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1'))
    assert AllChem.EmbedMolecule(mol, randomSeed=2026) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, str(output / 'benzene.sdf'))
    # A software probe outside the protein envelope, at the membrane center.
    # This is not an experimentally established glycophorin binding site.
    conformer = mol.GetConformer()
    positions = conformer.GetPositions()
    positions += [-14., -12., 0.]
    pose = []
    for i in range(6):
        x, y, z = positions[i]
        pose.append(f'ATOM  {i+1:5d}  C{i+1:<2d} BEN A   1    {x:8.3f}{y:8.3f}{z:8.3f}'
                    '  1.00  0.00     0.000 A')
    docking = output / 'docking'
    docking.mkdir()
    (docking / 'protein__S1__benzene_12345678.pdbqt').write_text(
        'MODEL 1\nREMARK VINA RESULT: -1 0 0\n' + '\n'.join(pose) + '\nENDMDL\n')
    transform = np.eye(4)
    transform[:3, 3] = translation
    config = {'system_type': 'membrane', 'receptor_reviewed': True, 'receptor_id': 'protein',
              'engine': 'vina', 'receptor_pdb': 'protein.pdb', 'docking_receptor_pdbqt': 'protein.pdbqt',
              'force_field': 'amber99sb-ildn', 'water_model': 'tip3p', 'box_shape': 'triclinic',
              'padding_nm': 1.4, 'salt_molar': .15, 'salt_basis': 'water-count',
              'temperature_k': 303, 'production_ns': .002, 'nvt_ps': 2, 'npt_ps': 2,
              'dt_ps': .002, 'seed': 2026, 'threads': 2,
              'ligands': {'benzene': {'sdf': 'benzene.sdf', 'net_charge': 0,
                                    'atom_map': {str(i+1): i for i in range(6)}}},
              'membrane': {'gro': 'bilayer.gro', 'topology_dir': 'lipid_parameters',
                           'compatible_force_field': 'amber99sb-ildn',
                           'include_files': ['lipid_globals.itp', 'POPC.itp'],
                           'molecules': [{'name': 'POPC', 'count': 128, 'atom_count': 134}],
                           'orientation_reviewed': True, 'transform': transform.tolist(),
                           'hydrophobic_z_nm': [center_z-1.595, center_z+1.595], 'water_pores': [],
                           'clash_distance_nm': .25,
                           'mdp_nonbonded': {'rvdw': 1.4, 'rcoulomb': 1.4,
                                            'vdw-modifier': 'Potential-shift', 'DispCorr': 'EnerPres'}}}
    (output / 'protocol.json').write_text(json.dumps(config, indent=2) + '\n')
    review = {'scope': 'Synthetic ligand / OPM glycophorin A / published POPC software smoke test only',
              'sources': {'slipids': 'https://zenodo.org/records/1149623',
                          'protein': 'https://opm-assets.storage.googleapis.com/pdb/1afo.pdb'},
              'archive_md5': ARCHIVE_MD5, 'opm_sha256': hashlib.sha256(opm.read_bytes()).hexdigest(),
              'protein': 'OPM 1AFO chains A/B, residues 66-101, standard residues; dummy membrane atoms removed; charged fragment termini, histidine states chosen by pdb2gmx and recorded in its log',
              'orientation': 'Retain OPM Z orientation and midplane; common translation only',
              'bilayer': 'POPC 303 K, 128 molecules, 134 atoms each; original XY lattice; lipid-only; Z solvent space extended to 11 nm',
              'parameter_changes': 'Atom-type names prefixed S16_; unrelated lipid types/terms excluded; all applicable POPC terms and numerical parameters retained',
              'assumptions': 'No claim of a physiological glycophorin-benzene binding pose, equilibrated inserted complex, or SERT validation',
              'files': {p.relative_to(output).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in output.rglob('*') if p.is_file()}}
    (output / 'review.json').write_text(json.dumps(review, indent=2) + '\n')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--slipids-archive', type=Path, required=True)
    parser.add_argument('--opm-pdb', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.slipids_archive, args.opm_pdb, args.output)
