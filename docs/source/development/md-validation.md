# MD validation and native tests

The branch includes regression tests, opt-in native tests and a recorded
acceptance report. Use regression tests for changes to selection, protocols,
system construction and restart behavior. Native tests additionally require the
[MD toolchain](../getting-started/md-installation.md).

## Regression suite

From the checkout, in an activated environment:

```bash
python -m pip install -e ".[dev,md]" openbabel-wheel
python -m pytest -q md-simulation/tests
```

Ordinary runs skip the opt-in native cases. To focus on membrane construction:

```bash
python -m pytest -q md-simulation/tests/test_membrane.py md-simulation/tests/test_systems.py
```

These cover whole-molecule lipid removal, retained leaflets, periodic pore-water
selection, water grouping, coordinate/topology ordering and box geometry.
Other MD tests cover scored-pose selection, input identity, chemistry preparation,
stage generation, failure isolation and continuation checks. This suite is
separate from the [docking and CLI tests](testing.md).

## Native soluble smoke test

With GROMACS, ACPYPE/AmberTools and Open Babel available:

```bash
ULTIDOCK_MD_NATIVE_TEST=1 python -m pytest -q -s md-simulation/tests/test_workflow.py
```

The opt-in case uses an artificial ethanol/villin complex with a heavy-atom
protein fixture derived from PDB 1VII. It runs the native toolchain through
minimization and 2 ps each of NVT, NPT and production. It checks execution and
topology consistency; it does not validate an experimental binding interaction.

## Native membrane fixture

The opt-in membrane test uses Slipids 2016 POPC coordinates/parameters and an
OPM-oriented glycophorin A dimer (1AFO). A synthetic benzene pose exercises
parameterization and the shared coordinate frame. It is a software fixture;
it is not a bundled SERT MD model or an experimentally supported binding pose.

Download the sources explicitly and prepare a new fixture directory:

```bash
mkdir -p md-simulation/workspace/membrane-sources
curl -fL https://zenodo.org/records/1149623/files/Slipids_2016.tar.gz \
  -o md-simulation/workspace/membrane-sources/Slipids_2016.tar.gz
curl -fL https://opm-assets.storage.googleapis.com/pdb/1afo.pdb \
  -o md-simulation/workspace/membrane-sources/1afo-opm.pdb
python md-simulation/tests/prepare_membrane_fixture.py \
  --slipids-archive md-simulation/workspace/membrane-sources/Slipids_2016.tar.gz \
  --opm-pdb md-simulation/workspace/membrane-sources/1afo-opm.pdb \
  --output md-simulation/workspace/membrane-fixture
ULTIDOCK_MD_NATIVE_TEST=1 \
ULTIDOCK_MD_MEMBRANE_FIXTURE="$PWD/md-simulation/workspace/membrane-fixture" \
  python -m pytest -q -s md-simulation/tests/test_membrane_native.py
```

The helper verifies pinned source checksums and records assumptions, transforms
and hashes in `review.json`. Pytest does not download the structures. A changed
upstream archive requires review before changing its expected checksum.

The native test exercises build, minimization and 2 ps each of NVT, NPT and
production. It checks transformed ligand coordinates, retained whole lipids,
solvent exclusion, ion counts, semi-isotropic coupling, finite energy samples
and absence of LINCS warnings. Source details, force-field compatibility and
fixture transformations are in the
{download}`complete membrane test guide <../../../md-simulation/MEMBRANE_TESTING.md>`.

## Recorded acceptance and its scope

The repository's {download}`acceptance report <../../../md-simulation/ACCEPTANCE.md>`
records a native run on **2026-10-05**, following a clean installation starting
from `gmx-beta` and subsequent fixes. Its results are historical evidence, not a
claim that the documentation build reran those simulations.

| Recorded check | Outcome |
| --- | --- |
| Ordinary SERT docking | Five successful pocket outputs, 45 database poses and site-aware analysis exports |
| Soluble T4 lysozyme L99A / benzene (181L) | Vina docking through preparation, EM, NVT, NPT and production |
| Invalid selected ligand | Reported preparation failure while the valid ligand completed |
| Interrupted production | Continued from a checkpoint with the unchanged TPR using `-cpi -append` |
| Modified immutable inputs | Rejected before preparation or execution |
| Full suite with native smoke enabled | 211 passed, 2 skipped in the recorded environment |
| POPC/glycophorin membrane fixture | Completed native stages with both leaflets retained and no warning bypass |

The recorded installation used Ubuntu 26.04.1, GROMACS 2025.4, ACPYPE 2026.9.4,
Open Babel and Meeko 0.8.0. These are tested versions in that report, rather than
additional minimum requirements. Generated evidence directories are ignored by
Git and may be absent from a fresh clone; use the source guides and tests to
reproduce software checks.

The soluble dynamics lasted 2 ps each for NVT/NPT and 4 ps for production; the
membrane dynamics lasted 2 ps per stage. They demonstrate native execution,
coordinate integrity and checkpoint continuation. Scientific equilibration,
SERT-specific membrane validation and production-scale performance remain
system-dependent. See [MD assumptions](../scientific-background/molecular-dynamics.md).
