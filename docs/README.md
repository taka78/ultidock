# Building and hosting the documentation

These sources describe `gmx-dev`, including optional GROMACS molecular dynamics.
The docking documentation base was brought from `release/v1.1.2` and reconciled
with this branch's actual examples, receptor handling and result formats.

This documentation follows the official [Read the Docs tutorial](https://docs.readthedocs.com/platform/stable/tutorial/)
and [tutorial template](https://github.com/readthedocs/tutorial-template), adapted
for Ultidock. The template's Makefile and Windows build script come from commit
`f594c40889a4d9949339bd4e0d64275fc99c22f5`.

## Layout

```text
.readthedocs.yaml          # Hosted build configuration (repository root)
docs/
├── Makefile              # Linux/macOS Sphinx build commands
├── make.bat              # Windows Sphinx build commands
├── requirements.txt      # Documentation dependencies
├── source/
│   ├── conf.py           # Sphinx configuration
│   ├── index.rst         # Site entry point and navigation
│   ├── getting-started/
│   ├── user-guide/
│   │   └── molecular-dynamics/
│   ├── scientific-background/
│   ├── tutorials/
│   ├── benchmarks/
│   ├── reference/
│   └── development/
└── build/                # Generated output, ignored by Git
```

The entry point uses reStructuredText, as in the template. The existing guide
pages retain Markdown through MyST, including their internal links and nested
navigation. The site uses the template's `sphinx_rtd_theme`. Documentation
dependencies remain pinned to the tested versions in `requirements.txt`.

The template's fictional Lumache package and API examples are replaced by
Ultidock's guides. Autodoc and package installation are unnecessary for these
pages: the documentation build does not import Ultidock or need native engines.
The maintainer notes `release_checklist.md` and `why_ultidock.md` remain outside
`source/`, so they are not included in the site.

## Local preview

From the repository root, in your activated Python environment:

```bash
python -m pip install -r docs/requirements.txt
make -C docs html SPHINXOPTS="-n -W --keep-going"
```

Open `docs/build/html/index.html`. To build without Make, use:

```bash
python -m sphinx -b html -n -W --keep-going docs/source docs/build/html
```

On Windows, from Command Prompt:

```bat
set SPHINXOPTS=-n -W --keep-going
docs\make.bat html
```

Edit pages under `docs/source/`; add navigation entries in `source/index.rst` or
the nested binding-site, docking-engine and molecular-dynamics indexes. Generated HTML should not be
committed. `make -C docs clean` removes generated documentation from `docs/build/`.

## Where root-guide topics live

| README.md / SETUP.md topic | Published documentation |
| --- | --- |
| One-command docking workflow | `user-guide/start-docking.md` |
| Hardware, OS, native packages, CUDA/OpenCL and Rusticl | `getting-started/requirements.md` |
| Python install choices, managed workspace and optional tools | `getting-started/installation.md` |
| Setup wizard, local/downloaded inputs, full run and artifacts | `getting-started/setup-and-run.md` |
| MD tools, AmberTools build and optional Conda environment | `getting-started/md-installation.md` |
| Docking-to-MD, soluble and membrane inputs, production and restart | `user-guide/molecular-dynamics/` |
| MD commands, protocol fields and failures | `reference/md-cli.md`, `reference/md-protocol.md`, `reference/md-troubleshooting.md` |
| MD assumptions, regression/native tests and acceptance evidence | `scientific-background/molecular-dynamics.md`, `development/md-validation.md` |
| CLI commands, MolGuard and script entry points | `reference/cli.md` |
| Generated config, flags, grid dials and concurrency | `reference/configuration.md` |
| Receptor and ligand handling | `user-guide/preparing-receptor.md`, `user-guide/preparing-ligands.md` |
| Site discovery and CaV-EMPS algorithm | `user-guide/binding-site-discovery/`, `scientific-background/cav-emps-methodology.md` |
| User-guide examples and high-throughput screens | `tutorials/examples.md`, `tutorials/virtual-screening.md` (both in User Guide navigation) |
| Benchmark commands and outputs | `benchmarks/index.md` |
| SQLite, poses, reports and visualization | `user-guide/results-reports.md` |
| Failures, repository layout, tests, citation and license | `reference/troubleshooting.md`, `development/` |

## MD source material

The MD pages use `README.md`, `SETUP.md`, `md-simulation/README.md`,
`md-simulation/ACCEPTANCE.md`, `md-simulation/MEMBRANE_TESTING.md`, example notes
and the executable workflow as their sources. Soluble and membrane guides
include/download the repository JSON templates directly, so the build needs
those files but no MD dependency installation. Acceptance reports are historical
software checks; the documentation build does not execute simulations.

## Read the Docs setup

1. Push the documentation changes, including `.readthedocs.yaml`, to your GitHub repository.
2. Sign in to Read the Docs with GitHub and import the existing Ultidock repository.
   The template is already integrated here; a separate tutorial repository is unnecessary.
3. Select `gmx-dev` for these branch-specific pages. Read the Docs reads
   `docs/source/conf.py` and installs `docs/requirements.txt` using the root YAML configuration.
4. Inspect the first build's logs, then use **View docs** after it succeeds.
5. In project settings, review pull-request previews and build-failure notifications.
   Activate release versions when their branches or tags contain the documentation setup.

The build fails on Sphinx warnings, matching the tutorial's validation guidance.
HTML hosting is configured; importing the repository into your Read the Docs
account is still required. PDF/EPUB hosting is optional and is not enabled here.
