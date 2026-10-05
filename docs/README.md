# Building and hosting the documentation

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
│   ├── scientific-background/
│   ├── tutorials/
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
the nested binding-site and docking-engine indexes. Generated HTML should not be
committed. `make -C docs clean` removes generated documentation from `docs/build/`.

## Read the Docs setup

1. Push the documentation changes, including `.readthedocs.yaml`, to your GitHub repository.
2. Sign in to Read the Docs with GitHub and import the existing Ultidock repository.
   The template is already integrated here; a separate tutorial repository is unnecessary.
3. Select the branch containing these files. Read the Docs reads
   `docs/source/conf.py` and installs `docs/requirements.txt` using the root YAML configuration.
4. Inspect the first build's logs, then use **View docs** after it succeeds.
5. In project settings, review pull-request previews and build-failure notifications.
   Activate release versions when their branches or tags contain the documentation setup.

The build fails on Sphinx warnings, matching the tutorial's validation guidance.
HTML hosting is configured; importing the repository into your Read the Docs
account is still required. PDF/EPUB hosting is optional and is not enabled here.
