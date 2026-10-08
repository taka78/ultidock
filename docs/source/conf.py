"""Sphinx configuration adapted from readthedocs/tutorial-template.

Build the documentation without importing the scientific runtime.
"""

from pathlib import Path
import re

# -- Project information

project = "Ultidock"
author = "Taha Turgut"
copyright = "2026, Taha Turgut"
release = re.search(r'^version = "([^"]+)"',
                    (Path(__file__).resolve().parents[2] / "pyproject.toml").read_text(), re.M).group(1)
version = release
release = f"{version} (gmx-dev)"

# -- General configuration

extensions = ["sphinx.ext.duration", "myst_parser"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
exclude_patterns = []
myst_enable_extensions = ["colon_fence", "deflist"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"
html_title = "Ultidock Documentation — gmx-dev"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": True}
html_show_sourcelink = True

# -- Options for EPUB output

epub_show_urls = "footnote"
