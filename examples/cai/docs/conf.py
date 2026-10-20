"""Sphinx configuration for the CompDam CAI model generator.

Build with ``uv run --group docs sphinx-build -b html docs docs/_build/html``.
The package is pure Python, so no Abaqus installation is needed.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "CompDam CAI model generator"
author = "CompDam_DGD contributors"
release = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]
napoleon_numpy_docstring = True
napoleon_google_docstring = False
autodoc_member_order = "bysource"
autodoc_mock_imports = []
html_theme = "alabaster"
