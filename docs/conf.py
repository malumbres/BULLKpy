from __future__ import annotations

import os
import sys
from datetime import date
from pathlib import Path

# --- make src/ importable ---
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

project = "BULLKpy"
author = "Marcos Malumbres"
copyright = f"{date.today().year}, {author}"

# If bullkpy is importable, try to pull version
try:
    import bullkpy  # noqa: F401
    version = getattr(bullkpy, "__version__", "0.0.0")
    release = version
except Exception:
    version = "0.0.0"
    release = version

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "myst_nb",
    "sphinx_copybutton",
    "sphinx_design",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/.ipynb_checkpoints",
]

# --- MyST / MyST-NB ---
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
    "tasklist",
    "attrs_inline",
]

nb_execution_mode = "off"  # IMPORTANT: fast & stable on RTD; you can set "auto" later
nb_execution_timeout = 300
nb_output_stderr = "remove"

# --- autosummary ---
autosummary_generate = True

# --- autodoc ---
autodoc_typehints = "description"  # move hints to description
autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# --- theme (Scanpy-like) ---
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_logo = None
html_favicon = None

html_theme_options = {
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navigation_with_keys": True,
    "show_toc_level": 2,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/malumbres/BULLKpy",
            "icon": "fa-brands fa-github",
        },
    ],
}

# --- intersphinx (optional but useful) ---
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", {}),
    "numpy": ("https://numpy.org/doc/stable", {}),
    "pandas": ("https://pandas.pydata.org/docs", {}),
    "anndata": ("https://anndata.readthedocs.io/en/stable", {}),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", {}),
}

# --- MyST parser file types ---
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Keep notebook UX clean
html_show_sourcelink = True
html_show_sphinx = False