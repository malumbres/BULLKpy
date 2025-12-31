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
    "myst_nb",  # this provides MyST markdown support
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
autosummary_imported_members = True

# --- autodoc ---
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": False,
}
autodoc_member_order = "bysource"
autodoc_typehints = "none"  
autodoc_typehints_format = "short"
python_use_unqualified_type_names = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_rtype = False

# --- theme (Scanpy-like) ---
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_logo = None
html_favicon = None

html_theme_options = {
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navigation_depth": 4,
    "navigation_with_keys": True,
    "secondary_sidebar_items": ["page-toc"],
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
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
}

# --- MyST parser file types ---
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
}

# Keep notebook UX clean
html_show_sourcelink = True
html_show_sphinx = False