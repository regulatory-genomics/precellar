# -- Path setup --------------------------------------------------------------

import re
import sys

import precellar
from precellar import utils
sys.modules['precellar.utils'] = utils

# -- Software version --------------------------------------------------------

# The short X.Y version (including .devXXXX, -devXXXX, rcX, b1 suffixes if present)
version = re.sub(r'(\d+\.\d+)\.\d+([-\.].*)?', r'\1\2', precellar.__version__)
version = re.sub(r'([-\.]dev\d+).*?$', r'\1', version)

# The full version, including alpha/beta/rc tags.
release = precellar.__version__

# pyData/Sphinx-Theme version switcher
if "dev" in version:
    switcher_version = "dev"
else:
    switcher_version = f"{version}"

print(f'Building documentation for precellar {release} (short version: {version}, switcher version: {switcher_version})')

# -- Project information -----------------------------------------------------

project = 'precellar'
copyright = '2024-2024, Regulatory Genomics Lab, Westlake University'
author = 'Kai Zhang'

# -- General configuration ---------------------------------------------------

suppress_warnings = ['ref.citation']
default_role = 'code'
add_function_parentheses = False

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "nbsphinx",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_plotly_directive",
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

myst_enable_extensions = [
    "amsmath",
    #"colon_fence",
    #"deflist",
    "dollarmath",
    #"fieldlist",
    #"html_admonition",
    #"html_image",
    #"linkify",
    #"replacements",
    #"smartquotes",
    #"strikethrough",
    #"substitution",
    #"tasklist",
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
}

smv_branch_whitelist = r'main'  # Include all branches

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_show_sphinx = False
html_show_sourcelink = False
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]

html_theme_options = {
    "logo": {
        "text": "precellar",
        "image_dark": "_static/logo-dark.svg",
        "alt_text": "precellar",
    },

    "github_url": "https://github.com/regulatory-genomics/precellar",
    "external_links": [
    ],
    "header_links_before_dropdown": 6,

    "navbar_center": ["version-switcher", "navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_align": "left",
    "show_version_warning_banner": switcher_version == "dev",

    "switcher": {
        "version_match": switcher_version,
        "json_url": "https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/docs/_static/versions.json",
    },
}