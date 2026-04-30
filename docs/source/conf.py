# Configuration file for the Sphinx documentation builder.

import os
import sys
import datetime

# Make radprocess importable for autodoc
sys.path.insert(0, os.path.abspath("../.."))

# -- Project information

project = "radprocess"
copyright = "{0}, Sacha Gavino".format(datetime.datetime.now().year)
author = "Sacha Gavino"
release = "0.1.0"

# -- General configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

# Napoleon settings (numpy-style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

templates_path = ["_templates"]

# -- Options for LaTeX output

latex_elements = {
    "sphinxsetup": "VerbatimColor={rgb}{0.95,0.95,0.95}",
}

# -- Options for HTML output

html_theme = "renku"

html_static_path = ["_static"]

# -- Options for EPUB output
epub_show_urls = "footnote"