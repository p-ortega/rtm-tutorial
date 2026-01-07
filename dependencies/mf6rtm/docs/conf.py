# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from datetime import datetime
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'mf6rtm'
copyright = f'{datetime.now().year}, Pablo Ortega'
# year = datetime.now().year
author = 'Pablo Ortega, Anthony Aufdenkampe and others'
# release = '0.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Sphinx extensions
extensions = [
    "sphinx.ext.autodoc",              # automatic API docs
    "sphinx.ext.napoleon",             # Google/NumPy docstrings
    "sphinx_autodoc_typehints",        # types from function annotations
    "myst_parser",                     # Markdown support
    "sphinx.ext.autosummary",
]

# Auto-generate summary tables for modules
autosummary_generate = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    "navigation_with_keys": True,   # optional: navigate with keyboard
}