# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'esrf_id10_surf'
copyright = '2025, Egor A. Bersenev'
author = 'Egor A. Bersenev'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
    'nbsphinx',
]

autodoc_member_order = 'groupwise'

autodoc_default_options = {
    'members': True,           # Document members (methods/attributes)
    'member-order': 'groupwise', # Enforce grouping by type
    'undoc-members': True,     # Include members without docstrings
    'show-inheritance': True,  # Show parent classes
}

autodoc_typehints = "description" # Moves type hints into the description instead of signature

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
