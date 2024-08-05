# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------


# Add top-level directory to path so python modules can be referred to
# using the normal relative path
import os
import sys
import importlib.util
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'skyCatalogs'
copyright = '2020-2024, LSSTDESC'
author = 'LSSTDESC'

# The full version, including alpha/beta/rc tags
# release = '1.7.0rc2'
# Use load_skycatalogs_version() to determine release dynamically

def load_skycatalogs_version():
    """Extract version of skyCatalogs without importing the whole module"""

    spec = importlib.util.spec_from_file_location(
        "skycatalogs_version",
        os.path.join(os.path.dirname(__file__), "..", "skycatalogs",
                     "_version.py"),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# -- General configuration ---------------------------------------------------

pygments_style = 'sphinx'
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

skycatalogs_version = load_skycatalogs_version()
release = skycatalogs_version.__version__

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Normally want this theme, but may need to comment out if module not available
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The sphinx_rtd_theme does not properly handle wrapping long lines in
# table cells when rendering to HTML due to a CSS issue (see
# https://github.com/readthedocs/sphinx_rtd_theme/issues/1505).  Until
# the issue is fixed upstream in sphinx_rtd_theme, we can simply
# override the CSS here.
rst_prolog = """
.. raw:: html

   <style>
   .wy-table-responsive table td,.wy-table-responsive table th{white-space:normal}
   </style>
"""
