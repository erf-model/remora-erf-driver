# Configuration file for the Sphinx documentation builder.
#
# This is deliberately minimal. It mirrors the ERF and REMORA sphinx
# configurations closely enough that pages can move between the three
# repositories, but it carries no options this repository does not use yet.
#
# Full reference: https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = 'ERF-REMORA Driver'
copyright = ' '
author = 'ERF-REMORA driver development team'

# The short X.Y version, and the full version including alpha/beta/rc tags.
version = '0.1'
release = '0.1'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'sphinx_rtd_theme',
]

master_doc = 'index'
source_suffix = '.rst'
language = 'en'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'

numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
htmlhelp_basename = 'ERFREMORAdriverdoc'

# -- Options for LaTeX output ------------------------------------------------

latex_documents = [
    (master_doc, 'ERFREMORAdriver.tex', 'ERF-REMORA Driver Documentation',
     author, 'manual'),
]
