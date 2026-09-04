# Configuration file for the Sphinx documentation builder.
#
# The extension set is chosen so that pages can actually move between this
# repository, ERF and REMORA. That is not automatic: ERF's build and
# getting-started pages use sphinx-design directives (`dropdown`, `tab-set`)
# and would fail here as "Unknown directive type" without sphinx_design, which
# is a hard error under the -W this tree is verified with.
#
# myst_parser is loaded because this repository's existing documentation is
# Markdown (README.md, INSTALL.md, Exec/README.md, Docs/*.md), so those pages
# can be included directly rather than converted to reStructuredText.
#
# Full reference: https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = 'ERF-REMORA Driver'
copyright = '2026, ERF-REMORA driver development team'
author = 'ERF-REMORA driver development team'

# Hardcoded because the repository has no tags to derive from. Castro and AMReX
# both compute this (get_version() / get_amrex_version()); switch to a
# git-describe once releases exist.
version = '0.1'
release = '0.1'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'myst_parser',
    'sphinx_design',
    'sphinx_copybutton',
    # Sphinx 6 dropped bundled jQuery; the RTD theme's search box and nav
    # collapse break silently without it. ERF loads it for the same reason.
    'sphinxcontrib.jquery',
    'sphinx_rtd_theme',
]

root_doc = 'index'
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
language = 'en'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'

numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}

# -- Cross-project references ------------------------------------------------

# Link into ERF and REMORA by object rather than by hardcoded URL, so anything
# already documented upstream (erf.use_coupled_sst, REMORA's bulk-flux source
# selectors) is referenced rather than duplicated here. Pattern taken from
# AMReX, which maps into amrex_tutorials and amrex-hydro the same way.
intersphinx_mapping = {
    'erf': ('https://erf.readthedocs.io/en/latest/', None),
    'remora': ('https://remora-ocean.readthedocs.io/en/latest/', None),
}
# Do not let a slow or unreachable upstream inventory stall or fail the build.
intersphinx_timeout = 10
intersphinx_disabled_reftypes = ['*']

# -- Link checking -----------------------------------------------------------

# myroms.org rejects linkcheck's requests; REMORA ignores it for the same
# reason. The driver docs cite ROMS when comparing remap behaviour to
# precedent, so the entry is needed here too.
linkcheck_ignore = [
    r'https://www\.myroms\.org/.*',
]

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
# The RTD theme clips wide tables; the driver.* parameter reference is one.
# Both AMReX and Castro ship an override for exactly this.
html_css_files = ['theme_overrides.css']
htmlhelp_basename = 'ERFREMORAdriverdoc'

# -- Options for LaTeX output ------------------------------------------------

latex_documents = [
    (root_doc, 'ERFREMORAdriver.tex', 'ERF-REMORA Driver Documentation',
     author, 'manual'),
]
