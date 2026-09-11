.. ERF-REMORA Driver documentation master file.

ERF-REMORA Driver
=================

This repository builds a coupled driver that advances ERF (atmosphere) and
REMORA (ocean) in one executable. The driver owns the coupling clock, the
per-variable exchange **lanes** and the caches behind them, and the pack/apply
calls that move state between the two solvers.

ERF and REMORA are documented in their own repositories. This documentation
covers only what the driver itself owns: the coupling contract, the
``driver.*`` runtime controls, the superbuild, the exchange diagnostics, and
the coupled test suite. Anything already documented upstream is referenced
there rather than restated here.

* ERF documentation: https://erf.readthedocs.io/en/latest/
* REMORA documentation: https://remora-ocean.readthedocs.io/en/latest/

.. note::

   Page migration is not complete. ``README.md``, ``INSTALL.md`` and
   ``Exec/README.md`` are not in the toctree below because they live at the
   repository root, outside this source tree, and a toctree entry that
   points above the source directory is unreadable and fails the ``-W``
   build. MyST is configured, so they can be added as-is once they are
   reachable from here.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   Glossary

Indices and tables
==================

* :ref:`search`
