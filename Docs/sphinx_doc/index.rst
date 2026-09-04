.. ERF-REMORA Driver documentation master file.

ERF-REMORA Driver
=================

This repository builds a coupled driver that advances ERF (atmosphere) and
REMORA (ocean) in one executable. The driver owns the coupling clock, the
interface caches, and the pack/apply calls used to exchange state between the
two solvers.

ERF and REMORA are documented in their own repositories. This documentation
covers only what the driver itself owns: the coupling contract, the
``driver.*`` runtime controls, the superbuild, the exchange diagnostics, and
the coupled test suite.

* ERF documentation: https://erf.readthedocs.io/en/latest/
* REMORA documentation: https://remora-ocean.readthedocs.io/en/latest/

.. toctree::
   :maxdepth: 1
   :caption: Contents

Indices and tables
==================

* :ref:`search`
