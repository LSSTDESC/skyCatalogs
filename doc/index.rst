.. skyCatalogs documentation master file, created by
   sphinx-quickstart on Fri Apr 12 13:47:18 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======================================================
skyCatalogs: provides API for simulated object catalogs
=======================================================
skyCatalogs both creates simulated object catalogs (from precursor catalogs)
and provides a uniform API to the catalogs it makes and to certain other
catalogs. These catalogs may be used as inputs to image simulation (originally
only for the Rubin Observatory; since also used in simulations for the
Roman Observatory) and as truth catalogs.

This documentation covers installation and usage instructions for the package.

The source code of *skyCatalogs* is hosted at
https://github.com/LSSTDESC/skyCatalogs.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage_read
   ops-rehearsals-3-4/galaxies
   ops-rehearsals-3-4/UW_stars
   ops-rehearsals-3-4/SSO
   roman-rubin-1.1.2/diffsky_galaxy
   roman-rubin-1.1.2/UW_stars
   roman-rubin-1.1.2/snana


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
