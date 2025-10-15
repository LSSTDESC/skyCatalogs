Installation Instructions
=========================

.. note::
   **Prerequisites**:

   To use **skyCatalogs>=1.7.0**

   all that is required is a reasonably current version of the  `LSST science pipelines <https://pipelines.lsst.io/>`_ .  See section :ref:`with-pipelines`
   
   For **skyCatalogs>=2.1.1**, having LSST Science Pipelines installed is optional; the only feature that is unavailable without it is the ability to use a Butler to include Gaia objects.

   If installing without the LSST Science Pipelines, go to the section
   :ref:`without-pipelines`

   To access trilegal objects one must install at least one additional package, possibly more if not using the LSST Science Pipelines.  See :ref:`trilegal`.

.. note::

   **Use with imSim**:

   If you intend to use skyCatalogs with imSim, you should follow the `imSim installation instructions <https://lsstdesc.org/imSim/install.html/>`_ . In that case, you need proceed no further with the instructions here; you're all set.  The remainder of these notes is largely extracted from the imSim instructions, omitting anything not required for skyCatalogs.

Installation with LSST science pipelines
------------------------------------------

.. _with-pipelines:

Installing LSST science pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are several methods of installation.  Only the simplest (using a prebuilt cvmfs version) is described here.  For other methods, see the imSim installation instructions.

If you are working at the USDF (Rubin Project computing) or at NERSC (DESC computing), perhaps the easiest way to setup and use *skyCatalogs* is to rely on the prebuilt versions of the pipelines contained in the cvmfs distribution which is installed there.  This solution is also appropriate for personal laptops and university/lab based computing systems if you are able to install the *cvmfs* system.

The `CernVM file system <https://cvmfs.readthedocs.io/>`_  (cvmfs) is a distributed read-only file system developed at CERN for reliable low-maintenance world-wide software distribution.  LSST-France distributes weekly builds of the Rubin science pipelines for both Linux and MacOS.  Details and installation instructions can  be found at `sw.lsst.eu <https://sw.lsst.eu/index.html>`_ .  The distribution includes conda and skyCatalogs dependencies from conda-forge along with the science pipelines.

Load and setup the science pipelines
++++++++++++++++++++++++++++++++++++

First you need to setup the science pipelines.  This involves sourcing a setup file and then using the Rubin *eups* commands to set them up.

.. note::

   Version  ``w_2024_20`` or later of the science pipelines is recommended. This will guarantee other dependencies of skyCatalogs, such as GalSim, are new enough.

   Also note: the cvmfs distribution is a read-only distribution.  This means you cannot add packages to the included conda environment and packages you install via *pip* will be installed in the user area.  If you need a *conda*  environment you will need to use a different installation method.

Source the appropriate setup script (note the -ext in the name) and then setup the distribution (if you are on MacOS use darwin-x86_64 instead of linux-x86_64).

.. code-block:: sh

   source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2024_20/loadLSST-ext.bash
   setup lsst_distrib


Install skyCatalogs itself
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the skyCatalogs package from GitHub.  Here we assume you are in
your installation directory SKYCATALOGS_HOME as described in the section :ref:`per-session` below.

.. code-block:: sh

   git clone https://github.com/LSSTDESC/skyCatalogs

at this point if you would only like to use *skyCatalogs* you can  ``pip install skyCatalog/`` however we instead suggest using the *eups* tool to simply setup the package for use without installing it. This will allow you to edit the package in place, use multiple versions, change branches etc. You should definitely do this if you plan to do any *skyCatalogs* development.

If you do not intend to do any development you may choose instead to clone the most recent release tag.  As of Aug., 2025, this is v.2.1.1

.. code-block:: sh

   git clone https://github.com/LSSTDESC/skyCatalogs.git --branch v.2.1.1

.. _trilegal

Accessing trilegal objects
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to access the trilegal objects you need to install the pystellibs package.  You can do something like this:

.. code-block :: sh

   git clone https://github.com/JoanneBogart/pystellibs.git
   cd pystellibs
   git checkout v1.0.0
   pip install  --user --no-deps --nobuild-isolation -e .
   cd ..

.. _per-session:

Per-session Setup 
~~~~~~~~~~~~~~~~~~

Here is a ``skycatalogs-setup.sh`` file you can use before each session

.. code-block:: sh

   source /cvmfs/...            # as above
   setup lsst_distrib

   export SKYCATALOGS_HOME=*PUT YOUR INSTALL DIRECTORY HERE*

   setup -k -r $SKYCATALOGS_HOME/skyCatalogs

You may want to add the per-session setup needed for the data files
(:ref:`install-data-files`) and, if applicable, for trilegal objects.

If you are using trilegal objects you also need to make sure your installation of pystellibs is accessible.  If not, and assuming PYTHONPATH is not null

.. code-block:: sh

   export PYTHONPATH=${SKYCATALOGS_HOME/pystellibs}:${PYTHONPATH}
   
Now go to section :ref:`install-data-files` below.

.. _without-pipelines:

Install without LSST science pipelines
--------------------------------------

Install skyCatalogs
~~~~~~~~~~~~~~~~~~~

All you need to do is pip install:

.. code-block:: sh
                
   pip install skyCatalogs

Per-session setup
~~~~~~~~~~~~~~~~~

Every session you will also need to define a ``SKYCATALOGS_HOME`` directory
where other needed files (see section :ref:`install-data-files` below) will go:

.. code-block:: sh
                
   export SKYCATALOGS_HOME=*PUT YOUR INSTALL DIRECTORY HERE*

.. _install-data-files:
   
Install needed data files
-------------------------

This step is necessary whether or not you are installing with LSST science pipelines.

Go to your `SKYCATALOGS_HOME` directory and download some needed data files (you will only need to do this once).

.. code-block:: sh

   mkdir -p rubin_sim_data/sims_sed_library
   curl https://s3df.slac.stanford.edu/groups/rubin/static/sim-data/rubin_sim_data/throughputs_2023_09_07.tgz | tar -C rubin_sim_data -xz
   curl https://s3df.slac.stanford.edu/groups/rubin/static/sim-data/sed_library/seds_170124.tar.gz  | tar -C rubin_sim_data/sims_sed_library -xz

The following exports must be done every session, not just when installing the data files.

.. code-block:: sh

   export RUBIN_SIM_DATA_DIR=$SKYCATALOGS_HOME/rubin_sim_data
   export SIMS_SED_LIBRARY_DIR=$SKYCATALOGS_HOME/rubin_sim_data/sims_sed_library

Using skyCatalogs
-----------------

You should now be able to import the code you need from the skyCatalogs package, e.g.

.. code-block:: python

   from skycatalogs.skyCatalogs import open_catalog
   from skycatalogs.utils.shapes import Disk

   skycatalog_root = "path_to/skycatalog_files"  # folder containing catalog
   config_file = "some_folder/skyCatalog.yaml"

   cat = open_catalog(config_file, skycatalog_root=skycatalog_root)

   # define disk at ra, dec = 45.0, -9.0 of radius 100 arcseconds
   disk = disk(45.0, -9.0, 100.0)

   # get galaxies and stars in the region
   objects = cat.get_objects_by_region(disk, obj_type_set={'galaxy', 'star'})
