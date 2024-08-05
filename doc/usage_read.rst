Reading Catalogs with skyCatalogs
=================================
Typically making use of skyCatalogs involves three steps:

1. "Open" the catalog. This is accomplished by reading a *config*: one or more
   yaml files which describes features of the catalog.  Config files will
   be described in detail below.
2. Select objects of specified types (e.g. galaxies or stars or both) within
   a region of the sky.
3. Obtain properties of the objects for use by your application.

The last two steps are often executed multiple times.

Opening a Catalog
-----------------
You open a catalog like this:

.. code-block:: python

   from skycatalogs.skyCatalogs import open_catalog

   skycatalog_root = "some_path/data_folder"
   config_file = "another_path/skyCatalog.yaml"

   cat = open_catalog(config_file, skycatalog_root=skycatalog_root)

Often the config file is in the same folder as the data but it need not be.

Structure of skyCatalogs Configs
================================
Configs contain a variety of information.   Some of it is required for the
skyCatalogs API to find and make sense of the binary data.  Other quantities
are primarily for human consumption. An up-to-date config (version 1.3.0
or later) is a mapping with only a few keys, currently ``skycatalog_root``,
``catalog_dir``, ``catalog_name``, ``schema_version`` and ``object_types``. The
first three are used to find the data. ``schema_version`` refers to the
version of the config structure.  The value of ``object_types``
is a mapping with a key for each object type having data in the catalog.
The bulk of the information may be found
under those keys. That information may be included literally in the
config file but it is preferable to isolate information pertaining to each
object type in its own yaml file and use the yaml **!include**
directive to incorporate it.

Using !include
--------------

.. code-block:: yaml

  catalog_dir:        .               # See Note 1.
  catalog_name:       skyCatalog      # default name for top-level yaml file
  object_types:
    galaxy: !include galaxy.yaml      # See Note 2.
    star: !include star.yaml
  skycatalog_root:    /some/path/often/absolute   # see Note 3.
  schema_version:     1.3.0

.. note::
   #. ``catalog_dir`` indicates where yaml and data files are to be found
      relative to skycatalog_root.  See also Note 3.
   #. The **!include** directive is not supported by yaml natively, but one can
      add an implementation.  This has been done for the skyCatalogs API.
   #. When the top-level yaml file is written, ``skycatalog_root`` is normally
      an absolute path. That file is in the ``catalog_dir`` subdirectory of
      ``skycatalog_root``. But, in case the catalog has been copied,
      ``skycatalog_root`` may be overridden at run time when invoking
      ``open_catalog``. ``catalog_dir`` is then interpreted to be relative to
      the new value of ``skycatalog_root``.

Description of an object type
-----------------------------
For the most part the configuration file fragment for an object type will
either be written by the code generating the binary data for that object
type or will be provided some other way by whoever provides the data.  It
is not normally something users of the skyCatalogs need be concerned with.
However, for those who are curious, such fragments usually include several
kinds of information:

* information allowing the API to find the data files for the object type
* provenance, which in turn may consist of some or all of the following:
  
  * inputs used in the creation of the data files
    
  * run options supplied to the creating program
    
  * schema version
    
  * version of the program and information about the state of the git repo containing the program source
    
* other information associated with the object type


Here is an example file for the ``star`` object type:
  
.. code-block:: yaml
   
  area_partition:
    nside: 32
    ordering: ring
    type: healpix
  data_file_type: parquet
  file_template: pointsource_(?P<healpix>\d+).parquet
  flux_file_template: pointsource_flux_(?P<healpix>\d+).parquet
  provenance:
    inputs:
      star_truth: /global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db
    run_options:
      catalog_dir: out_config
      catalog_name: skyCatalog
      config_path: null
      dc2: false
      flux_parallel: 16
      galaxy_magnitude_cut: 29.0
      galaxy_nside: 32
      galaxy_stride: 1000000
      galaxy_truth: null
      galaxy_type: cosmodc2
      include_roman_flux: false
      knots_magnitude_cut: 27.0
      log_level: DEBUG
      no_flux: true
      no_galaxies: true
      no_knots: false
      no_main: false
      no_pointsources: false
      options_file: local/out_config/star_main.yaml
      pixels:
      - 9556
      sed_subdir: galaxyTopHatSED
      skip_done: true
      skycatalog_root: null
      sso: false
      sso_sed: null
      sso_truth: null
      star_input_fmt: sqlite
    skyCatalogs_repo:
      git_branch: u/jrbogart/config_reorg
      git_hash: 6da4f9636cc63010480c1a1c086cbde8f6ca4dd4
      git_status:
      - UNCOMMITTED_FILES
      - UNTRACKED_FILES
    versioning:
      code_version: 1.7.0-rc4
      schema_version: 1.3.0
  sed_file_root_env_var: SIMS_SED_LIBRARY_DIR
  sed_model: file_nm
  file_nm:
    units: nm
  internal_extinction: None
                
