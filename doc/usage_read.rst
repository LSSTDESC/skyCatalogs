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
or later) is a mapping with only a few keys, currently `skycatalog_root`,
`catalog_dir`, `catalog_name`, `schema_version` and `object_types`. The
first three are used to find the data (though `skycatalog_root` may be
and often is overridden at run time). `schema_version` refers to the
version of the config structure.  The bulk of the information may be found
under `object_types`.
