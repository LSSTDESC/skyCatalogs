Reading Catalogs with skyCatalogs
=================================
Typically making use of skyCatalogs involves three steps:

1. "Open" the catalog. This is accomplished by reading a *config file*: a
   yaml file which describes features of the catalog.  Config files will
   be described in detail elsewhere.
2. Select objects (e.g. galaxies or stars) within a region of the sky.
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
