Creating New Catalogs
=====================
Catalogs - both the binary data files and yaml config - for some object types
can be created using the script `create_main.py` and
`create_flux.py`.  Creation of the main catalog must precede
creation of the flux catalog for the same sky region.
For other object types the data are
created by other means, but it is still necessary to make a suitable
config fragment so that the skyCatalogs API can access the data.


The script ``create_main.py`` and its options
---------------------------------------------
``create_main.py`` has a plethora of options, but few are of interest
for most invocations and even fewer are required.  The options determine

* object type for which files will be created (required)
* region of the sky for which files will be created.
* where the output goes
* inputs (usually defaulted, depending on object type)
* numerous options conditioning exactly how the output is to be created,
  most of which are usually defaulted

One particularly handy option, ``--options-file``, allows you to specify
everything else as key-value pairs in a yaml file.

Options table
+++++++++++++
Here is the complete list of options as they appear in an options file.
On the command line, prepend ``--`` and change all underscores to hyphens.

=====================  =========  ============  ===============================
Name                   Datatype   Default       Description
=====================  =========  ============  ===============================
object_type            string                   Required. One of {star, sso,
                                                cosmodc2_galaxy, diffsky_galaxy}
catalog_dir            string     "."           Location of catalog relative
                                                to skycatalog_root
                                                (see below)
catalog_name           string     "skyCatalog"  Name of top-level yaml config
                                                file
config_path            string     None          where to write config. If
                                                ``None``, same folder as data
dc2                    boolean    False         Use dc2 conventions
galaxy_magnitude_cut   float      29.0          Discard galaxies above cut.
                                                Ignored for non-galaxy
                                                object types
nside                  int        32            nside for healpixels
stride                 int        1_000_000     Max objects output per row
                                                group
truth                  string     None          Default depends on object_type
knots_magnitude_cut    float      27.0          Omit knots component from
                                                galaxies with i-mag above cut
log_level              string     "INFO"        Log level
no_knots               boolean    False         Omit knot component
options_file           string     None          Path to file where other
                                                options are set
pixels                 int list   [9556]        healpix pixels for which
                                                catalog will be created
skip_done              boolean    False         do not overwrite existing files
skycatalog_root        string     None          Abs. path. See catalog_dir
sso_sed                string     None          Path to file to SED to be
                                                used for all SSOs. Defaults
                                                to `solar_sed_thin.txt`,
                                                included in repo.
star_input_fmt         string     "sqlite"      Format of star truth
=====================  =========  ============  ===============================

The script ``create_flux.py`` and its options
-------------------------------------------
``create_flux.py`` has most of the same options as ``create_main.py``
plus a couple new ones. Only ``object_type`` is required but several
others are specified for most invocations.

Options table
+++++++++++++
Here is the complete list of options as they appear in an options file.
On the command line, prepend ``--`` and change all underscores to hyphens.

=====================  =========  ============  ===============================
Name                   Datatype   Default       Description
=====================  =========  ============  ===============================
object_type            string                   Required. One of {star, sso,
                                                cosmodc2_galaxy, diffsky_galaxy}
catalog_dir            string     "."           Location of catalog relative
                                                to skycatalog_root
                                                (see below)
catalog_name           string     "skyCatalog"  Name of top-level yaml config
                                                file
config_path            string     None          where to write config. If
                                                ``None``, same folder as data
flux_parallel          int        16            # processes to run in parallel
                                                when computing fluxes
iclude_roman_flux      boolean    False         If True calculate & store Roman
                                                as well as Rubin fluxes.
log_level              string     "INFO"        Log level
options_file           string     None          Path to file where other
                                                options are set
pixels                 int list   [9556]        healpix pixels for which
                                                catalog will be created
skip_done              boolean    False         do not overwrite existing files
skycatalog_root        string     None          Abs. path. See catalog_dir
sso_sed                string     None          Path to file to SED to be
                                                used for all SSOs. Defaults
                                                to `solar_sed_thin.txt`,
                                                included in repo.

Example options files
+++++++++++++++++++++
Create cosmodc2-style galaxies main file.  This file was one I wrote
primarily for testing the creation code. In order to speed things up, I
set ``galaxy_magnitude_cut`` down to 20.0. The output file is about 1000 times
smaller than with the default cut.

.. code-block:: yaml

   object_type:          cosmodc2_galaxy
   catalog_dir:          just_testing
   pixels:               [9683]
   galaxy_magnitude_cut: 20.0    # Default is 29.0
   log_level:            DEBUG   # Default is INFO

Suppose the file is called `main_galaxy.yaml`. It can be invoked as follows:

.. code-block:: sh

   python cosmodc2_galaxy --options-file main_galaxy.yaml

.. note::
   Since ``object_type`` is required it must be specified both on the
   command line and in the options file.

Create star flux file for a couple healpixels.
The default value of ``flux_parallel`` is rather conservative for Perlmutter,
so use something higher to make the process go faster.

.. code-block:: yaml

   catalog_dir:          just_testing
   pixels:               [9683, 9684]
   flux_parallel:        24      # Default is 16

.. note::
   The star main files for both healpixels must already exist in the output
   directory since they are input to the flux generation.
