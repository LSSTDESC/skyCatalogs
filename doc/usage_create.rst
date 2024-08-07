Creating New Catalogs
=====================
Catalogs - both the binary data files and yaml config - for some object types
can be created using the script `create_sc.py`. For others the data are
created by other means, but it is still necessary to make a suitable
config fragment so that the skyCatalogs API can access the data.  New object
types will require updates to the API.

The script ``create_sc.py`` and its options
-------------------------------------------
``create_sc.py`` has a plethora of options, but few are of interest for
most invocations and even fewer are required.  The options determine

* object type(s) for which files will be created
* region of the sky for which files will be created. 
* whether to create so-called *main* files (containing basic quantities
  such as ra, dec, object id, ..), flux files, or both
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
catalog_dir            string     "."           Location of catalog relative
                                                to skycatalog_root
catalog_name           string     "skyCatalog"  Name of top-level yaml config
                                                file 
config_path            string     None          where to write config. If
                                                ``None``, same folder as data
dc2                    boolean    False         Use dc2 conventions
flux_parallel          int        16            # of parallel processes to use
                                                when computing flux
galaxy_magnitude_cut   float      29.0          Discard galaxies above cut
galaxy_nside           int        32            nside for healpixels
galaxy_stride          int        1_000_000     Max galaxies output per row
                                                group
galaxy_truth           string     None          Default depends on galaxy_type
galaxy_type            string     "cosmodc2"    May be "cosmodc2" or "diffsky"
include_roman_flux     boolean    False         Output fluxes for Roman as
                                                well as for Rubin
knots_magnitude_cut    float      27.0          Omit knots component from
                                                galaxies with i-mag above cut
log_level              string     "INFO"        Log level
no_flux                boolean    False         Omit flux files
no_galaxies            boolean    False         Omit galaxies from output
no_knots               boolean    False         Omit knot component
no_main                boolean    False         Do not generate main files
no_pointsources        boolean    False         Do not generate star files
options_file           string     None          Path to file where other
                                                options are set
pixels                 int list   [9556]        healpix pixels for which
                                                catalog will be created
sed_subdir             string     N/A           obsolete
skip_done              boolean    False         do not overwrite existing files
skycatalog_root        string     None          Abs. path. See catalog_dir
sso                    boolean    False         Generate sso files
sso_sed                string     None          Path to file to SED to be
                                                used for all SSOs
sso_truth              string     None          Directory containing sso input
star_input_fmt         string     "sqlite"      Format of star truth
=====================  =========  ============  ===============================

Example options files
+++++++++++++++++++++
