## Version 2.1.1
- Minimize LSST stack dependency further to simulate Gaia objects without it
- Switch to `uv pip install` in CI
- Correct python requirement to be `>=3.9`

## Version 2.1.0
- Allow for non-imaging Roman bandpasses (i.e., grism and prism)
- Bugfix for PolygonalRegion.compute_mask
- Minimize LSST stack dependencies

## Version 2.0.1
- Bugfix for PolygonalRegion.compute_mask
- Pass object type for external catalog registration and loading

## Version 2.0.0

- Bugfix for get_config_value
- Bugfix for units
- Treat SSOs with zero length streaks as point sources
- Use thinned solar SEDs for SSOs.
- Refactor region handling
- Including provenance information when creating parquet catalog files
- Refactor of config file creation. Use internal code for yaml !INCLUDE directive handling
- Publishing documentation to https://lsstdesc.org/skyCatalogs/

## Version 1.7.0rc4

- Support for OR4 fields
- Use Galsim to normalize fluxes
- Allow Gaia data_dir to be an absolute path

## Version 1.7.0rc3

- Support alternate access mode (direct, not via Butler) to Gaia objects
- Use logger rather than try..catch for abnormal but not fatal conditions
- Support yaml !include directive in skyCatalogs config files
- Bug fixes to CI
- Improvements for SSO object type including option to simulate streaks
- Support for OR3

## Version 1.7.0rc2

- Prefilter PolygonalRegion using Disk rather than cuts on RA, Dec
- Support for diffsky galaxy object type
- Support for snana (SNe) object type
- Support parquet (as well as sqlite) input for star truth
- Remove redshift property from BaseObject

## Version 1.7.0rc1

- Refinements in Gaia stars support
- Allow galaxy output to be partitioned with healpix nside > 32
- Thin bandpasses
- Reorganization to aid pip install

## Version 1.6.0rc2

- Reorganize file structure, eliminating "desc" level

## Version 1.5.0rc2

- Remove reliance on SKYCATALOGS_DIR environment variable

## Version 1.5.0rc1

- Gaia star support
- Add and test; update existing tests

## Version 1.4.0rc5

Nothing substantive; increment version only

## Version 1.4.0rc4

Minor bug fies.

## Version 1.4.0rc1

Initial release.  Produces main and flux catalogs for DC2 galaxies, stars
and SNe
