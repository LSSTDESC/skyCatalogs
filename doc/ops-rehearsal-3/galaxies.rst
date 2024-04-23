+++++++++++++++++++++++++++++++++++++
Galaxy quantities for ops-rehearsal-3
+++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" (with information largely coming from the
cosmodc2 galaxy catalog) and a flux file. The latter contains
for each object only fluxes for lsst bands and the object id so it can be
joined with the main file. Both are parquet files. The main file can be
rather large (over 3 Gbytes) so, depending on your application, you may
want to iterate over row groups rather than reading it all in at once.

Galaxy main file
----------------

========================  ============   ==============  ========================
Name                      Datatype       Units           Description
========================  ============   ==============  ========================
galaxy_id                 int64          N/A             Unique object identifier
ra                        double         degrees         Object right ascension
dec                       double         degrees         Object declination
redshift                  double         N/A             Cosmological redshift
                                                         with line-of-sight motion
redshift_hubble           double         N/A             Cosmological redshift
peculiar_velocity         double         km/sec          Peculiar velocity
shear_1                   double         N/A             Shear (gamma) component 1,
                                                         treecorr/Galsim convention
shear_2                   double         N/A             Shear (gamma) component 2,
                                                         treecorr/Galsim convention
convergence               double         N/A             Convergence (kappa)
size_bulge_true           float          arcsec          Bulge half-light radius
                                                         (major axis) not lensed
size_minor_bulge_true     float          arcsec          Bulge half-light radius
                                                         (minor axis) not lensed
sersic_bulge              float          N/A             Sersic index of bulge
                                                         light profile
size_disk_true            float          arcsec          Disk half-light radius
                                                         (major axis) not lensed
size_disk_bulge_true      float          arcsec          Disk half-light radius
                                                         (minor axis) not lensed
sersic_disk               float          N/A             Sersic index of disk
                                                         light profile
ellipticity_1_disk_true   double         N/A             Ellipticity component 1
                                                         for disk, not lensed
ellipticity_2_disk_true   double         N/A             Ellipticity component 2
                                                         for disk, not lensed
ellipticity_1_bulge_true  double         N/A             Ellipticity component 1
                                                         for bulge, not lensed
ellipticity_2_bulge_true  double         N/A             Ellipticity component 2
                                                         for bulge, not lensed
sed_val_bulge             list(double)   4.4659e13 W/Hz  Integrated rest-frame AB
                                                         luminosities for tophat
                                                         filters (bulge component)
sed_val_disk              list(double)   4.4659e13 W/Hz  Like sed_val_bulge
bulge_magnorm             double         N/A             Normalization for SED
disk_magnorm              double         N/A             Normalization for SED
MW_rv                     float          N/A             Extinction parameter
MW_av                     float          N/A             Extinction parameter
                                                         from F19 dust map
sed_val_knots             list(double)   4.4659e13 W/Hz  Like sed_val_bulge
n_knots                   float          N/A             Number of knots
knots_magnorm             double         N/A             Normalization for SED
========================  ============   ==============  ========================

NOTE: For tophat definitions (start and width in units of angstroms) see yaml
config file associated with the data.

Galaxy flux file
----------------

=============   =========   ================  ================================
Name            Datatype    Units             Description
=============   =========   ================  ================================
galaxy_id       int64       N/A               Unique object identifier
lsst_flux_u     float       photons/sec/cm^2  Flux integrated over lsst u-band
lsst_flux_g     float       photons/sec/cm^2  Flux integrated over lsst g-band
lsst_flux_r     float       photons/sec/cm^2  Flux integrated over lsst r-band
lsst_flux_i     float       photons/sec/cm^2  Flux integrated over lsst i-band
lsst_flux_z     float       photons/sec/cm^2  Flux integrated over lsst z-band
lsst_flux_y     float       photons/sec/cm^2  Flux integrated over lsst y-band
=============   =========   ================  ================================
