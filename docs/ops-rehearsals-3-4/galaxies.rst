++++++++++++++++++++++++++++++++++++++++++++
Galaxy quantities for ops-rehearsals 3 and 4
++++++++++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" (with information largely coming from the
cosmodc2 galaxy catalog) and a flux file. The latter contains
for each object only fluxes for lsst bands and the object id so it can be
joined with the main file. Both are parquet files. The main file can be
rather large (over 3 Gbytes) so, depending on your application, you may
want to iterate over row groups rather than reading it all in at once.

Galaxy main file
----------------

=============================  ============   ==============  ========================
Name                           Datatype       Units           Description
=============================  ============   ==============  ========================
galaxy_id                      int64          N/A             Unique object identifier
ra                             double         degrees         Object right ascension
dec                            double         degrees         Object declination
redshift                       double         N/A             Cosmological redshift
                                                              with line-of-sight motion
redshiftHubble                 double         N/A             Cosmological redshift
peculiarVelocity               double         km/sec          Peculiar velocity
shear1                         double         N/A             Shear (gamma) component 1,
                                                              treecorr/Galsim convention
shear2                         double         N/A             Shear (gamma) component 2,
                                                              treecorr/Galsim convention
convergence                    double         N/A             Convergence (kappa)
spheroidHalfLightRadiusArcsec  float          arcsec          Bulge component half-light
                                                              radius
diskHalfLightRadiusArcsec      float          arcsec          Disk component half-light
                                                              radius
diskEllipticity1               double         N/A             Ellipticity component 1
                                                              for disk, not lensed
diskEllipticity2               double         N/A             Ellipticity component 2
                                                              for disk, not lensed
spheroidEllipticity1           double         N/A             Ellipticity component 1
                                                              for bulge, not lensed
spheroidEllipticity2           double         N/A             Ellipticity component 2
                                                              for bulge, not lensed
um_source_galaxy_obs_sm        float          solar mass      stellar mass
                                              (h=0.7)
MW_rv                          float          N/A             Extinction parameter
MW_av                          float          N/A             Extinction parameter
                                                              from F19 dust map
=============================  ============   ==============  ========================

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
