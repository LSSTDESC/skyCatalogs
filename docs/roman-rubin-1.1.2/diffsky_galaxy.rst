++++++++++++++++++++++++++++++++++++++++++++++++++
Galaxy quantities for Roman-Rubin joint simulation
++++++++++++++++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" (with information largely coming from the
diffsky galaxy catalog), a flux file and a SED file. The flux file contains for
each object only fluxes for lsst and Roman bands and the object id so it can be
joined with the main file. Both it and the main file are parquet files. The
SED file is in hdf5 format. The file names take the form galaxy_<hp>.parquet
(main file), galaxy_flux_<hp>.parquet (flux file) and galaxy_sed_<hp>.hdf5 for
SED files.  In all cases <hp> is the id of the healpixel covered by the file.

Galaxy main file
----------------

=============================  ========  ==========  ==========================
Name                           Datatype  Units       Description
=============================  ========  ==========  ==========================
galaxy_id                      int64     N/A         Unique object identifier
ra                             double    degrees     Object right ascension
dec                            double    degrees     Object declination
redshift                       double    N/A         Cosmological redshift
                                                     with line-of-sight motion
redshiftHubble                 double    N/A         Cosmological redshift
peculiarVelocity               double    km/sec      Peculiar velocity
shear1                         double    N/A         Shear (gamma) component 1,
                                                     treecorr/Galsim convention
shear2                         double    N/A         Shear (gamma) component 2,
                                                     treecorr/Galsim convention
convergence                    double    N/A         Convergence (kappa)
spheroidHalfLightRadiusArcsec  float     arcsec      Bulge component half-light
                                                     radius
diskHalfLightRadiusArcsec      float     arcsec      Disk component half-light
                                                     radius
diskEllipticity1               double    N/A         Ellipticity component 1
                                                     for disk, not lensed
diskEllipticity2               double    N/A         Ellipticity component 2
                                                     for disk, not lensed
spheroidEllipticity1           double    N/A         Ellipticity component 1
                                                     for bulge, not lensed
spheroidEllipticity2           double    N/A         Ellipticity component 2
                                                     for bulge, not lensed
um_source_galaxy_obs_sm        float     solar mass  Stellar mass
                                         (h=0.7)
MW_rv                          float     N/A         Extinction parameter
MW_av                          float     N/A         Extinction parameter
                                                     from F19 dust map
=============================  ========  ==========  ==========================



Galaxy flux file
----------------

===============  ========   ================  ====================================
Name             Datatype   Units             Description
===============  ========   ================  ====================================
galaxy_id        int64      N/A               Unique object identifier
lsst_flux_u      float      photons/sec/cm^2  Flux integrated over lsst u-band
lsst_flux_g      float      photons/sec/cm^2  Flux integrated over lsst g-band
lsst_flux_r      float      photons/sec/cm^2  Flux integrated over lsst r-band
lsst_flux_i      float      photons/sec/cm^2  Flux integrated over lsst i-band
lsst_flux_z      float      photons/sec/cm^2  Flux integrated over lsst z-band
lsst_flux_y      float      photons/sec/cm^2  Flux integrated over lsst y-band
roman_flux_W146  float      photons/sec/cm^2  Flux integrated over roman W146-band
roman_flux_R062  float      photons/sec/cm^2  Flux integrated over roman R062-band
roman_flux_Z087  float      photons/sec/cm^2  Flux integrated over roman Z087-band
roman_flux_Y106  float      photons/sec/cm^2  Flux integrated over roman Y106-band
roman_flux_J129  float      photons/sec/cm^2  Flux integrated over roman J129-band
roman_flux_H158  float      photons/sec/cm^2  Flux integrated over roman H158-band
roman_flux_F184  float      photons/sec/cm^2  Flux integrated over roman F184-band
roman_flux_K213  float      photons/sec/cm^2  Flux integrated over roman K213-band
===============  ========   ================  ====================================

Galaxy SED file
---------------

This file has two sections. One (**meta**) has an array **wave_list**, shape(312),
of wavelengths (angstroms), which is used for SEDs of all galaxies in the file's healpixel.
The other section, **galaxy**, contains groups of galaxies.  Each group has entries
for individual galaxies, keyed by galaxy id. Such an entry is an array, shape(3, 312),
consisting of SEDs for the three galaxy components (disk, bulge, knot). Units are fnu.
