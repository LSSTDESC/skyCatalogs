+++++++++++++++++++++++++++++++++++++
Star quantities for ops-rehearsal-3
+++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" (with information largely coming from the
UW star catalog) and a flux file. The latter contains for each object only
fluxes for lsst bands and the object id so it can be joined with the main
file. Both are parquet files.

Star main file
----------------
Note the same schema may be used for SNe. Some fields are unused for stars.

========================  ============   ==========  =========================
Name                      Datatype       Units       Description
========================  ============   ==========  =========================
object_type               string         N/A         In practice always 'star'
id                        string         N/A         Unique object identifier
ra                        double         degrees     Object right ascension
dec                       double         degrees     Object declination
host_galaxy_id            int64          N/A         Unused for stars
magnorm                   double         N/A         To be applied to SED
sed_filepath              string         N/A         Path to object's SED file
                                                     relative to env variable
                                                     SIMS_SED_LIBRARY_DIR
MW_rv                     float          N/A         Extinction parameter
MW_av                     float          N/A         Extinction parameter
                                                     from F19 dust map
mura                      double         milarcsec   RA proper motion
                                         per year
mudec                     double         milarcsec   Declin. proper motion
                                         per year
radial_velocity           double         km/sec      Radial velocity
parallax                  double         milarcsec   Parallax
variability_model         string         N/A         Unused for stars
salt2_params              int32          N/A         Unused for stars
is_variable               boolean        N/A         If true, the star has
                                                     sinusoidal variability
                                                     in magnitude
period                    double         days        period of variability
mag_amplitude             double         N/A         amplitude of magnitude
                                                     variability
phase                     double         radians     phase of variability
========================  ============   ==========  =========================

Star flux file
----------------

=============   =========   ================  ================================
Name            Datatype    Units             Description
=============   =========   ================  ================================
id              string      N/A               Unique object identifier
lsst_flux_u     float       photons/sec/cm^2  Flux integrated over lsst u-band
lsst_flux_g     float       photons/sec/cm^2  Flux integrated over lsst g-band
lsst_flux_r     float       photons/sec/cm^2  Flux integrated over lsst r-band
lsst_flux_i     float       photons/sec/cm^2  Flux integrated over lsst i-band
lsst_flux_z     float       photons/sec/cm^2  Flux integrated over lsst z-band
lsst_flux_y     float       photons/sec/cm^2  Flux integrated over lsst y-band
=============   =========   ================  ================================
