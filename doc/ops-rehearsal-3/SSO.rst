+++++++++++++++++++++++++++++++++++++
SSO quantities for ops-rehearsal-3
+++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" and a flux file. The latter contains
for each object only fluxes for lsst bands, the object id and the time (mjd)
of the observation. The latter two fields may be used to join with the
main file.

SSO main file
-------------

========================  ============   ==============  ===========================
Name                      Datatype       Units           Description
========================  ============   ==============  ===========================
id                        string         N/A             Unique object identifier
mjd                       double         days            Time of observation
ra                        double         degrees         Object right ascension
dec                       double         degrees         Object declination
trailed_source_mag        double         N/A             Mag. normalization for SED
ra_rate                   double         deg/day         RA rate of change * cos(dec)
dec_rate                  double         deg/day         Declination rate of change
========================  ============   ==============  ===========================

SSO flux file
----------------

=============   =========   ================  ================================
Name            Datatype    Units             Description
=============   =========   ================  ================================
id              string      N/A               Unique object identifier
mjd             double      days              Time of observation (Julian date)
lsst_flux_u     float       photons/sec/cm^2  Flux integrated over lsst u-band
lsst_flux_g     float       photons/sec/cm^2  Flux integrated over lsst g-band
lsst_flux_r     float       photons/sec/cm^2  Flux integrated over lsst r-band
lsst_flux_i     float       photons/sec/cm^2  Flux integrated over lsst i-band
lsst_flux_z     float       photons/sec/cm^2  Flux integrated over lsst z-band
lsst_flux_y     float       photons/sec/cm^2  Flux integrated over lsst y-band
=============   =========   ================  ================================
