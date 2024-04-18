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
object_type               string          -          In practice always 'star'
id                        string          -          Unique object identifier
ra                        double         degrees     object right ascension
dec                       double         degrees     object declination
host_galaxy_id            int64           -          Unused for stars
magnorm                   double          -          To be applied to SED
sed_filepath              string          -          Path to object's SED file
                                                     relative to env variable
						     SIMS_SED_LIBRARY_DIR
MW_rv                     float           -          Extinction parameter
MW_av                     float           -          Extinction parameter
                                                     from F19 dust map
mura                      double         milarcsec   RA proper motion
                                         per year
mudec                     double         milarcsec   Declin. propermotion
radial_velocity           double         km/sec      Radial velocity
parallax                  double         milarcsec   Parallax
variability_model         string          -          Unused for stars
salt2_params              int32           -          Unused for stars
is_variable               boolean         -          If true remaining
                                                     columns are meaningful
period                    double
mag_amplitude             double
phase                     double
========================  ============   ==========  =========================
