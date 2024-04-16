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

========================  ============   ==========  ========================
Name                      Datatype       Units       Description
========================  ============   ==========  ========================
galaxy_id                 int64          N/A         Unique object identifier
ra                        double         degrees     object right ascension
dec                       double         degrees     object declination
redshift                  double
redshift_hubble           double
peculiar_velocity         double
shear_1                   double
shear_2                   double
convergence               double
size_bulge_true           float
size_minor_bulge_true     float
sersic_bulge              float
size_disk_true            float
size_disk_bulge_true      float
sersic_disk               float
ellipticity_1_disk_true   double
ellipticity_2_disk_true   double
ellipticity_1_bulge_true  double
ellipticity_2_bulge_true  double
sed_val_bulge             list(double)
sed_val_disk              list(double)
bulge_magnorm             double
disk_magnorm              double
MW_rv                     float
MW_av                     float
sed_val_knots             list(double)
n_knots                   float
knots_magnorm             double
========================  ============   ==========  ========================

Galaxy flux file
----------------

=============   =========   ==============  ========================
Name            Datatype    Units           Description
=============   =========   ==============  ========================
galaxy_id       int64       N/A             Unique object identifier
lsst_flux_u     float
lsst_flux_g     float
lsst_flux_r     float
lsst_flux_i     float
lsst_flux_z     float
lsst_flux_y     float
=============   =========   ==============  ========================
