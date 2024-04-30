+++++++++++++++++++++++++++++++++++++
SNANA quantities for Roman-Rubin sims
+++++++++++++++++++++++++++++++++++++
Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" in parquet format, containing static
quantities associated with each source, and an hdf5 file for quantities
which vary with time or wavelength or both.


SNANA main file
----------------

========================  ============   ==========  ==========================
Name                      Datatype       Units       Description
========================  ============   ==========  ==========================
id                        int64          N/A         Unique object identifier
ra                        double         degrees     Object right ascension
dec                       double         degrees     Object declination
host_id                   int64          N/A         Host galaxy
gentype                   int32          N/A         Type of object simulated
model_name                string         N/A         Type name
start_mjd                 float          days        Earliest visibility
end_mjd                   float          days        Latest visibility
z_CMB                     float          N/A         redshift
mw_EBV                    float          N/A         MW extinction parameter
mw_extinction_applied     boolean        N/A         Here always False
AV                        float          N/A         For host galaxy, randomly
                                                     assigned
RV                        float          N/A         For host galaxy
v_pec                     float          km/sec      Peculiar velocity. Not
                                                     yet simulated so zero.
host_ra                   double         degrees     Host galaxy RA
host_dec                  double         degrees     Host galaxy decl.
host_mag_g                float          N/A         Host mag. in Rubin g-band
host_mag_i                float          N/A         Host mag. in Rubin i-band
host_mag_F                float          N/A         Host mag. in Roman F-band
host_sn_sep               float          arcsec      Separation from host 
peak_mjd                  float          days        Time of peak flux
                                                     (model-dependent)
peak_mag_g                float          N/A         g-band mag. at peak
peak_mag_i                float          N/A         i-band mag. at peak
peak_mag_F                float          N/A         F-band mag. at peak
lens_dmu_applied          boolean        N/A         False (distance modulus
                                                     lensing shift not applied)
model_param_names         list(string)   N/A         Model parameter names
model_param_values        list(string)   N/A         Model parameter values
MW_av                     float          N/A         MW extinction parameter
MW_rv                     float          N/A         MW extinction parameter
========================  ============   ==========  ==========================

SNANA varying quantities
------------------------

The top-level keys for each file are the ids for sources in that pixel, the
same as the ids appearing in the parquet file but here stored as strings.
Under each id are the following keys

=============  ===============  ================  ===============================
Name           Datatype         Units             Description
=============  ===============  ================  ===============================
mjd            float array      days              Time axis for flambda
lambda         double array     angstroms         Wavelength axis for flambda
flambda        float 2-d array  erg/sec/A/cm^2    Flux indexed by (mjd,lambda)
mag_<band>     float array      N/A               Magnitude for each mjd computed
                                                  using 5 angstrom bins.
magcor_<band>  float array      N/A               Correction to magnitude
                                                  computed using flambda
synmag_<band>  float array      N/A               Magnitude computed from flambda
=============  ===============  ================  ===============================

where <band> takes on values u, g, r, i, z, y (Rubin bandpasses) and F, H, J,
K, R, W, Y, Z (for Roman).

mjd values have uniform spacing of 1 day. lambda values have uniform spacing
of 100 angstroms.

synmag + magcor = mag
