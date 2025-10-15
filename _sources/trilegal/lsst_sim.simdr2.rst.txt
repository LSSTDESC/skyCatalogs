+++++++++++++++++++++++++++++++++++++++++++++++++++++
Quantities for TRILEGAL stars version lsst_sim.simdr2
+++++++++++++++++++++++++++++++++++++++++++++++++++++

These data are sourced from the [Astro Data Lab database](https://datalab.noirlab.edu/), in particular from the lsst_sim.simdr2 catalog.

Data are partitioned by (nside=32, ring ordering) healpixel. For each pixel
there is a so-called "main file" (with information largely coming from the
source catalog) and a flux file. The latter contains for each object only
fluxes for lsst bands and the object id so it can be joined with the main
file. Both are parquet files.

TRILEGAL main file
-------------------

This file contains a subset of the columns in the source table. It also has an
id column.  The ids are unique across the catalog (the source catalog has no
similar column). For a complete description of all columns in the source
catalog go to
https://datalab.noirlab.edu/data-explorer?showTable=lsst_sim.simdr2

========================  ============     ===================================
Name  (Data Lab name)     Datatype         Description
========================  ============     ===================================
id                        string           Unique object identifier
ra                        double           Right Ascension J2000 (deg)
dec                       double           Declination J2000 (deg)
av                        double           V-band Extinction
pmdec                     float            Proper motion in dec
pmracosd                  float            Proper motion in RA*cos(Dec)
vrad                      float            Radial velocity
mu0                       float            True distance modulus
evol_label (label)        int              Evolutionary phase
logT       (logte)        float            log(effective temperature) (deg K)
logg                      float            log(surface gravity)  (cgs)
logL       (logl)         float            log(luminosity)   (L_sun)
Z          (z)            float            Linear heavy element abundance
umag                      float            LSST u-band magnitude
gmag                      float            LSST g-band magnitude
rmag                      float            LSST r-band magnitude
imag                      float            LSST i-band magnitude
zmag                      float            LSST z-band magnitude
ymag                      float            LSST y-band magnitude
========================  ============     ===================================

TRILEGAL flux file
-------------------

Fluxes were calculated using the quantities from the main file and LSST
throughputs_2023_09_07.

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
