## Galaxy Schema

Taken from schema of parquet file produced by current code.
This draft only includes bulge and disk.  Some fields will be added for knots.  Most likely AGN will also go in the galaxy file.

| Fieldname | Type   | Remarks |
| --------- | ------ | ------- |
| galaxy_id | bigint |         |
| ra        | R8     |         |
| dec       | R8     |         |
| redshift  | R8     |         |
| shear_1   | R8     |         |
| shear_2   | R8     | treecor |
| convergence | R8 | |
| size_bulge_true | R4 | |
| size_minor_bulge_true | R4| |
| sersic_bulge | R4 | |
| size_disk_true | R4 | |
| size_minor_disk_true | R4| |
| sersic_disk | R4 | |
| position_angle | R8 | |
| sed_val_bulge | R4 * N | list of tophat flux values |
| sed_val_disk | R4 * N | list of tophat flux values |
| internalAv_bulge | R8 | |
| internalRv_bulge | R8 | |
| internalAv_disk | R8 | |
| internalRv_disk | R8 | |
| bulge_magnorm | R8 | |
| disk_magnorm | R8 | |
| MW_rv | R4 |  |
| MW_av_lsst_u | R4 | |
| MW_av_lsst_g | R4 | |
| MW_av_lsst_r | R4 | |
| MW_av_lsst_i | R4 | |
| MW_av_lsst_z | R4 | |
| MW_av_lsst_y | R4 | |

## Pointsource schema

Will include stars (variable or not) and SNe.  May also include AGN. This
draft assumes stars and SNe only.

| Field | Type   | Remarks |
| --------- | ------ | ------- |
| object type | string | 'star' or 'sne' |
| id          | bigint? | Current instance catalogs have strings for sn ids |
| ra | R8 | |
| dec | R8 | |
| host galaxy id | bigint | Only include column if we want it for SNe  |
| mag norm | R | |
| sed filepath | string | relative to $SIMS_SED_LIBRARY_DIR |
| MW_av_lsst_u | R4 | per-band extinction value for this ra, dec |
| MW_av_lsst_g | R4 | |
| MW_av_lsst_r | R4 | |
| MW_av_lsst_i | R4 | |
| MW_av_lsst_z | R4 | |
| MW_av_lsst_y | R4 | |

Several other fields which might have been included are not because either
* they are not relevant for point sources or
* they are constant for all objects in the file

In the latter case the value will be stored in the config.
