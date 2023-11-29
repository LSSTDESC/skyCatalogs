import pyarrow as pa
import logging

__all__ = ['make_galaxy_schema', 'make_galaxy_flux_schema',
           'make_pointsource_schema', 'make_star_flux_schema']


def _add_roman_fluxes(fields):
    fields += [pa.field('roman_flux_W146', pa.float32(), True),
               pa.field('roman_flux_R062', pa.float32(), True),
               pa.field('roman_flux_Z087', pa.float32(), True),
               pa.field('roman_flux_Y106', pa.float32(), True),
               pa.field('roman_flux_J129', pa.float32(), True),
               pa.field('roman_flux_H158', pa.float32(), True),
               pa.field('roman_flux_F184', pa.float32(), True),
               pa.field('roman_flux_K213', pa.float32(), True)]
    return fields


# This schema is not the same as the one taken from the data,
# probably because of the indexing in the schema derived from a pandas df.
def make_galaxy_schema(logname, sed_subdir=False, knots=True,
                       galaxy_type='cosmodc2'):
    logger = logging.getLogger(logname)  # maybe move this above if:
    if galaxy_type == 'cosmodc2':
        fields = [pa.field('galaxy_id', pa.int64()),
                  pa.field('ra', pa.float64(), True),
                  # metadata={"units" : "radians"}),
                  pa.field('dec', pa.float64(), True),
                  # metadata={"units" : "radians"}),
                  pa.field('redshift', pa.float64(), True),
                  pa.field('redshift_hubble', pa.float64(), True),
                  pa.field('peculiar_velocity', pa.float64(), True),
                  pa.field('shear_1', pa.float64(), True),
                  pa.field('shear_2', pa.float64(), True),
                  pa.field('convergence', pa.float64(), True),
                  pa.field('size_bulge_true', pa.float32(), True),
                  pa.field('size_minor_bulge_true', pa.float32(), True),
                  pa.field('sersic_bulge', pa.float32(), True),
                  pa.field('size_disk_true', pa.float32(), True),
                  pa.field('size_minor_disk_true', pa.float32(), True),
                  pa.field('sersic_disk', pa.float32(), True),

                  # Depending on value of --dc2-like option, value for
                  # ellipticity_2_true column will differ
                  pa.field('ellipticity_1_disk_true', pa.float64(), True),
                  pa.field('ellipticity_2_disk_true', pa.float64(), True),
                  pa.field('ellipticity_1_bulge_true', pa.float64(), True),
                  pa.field('ellipticity_2_bulge_true', pa.float64(), True),
                  pa.field('sed_val_bulge',
                           pa.list_(pa.float64()), True),
                  pa.field('sed_val_disk',
                           pa.list_(pa.float64()), True),
                  pa.field('bulge_magnorm', pa.float64(), True),
                  pa.field('disk_magnorm', pa.float64(), True),
                  pa.field('MW_rv', pa.float32(), True),
                  pa.field('MW_av', pa.float32(), True)]
        if knots:
            logger.debug("knots requested")
            fields.append(pa.field('sed_val_knots',
                                   pa.list_(pa.float64()), True))
            # For sizes API can alias to disk sizes
            # position angle, shears and convergence are all
            # galaxy-wide quantities.
            fields.append(pa.field('n_knots', pa.float32(), True))
            fields.append(pa.field('knots_magnorm', pa.float64(), True))

        if sed_subdir:
            fields.append(pa.field('bulge_sed_file_path', pa.string(), True))
            fields.append(pa.field('disk_sed_file_path', pa.string(), True))

    elif galaxy_type == 'diffsky':
        fields = [pa.field('galaxy_id', pa.int64()),
                  pa.field('ra', pa.float64(), True),
                  # metadata={"units" : "radians"}),
                  pa.field('dec', pa.float64(), True),
                  # metadata={"units" : "radians"}),
                  pa.field('redshift', pa.float64(), True),
                  pa.field('redshiftHubble', pa.float64(), True),
                  pa.field('peculiarVelocity', pa.float64(), True),
                  pa.field('shear1', pa.float64(), True),
                  pa.field('shear2', pa.float64(), True),
                  pa.field('convergence', pa.float64(), True),
                  pa.field('spheroidHalfLightRadiusArcsec', pa.float32(), True),
                  pa.field('diskHalfLightRadiusArcsec', pa.float32(), True),

                  # Not sure these are what we want
                  pa.field('diskEllipticity1', pa.float64(), True),
                  pa.field('diskEllipticity2', pa.float64(), True),
                  pa.field('spheroidEllipticity1', pa.float64(), True),
                  pa.field('spheroidEllipticity2', pa.float64(), True),
                  pa.field('um_source_galaxy_obs_sm', pa.float32(), True),
                  pa.field('MW_rv', pa.float32(), True),
                  pa.field('MW_av', pa.float32(), True)]

    debug_out = ''
    for f in fields:
        debug_out += f'{f.name}\n'
    logger.debug(debug_out)
    return pa.schema(fields)


def make_galaxy_flux_schema(logname, galaxy_type='cosmodc2',
                            include_roman_flux=False):
    '''
    Will make a separate parquet file with lsst flux for each band
    and galaxy id for joining with the main galaxy file
    '''
    logger = logging.getLogger(logname)
    logger.debug('Creating galaxy flux schema')
    fields = [pa.field('galaxy_id', pa.int64()),
              # should flux fields be named e.g. lsst_cmodel_flux_u?
              pa.field('lsst_flux_u', pa.float32(), True),
              pa.field('lsst_flux_g', pa.float32(), True),
              pa.field('lsst_flux_r', pa.float32(), True),
              pa.field('lsst_flux_i', pa.float32(), True),
              pa.field('lsst_flux_z', pa.float32(), True),
              pa.field('lsst_flux_y', pa.float32(), True)]
    if include_roman_flux:
        fields = _add_roman_fluxes(fields)
    return pa.schema(fields)


def make_star_flux_schema(logname, include_roman_flux=False):
    '''
    Will make a separate parquet file with lsst flux for each band
    and id for joining with the main star file
    '''
    logger = logging.getLogger(logname)
    logger.debug('Creating star flux schema')
    fields = [pa.field('id', pa.string()),
              pa.field('lsst_flux_u', pa.float32(), True),
              pa.field('lsst_flux_g', pa.float32(), True),
              pa.field('lsst_flux_r', pa.float32(), True),
              pa.field('lsst_flux_i', pa.float32(), True),
              pa.field('lsst_flux_z', pa.float32(), True),
              pa.field('lsst_flux_y', pa.float32(), True)]
    if include_roman_flux:
        fields = _add_roman_fluxes(fields)
    return pa.schema(fields)


def make_pointsource_schema():
    '''
    Ultimately should handle stars both static and variable, SN, and AGN
    For now add everything needed for SN and put in some additional
    star fields, but not structs for star variability models
    '''

    salt2_fields = [
        pa.field('z', pa.float64(), True),
        pa.field('t0', pa.float64(), True),
        pa.field('x0', pa.float64(), True),
        pa.field('x1', pa.float64(), True),
        pa.field('c', pa.float64(), True)]
    fields = [pa.field('object_type', pa.string(), False),
              pa.field('id', pa.string(), False),
              pa.field('ra', pa.float64(), False),
              pa.field('dec', pa.float64(), False),
              pa.field('host_galaxy_id', pa.int64(), True),
              pa.field('magnorm', pa.float64(), True),
              pa.field('sed_filepath', pa.string(), True),
              pa.field('MW_rv', pa.float32(), True),
              pa.field('MW_av', pa.float32(), True),
              pa.field('mura', pa.float64(), True),
              pa.field('mudec', pa.float64(), True),
              pa.field('radial_velocity', pa.float64(), True),
              pa.field('parallax', pa.float64(), True),
              pa.field('variability_model', pa.string(), True),
              pa.field('salt2_params', pa.struct(salt2_fields), True)
              ]
    return pa.schema(fields)


def make_pointsource_flux_schema(logname, include_roman_flux=False):
    '''
    Will make a separate parquet file with lsst flux for each band
    and id for joining with the main star file.
    For static sources mjd field could be -1. Or the field could be
    made nullable.
    '''
    logger = logging.getLogger(logname)
    logger.debug('Creating pointsource flux schema')
    fields = [pa.field('id', pa.string()),
              pa.field('lsst_flux_u', pa.float32(), True),
              pa.field('lsst_flux_g', pa.float32(), True),
              pa.field('lsst_flux_r', pa.float32(), True),
              pa.field('lsst_flux_i', pa.float32(), True),
              pa.field('lsst_flux_z', pa.float32(), True),
              pa.field('lsst_flux_y', pa.float32(), True),
              pa.field('mjd', pa.float64(), True)]
    if include_roman_flux:
        fields = _add_roman_fluxes(fields)
    return pa.schema(fields)
