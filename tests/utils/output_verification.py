import pyarrow.parquet as pq
import numpy as np
import logging

'''
Collect code which could be useful for verifying output, regression testing
'''
logger = logging.getLogger('compare_sky')

ch = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)

def compare_parquet(path1, path2, meta_only=False, columns=None,
                    row_limit=None, log_level='INFO'):
    '''
    Compare parquet files.
    Parameters
    ----------
    path1, path2      strings      Paths to files to be compared
    meta_only         boolean      If true only compare metadata and schema
    columns           list of column names  If None compare all columns
    '''

    logger.setLevel(log_level)

    logger.info(" ")
    logger.info(f"Comparing files {path1} and {path2}")
    # Compare schemas
    schema1 = pq.read_schema(path1)
    schema2 = pq.read_schema(path2)
    schema_match = (schema1 == schema2)
    if not schema_match:
        logger.warning("Schema mismatch")
        logger.warning(f"Schema 1: {schema1}")
        logger.warning(f"Schema 2: {schema2}")

    meta1 = pq.read_metadata(path1)
    meta2 = pq.read_metadata(path2)

    metadata_match = (meta1 == meta2)
    if metadata_match:
        logger.info("Metadata matches")
    else:
        logger.warning("Metadata mismatch")
        for k in ('created_by', 'num_columns', 'num_rows', 'num_row_groups',
                  'format_version', 'serialized_size'):
            v1 = meta1.__getattribute__(k)
            v2 = meta2.__getattribute__(k)
            if v1 != v2:
                logger.warning(f'for key {k} values are {v1} != {v2}')

    if not row_limit and not metadata_match:
        row_limit = min(meta1.num_rows, meta2.num_rows)
    if not columns:
        if not schema_match:
            logger.warning("Schema mismatch.  Cannot compare all columns")
            return
        else:
            columns = list(schema1.names)
    else:
        s1_names = set(schema1.names)
        s2_names = set(schema2.names)
        if not (set(columns)).issubset(s1_names.intersection(s2_names)):
            logger.error(f"Cannot compare columns {columns}; not in both schemas")
            return
    d1 = pq.read_table(path1, columns=columns).to_pandas()
    d2 = pq.read_table(path2, columns=columns).to_pandas()

    bad_cols = 0
    for c in columns:
        logger.debug(f"Working on column {c}")
        logger.debug(f"Type is {type(d1[c][0])}")
        if type(d1[c][0]) == np.ndarray:
            logger.debug("Skipping column with list-valued items")
            continue     # compare of lists can't be done this way
        if row_limit:
            d1_limited = d1[c][:row_limit]
            d2_limited = d2[c][:row_limit]
            cmp_c = d1_limited.compare(d2_limited)
        else:
            cmp_c = d1[c].compare(d2[c])
        if len(cmp_c) > 0:
            logger.warning(f"{len(cmp_c)} differences in column {c}")
            bad_cols += 1

    logger.info(f"Found {bad_cols} mismatched columns")

if __name__ == '__main__':
    import os
    sky_root = '/global/cscratch1/sd/jrbogart/desc/skycatalogs'
    fname = 'galaxy_9556.parquet'
    #fname = 'galaxy_flux_9556.parquet'
    p1 = os.path.join(sky_root, 'gal_flux_save', fname)
    p2 = os.path.join(sky_root, 'gal_flux', fname)
    p3 = os.path.join(sky_root, 'true_divide', fname)

    columns = ['galaxy_id', 'ra', 'dec']
    ##columns = ['galaxy_id', 'lsst_flux_u', 'lsst_flux_g']

    compare_parquet(p1, p2, log_level='DEBUG')
    compare_parquet(p1, p3, columns=columns)
    compare_parquet(p2, p3, columns=columns)
