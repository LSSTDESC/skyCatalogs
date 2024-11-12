import os
import sys
from multiprocessing import Process, Pipe
import sqlite3
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import healpy
import json
from .objects.base_object import LSST_BANDS
from .objects.trilegal_object import TrilegalConfigFragment
from .utils.config_utils import assemble_provenance
from .utils.config_utils import assemble_file_metadata

"""
Code for creating sky catalogs for trilegal stars
"""

__all__ = ['TrilegalMainCatalogCreator', 'TrilegalFluxCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 1000000
_DEFAULT_TRUTH_CATALOG = 'lsst_sim.simdr2'
_DEFAULT_START_EPOCH = 2000

def _do_trilegal_flux_chunk(send_conn, trilegal_collection, instrument_needed,
                            l_bnd, u_bnd):
    '''
    end_conn         output connection
    trilegal_collection  information from main file
    instrument_needed List of which calculations should be done
    l_bnd, u_bnd     demarcates slice to process

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    # Not implemented yet
    return out_dict


class TrilegalMainCatalogCreator:
    # Note dl is not part of desc-python or LSST Pipelines; it must be
    # separately installed
    from dl import queryClient as qc

    def __init__(self, catalog_creator, truth_catalog=_DEFAULT_TRUTH_CATALOG,
                 start_epoch=_DEFAULT_START_EPOCH):
        '''
        Parameters
        ----------
        catalog_creator  instance of MainCatalogCreator
        input_catalog    name of Trilegal catalog to be queried
        '''
        self._catalog_creator = catalog_creator
        self._truth_catalog = truth_catalog
        self._start_epoch = start_epoch
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._output_dir
        self._row_group_size = _DEFAULT_ROW_GROUP_SIZE

    @property
    def trilegal_truth(self):
        return self._truth_catalog

    def _create_main_schema(self, metadata_input=None,
                            metadata_key='provenance'):

        fields = [
            pa.field('id', pa.string()),
            pa.field('ra', pa.float64()),
            pa.field('dec', pa.float64()),
            pa.field('pmdec', pa.float64()), # proper motion in dec
            pa.field('pmracosd', pa.float64()), # proper motion in cos(ra)
            pa.field('vrad', pa.float64()), # radial velocity
            pa.field('mu0', pa.float64()), # true distance modulus
            #  pa.field('parallax', pa.float64()),
            ]
        # The truth catalog supplies only mu0, not parallax
        # from https://en.wikipedia.org/wiki/Distance_modulus
        # luminous distance in parsecs =
        #  (1 + mu0/5) ** 10
        #
        # Formula for distance using parallax:
        #     d = E-sun/tan(0.5*theta) where E-sun is one AU
        # and 0.5*theta is parallax.   Or, approximately for small angles
        #    d (in parsecs) ~= 1/parallax

        if metadata_input:
            metadata_bytes = json.dumps(metadata_input).encode('utf8')
            final_metadata = {metadata_key: metadata_bytes}
        else:
            final_metadata = None

        return pa.schema(fields, metadata=final_metadata)

    def _write_hp(self, hp, arrow_schema):
        from dl import queryClient as qc
        _NSIDE = 32
        _TRILEGAL_RING_NSIDE = 256

        def _next_level(pixel):
            return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

        # Find all subpixels contained in our pixel with nside 256,
        # which is what the Trilegal catalog uses
        current_nside _NSIDE
        pixels = [healpy.ring2nest(_NSIDE, hp)]

        while current_nside < _TRILEGAL_RING_NSIDE:
            pixels = [pix for p in pixels for pix in _next_level(p)]
            current_nside = current_nside * 2

        ring_pixels = [healpy.nest2ring(_TRILEGAL_RING_NSIDE, p) for p in pixels]
        in_pixels = ','.join(str(p) for p in ring_pixels)

        # Form query
        to_select = ['ra', 'dec', 'pmracosd', 'pmdec', 'vrad', 'mu0']
        q = 'select ' + ','.join(to_select)
        q += f' from {self._truth_catalog} where ring256 in ({in_pixels})'

        results = qc.query(adql=q, fmt='pandas')
        n_row = len(results['ra'])

        print(f'rows returned: {nrow}')

        # generate id
        id_prefix = f'{self._truth_catalog}_hp{hp}_'
        results['id'] = [f'{id_prefix}{n}' for n in range(n_row)]
