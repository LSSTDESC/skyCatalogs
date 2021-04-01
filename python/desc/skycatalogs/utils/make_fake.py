import io
import argparse

'''
Read in a bit of old instance catalogs (single healpixel; bulge, disk and
knots).  For each of the three:
1. Transpose; store in dataframe. Columns are
prefix (string 'object')                           ignore
id     (identifies galaxy)                         use
ra,dec                                             use
magnitude   (for band observed in this visit)        ??
sedfilepath                                         probably ignore
redshift                                            use
gamma1,gamma2,kappa                                 use
raOffset,decOffset                                  ignore
spatialmodel                                        goes in config
majorAxis,minorAxis,positionAngle,sindex            use (spatial_params)
internalExtinctionModel                             goes in config
internalAv,internalRv                               use
galacticExtinctionModel                             goes in config
galacticAv,galacticRv                               use

2. Add column for component type
3. Add a column for tophat values (make them up?)
4. Add empty column for rel. SED filepath(s)
5. Sort by galaxy id
6. write out as parquet

'''
