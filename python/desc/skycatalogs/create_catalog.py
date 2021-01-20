"""
Code to create a sky catalog for a particular object type
"""


class SkyCatalog(object):
    """
    A base class with derived classes for galaxies, static (w.r.t. coordinates)
    point sources, SSOs
    """

    def __init__(self, object_type, output_partition, output_type='parquet'):
        """
        Parameters
        ----------
        object_type
                     one of a list of known types, for now represented
                     as strings: ('galaxy', 'sn', 'star', 'agn',...)
        output_partition
                     e.g. healpix parameters
        output_type
                     one of a list of known types, e.g. ('parquet', 'csv',..)
        """
        self.object_type = object_type
        self.output_type = output_type


    def create_catalog(pixels):
        raise NotImplementedError

class GalaxySkyCatalog(SkyCatalog):

    def __init__(self, galaxy_truth, output_type='parquet'):
        """
        Parameters
        ----------
        galaxy_truth: object with methods to describe and read in
                      galaxy truth information
        output_type:  file (or in-memory?) type of output

        self.super('galaxy', output_type)
        self.galaxy_truth = galaxy_truth
        """


# May want a base truth class for this
class GalaxyTruth():
        """
        Responsible for reading from a source like CosmoDC2 catalog
        and making available to GalaxySkyCatalog object
        """
