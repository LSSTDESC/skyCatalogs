##from ._version import *
try:
    # For Python >= 3.8
    from importlib import metadata
except ImportError:
    # For Python < 3.8
    import importlib_metadata as metadata

__version__ = metadata.version("mydescpackage")
from .skyCatalogs import *
from .catalog_creator import *
