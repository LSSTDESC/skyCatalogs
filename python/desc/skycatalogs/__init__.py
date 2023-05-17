##from ._version import *
try:
    # For Python >= 3.8
    from importlib import metadata
except ImportError:
    # For Python < 3.8
    import importlib_metadata as metadata

try:
    __version__ = metadata.version("skyCatalogs")
except metadata.PackageNotFoundError:
    pass

from .skyCatalogs import *
from .catalog_creator import *
