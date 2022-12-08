import os
import re

from setuptools import setup, find_packages

with open('README.md') as file:
    long_description = file.read()

packages = find_packages(where="python")
print('packages = ', packages)

def all_files_from(dir, ext=''):
    """Quick function to get all files from directory and all subdirectories
    """
    files = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in filenames:
            if filename.endswith(ext) and not filename.startswith('.'):
                files.append(os.path.join(root, filename))
    return files

# Need to clarify if this belongs in the distribution or not
## shared_data = all_files_from('data')
## print('shared_data = ', shared_data)

## configs = [os.path.join('cfg', 'latest.yaml')]

# Read in the version from python/desc/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package

version_file = os.path.join('python', 'desc', 'skycatalogs', '_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    skycatalogs_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print('skyCatalogs version is %s'%(skycatalogs_version))

dist = setup(name="skyCatalogs",
             version=skycatalogs_version,
             maintainer="Joanne Bogart",                 # also author line?
             maintainer_email="jrb@slac.stanford.edu",
             license="BSD",
             description="Writes, reads catalogs input to LSST DESC simulations",
             long_description=long_description,
             package_dir={"": "python"},
             packages=find_packages(where="python"),
             ##package_data={"skycatalogs": shared_data + configs},
             url="https://github.com/LSSTDESC/skyCatalogs",
             classifiers=[
                 "License :: OSI Approved :: BSD License",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "Programming Language :: Python",
                 "Development Status :: 3 - Alpha",
                 ],
             install_requires=['numpy', 'healpy', 'astropy', 'pyarrow',
                               'pandas', 'galsim', 'dust_extinction']
             )
