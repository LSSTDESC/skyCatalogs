![CI](https://github.com/LSSTDESC/skyCatalogs/actions/workflows/ci.yml/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/LSSTDESC/skyCatalogs/branch/main/graph/badge.svg?branch=main)](https://codecov.io/gh/LSSTDESC/skyCatalogs)


# The skyCatalogs package

This package contains
* code to create sky catalogs from input like cosmoDC2 (for galaxies) and similar catalogs for other source types
* an API to access information in the catalogs and products, such as flux calculations, derived from the catalogs

Sky Catalogs are chiefly intended for use by simulation programs. The physical representation is meant to compactly represent information needed to compute fluxes at a given time, under specified conditions (e.g. band).

For more complete information on installation and use see
[https://lsstdesc.org/skyCatalogs/](https://lsstdesc.org/skyCatalogs/)

## Set-up and testing
From bash
```
$ source <skyCatalogs install directory>/setup/setup.sh
$ nosetests <skyCatalogs install directory>
```

## Demo

## People
* [Joanne Bogart](https://github.com/LSSTDESC/skyCatalogs/issues/new?body=@JoanneBogart) (SLAC)

## License, etc.

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/LSSTDESC/skyCatalogs/issues/new).
