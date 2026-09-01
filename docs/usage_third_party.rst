Third-party catalogs
====================
For simulations requiring only relatively small-scale catalogs of a new type,
the effort involved in creating and supporting that type in the production
repos ``skyCataglogs`` and ``skyCatalogs_creator``
may be prohibitive. There is an alternative, which is to write code satisfying

- the minimal interface needed by imSim's use of skyCatalogs
- requirements skyCatalogs has on non-native object types

The code may or may not depend on input files of some sort,
but in either case no changes to the skyCatalogs code will be needed.

While reading the following sections you can refer to a stripped-down example
of what's required in the
`tests <https://github.com/LSSTDESC/skyCatalogs/tree/main/tests>`_
directory.  See the files ``external_catalog.py`` for code and
``external_skycat.yaml`` for configuration.

Object and ObjectCollection
---------------------------
A source for skyCatalogs is represented by a class derived from ``BaseObject``.
(See the file `base_object.py <https://github.com/LSSTDESC/skyCatalogs/blob/main/skycatalogs/objects/base_object.py>`_.) The derived class must, at a minimum,
implement the routines ``get_observer_sed_component`` and
``get_gsobject_components``.  A realistic implementation of the former would
return a galsim SED, computed from information read in on the fly or when the
object was created.

To implement a third-party catalog you also need a class derived from
``ObjectCollection`` which is a container for the source objects. This class
must implement a static method ``load_collection`` which, given configuration
information (the ``sky_catalog`` argument), and a region of the sky, returns
an object collection.  In typical implementations, ``load_collection``
would use the configuration to find the file or files
containing data for the region and then use that to create the
collection.

In order to make skyCatalogs aware of your object type, the derived
collection class also needs to implement a static method ``register``
looking very much like the one in the example: just substitute appropriate
values for ``object_class`` and ``collection_class``.
Your module must also have
an implementation for ``register_objects`` like the example, substituting your
collection class name for "ExternalCollection".

Configuration
-------------
The example file ``external_skycat.yaml`` defines two object types to be
handled by the same code (indicated by the value of the ``module`` field),
but differentiated by the value of ``object_param``. (Alternatively one could
write different code for the two object types in different modules.)
The minimum configuration needed would omit the second object type and
the ``object_param`` field from the remaining object type description.
You are free to add other fields to the configuration for your object
type, for example where to find input files.
Natively-supported catalogs are often large
enough that they are partitioned into multiple files, in which case the
value for ``area_partition`` might be something like
``{type: healpix, order: ring, nside: 32}``, but for third-party catalogs
``None`` is usually the right value.

The value of ``module`` is
required.  It should be
your python implementation module, containing the classes derived from
``BaseObject`` and ``ObjectCollection`` as described above.

It is possible to run a simulation making use of other object types
along with yours,
including natively-supported ones, by describing all
of them in the same configuration file in the ``object_types`` section.
