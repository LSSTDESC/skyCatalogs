Third-party catalogs
====================
The required interface for using catalogs with imSim just entails
implementing a few methods in subclasses of the ``BaseObject`` and
``ObjectCollection`` classes.  Thus, it's fairly straightforward to
implement new object types. We've provided an interface to enable
those object types to be imported and configured via a third party package.
Isolating user code this way has a couple notable benefits:

- It avoids having to coordinate adding new object types to skyCatalogs itself
- Any additional software dependencies introduced by the new code
  won't be imposed on other skyCatalogs users.

While reading the following sections you can refer to a small
example of what's required in this (docs)
`docs <https://github.com/LSSTDESC/skyCatalogs/tree/main/docs>`_
directory.  See the file ``external_example.py`` for code (but don't
try to use it as it; it's not valid python).

Object and ObjectCollection
---------------------------
A source for skyCatalogs is represented by a class derived from ``BaseObject``.
(See the file `base_object.py <https://github.com/LSSTDESC/skyCatalogs/blob/main/skycatalogs/objects/base_object.py>`_.) The derived class must, at a minimum,
implement the routines ``get_observer_sed_component`` and
``get_gsobject_components``.

To implement a third-party catalog you also need a class derived from
``ObjectCollection`` which is a container for the source objects. This class
must implement a static method ``load_collection`` which, given configuration
information (the ``sky_catalog`` argument), and a region of the sky, returns
an object collection.  In typical implementations, ``load_collection``
would use the configuration to find a source (such as a file or files)
of data for the region and then use that to create the
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
The minimal configuration file needed for a third-party object type looks like
this:

.. code-block:: yaml

   catalog_dir: .
   catalog_name: my_external_catalog
   object_types:
     my_object_type:
       area_partition: None
       module: my_module
   skycatalog_root: .
   provenance:
     versioning:
       schema_version: 1.3.0

It's usually convenient to set ``catalog_dir`` and ``skycatalog_root`` to
current directory (that is, the one containing the configuration file) as
is done above.

Natively-supported catalogs are often large
enough that they are partitioned into multiple files, in which case the
value for ``area_partition`` might be something like
``{type: healpix, order: ring, nside: 32}``, but for third-party catalogs
``None`` is usually the right value.

The value of ``module``  should be the name of
your python implementation module, containing the classes derived from
``BaseObject`` and ``ObjectCollection`` as described above. It needs
to be accessible so that it can be imported.

The configuration used for the example code
is a little more complex.  It illustrates how two object types can
be handled by the same code (indicated by the value of the ``module`` field),
but differentiated by the value of the additional field ``object_param``.
It looks like this:

.. code-block:: yaml

   catalog_dir: .
   catalog_name: my_external_catalog
   object_types:
     object_type_1:
       area_partition: None
       module: external_example
       object_param: 111.0
     object_type_2:
       area_partition: None
       module: external_example
       object_param: 2.1
   skycatalog_root: .
   provenance:
     versioning:
       schema_version: 1.3.0

Just as the field ``object_param`` has been added for this example,
you are free to add whatever you want to the configuration for your object
type long as it's valid yaml syntax.  It will be up to your module to
interpret it properly.

It is possible to run a simulation making use of other object types
along with yours, including natively-supported ones, by describing all
of them in the same configuration file in the ``object_types`` section.
