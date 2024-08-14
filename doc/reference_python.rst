===========================
The ``skyCatalogs`` package
===========================

Reference documentation for core objects with the ``sktyCatalogs`` package

.. _skyCatalog_class:

The skyCatalog class
=====================

The ``skyCatalog`` class represents a catalog. Users gain access
to a catalog by calling its static method ``open_catalog``.

.. autoclass:: skycatalogs.skyCatalogs.SkyCatalog
   :members:

   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.open_catalog
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_objects_by_region
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_object_type_by_region
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_object_type_names
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_object_type_by_hp
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_hps_by_region
   .. automethod:: skycatalogs.skyCatalogs.SkyCatalog.get_hps_by_type

.. _baseobject_class:

There are several classes (one per supported object type) representing
individual objects as observed at a particular time, all subclassing
``BaseObject``.

.. autoclass:: skycatalogs.objects.base_object.BaseObject
   :members:

   .. automethod:: skycatalogs.objects.base_object.BaseObject.id
   .. automethod:: skycatalogs.objects.base_object.BaseObject.ra
   .. automethod:: skycatalogs.objects.base_object.BaseObject.dec
   .. automethod:: skycatalogs.objects.base_object.BaseObject.object_type
   .. automethod:: skycatalogs.objects.base_object.BaseObject.get_native_attribute
   .. automethod:: skycatalogs.objects.base_object.get_gsobject_components
   .. automethod:: skycatalogs.objects.base_object.get_observer_sed_component
   .. automethod:: skycatalogs.objects.base_object.get_observer_sed_components
   .. automethod:: skycatalogs.objects.base_object.get_total_observer_sed
   .. automethod:: skycatalogs.objects.base_object.get_flux
   .. automethod:: skycatalogs.objects.base_object.get_LSST_flux
   .. automethod:: skycatalogs.objects.base_object.get_LSST_fluxes
   .. automethod:: skycatalogs.objects.base_object.get_roman_flux
   .. automethod:: skycatalogs.objects.base_object.get_roman_fluxes

.. _objectlist_class:

The ``ObjectList`` class is a container for instances of ``BaseObject``.
It subclasses ``Sequence``.

.. autoclass:: skycatalogs.objects.base_object.ObjectList


