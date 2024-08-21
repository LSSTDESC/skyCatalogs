===========================
The ``skyCatalogs`` package
===========================

Reference documentation for core objects with the ``skyCatalogs`` package

.. _skyCatalog_class:

The skyCatalog class
=====================

The ``skyCatalog`` class represents a catalog. Users gain access
to a catalog by calling its static method ``open_catalog``.

..
   .. autofunction:: skycatalogs.skyCatalogs.SkyCatalog.open_catalog

.. autoclass:: skycatalogs.skyCatalogs.SkyCatalog
   :members:

.. _baseobject_class:

There are several classes (one per supported object type) representing
individual objects as observed at a particular time, all subclassing
``BaseObject``.

.. autoclass:: skycatalogs.objects.base_object.BaseObject
   :members:


.. _objectlist_class:

The ``ObjectList`` class is a container for instances of ``BaseObject``.
It subclasses ``Sequence``.

.. autoclass:: skycatalogs.objects.base_object.ObjectList
   :members:               


