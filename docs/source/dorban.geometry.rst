Geometry
========
The geometry subpackage holds the tools to construct :class:`~dorban.geometry.geometry.Geometry` objects,
including the boundary conditions at the boundary of the geometry.

.. inheritance-diagram::
   dorban.geometry.polybox.PolyBox
   dorban.geometry.cartesian.Cartesian
   dorban.geometry.triangular_tessellation.TriangularTessellation
   dorban.geometry.spherical.Spherical
   :parts: 1
   :top-classes: dorban.geometry.geometry.FiniteGeometry

There are three types of geometries:
    - Cartesian geometry
    - Triangular geometry
    - Spherical geometry

There are two geometry classes which represent the cartesian geometry, they are
the :class:`~dorban.geometry.polybox.PolyBox` class and the
:class:`~dorban.geometry.cartesian.Cartesian` class.

There is a unique class representing the triangular geometry which is the
:class:`~dorban.geometry.triangular_tessellation.TriangularTessellation` class.

The :class:`~dorban.geometry.spherical.Spherical` class represents an n-dimensional sphere which is homogenous
in the radial directions.


Abstract Base Geometry
----------------------

.. automodule:: dorban.geometry.geometry
   :members:
   :no-undoc-members:
   :show-inheritance:

Cartesian
---------

.. automodule:: dorban.geometry.cartesian
   :members:
   :no-undoc-members:
   :show-inheritance:

PolyBox
-------

.. automodule:: dorban.geometry.polybox
   :members:
   :no-undoc-members:
   :show-inheritance:

Triangular Tessellation
-----------------------

.. automodule:: dorban.geometry.triangular_tessellation
   :members:
   :no-undoc-members:
   :show-inheritance:



Spherical
---------

.. automodule:: dorban.geometry.spherical
   :members:
   :no-undoc-members:
   :show-inheritance:
