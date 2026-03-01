Boundary Conditions
===================
The boundary conditions are implemented by having cells in a geometry whose neihgbors aren't numbers of other cells.
Rather the neighbors are instances of a Boundary object.
The supported Boundary objects are presented below.


.. inheritance-diagram::
   dorban.geometry.boundary_conditions.ExtrapolationLength
   dorban.geometry.boundary_conditions.CurrentCondition
   dorban.geometry.boundary_conditions.IncomingCurrent
   dorban.geometry.boundary_conditions.KnownCurrent
   dorban.geometry.boundary_conditions.OrangePortal
   dorban.geometry.boundary_conditions.BluePortal
   dorban.geometry.boundary_conditions.Reflector
   dorban.geometry.boundary_conditions.Void
   dorban.geometry.boundary_conditions.ZeroFlux
   :parts: 1

.. automodule:: dorban.geometry.boundary_conditions
   :members:
   :no-undoc-members:
   :show-inheritance:
