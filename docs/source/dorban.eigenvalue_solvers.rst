Eigenvalue Solvers
==================

One of the most common application of the neutron diffusion approximation is to find the reactor's mulitplication factor and eigenflux.
Naturally this requires solving an eigenvalue problem. DORBAN offers some eigensolvers to solve this eigenproblem and they are contained in this package. There 3 solvers in this package: solver based on the power iteration algorithm, solver based on arpacks implementation of the arnoldi algrorithm and a tool to use all the solvers offered by the SLEPc package.

.. automodule:: dorban.eigenvalue_solvers
   :members:
   :no-undoc-members:
   :show-inheritance:

Arnoldi ARPACK
--------------

.. automodule:: dorban.eigenvalue_solvers.arnoldi_arpack
   :members:
   :no-undoc-members:
   :show-inheritance:

Power Iteration
---------------

.. automodule:: dorban.eigenvalue_solvers.power_iteration
   :members:
   :no-undoc-members:
   :show-inheritance:

SLEPc Methods
-------------

.. automodule:: dorban.eigenvalue_solvers.slepc_methods
   :members:
   :no-undoc-members:
   :show-inheritance:
