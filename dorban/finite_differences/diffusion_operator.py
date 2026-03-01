r"""
This module contains the function to construct diffusion operator
:math:`D\nabla\phi` as a sparse matrix
"""
import numpy as np
import scipy.sparse as spar

from dorban.current_calculator import CurrentCalculator
from dorban.geometry.boundary_conditions import Boundary
from dorban.geometry.geometry import FiniteGeometry


def diffusion_matrix(geometry: FiniteGeometry, E: int,
                     current_calc: CurrentCalculator) -> spar.csr_matrix:
    r"""
    This function will generate a sparse matrix representing the diffusion
    operator.



    Parameters
    ----------
    geometry: FiniteGeometry
     The geometry of the system
    E: int
     The number of energy groups
    current_calc: CurrentCalculator
     The current calculator of the system

    Returns
    -------
    spar.csr_matrix
     sparse matrix in csr format representing the linear operator flux goes to
     integral of :math:`D\nabla \phi` , :math:`\phi` being the flux.

     Warnings
     --------
     This function is one of the bottlenecks of the code, be careful to check
     the performance when editing it.
    """
    if geometry.dim == 0:
        return spar.csr_matrix((E, E))
    C = geometry.cells

    def generate_coefficents():
        for cell in range(C):
            adj = geometry.neighbors[cell]
            for j, neighbor in enumerate(adj):
                if not isinstance(neighbor, Boundary):
                    yield current_calc.compute_current_coefficients(cell,
                                                                    neighbor,
                                                                    j,
                                                                    geometry)
                else:
                    yield neighbor.compute_boundary_coefficient(geometry, E,
                                                                current_calc,
                                                                cell, j)

    row_indices = []
    col_indices = []
    values = []
    for d in generate_coefficents():
        row_indices.append(d[0])
        col_indices.append(d[1])
        values.append(d[2])
    return spar.coo_matrix(
        (np.hstack(values), (np.hstack(row_indices), np.hstack(col_indices))),
        shape=(E * C, E * C)).tocsr()
