"""
This module contains the CurrentCalculator protocol
"""

from typing import Protocol, Sequence, Tuple

import numpy as np

from dorban.geometry import FiniteGeometry


class CurrentCalculator(Protocol):
    """
    Such classes are responsible for calculating the current between 2 points
    in the system. They hold all the needed information like the diffusion
    coefficients and the discontinuity factors.
    """

    dc: np.array

    def compute_current_coefficients(
        self, cell: int, neighbor: int, face: int, geometry: FiniteGeometry
    ) -> Tuple[np.array, np.array, np.array]:
        """
        This method computes the coefficients in the diffusion matrix
        that represents the current between a cell and its neighbor

        Parameters
        ----------
        cell: int
            the number of the cell
        neighbor: int
            the number of the neighbor
        face: int
            the number of the face between cell and neighbor
        geometry: Geometry
            the geometry of the system

        Returns
        -------
        Tuple[np.array, np.array, np.array]
            tuple of the row indices in the diffusion matrix in which the coefficents
            appear, the column indices in the diffusion matrix in which the coefficents
            appear, the values of the diffusion matrix
        """
        raise NotImplementedError

    def mesh_refine(self, geo: FiniteGeometry, split: Sequence, refgeo: FiniteGeometry) -> "CurrentCalculator":
        """
        creates and returns a new CurrentCalculator which can calculate the
        current for the refined system.

        Parameters
        ----------
        geo: FiniteGeometry
            the geometry of the system
        split:
            split data whose structure is determined by the geometry of the system
        refgeo: FiniteGeometry
            the refined geometry
        """
        raise NotImplementedError

    def dc_on_face(self, cell: int, face: int, geometry: FiniteGeometry) -> float:
        """
        returns the diffusion coefficient across a given face

        Parameters
        ----------
        cell: int
         the index of the cell
        face: int
         the index of the face
        geometry: FiniteGeometry
         the geometry of the system
        """
        raise NotImplementedError
