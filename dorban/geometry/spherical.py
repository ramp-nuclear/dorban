"""
Contains the Spherical Geometry class, which can be used to model spherical
system for example the Godiva and Jezebel experiments. This geometry can't be
used to model a real reactor but can be used for toy models.
"""

from typing import Sequence

import numpy as np
from more_itertools import pairwise
from scipy.special import gamma

from dorban.geometry.boundary_conditions import Boundary
from dorban.geometry.geometry import FiniteGeometry


class Spherical(FiniteGeometry):
    r"""A class that represents a spherical geometry.

    Parameters
    ----------

    lengths: np.array
     Array of the lengths of the cells. 1 dimensional array.
    boundary_condition: Boundary
     The boundary condition around the sphere
    dim: int
     The dimension of the system

    Attributes
    ----------
    cells: int
     Number of cells
    volumes: np.array
     Array which contains at index `i` the volume of the cell numbered `i`.
    neighbors: List[List[int]]
     List that at index `i` has a list of the neighbors of the cell numbered `i`.

    """

    def __init__(self, lengths: np.array, boundary_condition: Boundary, dim: int):
        self.dim = dim
        self.lengths = lengths
        self.boundary_conditions = boundary_condition
        self.cells = len(lengths)
        self.volume_constant = np.power(np.pi, self.dim / 2) / gamma(self.dim / 2 + 1)
        self.surface_area_constant = 2 * np.power(np.pi, self.dim / 2) / gamma(self.dim / 2)
        self.neighbors = self.neighbors_calc()
        self.radii = np.hstack([[0], np.cumsum(self.lengths)])
        self.volumes = np.array(
            [
                self.volume_constant * (np.power(r1, self.dim) - np.power(r0, self.dim))
                for r0, r1 in pairwise(self.radii)
            ]
        )

    def direction(self, face: int) -> np.array:
        """
        Method to compute the direction in space of the unit vector from the
        center of a cell to the face.

        Parameters
        ----------
        face: int
         the number of the face

        Returns
        -------
        np.array
         Unit vector that represents the direction to the face
        """
        d = np.zeros(self.dim)
        d[0] = 2 * (face % 2) - 1
        return d

    def neighbors_calc(self) -> Sequence[Sequence[int]]:
        """
        Function to compute the neighbors of each cell

        Returns
        -------
        Sequence[Sequence[int]]
         List of lists, the list at index `i` are the neighbors of the cell
         numbered `i`.
        """
        return [[1]] + [[x - 1, x + 1] for x in range(1, self.cells - 1)] + [[self.cells - 2, self.boundary_conditions]]

    def surface_area(self, cell: int, face: int) -> float:
        """
        Function to compute the surface area of a face of a cell
        the faces are numbered according to the order: (left, right)

        Parameters
        ----------
        cell: int
         The cell number
        face: int
         The face number
        """
        return self.surface_area_constant * np.power(
            self.lengths[0] if cell == 0 else self.radii[cell + face], self.dim - 1
        )

    def distance_to_face(self, cell: int, face: int) -> float:
        """
        Function that computes the distance from the center of the cell to the face.
        The faces are numbered according to the order:
        (left, right, forward, back, up, down)
        face//2 is the number of the axis of the line from the center of the
        cell to the face

        Parameters
        ----------
        cell: int
         The cell number
        face: int
         The face number
        """
        return self.lengths[0] if cell == 0 else self.lengths[cell] / 2

    def opposite_face(self, face: int) -> int:
        """
        Method to compute the number of the face of the cell who is adjacent to
        a given cell along a face.

        Parameters
        ----------
        face: int
         The number of the face in a given cell along which
         the two cells are adjacent. An integer between 0 and 6

        Returns
        -------
        int
         The number of the face in the other cell, a number between 0 and 6
        """
        return face - (2 * (face % 2) - 1)

    def array_refine(self, array: np.array, split: np.array, **kwargs) -> np.array:
        """
        Refines an array according to the split parameters. For example this
        function is used in order to refine an array of densities of some
        Isotope.

        Parameters
        ----------
        array: np.array
         The array that has to be refined
        split: np.array
         Array of ints, of the system's size. Each of its values
         specifies into how many parts each cell has to be divided.

        Returns
        -------
        np.array
         The new array containing the same information but suitable to the
         refined geometry
        """
        return np.repeat(array, split, axis=0)

    def refine_mesh(self, split: np.array) -> "Spherical":
        """
        This function returns a new Spherical Geometry with refined mesh
        according to the split data.

        Parameters
        ----------
        split: np.array
         Array of ints, of the system's size. Each of its values
         specifies into how many parts each cell has to be divided.

        Returns
        -------
        Spherical
         new Geometry with a refined mesh
        """
        return Spherical(self.array_refine(self.lengths / split, split), self.boundary_conditions, self.dim)

    def uniform_split(self, subcells: int) -> Sequence[np.array]:
        """
        Method that creates a split suitable for the geometry, which splits
        each cell into a constant number of cells in each dimension.

        Parameters
        ----------
        subcells: int
         The number of sub cells to split each cell in each dimension.
         Overall the number of cells is multiplied by this parameter to the
         power of the dimension.

        Returns
        -------
        Sequence[np.array]
         Split data to be used with the Spherical geometry.
        """
        return np.full(self.cells, subcells)
