"""
This module contains a data class named Core the holds all the
information about the system. This is the main object of the code.

A Core is made out of isotopes, energy groups geometry and a current calculator.,
 - The isotopes are given as list of :class:`CrossSectionData <dorban.materials.CrossSectionData>` objects
 - The energy groups can be given either by a sequence of the average energy at each group or by just the number of energy groups.
 - The geometry is a :class:`FiniteGeometry <dorban.geometry.geometry.FiniteGeometry>` object
 - The current calculator is a :class:`CurrentCalculator <dorban.current_calculator.CurrentCalculator>` object
"""
from dataclasses import dataclass
from typing import Sequence, Union

import numpy as np

from dorban.current_calculator import CurrentCalculator
from dorban.geometry.geometry import FiniteGeometry
from dorban.materials import CrossSectionData

IsotopeData = Sequence[CrossSectionData]


@dataclass
class Core:
    """

    Attributes
    ----------
    isotopes: Sequence[CrossSectionData]
        Sequence of the data about the cross-sections in each cell
    E: Union[int, Sequence[float]]
        The number of energy groups, or a sequence of the energy levels
    geometry: FiniteGeometry
        The geometry of the system
    current_calc: CurrentCalculator
        Object that allows to compute the current between any 2
        adjacent cells in the system
    """
    isotopes: IsotopeData
    E: Union[int, Sequence]
    geometry: FiniteGeometry
    current_calc: CurrentCalculator

    def __post_init__(self):
        try:
            self.energy_range = self.E
            self.E = len(self.E)
        except TypeError:
            pass
        self.size = self.E * self.geometry.cells
        for mat in self.isotopes:
            if mat.E != self.E:
                raise ValueError(
                    f"The material {mat.name} has {mat.E} energy groups "
                    f"which doens't equal to {self.E}, "
                    f"the core's number of energy groups")


def mesh_refinement(system: Core, split: Sequence) -> Core:
    """
    This function generates a new Core with the same materials but a
    refined mesh. The refinement happens in each direction according to the
    splitting data. The splitting data specifies to how many cell
    to divide each cell in each direction. The exact type of the splitting data
    is determined by the geometry of the core to be refined.

    Parameters
    ----------
    system:Core
        the original core
    split:Sequence
     splitting data which is consistent with the geometry of the system

    Returns
    -------
    Core
        The refined core.
    """
    geo = system.geometry
    geometry = geo.refine_mesh(split)
    current_calc = system.current_calc.mesh_refine(geo, split, geometry)
    assemblies = geo.array_refine(np.arange(geo.cells), split).astype(int)
    isotopes = [system.isotopes[assembly] for assembly in assemblies]
    return Core(isotopes, system.energy_range, geometry, current_calc)


