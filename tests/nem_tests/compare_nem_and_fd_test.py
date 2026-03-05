import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis.strategies import floats

from dorban import Core
from dorban.finite_differences.finite_difference_current_calculator import (
    CurrentCalculatorFD,
)
from dorban.finite_differences.solve_finite_difference import FDSettings
from dorban.geometry.boundary_conditions import Reflector
from dorban.geometry.cartesian import Cartesian
from dorban.materials import Fissionable, Isotope
from dorban.nem.cmfd_current_calculator import CMFDCurrentCalculator
from dorban.nem.solve_nem import NEMSettings
from dorban.solve_equation import solve_k


@given(*[floats(1, 10, allow_nan=False, allow_infinity=False)] * 4)
@settings(deadline=None, max_examples=20)
def test_four_cells_rectangle(a, b, c, d):
    geometry = Cartesian([np.array([a, b]), np.array([c, d])], [Reflector()] * 4)
    mat1 = Fissionable(
        "Fuel",
        scatter=np.array([[0, 0], [0.02, 0]]),
        absorb=np.array([0.01, 0.08]),
        nusigmaf=np.array([0, 0.135]),
        chi=np.array([1, 0]),
        diffusion=np.array([1.5, 0.4]),
    )
    mat2 = Isotope(
        "Reflector",
        scatter=np.array([[0, 0], [0.04, 0]]),
        absorb=np.array([0, 0.01]),
        diffusion=np.array([2, 0.3]),
    )
    materials = [mat1, mat2, mat1, mat2]
    current_calc_fd = CurrentCalculatorFD.from_isotopes(materials)
    current_calc_cmfd = CMFDCurrentCalculator.from_isotopes(materials, 2)
    core_fd = Core(materials, 2, geometry, current_calc_fd)
    core_nem = Core(materials, 2, geometry, current_calc_cmfd)
    settings_fd = FDSettings(split=geometry.uniform_split(5 + 5 * int(max(a, b, c, d))))
    settings_nem = NEMSettings(split=geometry.uniform_split(3))
    fd_k, _ = solve_k(core_fd, settings_fd)
    nem_k, _ = solve_k(core_nem, settings_nem)
    assert nem_k - fd_k == pytest.approx(0, abs=1e-2)
