"""
This file will contain benchmarks for calculations using discontinuity factors
"""

import numpy as np
import pytest

from dorban import solve_k
from dorban.finite_differences.finite_difference_current_calculator import (
    DiscontinuityCurrentCalculator,
)
from dorban.geometry.boundary_conditions import Reflector
from dorban.geometry.cartesian import Cartesian
from dorban.materials import Fissionable
from dorban.settings import FDSettings
from dorban.system import Core


def test_partially_poisened():
    """
    This is the benchmark 1.1 from the appendix to Smith's paper about
    discontinuity factors
    """
    geo = Cartesian(
        lengths=[[21, 21, 21]], boundary_condition=[Reflector(), Reflector()]
    ).to_polybox()
    A = Fissionable(
        "A",
        scatter=np.array([[0, 0], [0.017, 0]]),
        absorb=np.array([0.009, 0.08]),
        nusigmaf=np.array([0.006, 0.104]),
        chi=np.array([1, 0]),
        diffusion=np.array([1.32, 0.383]),
    )
    B = Fissionable(
        "B",
        scatter=np.array([[0, 0], [0.017, 0]]),
        absorb=np.array([0.009, 0.09]),
        nusigmaf=np.array([0.006, 0.104]),
        chi=np.array([1, 0]),
        diffusion=np.array([1.32, 0.383]),
    )
    A_B = Fissionable(
        "AB",
        scatter=np.array([[0, 0], [0.017, 0]]),
        absorb=np.array([0.009, 0.084219]),
        nusigmaf=np.array([0.006, 0.104]),
        chi=np.array([1, 0]),
        diffusion=np.array([1.32, 0.383]),
    )
    isotopes = np.array([A, A_B, B])
    discontinuity = np.array(
        [[[1, 1], [1, 1]], [[1.036, 1.0856], [0.9508, 0.8919]], [[1, 1], [1, 1]]]
    )
    diff = np.vstack([mat.diffusion for mat in isotopes])
    calc = DiscontinuityCurrentCalculator(diff, discontinuity)
    system = Core(isotopes, 2, geo, calc)
    subcells = 1000
    split = geo.uniform_split(subcells)
    settings = FDSettings(split=split)
    k, flux = solve_k(system, settings)
    assert k == pytest.approx(1.04514, abs=1e-5)
