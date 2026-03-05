"""
This file will include tests parameters for running tests
"""

import numpy as np
import pytest
from numpy import array as ar

from dorban import Core, CurrentCalculatorFD
from dorban.geometry.boundary_conditions import Reflector
from dorban.geometry.polybox import Box, PolyBox
from dorban.materials import Fissionable, Isotope

scatterplut = np.array(
    [
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        [0.20, 0.00, 0.00, 0.00, 0.00, 0.00],
        [0.27, 0.18, 0.00, 0.00, 0.00, 0.00],
        [0.45, 0.50, 0.45, 0.00, 0.00, 0.00],
        [0.31, 0.35, 0.30, 6.29, 0.00, 0.00],
        [0.04, 0.05, 0.06, 0.05, 0.05, 0.00],
    ]
)

sigmaGammaPlut = np.array([0.03, 0.05, 0.07, 0.11, 0.17, 0.5])

nusigmafPlut = np.array((3.48 * 1.9, 3.09 * 1.95, 2.99 * 1.83, 2.93 * 1.7, 2.88 * 1.67, 2.86 * 2.05))

fissionPlut = np.array([1.9, 1.95, 1.83, 1.7, 1.67, 2.05])

chiPlut = np.array([0.204, 0.344, 0.168, 0.180, 0.090, 0.014])


@pytest.fixture()
def plut6():
    return Fissionable("Pu239-6groups", scatterplut, sigmaGammaPlut + fissionPlut, nusigmafPlut, chiPlut)


@pytest.fixture()
def PlutMatrix():
    return np.linalg.inv(PlutAbsorb) @ PlutFission


@pytest.fixture()
def PlutFission():
    return np.outer(chiPlut, nusigmafPlut)


@pytest.fixture()
def PlutAbsorb():
    return np.diag(sigmaGammaPlut + fissionPlut + np.sum(scatterplut, axis=0)) - scatterplut


@pytest.fixture()
def water():
    return Isotope("H20-1group", absorb=ar([1.983e-4]), scatter=ar([0]), transport=ar([2.476]))


@pytest.fixture()
def uranium():
    return Fissionable(
        "U235-lgroup",
        scatter=ar([0]),
        absorb=ar([8.983e-1]),
        nusigmaf=np.array([5.925e-3]),
        chi=np.array([1.0]),
        transport=np.array([2.531e-1]),
    )


@pytest.fixture()
def CAR_2dim():
    E = 2
    ref = Reflector()
    lengths = np.array(
        [0.47498, 0.34544, 1.87452, 1.87452, 1.87452, 1.87452, 1.87452, 1.87452, 1.87452, 0.34544, 0.47625, 0.47625]
    )
    box = Box([12, 12], boundary_conditions=[ref] * 4)
    geo = PolyBox(2, box.neighbors, [[lengths[i], lengths[11 - j]] for j in range(12) for i in range(12)])
    mat1 = Fissionable(
        "1",
        scatter=np.array([[0, 0], [1.069e-2, 0]]),
        absorb=np.array([8.983e-3, 5.892e-2]),
        nusigmaf=np.array([5.925e-3, 9.817e-2]),
        chi=np.array([1, 0]),
        transport=np.array([2.531e-1, 5.732e-1]),
        fission=np.array([2.281e-3, 4.038e-2]),
    )
    mat2 = Fissionable(
        "2",
        scatter=np.array([[0, 0], [1.095e-2, 0]]),
        absorb=np.array([8.726e-3, 5.174e-2]),
        nusigmaf=np.array([5.242e-3, 8.228e-2]),
        chi=np.array([1, 0]),
        transport=np.array([2.536e-1, 5.767e-1]),
        fission=np.array([2.003e-3, 3.385e-2]),
    )
    mat3 = Fissionable(
        "3",
        scatter=np.array([[0, 0], [1.112e-2, 0]]),
        absorb=np.array([8.587e-3, 4.717e-2]),
        nusigmaf=np.array([4.820e-3, 7.2e-2]),
        chi=np.array([1, 0]),
        transport=np.array([2.535e-1, 5.797e-1]),
        fission=np.array([1.830e-3, 2.962e-2]),
    )
    mat4 = Fissionable(
        "4",
        scatter=np.array([[0, 0], [1.113e-2, 0]]),
        absorb=np.array([8.48e-3, 4.14e-2]),
        nusigmaf=np.array([4.337e-3, 5.9e-2]),
        chi=np.array([1, 0]),
        transport=np.array([2.533e-1, 5.837e-1]),
        fission=np.array([1.632e-3, 2.428e-2]),
    )
    mat5 = Fissionable(
        "5",
        scatter=np.array([[0, 0], [1.016e-2, 0]]),
        absorb=np.array([9.593e-3, 1.626e-1]),
        nusigmaf=np.array([5.605e-3, 2.424e-2]),
        chi=np.array([1, 0]),
        transport=np.array([2.506e-1, 5.853e-1]),
        fission=np.array([2.155e-3, 9.968e-3]),
    )
    mat19 = Isotope(
        "19",
        scatter=np.array([[0, 0], [9.095e-3, 0]]),
        absorb=np.array([1.043e-3, 4.394e-3]),
        transport=np.array([2.172e-1, 4.748e-1]),
    )
    mat20 = Isotope(
        "20",
        scatter=np.array([[0, 0], [3.682e-2, 0]]),
        absorb=np.array([1.983e-4, 7.796e-3]),
        transport=np.array([2.476e-1, 1.123]),
    )
    # fmt: off
    composition = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
        , 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
        , 20, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20
        , 20, 19, 3, 2, 2, 2, 3, 3, 4, 19, 20, 20
        , 20, 19, 1, 1, 1, 5, 1, 2, 3, 19, 20, 20
        , 20, 19, 1, 1, 1, 1, 1, 1, 3, 19, 20, 20
        , 20, 19, 1, 5, 1, 1, 1, 5, 2, 19, 20, 20
        , 20, 19, 1, 1, 1, 1, 1, 1, 2, 19, 20, 20
        , 20, 19, 1, 1, 1, 5, 1, 1, 2, 19, 20, 20
        , 20, 19, 2, 1, 1, 1, 1, 1, 3, 19, 20, 20
        , 20, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20
        , 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
    # fmt: on
    names = {1: mat1, 2: mat2, 3: mat3, 4: mat4, 5: mat5, 19: mat19, 20: mat20}
    isotopes = [names[c] for c in composition]
    calc = CurrentCalculatorFD.from_isotopes(isotopes, transport=True)
    return Core(isotopes, E, geo, calc)
