import numpy as np

from dorban import Core
from dorban.cross_sections_operators import fission_matrix, sigma_a
from dorban.geometry.geometry import Point


def test_fission_0dim(plut6,PlutFission):
    isotopes = [plut6]
    system = Core(isotopes, 6, Point(), [1])
    assert np.allclose(fission_matrix(system).toarray(), PlutFission)


def test_absorbtion_0dim(plut6,PlutAbsorb):
    isotopes = [plut6]
    system = Core(isotopes, 6, Point(), [1])
    assert np.allclose(sigma_a(system).toarray(), PlutAbsorb)
