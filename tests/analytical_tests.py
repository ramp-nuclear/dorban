import math
from typing import Tuple

import numpy as np
import pytest
import scipy as sp
import scipy.optimize as op
from scipy.optimize import fsolve

from dorban.eigenvalue_solvers.slepc_methods import generalized_eigenvalue_slepc
from dorban.finite_differences.finite_difference_current_calculator import (
    CurrentCalculatorFD,
)
from dorban.finite_differences.solve_finite_difference import solve_k_diffusion
from dorban.geometry.boundary_conditions import Void
from dorban.geometry.cartesian import Cartesian
from dorban.materials import CrossSectionData, Fissionable
from dorban.system import Core


def diml_group1(absorb: float, fission: float, diff: float,
                length: float, cells: int) -> Tuple[float, np.array]:
    """
    measures the difference between the given k and the actual eigenvalue of
    the system. system should represent a 1 dimensional homogenous system
    with void boundary conditions
    Parameters
    ----------
    absorb the macroscopic cross-section of absorption
    fission the macroscopic cross-section of nusigmaf
    diff the diffusion coefficient
    lengths the length of the system
    """
    l = length / 2
    a = op.bisect(lambda x: x * math.tan(x) - l / (2 * diff), 0, 1.1) / l
    k = fission / (a ** 2 * diff + absorb)
    cells *= 2
    flux = np.cos([a * l * (-1 + 2 * i / cells) for i in
                   np.arange(cells)])[1::2]
    assert np.allclose(flux, flux[::-1])
    return k[0], flux


def fuel_in_water(fuel: CrossSectionData, water: CrossSectionData,
                  fuel_length: float,
                  water_length: float, diffusion: Tuple[float, float],
                  fuel_cells: int, water_cells: int) -> Tuple[float, np.array]:
    """
    Compute the k eigenvalue of 1 dimensional system of a fuel slab
     surrounded by 2 water slabs in 1 energy group

    Parameters
    ----------
    fuel - Fissionable Isotope with cross sections in 1 energy group
    water - non Fissionable Isotope with cross sections in 1 energy group
    fuel_length - half the thickness of the fuel slab
    water_length - the thickness of 1 water slab
    diffusion_coefficients - tuple of length 2 containing the diffusion
    coefficients of the fuel at 0 and of the water at 1,
    Returns
    -------
    The k eigenvalue of the problem and a dict whose keys are points on the x
    axis and its values are the values of flux at those points.
    """
    fuel_cells *= 2
    water_cells *= 2
    water_laplace_coefficient = -water.absorb[0] / diffusion[1]
    a_2 = sp.sqrt(water_laplace_coefficient)[0]
    length = fuel_length + water_length
    C = (np.cos(a_2 * length) + 2 * diffusion[1] * a_2 * np.sin(
        a_2 * length)) / \
        (2 * diffusion[1] * a_2 * np.cos(a_2 * length) - np.sin(a_2 * length))
    C = C[0]
    tan_a1 = fuel_length * (a_2 * diffusion[1] * np.sin(a_2 * fuel_length) -
                            C * a_2 * diffusion[1] * np.cos(
                a_2 * fuel_length)) / \
             (diffusion[0] * np.cos(a_2 * fuel_length) +
              C * diffusion[0] * np.sin(a_2 * fuel_length))
    a_1 = fsolve(lambda x: x * math.tan(x) - np.XXXX(tan_a1), np.array([0]),
                 xtol=1e-10)[0] \
          / fuel_length
    k = fuel.nusigmaf[0] / (fuel.absorb + diffusion[0] * a_1 ** 2)
    point = length * a_2
    ratio = (2 * diffusion[1] * a_2 * np.sin(point) - np.cos(point)) / (
            2 * diffusion[1] * a_2 * np.cos(point) + np.sin(point))
    fuel_flux = np.array(
        [np.cos((a_1 * fuel_length * (-1 + 2 * i / fuel_cells)))
         for i in range(fuel_cells)])
    water_flux = np.array([
        np.cos(a_2 * (fuel_length + i * water_length)) + ratio * np.sin(
            a_2 * (fuel_length + i * water_length))
        for i in np.arange(water_cells) / (water_cells - 1)])
    water_flux = np.array([a[0].XXXX for a in water_flux])
    constant = water_flux[0] / fuel_flux[0]
    fuel_flux *= constant
    flux = np.hstack([water_flux[-1::-2], fuel_flux[1::2], water_flux[1::2]])
    return float(k), flux


def test_hom_1dim_1group():
    cells = 10000
    length = 10

    fuel = Fissionable("1", scatter=np.array([0]),
                       absorb=np.array([8.983e-3]),
                       nusigmaf=np.array([5.925e-3]),
                       chi=np.array([1]),
                       transport=np.array([2.531e-1]),
                       fission=np.array([2.281e-3]))
    diff_coefficent = np.array([1 / (3 * fuel.transport)])
    isotope = fuel
    isotopes = np.array([fuel] * cells)
    geo = Cartesian([np.full((cells,), length / cells)],
                    (Void(), Void())).to_polybox()
    current = CurrentCalculatorFD(
        np.reshape(np.full((cells,), diff_coefficent), (cells, 1)))
    system = Core(isotopes, 1, geo, current)
    solver = generalized_eigenvalue_slepc
    k, flux = solve_k_diffusion(system, solver=solver)
    k_analytical, flux_analytical = diml_group1(
        isotope.absorb,
        isotope.nusigmaf, system.current_calc.dc[0],
        length,
        cells)
    flux_analytical = flux_analytical * np.sum(
        flux) / np.sum(
        flux_analytical)
    kerror = np.abs(k_analytical - k)
    fluxerror = np.linalg.norm(
        np.abs(flux - flux_analytical) / flux_analytical) / cells
    assert kerror < 1e-2
    assert fluxerror < 1e-2


def test_fuel_in_water(water, uranium):
    """
    This test compares two different solutions of the 1dim one-speed probelm of fuel surrounded by water from
    both sides with void boundary conditions after the water. The first solution is obtained from the finite
    differences kernel of DORBAN and the other is the analytic solution. Both the k eigenvalue and the
     corresponding eigenflux are obtained. The eigenvalues are compared as is and the fluxes are compared
     after normalization as the eigenflux isn't well-defined.
    """
    water_den = 1
    ur_den = 1
    diff_water = 1 / (3 * water.transport)
    diff_fuel = 1 / (3 * uranium.transport)
    water_length = 10
    fuel_length = 20
    cells_water = 3000
    cells_fuel = 6000
    diff_coefficient = np.hstack([np.full((cells_water,), diff_water),
                                  np.full((cells_fuel,), diff_fuel),
                                  np.full((cells_water,), diff_water)])
    den = {water: (water_den, 0, water_den), uranium: (0, ur_den, 0)}
    cells = (cells_water, cells_fuel, cells_water)
    help_isotopes = {mat: np.hstack(
        [np.full((c,), d) for c, d in zip(cells, densities)])
        for mat, densities in den.items()}
    isotopes = np.zeros_like(help_isotopes[water], dtype=object)
    for mat, d in help_isotopes.items():
        density = d[np.nonzero(d)[0]][0]
        isotopes[np.nonzero(d)] = mat * density
    clengths = (water_length, fuel_length, water_length)
    lengths = np.hstack([np.full((c,), length / c)
                         for c, length in zip(cells, clengths)])
    geo = Cartesian((lengths,), (Void(), Void())).to_polybox()
    system = Core(isotopes, 1, geo,
                  CurrentCalculatorFD(np.reshape(diff_coefficient,
                                                 (len(diff_coefficient), 1))))
    solver = generalized_eigenvalue_slepc
    k, flux = solve_k_diffusion(system, solver=solver,
                                lin_atol=1e-10, lin_rtol=1e-10,
                                atol_vector=1e-12)

    k_analytic, flux_analytic = fuel_in_water(uranium, water,
                                              fuel_length / 2,
                                              water_length,
                                              (diff_fuel,
                                               diff_water),
                                              cells_fuel,
                                              cells_water)
    flux_analytic = flux_analytic * np.sum(flux) / np.sum(flux_analytic)
    assert np.isclose(k, k_analytic, atol=1e-6)
    assert np.max(np.abs(flux - flux_analytic) / flux_analytic) < 1e-3


