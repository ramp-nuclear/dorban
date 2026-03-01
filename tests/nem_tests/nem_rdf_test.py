from dorban import Isotope, Fissionable, Cartesian, Core, CurrentCalculatorFD, Reflector
import numpy as np

from dorban.cross_sections_operators import _absorb as absorb, _nufission as nufission
from dorban.rdfs.nem_rdf import compute_boundary_fluxes_using_nem, compute_rdfs
import scipy.linalg as la


def one_dim_colorset_flux_boundary(assembly: Core, current: np.array,
                                   k: float) -> np.array:
    """
    Function used to compute the flux on the boundary induced by a given current
    from one side and reflector on the other, used for computing rdf of a 1
    dimensional 2x1 color set problem.
    Parameters
    ----------
    current - the current on the boundary
    k- the k of the color set problem

    Returns
    -------
    the flux on the boundary
    """
    d = assembly.current_calc.dc_on_face(0, 0, assembly.geometry)
    if len(d.shape) == 1:
        d = np.diag(d)
    mat = assembly.isotopes[0]
    a = absorb(mat) - 1 / k * nufission(mat)
    alpha = la.sqrtm(np.linalg.inv(d) @ a)
    l = assembly.geometry.distance_to_face(0, 0) * 2
    v = la.solve(-d @ alpha @ la.sinm(l * alpha), current)
    return la.cosm(l * alpha) @ v


def test_compare_nem_flux_on_boundary_with_analytic_flux_on_boundary_in_2x1_problem():
    ref = Isotope(name="ref", scatter=np.array(
        [[1.86347243e+00, 2.89718656e-02],
         [8.35668223e-04, 6.30553914e-01]]),
                  absorb=np.array([0.04193057, 0.00242824]),
                  diffusion=np.array([0.27593463, 1.17322397]))
    fuel = Fissionable(name="fuel", scatter=np.array(
        [[1.21111292, 0.01692339],
         [0.00179553, 0.51985346]]),
                       absorb=np.array([0.099985, 0.0100711]),
                       nusigmaf=np.array([0.13072259, 0.00678641]),
                       chi=np.array([0, 1]),
                       diffusion=np.array([0.4039059, 1.33992157]))
    current = np.array([0.01044089, 0.129204])
    d = 0.1
    geometry = Cartesian([[d, d]],
                         boundary_condition=[Reflector()] * 2).to_polybox()
    color_set = Core([fuel, ref], 2, geometry,
                     CurrentCalculatorFD.from_isotopes([fuel, ref]))
    k = 0.91436
    zero = np.zeros((1, 2))
    flux_1 = one_dim_colorset_flux_boundary(color_set, current, k)
    average_flux = np.array([52.23043424, 320.29340333])
    flux_2 = compute_boundary_fluxes_using_nem(current, zero, zero, d, d,
                                               d, d, zero, zero, fuel, k)[1]
    assert np.all(np.isclose(flux_1, flux_2, atol=1e-2))


def test_nem_rdfs_2x1_problem_regression(file_regression):
    ref = Isotope(name="ref", scatter=np.array(
        [[1.86347243e+00, 2.89718656e-02],
         [8.35668223e-04, 6.30553914e-01]]),
                  absorb=np.array([0.04193057, 0.00242824]),
                  diffusion=np.array([0.27593463, 1.17322397]))
    fuel = Fissionable(name="fuel", scatter=np.array(
        [[1.21111292, 0.01692339],
         [0.00179553, 0.51985346]]),
                       absorb=np.array([0.099985, 0.0100711]),
                       nusigmaf=np.array([0.13072259, 0.00678641]),
                       chi=np.array([0, 1]),
                       diffusion=np.array([0.4039059, 1.33992157]))
    current = np.array([0.01044089, 0.129204])
    d = 0.1
    geometry = Cartesian([[d, d], [1, 1]],
                         boundary_condition=[Reflector()] * 4).to_polybox()
    color_set = Core([fuel, ref, fuel, ref], 2, geometry,
                     CurrentCalculatorFD.from_isotopes([fuel, ref, fuel, ref]))
    k = 0.91436
    currents = np.zeros((4, 4, 2))
    currents[0, 1] = current
    currents[1, 0] = current
    currents[2, 1] = current
    currents[3, 0] = current
    rdfs = compute_rdfs(color_set, currents, k)
    file_regression.check(str(rdfs))
