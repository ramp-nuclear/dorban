from typing import Callable, Tuple

import numpy as np

from dorban.eigenvalue_solvers.slepc_methods import generalized_eigenvalue_slepc
from dorban.finite_differences import fd_cross_sections_operators as fd
from dorban.finite_differences.diffusion_operator import diffusion_matrix
from dorban.settings import FDSettings
from dorban.system import Core, mesh_refinement
from dorban.utils import collapse_flux


def solve_k_fd(system: Core,
               settings: FDSettings,
               solver: Callable = generalized_eigenvalue_slepc
               ) -> Tuple[float, np.array]:
    r"""
    function to solve the diffusion k problem using the finite differences method

    Parameters
    ----------
    system: Core
     the Core object that represents the core.
    settings: FDSettings
     the settings to use for the solution
    solver: Callable
     eigenvalue problem solver, it gets two matrices A,B and solves the linear problem
     :math:`A\phi=kB\phi`. Default to the slepc solver :func:`generalized_slepc <dorban.eigenvalue_solvers.slepc_methods.generalized_slepc>`.
     The solver may also except additional keyword arguments used to determine various tolerances.

    Returns
    -------
    Tuple[float,np.array]
     tuple of the k eigenvalue and the corresponding flux.
    """
    equation_type = solve_k_diffusion if settings.adjoint is False else solve_adjoint
    if settings.split is None:
        return equation_type(system, solver, **settings.kwargs)
    core = mesh_refinement(system, settings.split)
    k, fine_flux = equation_type(core, solver, **settings.kwargs)
    return k, collapse_flux(system.geometry,fine_flux,core.geometry.volumes,settings.split)


def solve_k_diffusion(system: Core, solver: Callable = generalized_eigenvalue_slepc,
                      **kwargs) -> Tuple[
    float, np.array]:
    """
    Solves the diffusion k eigenvalue problem of the system

    Parameters
    ----------
    system:Core
     class that contains all the data about the system
    solver:Callable
     the eigenvsolver to be used
    kwargs
     keyword arguments to be passed down to the eigenvalue solver

    Returns
    -------
    Tuple[float, np.array]
     tuple that contains the k eigenvalue and the corresponding direct flux
    """
    return solver(
        fd.cmfd_fission(system),
        M=(fd.cmfd_absorption(system) +
           diffusion_matrix(system.geometry, system.E, system.current_calc)),
        **kwargs)


def solve_adjoint(system: Core,solver: Callable = generalized_eigenvalue_slepc
                  , **kwargs) -> Tuple[float, np.array]:
    """
    Solves the adjoint of the diffusion k eigenvalue problem of the system

    Parameters
    ----------
    system: Core
        class that contains all the data about the system
    solver: Callable
        eigenvalue problem solver, default is the SLEPc jd solver.
    kwargs
        keyword arguments to be passed down to the eigenvalue solver

    Returns
    -------
    Tuple[float, np.array]
        the k eigenvalue and the corresponding adjoint flux
    """
    fission = fd.cmfd_fission(system)
    absorber = (diffusion_matrix(system.geometry, system.E,
                                 system.current_calc)
                + fd.cmfd_absorption(system))
    return solver(fission.transpose(),
                                  M=absorber.transpose(),
                                  **kwargs)