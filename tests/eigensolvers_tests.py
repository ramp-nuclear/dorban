"""
tests for the various eigensolvers
"""
import numpy as np
from scipy.sparse import csr_matrix

from dorban.eigenvalue_solvers.arnoldi_arpack import generalized_arnoldi
from dorban.eigenvalue_solvers.power_iteration import generalized_eigenvalue
from dorban.eigenvalue_solvers.slepc_methods import generalized_eigenvalue_slepc


def test_generalized_eigenvalue():
    A = np.diag([1., 3., 53., 59.,1.,1.,1.,1.,1.])
    M = np.diag([7., 6., 4., 12.,1.,1.,1.,1.,1.])
    k, v = generalized_eigenvalue(A, M, rtol_vector=1e-5, atol_vector=1e-10)
    k1, v1 = generalized_arnoldi(A, M)
    k2, v2 = generalized_eigenvalue_slepc(csr_matrix(A), csr_matrix(M))
    v /= np.sum(v)
    v1 /= np.sum(v1)
    v2 /= np.sum(v2)
    assert np.isclose(k, 13.25)
    assert np.isclose(k1.XXXX, 13.25)
    assert np.isclose(k2, 13.25)
    assert np.allclose(v1, v)
    assert np.allclose(v2,v)
