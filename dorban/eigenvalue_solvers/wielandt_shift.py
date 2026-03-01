"""
This module contains the Wielandt shift solver for the eigenvalue solver with the same implementation of the PARCS
code documentation
"""

from typing import Union

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as la


def wielandt_shift(
        A: Union[sparse.spmatrix, np.array, la.LinearOperator],
        M: Union[sparse.spmatrix, np.array, la.LinearOperator],
        delta: float,
        **kwargs) -> (complex, np.array):
    """
    solves the generalized eigenvalue problem Av=kMv.
    Finds the largest eigenvalue satisfying Av=kMv and its eigenvector
    Parameters
    ----------
    A - linear operator
    M - Invertible linear operator matrix
    k - guess for the top eigenvalue
    v - guess for the eigenvector
    kwargs- key word argument to be passed down to the generator
    Returns
    -------
    tuple whose first value is the eigenvalue and whose second value is the
    eigenvector
    """
    solver = GeneratorEigenvalue(A, M, delta=delta, **kwargs)
    for k, v in solver:
        continue
    return float(k), v


class GeneratorEigenvalue():
    """
    generator class that solves the generalized eigenvalue problem Av=kMv.
    each iteration is like applying M^-1A to a vector, this way after enough
    iterations we get the largest eigenvalue and its eigenvector.
    Parameters
    ----------
    A - sparse matrix
    M - Invertible sparse matrix
    k - guess for the top eigenvalue
    v - guess for the eigenvector
    tol_value - Absolute convergence tolerance for the eigenvalue
    rtol_vector - relative convergence tolerance for the eigenvector in l\inf norm
    atol_vector - absolute convergence tolerance for the eigenvector in l\inf norm

    """

    def __init__(self, A: Union[sparse.spmatrix, np.array],
                 M: Union[sparse.spmatrix, np.array], k: float = 1,
                 v: np.array = None, tol_value: float = 1e-5,
                 rtol_vector: float = 1e-4,
                 atol_vector: float = 1e-5, lin_rtol: float = 1e-8,
                 lin_atol=1e-8, max_iter=np.inf, delta=0.04):
        self.max_iter = max_iter
        self.lin_atol = lin_atol
        self.lin_rtol = lin_rtol
        self.iter = 0
        self.A = A
        self.M = M
        self.k = k
        self.v = v if v is not None else np.ones(A.shape[0])
        self.tol_value = tol_value
        self.atol_vector = atol_vector
        self.rtol_vector = rtol_vector
        self.delta = delta

    def __next__(self):
        k_old = self.k + self.delta
        u, info = \
            la.bicgstab(self.M - 1 / k_old * self.A, self.A @ self.v,
                        self.v * self.k,
                        atol=self.lin_atol,
                        rtol=self.lin_rtol)
        assert info == 0
        k_new = np.linalg.norm(u)
        u = u / k_new
        k = 1 / (1 / k_new + 1 / k_old)
        if not self.should_stop(k, u):
            self.k = k
            self.v = u
            return self.k, self.v
        else:
            raise StopIteration

    def should_stop(self, k, u):
        return (abs(k - self.k) < self.tol_value and np.allclose(u, self.v,
                                                                 rtol=self.rtol_vector,
                                                                 atol=self.atol_vector)) \
               or self.iter > self.max_iter

    def preconditioner(self, M: sparse.spmatrix):
        pre = la.spilu(M)

        def solve(x):
            return pre.solve(M @ x)

        return pre, la.LinearOperator(matvec=solve, shape=M.shape,
                                      dtype=np.float)

    def __iter__(self):
        return self
