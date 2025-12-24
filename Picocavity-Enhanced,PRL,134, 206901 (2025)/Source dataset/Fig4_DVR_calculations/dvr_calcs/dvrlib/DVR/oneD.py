#!/usr/bin/env python
"""
Created on Sat Jun  2 15:33:19 2018

@author: gt317

Discrete variable representation approach to solving the Schroedinger
equation in low-dimensional systems.

"""

from __future__ import division, print_function
import numpy as np
from dvrlib import special, units
from scipy import sparse
from scipy.sparse import linalg as SLA
from typing import Any, Optional

class OpenInterval(object):

    def __init__(self, PES: Any, lb: float, ub: float, dq: float, 
                 E_cut: Optional[float] = None, mass: Optional[float] = None):
            
        """Solve the one-dimensional Schroedinger equation using the
        discrete variable representation by Miller et al. for a potential
        on the open (–∞, +∞) interval. See DOI 10.1063/1.462100

        Args:
            PES (object): the potential energy surface
            lb (float): lower bound of the x-coordinate grid
            ub (float): upper bound of the x-coordinate grid
            dq (float): mass-weighted grid spacing, Δx⋅√m
            E_cut (float, optional): energy cut-off for the grid, default is None.
            mass (float, optional): particle mass (by default read from PES)

        Notes:
            All quantities are in the base units of the system `PES.UNITS`
        """

        self.PES = PES
        self.UNITS: units.SI = PES.UNITS
        self.hbar: float = self.UNITS.hbar
        if mass is None:
            self.m: float = self.PES.mass
        else:
            self.m: float = mass
        self.dq: float = dq   
        self.hh: float = dq/np.sqrt(self.m)
        self.E_cut: float = np.inf if E_cut is None else E_cut
        self.grid: np.ndarray
        self.idx: np.ndarray
        self.KE: np.ndarray
        self.grid, self.idx, self.KE = self.make_grid(lb, ub, self.hh)
        self.potgrid: np.ndarray
        self.mask: list[bool]
        self.potgrid, self.mask = self.make_potgrid()
        self.H: sparse.dok_array = self.make_hamiltonian()

    def make_grid(self, lb: float, ub: float, hh: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Construct the position grid and the corresponding kinetic energy matrix
        [Eq. (2.1) of Colbert and Miller JCP 96, 1982-1991 (1992), http://aip.scitation.org/doi/10.1063/1.462100]

        Args:
            lb (float): lower bound on the grid
            ub (float): upper bound on the grid
            hh (float): grid step size

        Returns:
            grid (ndarray): position grid
            idx (ndarray): index array
            KE (ndarray): kinetic energy matrix
        """
        grid = np.arange(lb, ub+hh, hh)
        idx = np.arange(len(grid), dtype=int)
        KE = self.ke_matrix(idx)
        return grid, idx, KE


    def ke_matrix(self, idx : np.ndarray) -> np.ndarray :
        """Eq. (2.1) of Colbert and Miller JCP 96, 1982-1991 (1992), http://aip.scitation.org/doi/10.1063/1.462100"""
        diff = idx[:,None] - idx[None,:]
        sgn = 1 - 2*np.mod(diff,2)
        diag = np.einsum('ii->i', diff)
        diag[:] = 1 # avoid division by zero
        ans = 2/diff**2
        diag = np.einsum('ii->i', ans)
        diag[:] = np.pi**2/3
        ans *= self.hbar**2 * sgn / (2*self.dq**2)
        return ans
            
    def make_potgrid(self) -> tuple[np.ndarray, list[bool]]:
        potgrid : list[float] = []
        mask : list[bool] = []
        for x in self.grid:
            V = self.PES.potential([x])
            potgrid.append(V)
            mask.append(V < self.E_cut)
        return np.asarray(potgrid), mask

    def make_hamiltonian(self) -> sparse.dok_array:
        N = np.count_nonzero(self.mask)
        H = sparse.dok_array((N, N), dtype=float)
        V_mask = self.potgrid[self.mask]
        for i, V in enumerate(V_mask):
            H[i, i] = V
        i_row = -1
        for row, _ in enumerate(self.grid):
            if not self.mask[row]:
                # This row has been masked out
                continue
            i_row += 1
            i_col = -1
            for col, _ in enumerate(self.grid):
                if not self.mask[col]:
                    # This column has been masked out
                    continue
                i_col += 1
                # all hurdles cleared
                if (row == col):
                    # Diagonal already has potential
                    H[i_row, i_col] = H[i_row, i_col] + self.KE[row, col]
                else:
                    H[i_row, i_col] = self.KE[row, col]
        return H
    
    def dkdxk(self, k=1):
        """Return the matrix elements of the operator d^k/dx^k (k = 1 or 2) 
        in the DVR sinc basis. See lecture notes by Gerrit Groenenboom on
        https://www.acmm.nl/molsim/han/2001/gcg2001.pdf
        """
        grid_diff = self.idx[:,None] - self.idx[None,:]
        return special.sinc(np.pi*grid_diff, nu=k) * (np.pi/self.hh)**k
    
    def solve(self, 
              k : int, 
              E0 : Optional[float] = None,
              which : str = 'LA',
              mode='normal',
              **kwargs
              ) -> tuple[np.ndarray, np.ndarray] :
        """Compute the numerical solutions of the Schroedinger equation.

        Args:
            k (int): number of eigenvalues to compute
            E0 (optional): energy around which to search for eigenvalues. Defaults to None.
            which (str, optional): what range of eigenvalues to explore. Defaults to 'LA' (largest algebraic).
            mode (str, optional): strategy for the shift-invert mode. Defaults to 'normal'.

        Returns:
            w (np.ndarray): 1d array of k eigenvalues
            v (np.ndarray): 2d array, v[:,i] contains the eigenvector corresponding to the i-th eigenvalue

        Notes:
            see documentation for `scipy.sparse.linalg.eigsh` for explanation of `which`, `mode` and `**kwargs`.
        """
        w, v = SLA.eigsh(self.H, k=k, sigma=E0, which=which, mode=mode, **kwargs)
        v /= np.sqrt(np.prod(self.hh))
        return w, v

    def integrate(self, integrand, *args, **kwargs):
        """Compute a numerical integral in accordance to the DVR quadrature rule
        """
        return np.sum(self.hh*integrand, axis=kwargs.get('axis', -1))
    
    def dvr_to_eig(self, dvr_me: np.ndarray, eig: np.ndarray) -> np.ndarray:
        """Convert a matrix element in the DVR basis to the eigenbasis. The matrix element
        must already have been "masked out", with entries corresponding
        to high-energy parts of space deleted from the input

        Args:
            dvr_me (np.ndarray): matrix element in DVR representation
            eig (np.ndarray): eigenvector array (eigenvectors stored as columns)

        Returns:
            eigbasis_me (np.ndarray) 
        """

        return np.linalg.multi_dot([eig.T, dvr_me, eig])
    

if __name__ == "__main__":
    
    from dvrlib.PES import polynomial
    from dvrlib.DVR.oneD import OpenInterval
    import matplotlib.pyplot as plt
    
    PES = polynomial.PES(coeffs=[0.0, 0.0, 0.0, 0.0, 0.25], UNITS="atomic")
    lb, ub = -7, 7
    dq = 0.025
    solver = OpenInterval(PES, lb, ub, dq, mass=PES.mass)
    evals, evecs = solver.solve(k=5, which='SA')
    print(solver.integrate(evecs[:,0]*evecs[:,0]))
    
    plt.plot(solver.grid, solver.potgrid, 'k-')
    for i,E in enumerate(evals):
        plt.plot(solver.grid, evecs[:,i] + E)
    plt.xlim(-4,4)
    plt.ylim(-1,8)
    plt.show()
