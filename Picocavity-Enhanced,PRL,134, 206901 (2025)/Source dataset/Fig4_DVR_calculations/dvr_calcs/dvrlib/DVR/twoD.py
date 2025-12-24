#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   twoD.py
@Time    :   2024/04/11 15:39:55
@Author  :   George Trenins
@Desc    :   Discrete variable representation approach to solving the Schroedinger
equation in two-dimensional systems.
'''


from __future__ import print_function, division, absolute_import
import numpy as np
import math
from scipy import sparse
from dvrlib import special
from dvrlib.DVR import oneD
from typing import Any, Optional, Union, Tuple


class OpenInterval(oneD.OpenInterval):

    def __init__(self, 
                 PES : Any,
                 lb : tuple[float, float],
                 ub : tuple[float, float],
                 dq : float,
                 E_cut : Optional[float] = None,
                 mass : Optional[Union[float, tuple[float, float]]] = None,
                 ):
        
        self.PES = PES
        self.UNITS = PES.UNITS
        self.hbar = self.UNITS.hbar
        if mass is None:
            self.mx = 1.0
            self.my = 1.0
        elif np.isscalar(mass):
            self.mx = mass
            self.my = mass
        else:
            self.mx, self.my = mass
        self.dq = dq    # mass-weighted grid spacing
        self.hh = [dq / math.sqrt(m) for m in (self.mx, self.my)]
        self.E_cut = np.inf if E_cut is None else E_cut
        self.grid = []
        self.idx = []
        self.KE = []
        for l, u, hh in zip(lb, ub, self.hh):
            grid, idx, KE = self.make_grid(l, u, hh)
            self.grid.append(grid)
            self.idx.append(idx)
            self.KE.append(KE)
        self.potgrid, self.mask = self.make_potgrid()
        self.H = self.make_hamiltonian()

    def make_potgrid(self) -> tuple[np.ndarray, list[bool]]:
        potgrid : list[float] = []
        mask : list[bool] = []
        for x in self.grid[0]:
            for y in self.grid[1]:
                V = self.PES.potential([x,y])
                potgrid.append(V)
                mask.append(V < self.E_cut)
        return np.asarray(potgrid), mask
    
    
    @staticmethod
    def reconstruct_function(
        arr: np.ndarray,
        shape: tuple[int],
        mask: np.ndarray):
        """Reverse the operation arr = np.ravel(ans[mask]) where ans.shape = shape
        """
        ans = np.zeros(shape)
        flat = np.reshape(ans, -1)
        j = 0
        for k in range(np.prod(shape)):
            if mask[k]:
                flat[k] = arr[j]
                j += 1
        return ans

    def ke_matrix(self, idx : np.ndarray) -> np.ndarray :
        diff = idx[:,None] - idx[None,:]
        sgn = 1 - 2*np.mod(diff,2)
        diag = np.einsum('ii->i', diff)
        diag[:] = 1 # avoid division by zero
        ans = 2/diff**2
        diag = np.einsum('ii->i', ans)
        diag[:] = np.pi**2/3
        ans *= self.hbar**2 * sgn / (2*self.dq**2)
        return ans
    
    @staticmethod
    def orthogonal_build(mask, nx, ny, block_x, block_y, dtype=float):
        N = np.count_nonzero(mask)
        mat = sparse.dok_array((N, N), dtype=dtype)
        row = 0
        i_outer = 0
        for row_x in range(nx):
            for row_y in range(ny):
                i_outer += 1
                if not mask[i_outer-1]: continue
                col = 0
                i_inner = 0
                for col_x in range(nx):
                    for col_y in range(ny):
                        i_inner += 1
                        if not mask[i_inner-1]: continue
                        if row_x == col_x:
                            mat[row, col] += block_y[row_y, col_y]
                        if row_y == col_y:
                            mat[row, col] += block_x[row_x, col_x]
                        col += 1
                row += 1
        return mat
    
    def make_hamiltonian(self):
        # Get KE part
        H = self.orthogonal_build(
            self.mask,
            len(self.grid[0]),
            len(self.grid[1]),
            self.KE[0],
            self.KE[1]
            )
        # Add PE
        V_mask = self.potgrid[self.mask]
        for i, V in enumerate(V_mask):
            H[i, i] += V
        return H
    
    def dkdxk(self, k=1):
        """Return the matrix elements of the operator d^k/dx^k (k = 1 or 2) 
        in the DVR basis.
        """
        ans = []
        for idx, hh in zip(self.idx, self.hh):
            grid_diff = idx[:,None] - idx[None,:]
            ans.append(
                special.sinc(np.pi*grid_diff, nu=k) * (np.pi/hh)**k
            )
        return ans
    
    def integrate(self, integrand, *args, **kwargs):
        raise NotImplementedError