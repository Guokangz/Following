#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   onedim_model.py
@Time    :   2025/05/14 10:10:02
@Author  :   George Trenins
@Desc    :   Solve the one-dimensional Schrödinger equation for the vibrationally adiabatic potential Eq. (S4) of the SM.
'''


from __future__ import print_function, division, absolute_import
from dvrlib.PES._base import BasePES, check_ndims
from typing import Optional, Union
from scipy.interpolate import PPoly, CubicSpline
from dvrlib import units
from dvrlib.utils import add_boundary_knots
from pathlib import Path
import argparse
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).parent))
from constants import freq_scaling, mH2, mD2, UNITS


class Model1DPES(BasePES):
    ndim = 1

    def __init__(self, spl: PPoly, mass: Optional[float] = 1.0, UNITS: Optional[str] = "atomic") -> None:
        """Model one-dimensional PES based on a spline interpolation of numerical data.

        Args:
            spl (PPoly): _description_
            mass (Optional[float], optional): _description_. Defaults to 1.0.
            UNITS (Optional[str], optional): _description_. Defaults to "atomic".
        """
        self.mass: float = mass
        self._spline = spl
        self.UNITS: units.SI = getattr(units, UNITS)()

    @classmethod
    def from_csv(
        cls,
        mass: float,
        transpot: Union[str, Path],
        vibfreq: Optional[Union[str, Path]] = None,
        UNITS: Optional[str] = "atomic",
        bc_type: Optional[str] = "natural",
        n: Optional[int] = 0
        ) -> "Model1DPES":
        """Create a model 1D PES from numerical data in CSV format.

        Args:
            mass (float): mass of the molecule (either H2 or D2)
            transpot (str): Path to a two-column CSV file storing the translational potential energy as a function of the translational coordinate. The positions must be in angstrom and the energies in eV.
            vibfreq (str): Path to a two-column CSV file storing the vibrational frequencies of H2 as a function of the translational coordinate. The positions must be in angstrom and the frequencies in wave numbers. If not specified, the bare translational potential is used to initialise the model.
            UNITS (Optional[str], optional): name of unit system. Defaults to "atomic".
            bc_type (Optional[str], optional): boundary conditions for cubic spline interpolation. Defaults to "natural".
            n (Optional[int], optional): number of quanta in the vibrational mode [υ in Eq. (S3)]. Defaults to 0.

        Returns:
            Model1DPES
        """
        u: units.SI = getattr(units, UNITS)()
        # Translational data
        potpos, pot = np.loadtxt(transpot, unpack=True)
        potpos *= u.str2base('1 aa')
        pot *= u.str2base('1 eV')
        V_spl = CubicSpline(potpos, pot, bc_type=bc_type)
        add_boundary_knots(V_spl)
        if vibfreq is None:
            pot_spl = V_spl
        else:
            # Vibrational data
            freqpos, freq = np.loadtxt(vibfreq, unpack=True)
            # convert positions
            freqpos *= u.str2base('1 aa')
            # convert wavenumbers to frequency units
            freq *= u.wn2omega(1.0) * freq_scaling / np.sqrt(mass/mH2)
            w_spl = CubicSpline(freqpos, freq, bc_type=bc_type)
            add_boundary_knots(w_spl)
            # total energy for specified number of vibrational quanta
            W = pot + u.hbar * w_spl(potpos) * (n + 0.5)
            pot_spl = CubicSpline(potpos, W, bc_type=bc_type)
            add_boundary_knots(pot_spl)
        return Model1DPES(pot_spl, mass, UNITS=UNITS)

    @check_ndims
    def potential(self, x):
        return self._spline(x[...,0])
    
def get_input_paths(dtip: float, density: str) -> tuple[Path, Path]:
    datadir = Path(__file__).parent.resolve().parent / "data"
    density_map = {"low": "006ML", "high": "069ML"}
    try:
        subdir = density_map[density]
    except KeyError:
        raise ValueError(f"Unknown option {density = }.")
    datadir = datadir / subdir
    transpot = datadir / f"{dtip:.1f}A-energy-shifted.dat"
    if not transpot.is_file():
        raise ValueError((f"Translational potential data for surface-tip separation {dtip = } unavailable."))
    vibfreq = datadir/f"{dtip:.1f}A-freq-shifted.dat"
    if not vibfreq.is_file():
        raise ValueError((f"Vibrational frequency data for surface-tip separation {dtip = } unavailable."))
    return transpot, vibfreq
    
    
def main(dtip: float, isotope: str, density: str, dvr_json: str) -> None:
    import json
    from dvrlib.DVR.oneD import OpenInterval as DVR1D
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    # compute ZPE and mean surf-mol distance
    transpot, vibfreq = get_input_paths(dtip, density)
    mass: float = mH2 if isotope == "H" else mD2
    tPES: Model1DPES = Model1DPES.from_csv(mass, transpot, UNITS=UNITS, bc_type="natural")
    # Find the relevant (meta)stable minimum
    fun = lambda x: tPES._spline(x, nu=0)
    grad = lambda x: tPES._spline(x, nu=1)
    x0 = tPES.UNITS.str2base('4.0 aa')
    res = minimize(fun, x0, jac=grad, method="Newton-CG")
    xmin = res.x
    Vmin = fun(xmin)
    # Shift the PES in place
    tPES._spline.c[3, :] -= Vmin
    # Solve the DVR
    with open(dvr_json, 'r') as f:
        dvr_kwargs = json.load(f)
    solver = DVR1D(tPES, mass=mass, **dvr_kwargs)
    E, v = solver.solve(4, E0=0.0)
    # Determine the mean surface-molecule distance for all states:
    xgrid = np.copy(solver.grid)
    xmean = solver.integrate(v.T**2 * xgrid)
    # At low converage for 9.4 A, the relevant minimum is metastable (Fig. S8)
    # This selects the correct eigenvalue index
    for i, x in enumerate(xmean):
        n_eig = i
        if x < tPES.UNITS.str2base('5.0 aa'):
            break
    print(f"ZPE = {E[n_eig]/tPES.UNITS.str2base('1 eV')*1000:.1f}", end="  ")
    print(f"<x> = {xmean[n_eig]:.2f}",  end="  ")
    # Now compute the vibrational transitions
    
    adPES0 = Model1DPES.from_csv(mass, transpot, vibfreq=vibfreq,
                                 UNITS=UNITS, bc_type="natural", n=0)
    solver0 = DVR1D(adPES0, mass=mass, **dvr_kwargs)
    E0, v0 = solver0.solve(4, E0=0.0)
    adPES1 = Model1DPES.from_csv(mass, transpot, vibfreq=vibfreq,
                                 UNITS=UNITS, bc_type="natural", n=1)
    solver1 = DVR1D(adPES1, mass=mass, **dvr_kwargs)
    E1, v1 = solver1.solve(4, E0=0.0)

    # fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True)
    # scale=0.01
    # ax0.plot(solver0.grid, solver0.potgrid, 'k-', lw=2)
    # for e0, psi in zip(E0, v0.T):
    #     ax0.plot(solver0.grid, psi*scale+e0)
    # amplitude = E0[-1] - E0[0]
    # ax0.set_ylim([E0[0] - 0.75*amplitude, E0[-1] + 0.75*amplitude])
    # ax1.plot(solver1.grid, solver1.potgrid, 'k-', lw=2)
    # for e1, psi in zip(E1, v1.T):
    #     ax1.plot(solver1.grid, psi*scale+e1)
    # amplitude = E1[-1] - E1[0]
    # ax1.set_ylim([E1[0] - 0.75*amplitude, E1[-1] + 0.75*amplitude])
    # plt.show()

    # Indices of lowest-energy states localized in the relevant wells (n0, n1)
    idx_dict = {
        'H' : {
            'low' : {'11.4' : (0,0),
                     '9.4' : (1,0),
                     '7.9' : (0,0),
                     '7.4' : (0,0)},
            'high': {'11.4' : (0,0),
                     '9.4' : (0,0),
                     '7.9' : (0,0),
                     '7.4' : (0,0)},
            },
        'D' : {
            'low' : {'11.4' : (0,0),
                     '9.4' : (1,1),
                     '7.9' : (0,0),
                     '7.4' : (0,0)},
            'high': {'11.4' : (0,0),
                     '9.4' : (0,0),
                     '7.9' : (0,0),
                     '7.4' : (0,0)},
            },
    }
    n0, n1 = idx_dict[isotope][density][f"{dtip:.1f}"]
    print(f"ν (cm−1 ) = {tPES.UNITS.energy2wn(E1[n1]-E0[n0]):.2f}")


    

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the zero-point energy (relative to the bottom of the well) and the mean molecule-surface distance for the translational part of the PES + the vibrational transition frequency for the 1D vibrationally adiabatic PES.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dtip", type=float, help="tip-surface separation in angstrom", choices=[7.4, 7.9, 9.4, 11.4])
    parser.add_argument("dvr", type=str, help="path to JSON file with arguments for initialising dvrlib.DVR.oneD.OpenInterval, *except* for `PES` and `mass`")
    parser.add_argument("--isotope", type=str, default="H", choices=["H", "D"])
    parser.add_argument("--density", type=str, default="low", choices=["low", "high"])
    args = parser.parse_args()
    main(args.dtip, args.isotope, args.density, args.dvr)


