#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   onedim_model.py
@Time    :   2025/05/14 10:10:02
@Author  :   George Trenins
@Desc    :   Solve the two-dimensional Schrödinger equation for the 2D potential Eq. (S1) of the SM.
'''

from __future__ import print_function, division, absolute_import
from dvrlib.PES._base import BasePES, check_ndims
from typing import Optional, Union
from scipy.interpolate import PPoly, CubicSpline
from dvrlib import units
from dvrlib.utils import add_boundary_knots
from dvrlib.DVR.twoD import OpenInterval as DVR2D
from pathlib import Path
import argparse
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).parent))
from constants import freq_scaling, mH2, muH, mD2, muD, UNITS


class Model2DPES(BasePES):
    ndim = 2
    def __init__(self, tPES: PPoly, kspring: PPoly, mass: Optional[tuple[float,float]] = (1.0, 1.0), UNITS : Optional[str] = "atomic") -> None:
        # translational part 
        self.tPES = tPES
        # spring constant (depends on translational coordinate)
        self.kspring = kspring
        self.UNITS = getattr(units, UNITS)()
        self.mass = np.asarray(mass)

    @classmethod
    def from_csv(
        cls,
        mass: tuple[float,float],
        transpot: Union[str, Path],
        vibfreq: Optional[Union[str, Path]] = None,
        UNITS: Optional[str] = "atomic",
        bc_type: Optional[str] = "natural"
        ) -> "Model2DPES":
        """Create a model 2D PES from numerical data in CSV format.

        Args:
            mass (float, float): mass of the molecule and reduced vibrational mass (either H2 or D2)
            transpot (str): Path to a two-column CSV file storing the translational potential energy as a function of the translational coordinate. The positions must be in angstrom and the energies in eV.
            vibfreq (str): Path to a two-column CSV file storing the vibrational frequencies of H2 as a function of the translational coordinate. The positions must be in angstrom and the frequencies in wave numbers. If not specified, the bare translational potential is used to initialise the model.
            UNITS (Optional[str], optional): name of unit system. Defaults to "atomic".
            bc_type (Optional[str], optional): boundary conditions for cubic spline interpolation. Defaults to "natural".

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
        # Vibrational data
        freqpos, freq = np.loadtxt(vibfreq, unpack=True)
        # convert positions
        freqpos *= u.str2base('1 aa')
        # convert wavenumbers to frequency units
        freq *= u.wn2omega(1.0) * freq_scaling / np.sqrt(mass[0]/mH2)
        mu = mass[1]
        k = mu*freq**2
        k_spl = CubicSpline(freqpos, k, bc_type=bc_type)
        add_boundary_knots(k_spl)
        return Model2DPES(V_spl, k_spl, mass=mass, UNITS=UNITS)

    @check_ndims
    def potential(self, x):
        translational = self.tPES(x[...,0])
        k = self.kspring(x[...,0])
        harmonic = k/2 * x[...,1]**2
        return translational + harmonic


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

def plot_eigstates(res):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    fig, axarr = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(15,7))
    for ax, w0 in zip( np.ravel(axarr[:,:2]), res['W0'].T ):
        w0 = DVR2D.reconstruct_function(w0, res['X'].shape, res['mask'])
        ax.imshow(w0.T, extent=[res['X'].min(), res['X'].max(), res['Y'].min(), res['Y'].max()], 
                  origin='lower', cmap=cm.RdYlGn, vmax=np.abs(w0).max(), vmin=-np.abs(w0).max(), aspect=4)
        CS = ax.contour(res['X'], res['Y'], res['pot'], levels=np.linspace(0,0.2,6))
        ax.clabel(CS, inline=True, fontsize=10)
    for ax, w1 in zip( np.ravel(axarr[:,2:]), res['W1'].T ):
        w1 = DVR2D.reconstruct_function(w1, res['X'].shape, res['mask'])
        ax.imshow(w1.T, extent=[res['X'].min(), res['X'].max(), res['Y'].min(), res['Y'].max()], 
                  origin='lower', cmap=cm.RdYlGn, vmax=np.abs(w1).max(), vmin=-np.abs(w1).max(), aspect=4)
        CS = ax.contour(res['X'], res['Y'], res['pot'], levels=np.linspace(0,0.2,6))
        ax.clabel(CS, inline=True, fontsize=10)
    plt.show()


def main(dtip: float, isotope: str, density: str, dvr_json: str) -> None:
    import json
    transpot, vibfreq = get_input_paths(dtip, density)
    mass: float = (mH2, muH) if isotope == "H" else (mD2, muD)
    PES: Model2DPES = Model2DPES.from_csv(mass, transpot, vibfreq, UNITS=UNITS, bc_type="natural")
    u: units.SI = PES.UNITS
    omega = np.sqrt(PES.kspring(u.str2base("3.5 aa"))/mass[1])
    with open(dvr_json, 'r') as f:
        dvr_kwargs = json.load(f)
    solver = DVR2D(PES, mass=mass, **dvr_kwargs)
    show = False
    if show:
        k=4
    else:
        k=1
    E0, W0 = solver.solve(k=k, E0=u.hbar*omega/2, which='LA', mode='normal')
    E1, W1 = solver.solve(k=k, E0=3*u.hbar*omega/2, which='LA', mode='normal')
    x, y = solver.grid
    X, Y = np.meshgrid(x, y, indexing='ij')
    mask = np.ravel(solver.potgrid) < 2.0
    pot = np.where(mask, solver.potgrid, None).reshape((len(x), len(y)))
    res = dict(mass=mass, E0=E0, W0=W0, E1=E1, W1=W1, X=X, Y=Y, mask=solver.mask, pot=pot)
    if show:
        plot_eigstates(res)
    idx_dict = {
        'H' : {
            'high': {'9.4' : (0,0),
                     '7.4' : (0,0)},
            },
        'D' : {
            'high': {'9.4' : (0,0),
                     '7.4' : (0,0)},
            },
    }
    n0, n1 = idx_dict[isotope][density][f"{dtip:.1f}"]
    print(f"ν (cm−1 ) = {PES.UNITS.energy2wn(E1[n1]-E0[n0]):.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the zero-point energy (relative to the bottom of the well) and the mean molecule-surface distance for the translational part of the PES + the vibrational transition frequency for the 1D vibrationally adiabatic PES.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dtip", type=float, help="tip-surface separation in angstrom", choices=[7.4, 9.4])
    parser.add_argument("dvr", type=str, help="path to JSON file with arguments for initialising dvrlib.DVR.oneD.OpenInterval, *except* for `PES` and `mass`")
    parser.add_argument("--isotope", type=str, default="H", choices=["H", "D"])
    parser.add_argument("--density", type=str, default="high", choices=["high"])
    args = parser.parse_args()
    main(args.dtip, args.isotope, args.density, args.dvr)


