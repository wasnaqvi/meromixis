#!/usr/bin/env python
"""Regenerate the thermochemical constants hard-coded in the package.

Needs Cantera (``pip install -e ".[thermo]"``).  Prints

* the cp(T) table in ``meromixis.constants.CP_MOLAR`` and
* the Delta_G(T) fit coefficients in ``meromixis.chemistry.GIBBS_FITS``

from the NASA polynomials shipped with GRI-Mech 3.0.  Paste the output back into the
source if you change the temperature grid or the species list.
"""
import cantera as ct
import numpy as np

gas = ct.Solution("gri30.yaml")
T_GRID = [200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000]


def g0(species: str, T: float) -> float:
    s = gas.species(species).thermo
    return (s.h(T) - T * s.s(T)) / 1000.0  # J mol-1


print("CP_T_GRID =", tuple(float(t) for t in T_GRID))
for sp in ["H2", "CH4", "NH3", "H2O", "N2", "CO"]:
    th = gas.species(sp).thermo
    print(f'    "{sp}": {tuple(round(th.cp(T) / 1000.0, 3) for T in T_GRID)},')

REACTIONS = {
    "CH4+H2O=CO+3H2": ({"CO": 1, "H2": 3}, {"CH4": 1, "H2O": 1}),
    "N2+3H2=2NH3": ({"NH3": 2}, {"N2": 1, "H2": 3}),
}
Ts = np.linspace(300.0, 3000.0, 271)
basis = np.vstack([np.ones_like(Ts), Ts, Ts * np.log(Ts), Ts**2]).T
for name, (prod, reac) in REACTIONS.items():
    dg = np.array([sum(n * g0(s, T) for s, n in prod.items()) - sum(n * g0(s, T) for s, n in reac.items()) for T in Ts])
    coef, *_ = np.linalg.lstsq(basis, dg, rcond=None)
    print(f'    "{name}": {tuple(float(c) for c in coef)},   # max |error| = {np.abs(basis @ coef - dg).max():.0f} J/mol')
