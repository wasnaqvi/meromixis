"""Thermochemical equilibrium and quenching for the two diagnostic couples.

    N2  + 3 H2  <->  2 NH3          (the ammonia "factory")
    CH4 + H2O   <->  CO + 3 H2      (the CO discriminator)

Equilibrium constants come from fits to Delta_G(T) built from NASA polynomials (GRI-Mech
3.0 thermo, evaluated with Cantera; see ``scripts/dev/make_thermo_tables.py``).  The
fits reproduce the tabulated Delta_G to better than 1 kJ/mol over 300-3000 K, e.g.
Kp(steam reforming, 1000 K) = 26.5 bar^2.

Chemical timescales follow the Zahnle & Marley (2014) quench approximations.
VERIFY: the prefactors/exponents below were transcribed from memory of that paper
(their eqs. for t_CO and t_NH3); check them against the published text before any
number derived from them goes into a manuscript (docs/VERIFY.md, item C1).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .column import DRY, MOIST, Column
from .constants import R_GAS

# Delta_G(T) = a + b T + c T ln T + d T^2   [J mol-1], standard state 1 bar
GIBBS_FITS = {
    "CH4+H2O=CO+3H2": (2.00501780e05, 6.96948680e00, -3.57927883e01, 1.25312219e-02),
    "N2+3H2=2NH3": (-8.70396700e04, 1.74447710e-01, 3.19361608e01, -9.95345957e-03),
}


def delta_g(reaction: str, T):
    a, b, c, d = GIBBS_FITS[reaction]
    T = np.asarray(T, float)
    return a + b * T + c * T * np.log(T) + d * T**2


def k_eq(reaction: str, T):
    """Equilibrium constant in bar units (Delta_n = +2 for CO formation, -2 for NH3)."""
    return np.exp(-delta_g(reaction, T) / (R_GAS * np.asarray(T, float)))


def nh3_equilibrium(T, p, x_n, x_h2):
    """Equilibrium NH3 mole fraction given total N (as atoms) mole fraction ``x_n``.

    Solves x_NH3^2 = Kp p^2 x_H2^3 x_N2 with x_N = x_NH3 + 2 x_N2 (p in Pa).
    """
    p_bar = np.asarray(p, float) / 1.0e5
    b = 0.5 * k_eq("N2+3H2=2NH3", T) * p_bar**2 * x_h2**3
    # stable form of (-b + sqrt(b^2 + 4 b x_n)) / 2
    return 2.0 * b * x_n / (b + np.sqrt(b * b + 4.0 * b * x_n))


def co_equilibrium(T, p, x_c, x_o, x_h2):
    """Equilibrium CO mole fraction (CO2 neglected) given total C and O mole fractions."""
    p_bar = np.asarray(p, float) / 1.0e5
    r = k_eq("CH4+H2O=CO+3H2", T) / (p_bar**2 * x_h2**3)
    bb = r * (x_c + x_o) + 1.0
    disc = np.sqrt(np.maximum(bb * bb - 4.0 * r * r * x_c * x_o, 0.0))
    return 2.0 * r * x_c * x_o / (bb + disc)


def t_chem_nh3(T, p):
    """NH3 <-> N2 interconversion time, s (Zahnle & Marley 2014).  p in Pa.  VERIFY C1."""
    return 1.0e-7 / (np.asarray(p, float) / 1.0e5) * np.exp(52000.0 / np.asarray(T, float))


def t_chem_co(T, p, metallicity=1.0):
    """CO <-> CH4 interconversion time, s (Zahnle & Marley 2014).  p in Pa.  VERIFY C1."""
    return 1.5e-6 / (np.asarray(p, float) / 1.0e5) * metallicity**-0.7 * np.exp(42000.0 / np.asarray(T, float))


def kzz_mixing_length(col: Column):
    """Free-convection estimate K = w H, w = (R F / (mu rho cp))^(1/3), with F = sigma Tint^4.

    Only meaningful in convective zones.  m2 s-1.
    """
    from .constants import SIGMA_SB

    g = col.params.gravity
    flux = SIGMA_SB * col.params.tint**4
    rho = col.p * col.mu / (R_GAS * col.T)
    cp = col.comp.cp_molar_mix(col.T, col.x) / col.mu
    h = R_GAS * col.T / (col.mu * g)
    w = (R_GAS / col.mu * flux / (rho * cp)) ** (1.0 / 3.0)
    return w * h


@dataclass(frozen=True)
class Quench:
    species: str
    p: float  # Pa
    T: float  # K
    x: float  # quenched mole fraction
    t_mix: float  # s
    in_domain: bool  # False if the quench level lies below the model bottom


def quench(col: Column, species: str, k_deep: float | None = None) -> Quench:
    """Quench level and abundance of 'NH3' or 'CO' in the deep convective zone.

    ``k_deep`` (m2 s-1) is the mixing in the deep convective zone; default is the
    mixing-length estimate.  The mixing time is H^2 / K.
    """
    comp = col.comp
    xd = comp.x_deep
    h = R_GAS * col.T / (col.mu * col.params.gravity)
    k = kzz_mixing_length(col) if k_deep is None else np.full_like(col.T, k_deep)
    t_mix = h**2 / k
    if species == "NH3":
        t_ch = t_chem_nh3(col.T, col.p)
        x_eq = nh3_equilibrium(col.T, col.p, xd["NH3"], xd["H2"])
    elif species == "CO":
        t_ch = t_chem_co(col.T, col.p, comp.metallicity)
        x_eq = co_equilibrium(col.T, col.p, xd["CH4"], xd["H2O"], xd["H2"])
    else:
        raise ValueError("species must be 'NH3' or 'CO'")
    # search only the deep convective zone (contiguous convective levels above the bottom)
    conv = (col.regime == DRY) | (col.regime == MOIST)
    i_bot = len(col.p) - 1
    i = i_bot
    while i > 0 and conv[i - 1]:
        i -= 1
    sl = slice(i, i_bot + 1)
    d = np.log(t_ch[sl]) - np.log(t_mix[sl])  # > 0: chemistry slower than mixing (quenched)
    lp = np.log(col.p[sl])
    # Drop non-finite levels before searching for the crossing.  They occur when the
    # deep moist adiabat reaches the water critical point inside the domain, and both
    # guards below would silently fall through on NaN (every NaN comparison is False),
    # leaving the crossing search to index an empty array.  A column that had to be
    # truncated cannot be called converged, so in_domain is forced False.
    good = np.isfinite(d) & np.isfinite(lp) & np.isfinite(col.T[sl]) & np.isfinite(x_eq[sl])
    if not good.any():
        return Quench(species, float("nan"), float("nan"), float("nan"), float("nan"), False)
    truncated = not bool(good.all())
    d, lp = d[good], lp[good]
    t_sl, x_sl, tm_sl = col.T[sl][good], x_eq[sl][good], t_mix[sl][good]
    if d[-1] > 0:  # still quenched at the deepest usable level
        return Quench(species, float(np.exp(lp[-1])), float(t_sl[-1]), float(x_sl[-1]), float(tm_sl[-1]), False)
    if d[0] <= 0:  # chemistry fast throughout the convective zone: quench at its top
        return Quench(species, float(np.exp(lp[0])), float(t_sl[0]), float(x_sl[0]), float(tm_sl[0]), not truncated)
    j = int(np.nonzero(d <= 0)[0][0])  # first (shallowest) equilibrated level
    w = d[j - 1] / (d[j - 1] - d[j])
    lpq = lp[j - 1] + w * (lp[j] - lp[j - 1])
    tq = np.interp(lpq, lp, t_sl)
    xq = np.exp(np.interp(lpq, lp, np.log(x_sl)))
    return Quench(species, float(np.exp(lpq)), float(tq), float(xq), float(np.interp(lpq, lp, tm_sl)), not truncated)
