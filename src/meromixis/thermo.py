"""Water saturation thermodynamics and adiabatic gradients for a non-dilute H2-H2O gas.

Two saturation curves are provided:

* ``"wp"`` (default): Murphy & Koop (2005) over ice below the triple point and the
  Wagner & Pruss (2002, IAPWS-95) auxiliary equation over liquid up to the critical
  point.
* ``"cc"``: Clausius-Clapeyron with constant latent heat, the functional form the
  Generic PCM uses for generic condensable species
  (``psat = Pref * exp(-dH/R * (1/T - 1/Tref))``).  Use it when comparing against PCM
  output so the emulator and the GCM agree on what "saturated" means.

Everything downstream needs only ``psat(T)`` and ``gamma(T) = dln psat / dln T``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .composition import Composition, cp_molar
from .constants import (
    MOLAR_MASS,
    P_CRIT_H2O,
    P_TRIPLE_H2O,
    R_GAS,
    T_CRIT_H2O,
    T_TRIPLE_H2O,
)

M_V = MOLAR_MASS["H2O"]
_WP = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)


def _lnpsat_wp_scalar(T: float) -> float:
    """ln psat [Pa]: ice (Murphy & Koop 2005) below 273.16 K, liquid (Wagner & Pruss) above."""
    if T >= T_CRIT_H2O:
        return math.inf
    if T < T_TRIPLE_H2O:
        return 9.550426 - 5723.265 / T + 3.53068 * math.log(T) - 0.00728332 * T
    th = 1.0 - T / T_CRIT_H2O
    s = (
        _WP[0] * th
        + _WP[1] * th**1.5
        + _WP[2] * th**3
        + _WP[3] * th**3.5
        + _WP[4] * th**4
        + _WP[5] * th**7.5
    )
    return math.log(P_CRIT_H2O) + (T_CRIT_H2O / T) * s


@dataclass(frozen=True)
class Saturation:
    """Saturation curve of the condensable (water).

    kind : "wp" or "cc"
    latent_heat, t_ref, p_ref : only used for kind="cc" (J kg-1, K, Pa).
    """

    kind: str = "wp"
    latent_heat: float = 2.5e6
    t_ref: float = T_TRIPLE_H2O
    p_ref: float = P_TRIPLE_H2O

    def __post_init__(self):
        if self.kind not in ("wp", "cc"):
            raise ValueError("Saturation.kind must be 'wp' or 'cc'")

    # -- scalar fast paths (used inside the column integrator) -------------------
    def lnpsat_s(self, T: float) -> float:
        if self.kind == "cc":
            return math.log(self.p_ref) - (self.latent_heat * M_V / R_GAS) * (1.0 / T - 1.0 / self.t_ref)
        return _lnpsat_wp_scalar(T)

    def psat_s(self, T: float) -> float:
        v = self.lnpsat_s(T)
        return math.inf if v == math.inf else math.exp(v)

    def gamma_s(self, T: float) -> float:
        """dln psat / dln T.  Returns 0 above the critical point (no condensation)."""
        if self.kind == "cc":
            return self.latent_heat * M_V / (R_GAS * T)
        if T >= T_CRIT_H2O - 0.05:
            return 0.0
        h = 1e-4 * T
        return (self.lnpsat_s(T + h) - self.lnpsat_s(T - h)) / (math.log(T + h) - math.log(T - h))

    # -- array API ---------------------------------------------------------------
    def lnpsat(self, T):
        """ln psat [Pa] for arrays; +inf above the critical point."""
        T = np.asarray(T, dtype=float)
        if self.kind == "cc":
            return np.log(self.p_ref) - (self.latent_heat * M_V / R_GAS) * (1.0 / T - 1.0 / self.t_ref)
        tc = np.clip(T, 1.0, T_CRIT_H2O)
        ice = 9.550426 - 5723.265 / tc + 3.53068 * np.log(tc) - 0.00728332 * tc
        th = 1.0 - tc / T_CRIT_H2O
        liq = np.log(P_CRIT_H2O) + (T_CRIT_H2O / tc) * (
            _WP[0] * th + _WP[1] * th**1.5 + _WP[2] * th**3 + _WP[3] * th**3.5 + _WP[4] * th**4 + _WP[5] * th**7.5
        )
        out = np.where(T < T_TRIPLE_H2O, ice, liq)
        return np.where(T >= T_CRIT_H2O, np.inf, out)

    def psat(self, T):
        return np.exp(self.lnpsat(T))

    def gamma(self, T):
        """dln psat / dln T for arrays; 0 above the critical point."""
        T = np.asarray(T, dtype=float)
        if self.kind == "cc":
            return self.latent_heat * M_V / (R_GAS * T)
        ok = T < T_CRIT_H2O - 0.05
        ts = np.where(ok, T, 300.0)
        h = 1e-4 * ts
        g = (self.lnpsat(ts + h) - self.lnpsat(ts - h)) / (np.log(ts + h) - np.log(ts - h))
        return np.where(ok, g, 0.0)

    def latent_heat_cc(self, T):
        """Latent heat implied by Clausius-Clapeyron for an ideal vapour, J kg-1."""
        return self.gamma(T) * R_GAS * np.asarray(T, dtype=float) / M_V

    def x_sat(self, T, p):
        """Saturation mole fraction of water, capped at 1."""
        return np.minimum(1.0, self.psat(T) / np.asarray(p, dtype=float))


# ---- adiabatic gradients (dln T / dln p) -------------------------------------------

def dry_gradient(T: float, x_h2o: float, comp: Composition) -> float:
    """Adiabatic gradient of an unsaturated dry-gas + vapour mixture."""
    return R_GAS / float(comp.cp_molar_mix(T, x_h2o))


def moist_gradient(T: float, p: float, comp: Composition, sat: Saturation) -> float:
    """Saturated pseudo-adiabatic gradient dlnT/dlnp, valid for any vapour fraction.

    Derived from the first law for an ideal dry gas plus saturated ideal vapour with
    the condensate removed (Pierrehumbert 2010, eq. 2.33), converted from the dry
    partial pressure to total pressure.  Limits: R/cp_dry when psat << p, and 1/gamma
    (the vapour-pressure curve) when the dry gas vanishes.
    """
    ps = sat.psat_s(T)
    gam = sat.gamma_s(T)
    if gam <= 0.0:  # supercritical: no condensation
        return dry_gradient(T, 0.0, comp)
    if ps >= p:
        return 1.0 / gam
    pa = p - ps
    mu_d = comp.mu_dry
    r_a = R_GAS / mu_d
    r_c = R_GAS / M_V
    cp_a = float(comp.cp_dry_molar(T)) / mu_d
    cp_c = float(cp_molar("H2O", T)) / M_V
    latent = gam * r_c * T
    r = (M_V / mu_d) * ps / pa  # mass mixing ratio, vapour / dry gas
    grad_a = (r_a / cp_a) * (1.0 + latent * r / (r_a * T)) / (
        1.0 + r * (cp_c / cp_a + (latent / (cp_a * T)) * (latent / (r_c * T) - 1.0))
    )
    return grad_a / (pa / p + (ps / p) * gam * grad_a)
