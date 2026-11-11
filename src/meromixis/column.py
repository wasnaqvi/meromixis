"""A fast 1-D radiative-convective emulator with condensation-inhibited convection.

Purpose: scoping, unit-testing the diagnostics, and a 1-D baseline to hold the 3-D
Generic PCM results against.  It is *not* a substitute for the PCM's correlated-k
column model (``rcm1d``): radiation is double-grey (Guillot 2010) and clouds have no
radiative effect.  Calibrate ``kappa_th``/``kappa_pexp``/``gamma_v`` against ``rcm1d``
with :func:`calibrate_grey` before trusting any number quantitatively.

Structure, integrating downward from the top:

* radiative where the radiative gradient is below the local adiabat;
* moist-convective (saturated pseudo-adiabat) where water condenses and the
  radiative gradient is super-adiabatic;
* **inhibited**: saturated, super-adiabatic, but ``q_sat > q_cri``; convection is shut
  off and the layer follows the (steep) radiative gradient (Leconte+2017);
* dry-convective below the cloud base.

``teq`` is the *effective* equilibrium temperature, i.e. it already includes the Bond
albedo: ``teq = (1 - A_B)**0.25 * Teq(A_B = 0)``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from functools import cached_property

import numpy as np

from .coldtrap import ColdTrap, cold_trap, water_profile
from .composition import Composition, x_to_q
from .constants import SIGMA_SB
from .inhibition import (
    StableLayer,
    altitude,
    brunt_vaisala_sq,
    find_stable_layer,
    number_density,
)
from .thermo import Saturation, dry_gradient, moist_gradient

RADIATIVE, MOIST, DRY, INHIBITED = 0, 1, 2, 3
REGIME_NAMES = {RADIATIVE: "radiative", MOIST: "moist-convective", DRY: "dry-convective", INHIBITED: "inhibited"}
_SQRT3 = math.sqrt(3.0)


@dataclass(frozen=True)
class ColumnParams:
    gravity: float = 12.43  # m s-2 (K2-18 b: 8.63 M_E, 2.61 R_E)
    teq: float = 255.0  # K, effective (albedo included)
    tint: float = 50.0  # K, internal temperature
    metallicity: float = 100.0  # x solar
    kappa_th: float = 1.0e-2  # m2 kg-1 thermal opacity at p_ref
    kappa_pexp: float = 1.0  # kappa_th ~ p**kappa_pexp (1 = collision-induced)
    p_ref: float = 1.0e5  # Pa
    gamma_v: float = 0.4  # visible / thermal opacity ratio
    flux_factor: float = 1.0  # local irradiation / global mean (4 at the substellar point)
    inhibition: bool = True
    p_top: float = 1.0  # Pa
    p_bot: float = 1.0e8  # Pa
    n_lev: int = 400
    n_sub: int = 4
    saturation: Saturation = field(default_factory=Saturation)

    @property
    def composition(self) -> Composition:
        return Composition(self.metallicity)


@dataclass
class Column:
    params: ColumnParams
    p: np.ndarray  # Pa, top -> bottom
    T: np.ndarray  # K
    regime: np.ndarray  # int codes
    grad: np.ndarray  # dlnT/dlnp actually followed

    @cached_property
    def comp(self) -> Composition:
        return self.params.composition

    @cached_property
    def x(self) -> np.ndarray:
        """Water mole fraction (rain-out profile)."""
        return water_profile(self.p, self.T, self.comp.x_h2o_deep, self.params.saturation)

    @cached_property
    def q(self) -> np.ndarray:
        return x_to_q(self.x, self.comp.mu_dry)

    @cached_property
    def mu(self) -> np.ndarray:
        return self.comp.mu(self.x)

    @cached_property
    def z(self) -> np.ndarray:
        """Altitude above the bottom level, m."""
        return altitude(self.p, self.T, self.x, self.comp, self.params.gravity)

    @cached_property
    def n(self) -> np.ndarray:
        return number_density(self.p, self.T)

    @cached_property
    def n2(self) -> np.ndarray:
        return brunt_vaisala_sq(self.p, self.T, self.x, self.comp, self.params.gravity)

    @cached_property
    def cold_trap(self) -> ColdTrap:
        return cold_trap(self.p, self.T, self.comp.x_h2o_deep, self.params.saturation, p_min=self.params.p_top * 3)

    @cached_property
    def stable_layer(self) -> StableLayer:
        """Layer where convection is actually suppressed (regime == INHIBITED)."""
        idx = np.nonzero(self.regime == INHIBITED)[0]
        geo = find_stable_layer(
            self.p, self.T, self.x, self.comp, self.params.gravity, self.params.saturation, self.comp.x_h2o_deep
        )
        if idx.size == 0:
            return replace(geo, exists=False)
        return geo

    def at(self, p: float, what: str = "T") -> float:
        """Log-pressure interpolation of a profile attribute at pressure p."""
        return float(np.interp(math.log(p), np.log(self.p), getattr(self, what)))

    def to_dataset(self):
        import xarray as xr

        ds = xr.Dataset(
            {
                "T": ("p", self.T, {"units": "K"}),
                "x_h2o": ("p", self.x, {"long_name": "water mole fraction"}),
                "q_h2o": ("p", self.q, {"long_name": "water mass fraction"}),
                "regime": ("p", self.regime, {"flag_meanings": " ".join(REGIME_NAMES.values())}),
                "z": ("p", self.z, {"units": "m"}),
                "N2": ("p", self.n2, {"units": "s-2"}),
            },
            coords={"p": ("p", self.p, {"units": "Pa"})},
        )
        ds.attrs.update({k: v for k, v in self.params.__dict__.items() if isinstance(v, (int, float, bool))})
        return ds


def grey_temperature(tau, params: ColumnParams):
    """Guillot (2010) eq. 29 radiative-equilibrium temperature at thermal optical depth tau."""
    g = params.gamma_v
    irr = params.flux_factor * params.teq**4
    t4 = 0.75 * params.tint**4 * (2.0 / 3.0 + tau) + 0.75 * irr * (
        2.0 / 3.0 + 1.0 / (g * _SQRT3) + (g / _SQRT3 - 1.0 / (g * _SQRT3)) * np.exp(-g * tau * _SQRT3)
    )
    return t4**0.25


def optical_depth(p, params: ColumnParams):
    a = params.kappa_pexp
    return params.kappa_th * params.p_ref / (params.gravity * (1.0 + a)) * (np.asarray(p) / params.p_ref) ** (1.0 + a)


def solve_column(params: ColumnParams | None = None, **overrides) -> Column:
    """Integrate T(p) from the top of the model down to ``p_bot``."""
    params = replace(params or ColumnParams(), **overrides)
    comp = params.composition
    sat = params.saturation
    x_deep = comp.x_h2o_deep
    varpi = comp.varpi
    mu_d = comp.mu_dry
    from .constants import MOLAR_MASS

    m_v = MOLAR_MASS["H2O"]
    a1 = 1.0 + params.kappa_pexp
    tau_ref = params.kappa_th * params.p_ref / (params.gravity * a1)
    tint4 = params.tint**4
    irr = params.flux_factor * params.teq**4
    gv = params.gamma_v
    inhibit = params.inhibition

    def gradient(lnp: float, lnT: float):
        p = math.exp(lnp)
        T = math.exp(lnT)
        tau = tau_ref * (p / params.p_ref) ** a1
        dt4 = 0.75 * (tint4 + irr * (1.0 - gv * gv) * math.exp(-gv * _SQRT3 * tau))
        g_rad = a1 * tau * dt4 / (4.0 * T**4)
        ps = sat.psat_s(T)
        if ps >= x_deep * p:  # below cloud base (or supercritical): unsaturated
            g_ad = dry_gradient(T, x_deep, comp)
            return (g_ad, DRY) if g_rad > g_ad else (g_rad, RADIATIVE)
        g_m = moist_gradient(T, p, comp, sat)
        if g_rad <= g_m:
            return g_rad, RADIATIVE
        if inhibit:
            xs = ps / p
            q_s = xs * m_v / (xs * m_v + (1.0 - xs) * mu_d)
            gam = sat.gamma_s(T)
            if gam > 0.0 and q_s * varpi * gam > 1.0:
                return g_rad, INHIBITED
        return g_m, MOIST

    lnp = np.linspace(math.log(params.p_top), math.log(params.p_bot), params.n_lev)
    T = np.empty(params.n_lev)
    regime = np.empty(params.n_lev, dtype=int)
    grad = np.empty(params.n_lev)
    y = math.log(float(grey_temperature(optical_depth(params.p_top, params), params)))
    T[0] = math.exp(y)
    grad[0], regime[0] = gradient(lnp[0], y)
    h = (lnp[1] - lnp[0]) / params.n_sub
    for i in range(1, params.n_lev):
        x0 = lnp[i - 1]
        for j in range(params.n_sub):
            xj = x0 + j * h
            k1 = gradient(xj, y)[0]
            k2 = gradient(xj + 0.5 * h, y + 0.5 * h * k1)[0]
            k3 = gradient(xj + 0.5 * h, y + 0.5 * h * k2)[0]
            k4 = gradient(xj + h, y + h * k3)[0]
            y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        T[i] = math.exp(y)
        grad[i], regime[i] = gradient(lnp[i], y)
    return Column(params, np.exp(lnp), T, regime, grad)


def intheat_flux(tint: float) -> float:
    """Internal heat flux (W m-2) for the PCM ``intheat`` key."""
    return SIGMA_SB * tint**4


def calibrate_grey(p_ref_profile, t_ref_profile, base: ColumnParams, p_range=(1.0e2, 1.0e6)):
    """Fit (kappa_th, gamma_v) so the emulator matches a reference T(p), e.g. from ``rcm1d``.

    Minimises the rms difference in ln T between ``p_range`` limits.  Returns the
    calibrated :class:`ColumnParams` and the rms residual in K.
    """
    from scipy.optimize import minimize

    p_ref_profile = np.asarray(p_ref_profile, float)
    t_ref_profile = np.asarray(t_ref_profile, float)
    sel = (p_ref_profile >= p_range[0]) & (p_ref_profile <= p_range[1])
    order = np.argsort(p_ref_profile[sel])
    lp, lt = np.log(p_ref_profile[sel][order]), np.log(t_ref_profile[sel][order])

    def cost(v):
        col = solve_column(base, kappa_th=10.0 ** v[0], gamma_v=10.0 ** v[1], n_lev=200, n_sub=2)
        return float(np.mean((np.interp(lp, np.log(col.p), np.log(col.T)) - lt) ** 2))

    res = minimize(cost, [math.log10(base.kappa_th), math.log10(base.gamma_v)], method="Nelder-Mead",
                   options={"xatol": 1e-2, "fatol": 1e-7, "maxiter": 200})
    best = replace(base, kappa_th=10.0 ** res.x[0], gamma_v=10.0 ** res.x[1])
    col = solve_column(best)
    rms = float(np.sqrt(np.mean((np.interp(lp, np.log(col.p), col.T) - np.exp(lt)) ** 2)))
    return best, rms
