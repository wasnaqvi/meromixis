"""Condensation-inhibited convection (Guillot 1995; Leconte et al. 2017).

In a hydrogen-rich gas the condensable is heavier than the background, so a warm
saturated parcel carries more heavy vapour.  Moist convection shuts off where

    q_v  >  q_cri = 1 / (varpi * gamma),   varpi = 1 - mu_d/M_v,  gamma = dln psat/dln T,

with q_v the vapour mass fraction of the gas.  For constant latent heat this is the
familiar ``q_cri = R T / ((M_v - mu_d) L)``, the same expression the Generic PCM
evaluates in ``moistadj_generic`` when ``moist_convection_inhibition = .true.``.

The functions here are profile-agnostic: they take arrays of p, T and water mole
fraction, so the same code runs on the 1-D emulator and on GCM columns.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .composition import Composition, q_to_x, x_to_q
from .constants import K_B, MOLAR_MASS, R_GAS
from .thermo import Saturation

M_V = MOLAR_MASS["H2O"]


def q_crit(T, comp: Composition, sat: Saturation | None = None):
    """Critical vapour mass fraction above which moist convection is inhibited."""
    sat = sat or Saturation()
    gam = sat.gamma(T)
    safe = np.where(gam > 0, gam, 1.0)
    return np.where(gam > 0, 1.0 / (comp.varpi * safe), np.inf)


def x_crit(T, comp: Composition, sat: Saturation | None = None):
    """The same threshold expressed as a water mole fraction (clipped to 1)."""
    qc = np.minimum(q_crit(T, comp, sat), 1.0)
    return q_to_x(qc, comp.mu_dry)


def metallicity_crit(T, sat: Saturation | None = None, lo=0.1, hi=400.0) -> float:
    """Enrichment (x solar) at which the deep water abundance equals q_cri at temperature T."""
    from scipy.optimize import brentq

    def f(logz):
        c = Composition(10.0**logz)
        return np.log(c.q_h2o_deep) - np.log(float(q_crit(T, c, sat)))

    a, b = np.log10(lo), np.log10(hi)
    if f(a) > 0:
        return lo
    if f(b) < 0:
        return np.inf
    return 10.0 ** brentq(f, a, b)


def scale_height(T, x_h2o, comp: Composition, gravity: float):
    """Pressure scale height R T / (mu g), m."""
    return R_GAS * np.asarray(T) / (comp.mu(np.asarray(x_h2o)) * gravity)


def altitude(p, T, x_h2o, comp: Composition, gravity: float):
    """Hydrostatic altitude (m) above the highest-pressure level; p in any order."""
    p = np.asarray(p, float)
    order = np.argsort(-p)
    h = scale_height(np.asarray(T)[order], np.asarray(x_h2o)[order], comp, gravity)
    lnp = np.log(p[order])
    dz = -0.5 * (h[1:] + h[:-1]) * np.diff(lnp)
    # A non-finite layer (a column whose base is supercritical for water) must not
    # poison the cumulative sum for every level above it, the way a plain cumsum
    # would.  Skip its thickness and mark only that level unusable; every quantity
    # built from z uses differences over a selected range, and those ranges exclude
    # the offending base levels.
    finite = np.isfinite(dz)
    z_sorted = np.concatenate([[0.0], np.cumsum(np.where(finite, dz, 0.0))])
    z_sorted[1:][~finite] = np.nan
    z = np.empty_like(z_sorted)
    z[order] = z_sorted
    return z


def number_density(p, T):
    """Total number density, m-3."""
    return np.asarray(p, float) / (K_B * np.asarray(T, float))


def brunt_vaisala_sq(p, T, x_h2o, comp: Composition, gravity: float):
    """Ledoux N^2 (s-2) for unsaturated displacements, including the mu gradient.

    N^2 = (g / H) * (nabla_ad - nabla_T + nabla_mu), with nabla = d ln / d ln p.
    Positive means stable.  In a saturated layer the relevant criterion is the moist
    one (``q > q_cri``); this dry measure is what gravity waves and the large-scale
    overturning feel.
    """
    p = np.asarray(p, float)
    T = np.asarray(T, float)
    x = np.asarray(x_h2o, float)
    lnp = np.log(p)
    mu = comp.mu(x)
    nabla_t = np.gradient(np.log(T), lnp)
    nabla_mu = np.gradient(np.log(mu), lnp)
    nabla_ad = R_GAS / comp.cp_molar_mix(T, x)
    h = R_GAS * T / (mu * gravity)
    return (gravity / h) * (nabla_ad - nabla_t + nabla_mu)


@dataclass(frozen=True)
class StableLayer:
    """The compositionally stabilised layer under the water cloud deck."""

    exists: bool
    p_top: float = np.nan  # Pa, where q falls to q_cri
    p_base: float = np.nan  # Pa, cloud base (q reaches its deep value)
    t_top: float = np.nan
    t_base: float = np.nan
    thickness: float = np.nan  # m
    n_mid: float = np.nan  # number density at the log-pressure midpoint, m-3
    q_deep: float = np.nan
    q_crit_base: float = np.nan

    @property
    def delta_t(self) -> float:
        return self.t_base - self.t_top

    @property
    def supercriticality(self) -> float:
        """q_deep / q_cri at cloud base; > 1 means a layer can exist."""
        return self.q_deep / self.q_crit_base


def find_stable_layer(
    p,
    T,
    x_h2o,
    comp: Composition,
    gravity: float,
    sat: Saturation | None = None,
    x_deep: float | None = None,
    tol: float = 0.02,
) -> StableLayer:
    """Locate the inhibited layer in a column.

    Following Leconte+2017 the layer spans from the cloud base (water reaches its
    deep value) up to the level where the vapour mass fraction drops below q_cri(T).
    ``x_deep`` defaults to the value at the highest pressure.  ``tol`` is the
    fractional departure from ``x_deep`` that counts as "above cloud base".
    """
    sat = sat or Saturation()
    p = np.asarray(p, float)
    order = np.argsort(-p)  # bottom first
    pb, tb, xb = p[order], np.asarray(T, float)[order], np.asarray(x_h2o, float)[order]
    if x_deep is None:
        x_deep = float(xb[0])
    q = x_to_q(xb, comp.mu_dry)
    qc = q_crit(tb, comp, sat)
    q_deep = float(x_to_q(x_deep, comp.mu_dry))

    above = np.nonzero(xb < x_deep * (1.0 - tol))[0]
    if above.size == 0:  # water never condenses in the domain
        return StableLayer(False, q_deep=q_deep)
    i_base = max(above[0] - 1, 0)
    qc_base = float(qc[i_base])
    if q_deep <= qc_base:
        return StableLayer(
            False, q_deep=q_deep, q_crit_base=qc_base, p_base=float(pb[i_base]), t_base=float(tb[i_base])
        )
    sub = np.nonzero(q[i_base:] < qc[i_base:])[0]
    i_top = i_base + (sub[0] if sub.size else len(pb) - 1 - i_base)
    i_top = max(i_top, i_base + 1)
    seg = slice(i_base, i_top + 1)
    z = altitude(pb[seg], tb[seg], xb[seg], comp, gravity)
    p_mid = np.sqrt(pb[i_base] * pb[i_top])
    t_mid = np.interp(np.log(p_mid), np.log(pb[seg])[::-1], tb[seg][::-1])
    return StableLayer(
        True,
        p_top=float(pb[i_top]),
        p_base=float(pb[i_base]),
        t_top=float(tb[i_top]),
        t_base=float(tb[i_base]),
        thickness=float(z.max() - z.min()),
        n_mid=float(number_density(p_mid, t_mid)),
        q_deep=q_deep,
        q_crit_base=qc_base,
    )
