"""The ammonia barrier: is supply through the stable layer slower than photochemical loss?

A photochemically fragile gas (NH3) is made deep and destroyed aloft.  In steady state
the upward flux through every level equals the column loss above it.  The largest flux
diffusion can deliver from a reservoir of mixing ratio f_deep is

    Phi_max = f_deep / R,      R = integral dz / (K n)      (resistances in series),

so a thin layer of small K under the clouds can throttle the whole column.  If the
net (irreversible) photolysis loss  L = eta * Phi_UV  exceeds Phi_max, the upper
atmosphere is starved and the gas collapses until it becomes optically thin in the UV.

This module provides
* the K(p) profile implied by a column's regimes (:func:`kzz_profile`);
* the series resistance, maximum supply and the critical layer diffusivity
  (:func:`resistance`, :func:`max_supply`, :func:`k_crit`): the "kill test";
* a sharp-front steady-state profile for a self-shielded, photolysed gas
  (:func:`steady_state`) and a linear diffusion-loss solver for prescribed sinks
  (:func:`solve_diffusion_loss`, the 1-D analogue of a relaxation tracer in the GCM);
* a budget that separates the *thermal* effect of inhibition (hotter interior, lower
  quenched NH3; Leconte+2024) from the *transport* effect (:func:`ammonia_budget`).

Units are SI: K in m2 s-1 (1 m2 s-1 = 1e4 cm2 s-1), fluxes in m-2 s-1.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded

from .chemistry import Quench, kzz_mixing_length, quench
from .column import DRY, INHIBITED, MOIST, RADIATIVE, Column, solve_column


@dataclass(frozen=True)
class Mixing:
    """Eddy diffusivities by dynamical regime, m2 s-1.

    k_layer : inside the inhibited layer.  Leconte+2024's convection-resolving runs
        give ~0.1-10 m2 s-1 there (as quoted in the project primer; VERIFY L1).
    k_rad   : in radiative zones.
    k_conv  : in convective zones; ``None`` uses the mixing-length estimate.
    """

    k_layer: float = 1.0
    k_rad: float = 100.0
    k_conv: float | None = None


@dataclass(frozen=True)
class Photolysis:
    """Self-shielded photolysis sink.

    photon_flux : global-mean flux of photons the species can absorb, m-2 s-1.
        Default 1e17 (= 1e13 cm-2 s-1) is the hand estimate for NH3 (< 230 nm) from the
        explainer, *not* a measured spectrum; integrate a real one with
        :func:`photon_flux_from_spectrum`.
    eta : fraction of photolysis events that are irreversible (NH2 -> N2H4 -> N2
        rather than NH2 + H2 -> NH3).  The central unknown of the kill test.
    sigma : effective absorption cross-section, m2.
    p_floor : deepest pressure UV can reach whatever the species does (Rayleigh
        scattering / haze), Pa.
    """

    photon_flux: float = 1.0e17
    eta: float = 0.3
    sigma: float = 3.0e-22
    p_floor: float = 3.0e4

    @property
    def net_loss_max(self) -> float:
        return self.eta * self.photon_flux


NH3_PHOTOLYSIS = Photolysis()
CH4_PHOTOLYSIS = Photolysis(photon_flux=1.0e16, eta=0.5, sigma=1.8e-21, p_floor=1.0e2)


def photon_flux_from_spectrum(wavelength_nm, flux_w_m2_nm, band=(160.0, 230.0), geometry: float = 0.25) -> float:
    """Photon flux (m-2 s-1) in ``band`` from a spectrum at the planet's orbit.

    ``geometry`` = 1/4 converts the substellar flux to a global mean.
    """
    from .constants import N_A  # noqa: F401  (kept for symmetry with other helpers)

    h, c = 6.62607015e-34, 2.99792458e8
    wl = np.asarray(wavelength_nm, float)
    fl = np.asarray(flux_w_m2_nm, float)
    sel = (wl >= band[0]) & (wl <= band[1])
    photons = fl[sel] * (wl[sel] * 1e-9) / (h * c)
    return float(geometry * np.trapezoid(photons, wl[sel]))


# ---- mixing profile and resistances -------------------------------------------------

def layer_mask(col: Column) -> np.ndarray:
    """Levels where compositional stratification suppresses mixing.

    Any level flagged INHIBITED, plus the whole span from cloud base to the q = q_cri
    level when a stable layer exists (that span stays Ledoux-stable where it is radiative).
    """
    mask = col.regime == INHIBITED
    layer = col.stable_layer
    if layer.exists:
        mask = mask | ((col.p >= layer.p_top) & (col.p <= layer.p_base))
    return mask


def kzz_profile(col: Column, mix: Mixing = Mixing()) -> np.ndarray:
    """K(p) in m2 s-1, same ordering as ``col.p``."""
    k = np.empty_like(col.T)
    k_conv = kzz_mixing_length(col) if mix.k_conv is None else np.full_like(col.T, mix.k_conv)
    conv = (col.regime == DRY) | (col.regime == MOIST)
    k[conv] = k_conv[conv]
    k[col.regime == RADIATIVE] = mix.k_rad
    k[layer_mask(col)] = mix.k_layer
    return k


def resistance(col: Column, k: np.ndarray, p_from: float, p_to: float) -> float:
    """Series resistance integral dz / (K n) between two pressures, m2 s."""
    lo, hi = sorted((p_from, p_to))
    sel = (col.p >= lo) & (col.p <= hi)
    z = col.z[sel]
    order = np.argsort(z)
    return float(np.trapezoid(1.0 / (k[sel][order] * col.n[sel][order]), z[order]))


def max_supply(col: Column, k: np.ndarray, f_deep: float, p_from: float, p_to: float) -> float:
    """Largest diffusive flux (m-2 s-1) a reservoir of mixing ratio f_deep at p_from can deliver to p_to."""
    return f_deep / resistance(col, k, p_from, p_to)


def k_crit(net_loss: float, thickness: float, n_layer: float, f_deep: float) -> float:
    """Layer diffusivity below which supply cannot balance ``net_loss``: L dz / (n f)."""
    return net_loss * thickness / (n_layer * f_deep)


# ---- steady-state solver -------------------------------------------------------------

def solve_diffusion_loss(z, kn, n, k_loss, f_bottom: float) -> np.ndarray:
    """Solve d/dz( K n df/dz ) = k_loss n f with f(z0) = f_bottom and zero flux at the top.

    z ascending (m); kn = K*n and n at the nodes; k_loss in s-1.  Finite volumes with
    harmonic-mean interface conductances, so sharp K contrasts are handled exactly.
    """
    z = np.asarray(z, float)
    m = len(z)
    dzc = np.diff(z)
    cond = 2.0 * kn[:-1] * kn[1:] / (kn[:-1] + kn[1:]) / dzc  # interface conductance
    vol = np.empty(m)
    vol[1:-1] = 0.5 * (dzc[:-1] + dzc[1:])
    vol[0], vol[-1] = 0.5 * dzc[0], 0.5 * dzc[-1]
    ab = np.zeros((3, m))
    rhs = np.zeros(m)
    ab[1, 0], rhs[0] = 1.0, f_bottom
    for i in range(1, m):
        up = cond[i] if i < m - 1 else 0.0
        ab[1, i] = -(cond[i - 1] + up) - k_loss[i] * n[i] * vol[i]
        ab[2, i - 1] = cond[i - 1]  # sub-diagonal: coefficient of f[i-1] in row i
        if i < m - 1:
            ab[0, i + 1] = up  # super-diagonal: coefficient of f[i+1] in row i
    return solve_banded((1, 1), ab, rhs)


@dataclass
class TracerSolution:
    p: np.ndarray  # Pa, top -> bottom (the column from the top down to p_source)
    f: np.ndarray  # mixing ratio
    f_deep: float
    column_loss: float  # m-2 s-1
    supply_max: float  # m-2 s-1, largest flux deliverable to the UV floor
    p_front: float  # Pa, level above which the gas is photolysed away
    supply_limited: bool  # True: starved (the barrier wins); False: photon-limited

    def mean(self, p_lo: float = 1.0e2, p_hi: float = 1.0e4, floor: float = 1.0e-12) -> float:
        """Log-mean mixing ratio between p_lo and p_hi (default 1-100 mbar).

        Values are floored at ``floor``: above the front the sharp-front tail is not a
        prediction (recombination chemistry sets the real floor), only "gone".
        """
        sel = (self.p >= p_lo) & (self.p <= p_hi)
        return float(np.exp(np.mean(np.log(np.maximum(self.f[sel], floor)))))

    @property
    def depletion(self) -> float:
        """Band-mean mixing ratio relative to the deep value (floored; see :meth:`mean`)."""
        return self.mean() / self.f_deep

    def hidden(self, p_obs: float = 1.0e3) -> bool:
        """True if the photolysis front lies deeper than ``p_obs`` (gas absent where spectra probe)."""
        return self.p_front > p_obs


def steady_state(
    col: Column,
    f_deep: float,
    p_source: float,
    mix: Mixing = Mixing(),
    phot: Photolysis = NH3_PHOTOLYSIS,
) -> TracerSolution:
    """Steady profile of a gas held at ``f_deep`` at ``p_source`` and photolysed aloft.

    Sharp-front solution.  The gas is opaque in its own UV band, so photons are
    absorbed in a thin front and there is no loss below it; the flux is then uniform
    below the front and f falls linearly in the resistance coordinate
    rho(p) = integral dz/(K n) measured upward from the source:

    * photon-limited: the front sits where rho = f_deep / (eta Phi_UV);
    * starved: if that level lies deeper than ``p_floor`` (where UV cannot reach), the
      gas collapses at the floor to the optically thin value that balances supply.

    Above the front the gas decays over the reaction-diffusion length sqrt(K / J).
    """
    k_all = kzz_profile(col, mix)
    sel = col.p <= p_source
    order = np.argsort(col.z[sel])  # bottom first
    z, p, n, k = col.z[sel][order], col.p[sel][order], col.n[sel][order], k_all[sel][order]
    rho = np.concatenate([[0.0], np.cumsum(0.5 * (1.0 / (k[1:] * n[1:]) + 1.0 / (k[:-1] * n[:-1])) * np.diff(z))])
    loss_max = phot.net_loss_max
    j0 = phot.eta * phot.sigma * phot.photon_flux  # optically thin loss frequency, s-1
    ell = np.sqrt(k / j0)  # reaction-diffusion length, m
    i_floor = int(np.argmax(p <= phot.p_floor)) if (p <= phot.p_floor).any() else len(p) - 1
    rho_floor = max(rho[i_floor], 1e-300)
    supply_max = f_deep / rho_floor
    rho_front = f_deep / loss_max

    if rho_front > rho_floor:  # photon-limited
        starved = False
        if rho_front >= rho[-1]:  # front above the model top
            i_f, f = len(p) - 1, f_deep * (1.0 - rho / rho_front)
            loss = f_deep / max(rho[-1], 1e-300) if rho[-1] > 0 else loss_max
            loss = min(loss, loss_max)
        else:
            i_f = int(np.searchsorted(rho, rho_front))
            f = np.maximum(f_deep * (1.0 - rho / rho_front), 0.0)
            loss = loss_max
        f_edge = min(1.0 / (phot.sigma * n[i_f] * ell[i_f]), f_deep)
    else:  # starved: collapse at the UV floor
        from scipy.optimize import brentq

        starved, i_f = True, i_floor
        thin = phot.sigma * n[i_f] * min(ell[i_f], 3.0e5)

        def balance(fx):
            return (f_deep - fx) / rho_floor - loss_max * (-np.expm1(-thin * fx))

        f_edge = brentq(balance, 0.0, f_deep, xtol=1e-30, rtol=1e-10)
        f = f_deep - (f_deep - f_edge) * np.minimum(rho / rho_floor, 1.0)
        loss = (f_deep - f_edge) / rho_floor
    tail = f_edge * np.exp(-np.clip((z - z[i_f]) / ell[i_f], 0.0, 700.0))
    f = np.where(np.arange(len(p)) >= i_f, np.maximum(np.minimum(f, f_edge), tail) if starved else np.maximum(f, tail), f)
    return TracerSolution(p[::-1], f[::-1], f_deep, float(loss), float(supply_max), float(p[i_f]), starved)


def solve_selfshielded(
    col: Column,
    f_deep: float,
    p_source: float,
    mix: Mixing = Mixing(),
    phot: Photolysis = NH3_PHOTOLYSIS,
    n_iter: int = 4000,
    relax: float = 0.15,
    tol: float = 1.0e-7,
):
    """Numerical steady state with a *self-shielded* photolysis sink.

    Solves the equation :func:`steady_state` approximates analytically:

        d/dz ( K n df/dz ) = J(z) n f,
        J(z) = eta sigma Phi_UV exp(-tau(z)),
        tau(z) = sigma * integral_z^top n f dz',

    by fixed-point iteration on J with under-relaxation.  The sink is nonlinear
    because the gas shields itself, which is exactly the feature the sharp-front
    solution encodes; solving it directly is therefore a like-for-like check on
    ``steady_state``, where a sink fixed at the optically thin rate is not.

    Returns ``(p, f, info)`` with p top->bottom to match
    :class:`TracerSolution`, and ``info`` carrying the iteration diagnostics.
    """
    k_all = kzz_profile(col, mix)
    sel = col.p <= p_source
    order = np.argsort(col.z[sel])  # ascending z, i.e. bottom first
    z, p = col.z[sel][order], col.p[sel][order]
    n, k = col.n[sel][order], k_all[sel][order]
    good = np.isfinite(z) & np.isfinite(n) & np.isfinite(k) & (n > 0) & (k > 0)
    z, p, n, k = z[good], p[good], n[good], k[good]

    j0 = phot.eta * phot.sigma * phot.photon_flux  # optically thin frequency, s-1
    reachable = p <= phot.p_floor  # UV cannot get below this whatever the gas does
    f = np.full_like(z, f_deep)
    f[0] = f_deep
    last = None
    for it in range(n_iter):
        # column of absorber above each level (trapezoid, from the top down)
        nf = n * f
        seg = 0.5 * (nf[1:] + nf[:-1]) * np.diff(z)
        col_above = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
        tau = phot.sigma * col_above
        k_loss = np.where(reachable, j0 * np.exp(-np.clip(tau, 0.0, 700.0)), 0.0)
        f_new = solve_diffusion_loss(z, k * n, n, k_loss, f_deep)
        f_new = np.clip(f_new, 0.0, f_deep)
        f = (1.0 - relax) * f + relax * f_new
        if last is not None:
            change = float(np.max(np.abs(f - last)) / max(f_deep, 1e-300))
            if change < tol:
                break
        last = f.copy()
    nf = n * f
    seg = 0.5 * (nf[1:] + nf[:-1]) * np.diff(z)
    col_above = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    tau = phot.sigma * col_above
    k_loss = np.where(reachable, j0 * np.exp(-np.clip(tau, 0.0, 700.0)), 0.0)
    loss = float(np.trapezoid(k_loss * n * f, z))
    # front = deepest level at which the gas has lost half its deep value
    below_half = np.nonzero(f <= 0.5 * f_deep)[0]
    p_front = float(p[below_half[0]]) if below_half.size else float(p[-1])
    info = {
        "iterations": it + 1,
        "converged": bool(last is not None and change < tol),
        "max_change": float(change) if last is not None else np.nan,
        "column_loss": loss,
        "p_front": p_front,
        "tau_at_front": float(tau[below_half[0]]) if below_half.size else float(tau[-1]),
    }
    return p[::-1], f[::-1], info


def band_mean(p, f, p_lo: float = 1.0e2, p_hi: float = 1.0e4, floor: float = 1.0e-12) -> float:
    """Log-mean mixing ratio between two pressures, floored (cf. TracerSolution.mean)."""
    sel = (p >= p_lo) & (p <= p_hi)
    if not sel.any():
        return float("nan")
    return float(np.exp(np.mean(np.log(np.maximum(f[sel], floor)))))


# ---- the kill test -------------------------------------------------------------------

def k_layer_crit(
    col: Column,
    f_deep: float,
    p_source: float,
    mix: Mixing = Mixing(),
    phot: Photolysis = NH3_PHOTOLYSIS,
    p_obs: float = 1.0e3,
) -> float:
    """Layer diffusivity (m2 s-1) below which the photolysis front sinks beneath ``p_obs``.

    Series-resistance form of the kill test: the gas is hidden at ``p_obs`` when
    R_layer(k) + R_rest >= f_deep / (eta Phi_UV).  Returns ``inf`` if the rest of the
    column already hides it (the barrier is not needed) and ``nan`` if there is no layer.
    """
    layer = col.stable_layer
    if not layer.exists:
        return np.nan
    k_all = kzz_profile(col, mix)
    sel = (col.p <= p_source) & (col.p >= p_obs)
    order = np.argsort(col.z[sel])
    z, n, k = col.z[sel][order], col.n[sel][order], k_all[sel][order]
    inside = layer_mask(col)[sel][order]
    w = np.gradient(z)  # node weights
    geom_layer = float(np.sum((w / n)[inside]))  # integral dz / n over the layer
    r_rest = float(np.sum((w / (k * n))[~inside]))
    budget = f_deep / phot.net_loss_max - r_rest
    return np.inf if budget <= 0 else geom_layer / budget


def depletion_time(col: Column, f: float, p_level: float, net_loss: float) -> float:
    """Time (s) for ``net_loss`` to remove the gas column above ``p_level``."""
    above = col.p < p_level
    column = float(np.trapezoid(col.n[above][::-1], col.z[above][::-1]))
    return column * f / net_loss


@dataclass
class AmmoniaBudget:
    """Everything needed to judge the barrier for one planet."""

    quench_ref: Quench  # NH3, no inhibition
    quench_inh: Quench  # NH3, with inhibition
    sol_ref: TracerSolution
    sol_inh: TracerSolution
    k_crit_simple: float  # L dz / (n f): the back-of-envelope form, m2 s-1
    k_crit: float  # series-resistance form, m2 s-1
    k_layer: float  # the value assumed in ``mix``
    co_ref: Quench
    co_inh: Quench
    co_fill_time: float  # s, time for quenched CO to fill the column above the layer
    nh3_depletion_time: float  # s, time for photolysis to empty NH3 above the layer

    @property
    def thermal_factor(self) -> float:
        """Change in quenched NH3 caused by the hotter interior alone (Leconte+2024 effect)."""
        return self.quench_inh.x / self.quench_ref.x

    @property
    def co_enhancement(self) -> float:
        return self.co_inh.x / self.co_ref.x

    @property
    def hidden_by_barrier(self) -> bool:
        """NH3 hidden with the layer present but visible in the freely mixing reference."""
        return self.sol_inh.hidden() and not self.sol_ref.hidden()

    @property
    def margin(self) -> float:
        """k_crit / k_layer.  > 1: the barrier wins; < 1: the layer leaks too fast."""
        return self.k_crit / self.k_layer

    def verdict(self) -> str:
        if not np.isfinite(self.k_crit_simple):
            return "no stable layer: nothing to test"
        if self.sol_ref.hidden():
            return "NH3 is photolysed below the observable level even without a barrier"
        if self.sol_inh.hidden():
            return f"barrier hides NH3 (k_crit / k_layer = {self.margin:.2g})"
        return f"layer leaks too fast: NH3 stays visible (k_crit / k_layer = {self.margin:.2g})"


def ammonia_budget(
    params,
    mix: Mixing = Mixing(),
    phot: Photolysis = NH3_PHOTOLYSIS,
    col_inh: Column | None = None,
    col_ref: Column | None = None,
    p_obs: float = 1.0e3,
) -> AmmoniaBudget:
    """Run the inhibited and reference columns and separate thermal from transport effects."""
    col_inh = col_inh or solve_column(params, inhibition=True)
    col_ref = col_ref or solve_column(params, inhibition=False)
    q_ref, q_inh = quench(col_ref, "NH3", mix.k_conv), quench(col_inh, "NH3", mix.k_conv)
    sol_ref = steady_state(col_ref, q_ref.x, q_ref.p, mix, phot)
    sol_inh = steady_state(col_inh, q_inh.x, q_inh.p, mix, phot)
    layer = col_inh.stable_layer
    co_ref, co_inh = quench(col_ref, "CO", mix.k_conv), quench(col_inh, "CO", mix.k_conv)
    if layer.exists:
        kc_simple = k_crit(phot.net_loss_max, layer.thickness, layer.n_mid, q_inh.x)
        kc = k_layer_crit(col_inh, q_inh.x, q_inh.p, mix, phot, p_obs)
        k_all = kzz_profile(col_inh, mix)
        above = col_inh.p < layer.p_top
        n_col_above = float(np.trapezoid(col_inh.n[above][::-1], col_inh.z[above][::-1]))
        fill = n_col_above * resistance(col_inh, k_all, min(co_inh.p, col_inh.p[-1]), layer.p_top)
        t_dep = depletion_time(col_inh, q_inh.x, layer.p_top, phot.net_loss_max)
    else:
        kc_simple = kc = np.nan
        fill = t_dep = np.nan
    return AmmoniaBudget(q_ref, q_inh, sol_ref, sol_inh, kc_simple, kc, mix.k_layer, co_ref, co_inh, fill, t_dep)
