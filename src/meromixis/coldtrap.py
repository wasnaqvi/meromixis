"""Cold-trap diagnostics.

The water mole fraction that survives to the upper atmosphere is the smallest
saturation mole fraction met on the way up (rain-out, no supersaturation, no
re-evaporation aloft):

    x(p) = min( x_deep ,  min_{p' >= p} psat(T(p')) / p' ).

The *cold-trap transmission*  T_ct = x_above / x_deep  runs from ~0 (water hidden) to
1 (steam atmosphere).  "The cold-trap transition" is where T_ct crosses 1/2 as the
forcing increases; its width is the forcing interval between T_ct = 0.1 and 0.9.

All functions accept profiles in either vertical order and are reused unchanged on
GCM columns (see ``diagnostics3d``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .thermo import Saturation


@dataclass(frozen=True)
class ColdTrap:
    exists: bool  # False if water never saturates below p_min (no trap)
    p: float  # Pa, cold-trap level (hygropause)
    T: float  # K
    x: float  # water mole fraction above the trap
    x_deep: float

    @property
    def transmission(self) -> float:
        return self.x / self.x_deep


def water_profile(p, T, x_deep: float, sat: Saturation | None = None):
    """Rain-out water mole-fraction profile, returned in the input order."""
    sat = sat or Saturation()
    p = np.asarray(p, float)
    order = np.argsort(-p)  # bottom first
    xs = sat.x_sat(np.asarray(T, float)[order], p[order])
    x_sorted = np.minimum.accumulate(np.minimum(xs, x_deep))
    x = np.empty_like(x_sorted)
    x[order] = x_sorted
    return x


def cold_trap(p, T, x_deep: float, sat: Saturation | None = None, p_min: float = 0.1) -> ColdTrap:
    """Cold trap of a single column.  Levels with p < p_min (Pa) are ignored."""
    out = cold_trap_nd(np.asarray(p, float), np.asarray(T, float), x_deep, axis=0, sat=sat, p_min=p_min)
    return ColdTrap(bool(out["exists"]), float(out["p"]), float(out["T"]), float(out["x"]), float(x_deep))


def cold_trap_nd(p, T, x_deep, axis: int, sat: Saturation | None = None, p_min: float = 0.1) -> dict:
    """Vectorised cold trap along ``axis`` of arrays p, T (same shape).

    ``x_deep`` may be a scalar or an array broadcastable to the shape with ``axis``
    removed.  Returns a dict of arrays: exists, p, T, x, transmission.
    """
    sat = sat or Saturation()
    p = np.asarray(p, float)
    T = np.asarray(T, float)
    xs = np.where(p >= p_min, sat.x_sat(T, p), np.inf)
    idx = np.expand_dims(np.argmin(xs, axis=axis), axis)
    xs_min = np.take_along_axis(xs, idx, axis).squeeze(axis)
    x_deep = np.broadcast_to(np.asarray(x_deep, float), xs_min.shape)
    exists = xs_min < x_deep
    x = np.minimum(xs_min, x_deep)
    return {
        "exists": exists,
        "p": np.where(exists, np.take_along_axis(p, idx, axis).squeeze(axis), np.nan),
        "T": np.where(exists, np.take_along_axis(T, idx, axis).squeeze(axis), np.nan),
        "x": x,
        "transmission": x / x_deep,
    }


def skin_temperature(teq: float, gamma_v: float = 0.0, tint: float = 0.0) -> float:
    """Top-of-atmosphere temperature of the grey model (reduces to Teq / 2^(1/4))."""
    return (0.5 * tint**4 + teq**4 * (0.5 + 0.25 * np.sqrt(3.0) * gamma_v)) ** 0.25
