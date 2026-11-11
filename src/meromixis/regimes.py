"""Regime maps: one grid, two boundaries.

Every point of the (forcing, enrichment) plane is classified by the same column:

* ``steam``            water never saturates; no cold trap, no lid;
* ``cold-trapped``     water rains out, but q_deep < q_cri: clouds over a convecting interior;
* ``lidded``           water rains out *and* the cloud base is stabilised (q_deep > q_cri).

Result I lives on the boundary between ``steam`` and the other two; Result II lives
inside ``lidded``.  The same functions are evaluated on GCM output by replacing
:func:`emulate_point` with the column statistics of a run (see ``diagnostics3d``).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .column import ColumnParams, solve_column

STEAM, COLD_TRAPPED, LIDDED = 0, 1, 2
REGIME_LABELS = {STEAM: "steam", COLD_TRAPPED: "cold-trapped", LIDDED: "lidded"}


def emulate_point(params: ColumnParams) -> dict:
    """Scalar summary of one emulator column."""
    col = solve_column(params)
    ct, layer = col.cold_trap, col.stable_layer
    regime = STEAM if not ct.exists else (LIDDED if layer.exists else COLD_TRAPPED)
    return {
        "x_deep": col.comp.x_h2o_deep,
        "x_above": ct.x,
        "transmission": ct.transmission,
        "p_coldtrap": ct.p,
        "t_coldtrap": ct.T,
        "layer": float(layer.exists),
        "p_layer_top": layer.p_top if layer.exists else np.nan,
        "p_layer_base": layer.p_base if layer.exists else np.nan,
        "layer_dt": layer.delta_t if layer.exists else np.nan,
        "layer_dz": layer.thickness if layer.exists else np.nan,
        "supercriticality": layer.supercriticality,
        "t_deep": float(col.T[-1]),
        "regime": regime,
    }


def regime_map(teqs, metallicities, base: ColumnParams | None = None):
    """Emulator summary on a (metallicity, teq) grid, as an ``xarray.Dataset``."""
    import xarray as xr

    base = base or ColumnParams()
    teqs = np.asarray(teqs, float)
    mets = np.asarray(metallicities, float)
    first = emulate_point(replace(base, teq=float(teqs[0]), metallicity=float(mets[0])))
    out = {k: np.full((mets.size, teqs.size), np.nan) for k in first}
    for i, z in enumerate(mets):
        for j, t in enumerate(teqs):
            for k, v in emulate_point(replace(base, teq=float(t), metallicity=float(z))).items():
                out[k][i, j] = v
    ds = xr.Dataset(
        {k: (("metallicity", "teq"), v) for k, v in out.items()},
        coords={"metallicity": mets, "teq": teqs},
    )
    ds["regime"].attrs["flag_meanings"] = " ".join(REGIME_LABELS.values())
    ds.attrs["model"] = "meromixis grey emulator (scoping only; not PCM output)"
    return ds


def transition_teq(metallicity: float, base: ColumnParams | None = None, level: float = 0.5,
                   bracket=(150.0, 800.0), tol: float = 0.25) -> float:
    """Effective Teq at which the cold-trap transmission crosses ``level`` (bisection)."""
    base = base or ColumnParams()

    def tr(t):
        return solve_column(replace(base, teq=float(t), metallicity=metallicity)).cold_trap.transmission

    lo, hi = bracket
    if tr(lo) >= level:
        return lo
    if tr(hi) < level:
        return np.nan
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if tr(mid) >= level:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def transition_width(metallicity: float, base: ColumnParams | None = None) -> tuple[float, float, float]:
    """(T at transmission 0.1, 0.5, 0.9): the location and width of the transition."""
    return tuple(transition_teq(metallicity, base, lv) for lv in (0.1, 0.5, 0.9))


def transition_from_samples(teq, transmission, levels=(0.1, 0.5, 0.9)) -> tuple[float, ...]:
    """Forcing at which a sampled transmission curve crosses each level.

    Works on any set of runs: emulator columns, ``rcm1d`` twins, or the terminator
    ``transmission_actual`` of a set of GCM runs against their *realised* equilibrium
    temperature.  Interpolates log10(transmission) linearly in teq; returns NaN for a
    level the samples do not bracket.  The width of the transition is ``out[-1] - out[0]``.
    """
    teq = np.asarray(teq, float)
    order = np.argsort(teq)
    t, y = teq[order], np.log10(np.clip(np.asarray(transmission, float)[order], 1e-300, None))
    out = []
    for lv in levels:
        d = y - np.log10(lv)
        idx = np.nonzero((d[:-1] < 0) & (d[1:] >= 0))[0]
        if idx.size == 0:
            out.append(np.nan)
            continue
        i = idx[0]
        out.append(float(t[i] + (t[i + 1] - t[i]) * (-d[i]) / (d[i + 1] - d[i])))
    return tuple(out)
