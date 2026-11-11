"""Synthetic PCM-like output, so the whole analysis chain can run before any GCM does.

:func:`make_synthetic_run` writes fields with the *raw* Generic PCM names and layout
(``Time, altitude, latitude, longitude``; ``temp, ps, aps, bps, h2o_vap, w ...``) built from
emulator columns with a day-night forcing contrast.  It is a plumbing test, not a model:
winds are zero and there is no horizontal transport.

:func:`make_mixing_test` builds a standardised dataset whose resolved eddy flux
corresponds to a known diffusivity, to validate :func:`diagnostics3d.kzz_flux_gradient`.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import xarray as xr

from .column import ColumnParams, solve_column
from .composition import x_to_q
from .constants import R_GAS, SIGMA_SB


def _grid(nlon, nlat):
    lon = np.linspace(-180.0, 180.0, nlon, endpoint=False)
    lat = np.linspace(90.0, -90.0, nlat + 2)[1:-1]  # LMDZ orders north -> south
    return lon, lat


def make_synthetic_run(params: ColumnParams | None = None, nlon: int = 32, nlat: int = 24, nlev: int = 40,
                       ntime: int = 3, contrast: float = 0.3, p_bot: float = 8.0e6, p_top: float = 20.0,
                       n_columns: int = 7, seed: int = 0) -> xr.Dataset:
    """Raw-format synthetic ``diagfi`` built from emulator columns.

    ``contrast`` in [0, 1] scales the day-night forcing difference: 0 is a horizontally
    uniform (perfect weak-temperature-gradient) planet, 1 is no redistribution at all.
    """
    params = params or ColumnParams()
    rng = np.random.default_rng(seed)
    lon, lat = _grid(nlon, nlat)
    mu = np.cos(np.deg2rad(lat))[:, None] * np.cos(np.deg2rad(lon))[None, :]
    ff = 1.0 + contrast * (4.0 * np.maximum(mu, 0.0) - 1.0)  # local / global-mean forcing; area mean = 1
    p_lev = np.geomspace(p_bot, p_top, nlev)  # bottom first
    ffs = np.linspace(ff.min(), ff.max(), n_columns) if contrast > 0 else np.array([1.0, 1.0 + 1e-9])
    t_tab, x_tab = [], []
    for f in ffs:
        col = solve_column(replace(params, flux_factor=float(f)))
        t_tab.append(np.interp(np.log(p_lev), np.log(col.p), col.T))
        x_tab.append(np.exp(np.interp(np.log(p_lev), np.log(col.p), np.log(col.x))))
    t_tab, x_tab = np.array(t_tab), np.array(x_tab)
    idx = np.clip(np.searchsorted(ffs, ff) - 1, 0, len(ffs) - 2)
    w = (ff - ffs[idx]) / (ffs[idx + 1] - ffs[idx])
    temp = (t_tab[idx] * (1 - w)[..., None] + t_tab[idx + 1] * w[..., None]).transpose(2, 0, 1)
    lnx = (np.log(x_tab[idx]) * (1 - w)[..., None] + np.log(x_tab[idx + 1]) * w[..., None]).transpose(2, 0, 1)
    q = x_to_q(np.exp(lnx), params.composition.mu_dry)
    shape = (ntime, nlev, nlat, nlon)
    temp4 = np.broadcast_to(temp, shape) + rng.normal(0.0, 0.05, shape)
    q4 = np.broadcast_to(q, shape).copy()
    ps = np.full((ntime, nlat, nlon), p_bot * 1.02)
    olr = SIGMA_SB * (params.teq**4 + params.tint**4) * np.ones((ntime, nlat, nlon))
    asr = 4.0 * SIGMA_SB * params.teq**4 * np.broadcast_to(np.maximum(mu, 0.0), (ntime, nlat, nlon))
    dims = ("Time", "altitude", "latitude", "longitude")
    ds = xr.Dataset(
        {
            "temp": (dims, temp4),
            "h2o_vap": (dims, q4),
            "h2o_ice": (dims, np.zeros(shape)),
            "u": (dims, np.zeros(shape)),
            "v": (dims, np.zeros(shape)),
            "w": (dims, rng.normal(0.0, 1e-4, shape)),
            "ps": (dims[:1] + dims[2:], ps),
            "ASR": (dims[:1] + dims[2:], asr),
            "OLR": (dims[:1] + dims[2:], olr),
            "aps": ("altitude", np.zeros(nlev)),
            "bps": ("altitude", p_lev / ps.flat[0]),
        },
        coords={
            "Time": np.arange(ntime, dtype=float) * 100.0,
            "altitude": -np.log10(p_lev / p_bot),
            "latitude": lat,
            "longitude": lon,
        },
    )
    ds.attrs["title"] = "meromixis synthetic diagfi (emulator columns; NOT a GCM run)"
    return ds


def make_mixing_test(k_true: float = 5.0, nlon: int = 48, nlat: int = 32, nlev: int = 30, ntime: int = 6,
                     temperature: float = 300.0, gravity: float = 12.43, metallicity: float = 1.0,
                     w_rms: float = 0.05, seed: int = 1) -> xr.Dataset:
    """Standardised dataset whose eddy tracer flux equals -rho K_true dq/dz.

    An isothermal atmosphere with a tracer decaying upward; vertical-velocity eddies w'
    (zero area mean on every level) displace it by a mixing length l = K_true / w_rms, so
    that q' = -l (w'/w_rms) dq/dz and <w'q'> = -K_true dq/dz.
    """
    from .composition import Composition

    rng = np.random.default_rng(seed)
    lon, lat = _grid(nlon, nlat)
    comp = Composition(metallicity)
    h = R_GAS * temperature / (comp.mu(0.0) * gravity)
    p_lev = np.geomspace(1.0e6, 1.0e2, nlev)
    z = -h * np.log(p_lev / p_lev[0])
    q_bar = 1.0e-4 * np.exp(-z / (2.0 * h))
    dq_dz = -q_bar / (2.0 * h)
    shape = (ntime, nlev, nlat, nlon)
    aw = np.cos(np.deg2rad(lat))[:, None] * np.ones(nlon)[None, :]
    aw = aw / aw.sum()
    # one horizontal pattern per time, coherent in the vertical (so interpolation is exact)
    pattern = rng.normal(0.0, 1.0, (ntime, 1, nlat, nlon))
    pattern -= (pattern * aw).sum(axis=(2, 3), keepdims=True)  # zero area mean
    pattern /= np.sqrt((pattern**2 * aw).sum(axis=(2, 3)).mean())  # unit area/time-mean variance
    pattern = np.broadcast_to(pattern, shape)
    w = w_rms * pattern
    q = q_bar[None, :, None, None] - (k_true / w_rms) * pattern * dq_dz[None, :, None, None]
    dims = ("time", "lev", "lat", "lon")
    ds = xr.Dataset(
        {
            "T": (dims, np.full(shape, temperature)),
            "p": (dims, np.broadcast_to(p_lev[None, :, None, None], shape).copy()),
            "q_vap": (dims, np.zeros(shape)),
            "chi": (dims, q),
            "w": (dims, w),
        },
        coords={"time": np.arange(ntime, dtype=float), "lev": np.arange(nlev), "lat": lat, "lon": lon},
        attrs={"gravity": gravity, "metallicity": metallicity, "lon_substellar": 0.0},
    )
    return ds
