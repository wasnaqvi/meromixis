"""Diagnostics for 3-D runs, built on the same column functions as the 1-D emulator.

Input is a standardised dataset (see :mod:`meromixis.pcm.io`): dims (time, lev, lat,
lon) with ``lev`` bottom-first; variables T, p, q_vap and optionally w.

Result I (cold-trap transition)
    :func:`coldtrap_fields`, :func:`transmission_summary`
Result II (ammonia barrier)
    :func:`inhibition_fields`, :func:`ledoux_n2`, :func:`kzz_flux_gradient`
Housekeeping
    :func:`convergence`
"""
from __future__ import annotations

import numpy as np
import xarray as xr

from .coldtrap import cold_trap_nd
from .composition import Composition, q_to_x
from .constants import R_GAS
from .inhibition import q_crit
from .thermo import Saturation

LEV_AXIS = 1  # (time, lev, lat, lon)


def _comp(ds) -> Composition:
    return Composition(float(ds.attrs["metallicity"]))


# ---- geometry ---------------------------------------------------------------------------

def area_weights(ds) -> xr.DataArray:
    w = np.cos(np.deg2rad(ds["lat"]))
    return (w / w.sum()).rename("area_weight")


def stellar_mu(ds) -> xr.DataArray:
    """Cosine of the stellar zenith angle for a synchronous rotator, dims (lat, lon)."""
    lon0 = float(ds.attrs.get("lon_substellar", 0.0))
    return (np.cos(np.deg2rad(ds["lat"])) * np.cos(np.deg2rad(ds["lon"] - lon0))).rename("mu")


def region_masks(ds, limb_half_width: float = 10.0) -> dict:
    """Boolean (lat, lon) masks.

    ``limb_east`` is the terminator 90 deg east of the substellar point (the evening
    limb for prograde rotation with eastward equatorial flow), ``limb_west`` the other one.
    """
    mu = stellar_mu(ds)
    lon0 = float(ds.attrs.get("lon_substellar", 0.0))
    east = np.sin(np.deg2rad(ds["lon"] - lon0)) > 0
    limb = np.abs(mu) < np.sin(np.deg2rad(limb_half_width))
    ones = xr.ones_like(mu, dtype=bool)
    return {
        "global": ones,
        "day": mu > 0,
        "night": mu <= 0,
        "substellar": mu > np.cos(np.deg2rad(30.0)),
        "antistellar": mu < -np.cos(np.deg2rad(30.0)),
        "terminator": limb,
        "limb_east": limb & east,
        "limb_west": limb & ~east,
    }


def wmean(da: xr.DataArray, ds, mask: xr.DataArray | None = None) -> xr.DataArray:
    """Area-weighted mean over lat and lon (NaNs and masked cells excluded)."""
    w = area_weights(ds) * xr.ones_like(ds["lon"], dtype=float)
    if mask is not None:
        w = w.where(mask, 0.0)
    return da.weighted(w).mean(("lat", "lon"))


# ---- vertical interpolation ---------------------------------------------------------------

def interp_logp(p: np.ndarray, f: np.ndarray, p_levels, axis: int = LEV_AXIS) -> np.ndarray:
    """Interpolate ``f`` to ``p_levels`` linearly in ln p along ``axis`` (p decreasing along it).

    Targets outside a column's pressure range give NaN.
    """
    p = np.moveaxis(np.asarray(p, float), axis, 0)
    f = np.moveaxis(np.asarray(f, float), axis, 0)
    lnp = np.log(p)
    out = np.full((len(p_levels),) + p.shape[1:], np.nan)
    for k, pt in enumerate(np.asarray(p_levels, float)):
        n_above = (p > pt).sum(axis=0)  # number of levels with higher pressure
        ok = (n_above >= 1) & (n_above <= p.shape[0] - 1)
        i0 = np.clip(n_above - 1, 0, p.shape[0] - 2)[None]
        l0 = np.take_along_axis(lnp, i0, 0)[0]
        l1 = np.take_along_axis(lnp, i0 + 1, 0)[0]
        f0 = np.take_along_axis(f, i0, 0)[0]
        f1 = np.take_along_axis(f, i0 + 1, 0)[0]
        w = (np.log(pt) - l0) / (l1 - l0)
        out[k] = np.where(ok, f0 + w * (f1 - f0), np.nan)
    return np.moveaxis(out, 0, axis)


def to_pressure(ds, name: str, p_levels) -> xr.DataArray:
    """A variable on pressure levels, dims (time, plev, lat, lon)."""
    p_levels = np.atleast_1d(np.asarray(p_levels, float))
    data = interp_logp(ds["p"].values, ds[name].values, p_levels)
    return xr.DataArray(
        data, dims=("time", "plev", "lat", "lon"),
        coords={"time": ds["time"], "plev": p_levels, "lat": ds["lat"], "lon": ds["lon"]}, name=name,
    )


def water_mole_fraction(ds) -> xr.DataArray:
    return xr.apply_ufunc(q_to_x, ds["q_vap"], _comp(ds).mu_dry).rename("x_h2o")


# ---- Result I: the cold trap --------------------------------------------------------------

def coldtrap_fields(ds, sat: Saturation | None = None, p_obs=(1.0e1, 1.0e3), p_min: float | None = None) -> xr.Dataset:
    """Column-by-column cold trap, dims (time, lat, lon).

    Two transmissions are returned because they answer different questions:

    ``transmission_sat``     min(1, x_sat,min / x_deep) from each column's own T(p): what a
                             1-D model would say given the 3-D temperatures.
    ``transmission_actual``  x(p_obs) / x_deep from the transported vapour itself (log-mean
                             over ``p_obs``): what the GCM actually delivers aloft.

    Their ratio measures horizontal transport, sub-saturation and re-evaporation, the
    intrinsically 3-D part of the answer.
    """
    sat = sat or Saturation()
    x = water_mole_fraction(ds)
    x_deep = x.isel(lev=0)
    p, t = ds["p"].values, ds["T"].values
    p_min = float(np.nanmin(p)) * 2.0 if p_min is None else p_min
    ct = cold_trap_nd(p, t, x_deep.values, axis=LEV_AXIS, sat=sat, p_min=p_min)
    dims, coords = ("time", "lat", "lon"), {k: ds[k] for k in ("time", "lat", "lon")}
    out = xr.Dataset(
        {
            "p_coldtrap": (dims, ct["p"]),
            "t_coldtrap": (dims, ct["T"]),
            "x_sat_limit": (dims, ct["x"]),
            "transmission_sat": (dims, ct["transmission"]),
            "x_deep": x_deep.drop_vars("lev", errors="ignore"),
        },
        coords=coords,
    )
    plevs = np.geomspace(p_obs[0], p_obs[1], 5)
    lx = np.log(np.maximum(interp_logp(p, x.values, plevs), 1e-300))
    x_obs = np.exp(np.nanmean(lx, axis=LEV_AXIS))
    out["x_obs"] = (dims, x_obs)
    out["transmission_actual"] = out["x_obs"] / out["x_deep"]
    out.attrs["p_obs"] = tuple(float(v) for v in p_obs)
    return out


def transmission_summary(ds, sat: Saturation | None = None, p_obs=(1.0e1, 1.0e3)) -> xr.Dataset:
    """Region-by-region cold-trap statistics, time-averaged: the numbers behind Result I.

    Means are geometric (transit depths respond to log abundance).
    """
    ct = coldtrap_fields(ds, sat, p_obs)
    masks = region_masks(ds)
    rows = {}
    for name, mask in masks.items():
        rows[name] = {
            "log10_x_obs": float(wmean(np.log10(ct["x_obs"]), ds, mask).mean("time")),
            "log10_transmission_actual": float(wmean(np.log10(ct["transmission_actual"]), ds, mask).mean("time")),
            "log10_transmission_sat": float(wmean(np.log10(ct["transmission_sat"]), ds, mask).mean("time")),
            "t_coldtrap": float(wmean(ct["t_coldtrap"], ds, mask).mean("time")),
            "p_coldtrap": float(np.exp(wmean(np.log(ct["p_coldtrap"]), ds, mask).mean("time"))),
        }
    regions = list(rows)
    out = xr.Dataset(
        {k: ("region", [rows[r][k] for r in regions]) for k in rows[regions[0]]}, coords={"region": regions}
    )
    out.attrs["limb_asymmetry_dex"] = float(
        out["log10_x_obs"].sel(region="limb_east") - out["log10_x_obs"].sel(region="limb_west")
    )
    return out


# ---- Result II: the stable layer and its permeability -------------------------------------

def water_gradient(ds) -> xr.DataArray:
    """``d ln q_vap / d ln p``: positive where the vapour increases downward.

    This is the gate that bounds the inhibited layer from above.  A compositionally
    stabilised layer requires a mean-molecular-weight gradient, and in this problem the
    mu gradient *is* the water gradient.  Above the cold trap the vapour is a well-mixed
    tracer with no sources or sinks, so ``q`` is constant, the mu gradient vanishes, and
    nothing there is compositionally stabilised however large ``q`` happens to be.

    Without this gate both stability tests run away to the model top in a steam
    atmosphere: ``q > q_cri`` is satisfied at every level once the well-mixed vapour
    exceeds a few percent, and the radiative stratosphere above the trap is stable on
    its own, so ``_contiguous_from_base`` never finds a level that fails and never
    stops.  Verified on the first 3-D runs, where an ungated mask flagged 100% of the
    area from the cloud base to 0.24 Pa.
    """
    q = np.maximum(ds["q_vap"].values, 1e-300)
    lnp = np.log(ds["p"].values)
    g = np.gradient(np.log(q), axis=LEV_AXIS) / np.gradient(lnp, axis=LEV_AXIS)
    return xr.DataArray(g, dims=ds["T"].dims, coords=ds["T"].coords, name="nabla_q")


def _contiguous_from_base(mask: xr.DataArray, above_base: xr.DataArray) -> xr.DataArray:
    """Keep only the part of ``mask`` connected upward to the cloud base.

    The lid is the layer sitting on the cloud deck (Leconte et al. 2024: "a stable
    layer forms between the level where q=q_cri and the dry troposphere below").  A
    bare stability test would also flag the radiative stratosphere, which is stable
    for reasons that have nothing to do with the water gradient, so we walk upward
    from the cloud base and stop at the first level that fails.
    """
    a = np.asarray(mask.values, bool)
    ab = np.asarray(above_base.values, bool)
    nlev = a.shape[LEV_AXIS]
    idx = np.arange(nlev).reshape((1, nlev) + (1,) * (a.ndim - 2))
    i_base = np.argmax(ab, axis=LEV_AXIS)  # first level above cloud base, from the bottom
    has_base = ab.any(axis=LEV_AXIS)
    i_base_b = np.expand_dims(i_base, LEV_AXIS)
    below = idx < i_base_b
    run = np.cumprod(np.where(below, True, a), axis=LEV_AXIS).astype(bool)
    out = run & ~below & np.expand_dims(has_base, LEV_AXIS)
    return xr.DataArray(out, dims=mask.dims, coords=mask.coords)


def inhibition_fields(ds, sat: Saturation | None = None, tol: float = 0.02,
                      criterion: str = "moist", grad_tol: float = 1.0e-3) -> xr.Dataset:
    """Where the column is compositionally lidded.

    ``criterion`` selects the stability test:

    ``"moist"`` (default)
        q > q_cri(T) above the cloud base.  This is the test Leconte et al. (2024)
        apply to the *actual* vapour -- "a stable layer forms between the level where
        q=q_cri and the dry troposphere below" -- and the one the 1-D emulator uses,
        where the rain-out profile is saturated above cloud base so actual and
        saturation vapour coincide.
    ``"ledoux"``
        Ledoux N^2 > 0 including the mu gradient, kept contiguous from the cloud base.
        Leconte et al. (2024) find the layer "can be largely subsaturated and remain
        stable", because "the medium does not need to be saturated for the vapor
        gradient to have a stabilizing effect"; where the layer is sub-saturated it is
        the Ledoux criterion on the actual profile, not the moist one, that applies
        (their Sect. 4.1: "one needs to use the Ledoux (1947) criterion instead").
    ``"union"``
        Either test, contiguous from the cloud base.  This is the recommended mask for
        a GCM run, because a real layer is moist-inhibited near its saturated top and
        Ledoux-stable through its sub-saturated interior.

    All three masks are always returned (``lid_moist``, ``lid_ledoux``, ``lid_union``)
    together with ``lid_agreement``, the area fraction on which moist and union agree,
    so the sensitivity of any derived number to the choice is visible rather than
    implicit.  ``lid`` and the top/base pressures follow ``criterion``.
    """
    if criterion not in ("moist", "ledoux", "union"):
        raise ValueError("criterion must be 'moist', 'ledoux' or 'union'")
    sat = sat or Saturation()
    comp = _comp(ds)
    x = water_mole_fraction(ds)
    qc = xr.apply_ufunc(lambda t: q_crit(t, comp, sat), ds["T"])
    above_base = x < x.isel(lev=0) * (1.0 - tol)
    near_base = above_base | above_base.shift(lev=-1, fill_value=False)
    # Bound the layer above by the cold trap: no water gradient, no compositional
    # stabilisation.  See :func:`water_gradient` for why this is not optional.
    # The one exemption is the cloud-base level itself (the last level not yet above
    # base, whose neighbour above is), where a one-sided gradient can read flat.
    # ``above_base.shift(lev=-1)`` alone is true at every level below the top, which
    # switched the gate off and let the lid run to the model top again.
    cloud_base = ~above_base & above_base.shift(lev=-1, fill_value=False)
    gradient = (water_gradient(ds) > grad_tol) | cloud_base
    window = near_base & gradient

    lid_moist = (ds["q_vap"] > qc) & window
    lid_ledoux = _contiguous_from_base((ledoux_n2(ds) > 0) & window, above_base)
    lid_union = _contiguous_from_base(lid_moist | lid_ledoux, above_base)
    lid = {"moist": lid_moist, "ledoux": lid_ledoux, "union": lid_union}[criterion]

    p = ds["p"]
    out = xr.Dataset(
        {
            "lid": lid,
            "lid_moist": lid_moist,
            "lid_ledoux": lid_ledoux,
            "lid_union": lid_union,
            "p_lid_base": p.where(lid).max("lev"),
            "p_lid_top": p.where(lid).min("lev"),
            "supercriticality": (ds["q_vap"].isel(lev=0) / qc.where(near_base).max("lev")).drop_vars("lev", errors="ignore"),
        }
    )
    out["has_lid"] = lid.any("lev")
    out["lid_area_fraction"] = wmean(out["has_lid"].astype(float), ds)
    for name, m in (("moist", lid_moist), ("ledoux", lid_ledoux), ("union", lid_union)):
        out[f"lid_area_fraction_{name}"] = wmean(m.any("lev").astype(float), ds)
    agree = (lid_moist.any("lev") == lid_union.any("lev")).astype(float)
    out["lid_agreement"] = wmean(agree, ds)
    out.attrs["criterion"] = criterion
    out.attrs["grad_tol"] = grad_tol
    return out


def ledoux_n2(ds) -> xr.DataArray:
    """Brunt-Vaisala frequency squared including the water mu-gradient, dims as T."""
    comp, g = _comp(ds), float(ds.attrs["gravity"])
    x = water_mole_fraction(ds).values
    t, p = ds["T"].values, ds["p"].values
    mu = comp.mu(x)
    lnp = np.log(p)
    d = lambda a: np.gradient(a, axis=LEV_AXIS) / np.gradient(lnp, axis=LEV_AXIS)  # noqa: E731
    nabla_ad = R_GAS / comp.cp_molar_mix(t, x)
    n2 = (g * g * mu / (R_GAS * t)) * (nabla_ad - d(np.log(t)) + d(np.log(mu)))
    return xr.DataArray(n2, dims=ds["T"].dims, coords=ds["T"].coords, name="N2", attrs={"units": "s-2"})


def kzz_flux_gradient(ds, tracer: str = "q_vap", p_levels=None, mask: xr.DataArray | None = None) -> xr.Dataset:
    """Effective vertical diffusivity of the *resolved* flow from the flux-gradient relation.

    On each pressure level, with <.> the area- and time-mean (optionally within ``mask``):

        F     = < rho w q > - < rho w > < rho q > / < rho >       eddy + overturning flux
        K_eff = - F / ( <rho> d q_bar / dz )

    This is the large-scale-circulation leak through the stable layer: one of the three
    leak paths of the barrier, and the one only a global model can measure.  It excludes
    sub-grid mixing (convective adjustment, turbulent diffusion), so it is a *lower bound*
    on the true permeability.  K is NaN where the mean gradient is too flat to define it.
    """
    if "w" not in ds:
        raise KeyError("vertical velocity 'w' is required (add it to diagfi.def)")
    comp, g = _comp(ds), float(ds.attrs["gravity"])
    if p_levels is None:
        pm = ds["p"].mean(("time", "lat", "lon")).values
        p_levels = np.geomspace(pm.max() * 0.9, pm.min() * 1.1, 2 * ds.sizes["lev"])
    p_levels = np.asarray(p_levels, float)
    w = to_pressure(ds, "w", p_levels)
    q = to_pressure(ds, tracer, p_levels)
    t = to_pressure(ds, "T", p_levels)
    qv = q if tracer == "q_vap" else to_pressure(ds, "q_vap", p_levels)
    mu = comp.mu(xr.apply_ufunc(q_to_x, qv, comp.mu_dry))
    rho = q["plev"] * mu / (R_GAS * t)
    aw = area_weights(ds) * xr.ones_like(ds["lon"], dtype=float)
    if mask is not None:
        aw = aw.where(mask, 0.0)
    valid = np.isfinite(w) & np.isfinite(q) & np.isfinite(t)

    def mean(a):
        return a.where(valid).weighted(aw).mean(("time", "lat", "lon"))

    rho_m = mean(rho)
    q_bar = mean(rho * q) / rho_m
    flux = mean(rho * w * q) - mean(rho * w) * q_bar
    h = R_GAS * mean(t) / (mean(mu) * g)
    lnp = np.log(p_levels)
    if bool((q_bar > 0).all()):  # differentiate ln q: exact for power laws, well behaved at the ends
        dq_dlnp = q_bar * np.gradient(np.log(q_bar.values), lnp)
    else:
        dq_dlnp = xr.DataArray(np.gradient(q_bar.values, lnp, edge_order=2), dims="plev", coords=q_bar.coords)
    dq_dz = -dq_dlnp / h  # d/dz = -(1/H) d/dln p
    k = -flux / (rho_m * dq_dz)
    flat = np.abs(dq_dz) * h < 1e-6 * np.abs(q_bar)
    return xr.Dataset({"kzz": k.where(~flat), "flux": flux, "dq_dz": dq_dz, "q_bar": q_bar, "rho": rho_m, "scale_height": h})


def layer_permeability(ds, tracer: str = "q_vap", criterion: str = "moist") -> dict:
    """K_eff of the resolved flow inside the lid: the number that enters the kill test.

    ``criterion`` is passed to :func:`inhibition_fields`.  The returned dict carries
    it, and the agreement between the moist and union masks, because K_layer is
    sensitive to how the lid is delimited: on synthetic output at 100x solar it moves
    by a factor of ~285 between a fully saturated field and one at 20% relative
    humidity, purely through the mask thinning out.
    """
    inh = inhibition_fields(ds, criterion=criterion)
    if not bool(inh["has_lid"].any()):
        return {"has_lid": False, "kzz_layer": np.nan, "p_top": np.nan, "p_base": np.nan,
                "criterion": criterion,
                "lid_agreement": float(inh["lid_agreement"].mean())}
    p_top = float(np.exp(wmean(np.log(inh["p_lid_top"]), ds).mean("time")))
    p_base = float(np.exp(wmean(np.log(inh["p_lid_base"]), ds).mean("time")))
    plevs = np.geomspace(p_base, p_top, 9)[1:-1]
    k = kzz_flux_gradient(ds, tracer, plevs)["kzz"]
    # resistances add in series -> harmonic mean across the layer
    k_pos = k.where(k > 0)
    k_layer = float(1.0 / (1.0 / k_pos).mean()) if bool(k_pos.notnull().any()) else np.nan
    return {"has_lid": True, "kzz_layer": k_layer, "p_top": p_top, "p_base": p_base,
            "lid_area_fraction": float(inh["lid_area_fraction"].mean()),
            "criterion": criterion,
            "lid_area_fraction_moist": float(inh["lid_area_fraction_moist"].mean()),
            "lid_area_fraction_union": float(inh["lid_area_fraction_union"].mean()),
            "lid_agreement": float(inh["lid_agreement"].mean())}


# ---- housekeeping -------------------------------------------------------------------------

def global_profile(ds, name: str, p_levels) -> xr.DataArray:
    """Global-mean, time-mean profile of a variable on pressure levels."""
    return wmean(to_pressure(ds, name, p_levels), ds).mean("time")


def convergence(ds, p_deep: float = 1.0e6, intheat: float = 0.0) -> dict:
    """Spin-up monitors: TOA imbalance and drift of the deep temperature.

    ``t_deep_drift`` is in K per unit of the file's time axis (PCM: planetary days).
    """
    out = {}
    if "asr" in ds and "olr" in ds:
        asr, olr = float(wmean(ds["asr"], ds).mean()), float(wmean(ds["olr"], ds).mean())
        out.update(asr=asr, olr=olr, toa_imbalance=asr + intheat - olr)
    td = wmean(to_pressure(ds, "T", [p_deep]).isel(plev=0), ds)
    out["t_deep"] = float(td.isel(time=-1))
    if ds.sizes["time"] > 2:
        tt = np.asarray(ds["time"].values, float)
        out["t_deep_drift"] = float(np.polyfit(tt - tt[0], td.values, 1)[0])
    return out
