"""Read Generic PCM output (``diagfi*.nc``) into the layout the diagnostics expect.

Standard layout
---------------
dims      : time, lev (index 0 = bottom, as in the PCM), lat, lon
variables : T [K], p [Pa] (4-D, mid-layer), q_vap [kg/kg]; optional q_ice, u, v, w [m/s, up],
            asr, olr [W m-2]
attrs     : gravity [m s-2], metallicity [x solar], lon_substellar [deg]

``VARMAP`` maps PCM names to standard names.  The defaults are the names the generic
physics writes with ``writediagfi`` for the legacy water cycle; the generic-condensable
route names its tracers after the species in ``traceur.def`` (pass ``varmap={"q_vap":
"<name>_vap", ...}``).  VERIFY P4: confirm names with ``ncdump -h diagfi.nc`` on your
first run and adjust here once.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

VARMAP = {
    "T": "temp",
    "ps": "ps",
    "p": "p",
    "q_vap": "h2o_vap",
    "q_ice": "h2o_ice",
    "u": "u",
    "v": "v",
    "w": "w",
    "asr": "ASR",
    "olr": "OLR",
    "aps": "aps",
    "bps": "bps",
}
DIMMAP = {"Time": "time", "altitude": "lev", "latitude": "lat", "longitude": "lon"}


def standardize(raw: xr.Dataset, gravity: float, metallicity: float, lon_substellar: float = 0.0,
                varmap: dict | None = None, w_is_omega: bool = False) -> xr.Dataset:
    """Rename, rebuild pressure and attach the metadata the diagnostics need."""
    vm = {**VARMAP, **(varmap or {})}
    raw = raw.rename({k: v for k, v in DIMMAP.items() if k in raw.dims or k in raw.coords})
    out = xr.Dataset(coords={k: raw[k] for k in ("time", "lev", "lat", "lon") if k in raw.coords})
    for std in ("T", "q_vap", "q_ice", "u", "v", "w", "asr", "olr"):
        if vm[std] in raw:
            out[std] = raw[vm[std]]
    if "T" not in out or "q_vap" not in out:
        raise KeyError(f"need at least temperature ('{vm['T']}') and vapour ('{vm['q_vap']}') in the file")
    if vm["p"] in raw:
        out["p"] = raw[vm["p"]]
    elif vm["aps"] in raw and vm["bps"] in raw and vm["ps"] in raw:
        aps = raw[vm["aps"]].rename({raw[vm["aps"]].dims[0]: "lev"})
        bps = raw[vm["bps"]].rename({raw[vm["bps"]].dims[0]: "lev"})
        out["p"] = (aps + bps * raw[vm["ps"]]).transpose("time", "lev", "lat", "lon")
    else:
        raise KeyError("cannot build pressure: need 'p' or ('aps', 'bps', 'ps')")
    out = out.transpose("time", "lev", "lat", "lon", ...)
    if w_is_omega and "w" in out:  # Pa/s -> m/s (hydrostatic, mu taken as dry)
        from ..composition import Composition
        from ..constants import R_GAS

        h = R_GAS * out["T"] / (Composition(metallicity).mu_dry * gravity)
        out["w"] = -out["w"] * h / out["p"]
    if float(out["p"].isel(lev=0).mean()) < float(out["p"].isel(lev=-1).mean()):
        out = out.isel(lev=slice(None, None, -1))  # enforce bottom-first
    out.attrs.update(gravity=float(gravity), metallicity=float(metallicity), lon_substellar=float(lon_substellar))
    return out


def open_run(paths, gravity: float, metallicity: float, lon_substellar: float = 0.0,
             varmap: dict | None = None, w_is_omega: bool = False, last_n: int | None = None) -> xr.Dataset:
    """Open one or several ``diagfi`` segments (concatenated along time) and standardise them.

    ``last_n`` keeps only the final n time samples (the equilibrated tail).
    """
    if isinstance(paths, (str, Path)):
        paths = sorted(Path().glob(str(paths))) if any(c in str(paths) for c in "*?[") else [Path(paths)]
    raws = [xr.open_dataset(p, decode_times=False) for p in paths]
    tdim = "Time" if "Time" in raws[0].dims else "time"
    raw = xr.concat(raws, dim=tdim, data_vars="minimal", coords="minimal", compat="override") if len(raws) > 1 else raws[0]
    ds = standardize(raw, gravity, metallicity, lon_substellar, varmap, w_is_omega)
    return ds.isel(time=slice(-last_n, None)) if last_n else ds


def read_rcm1d_profile(path) -> tuple[np.ndarray, np.ndarray]:
    """(p [Pa], T [K]) from a 1-D ``diagfi.nc`` (last time step), for :func:`calibrate_grey`."""
    raw = xr.open_dataset(path, decode_times=False)
    t = raw[VARMAP["T"]].isel(Time=-1).squeeze().values
    if VARMAP["p"] in raw:
        p = raw[VARMAP["p"]].isel(Time=-1).squeeze().values
    else:
        p = (raw[VARMAP["aps"]] + raw[VARMAP["bps"]] * raw[VARMAP["ps"]].isel(Time=-1).squeeze()).values
    return np.asarray(p, float), np.asarray(t, float)
