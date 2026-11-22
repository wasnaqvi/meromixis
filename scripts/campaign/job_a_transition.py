"""Job A: the cold-trap transition under radiative uncertainty.

For one slice of the nuisance ensemble, measure for every metallicity
  T*(0.1), T*(0.5), T*(0.9)   -- cold-trap transmission crossings, and
  whether a lid exists at the transition and on the cold side,
then write one NetCDF per array task.

Scientific point: the README quotes T* ~ 245 K at 3x solar and ~282 K at 100x from
a single uncalibrated grey column.  The ensemble turns those into distributions, so
the GCM forcing nodes can be placed to straddle the *uncertainty*, not a point.

Usage: python job_a_transition.py <task_id> <n_tasks> <outdir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nuisance import (  # noqa: E402
    METALLICITIES, NUISANCE_NAMES, draws, params, nonfinite_profile)


def one_draw(row: np.ndarray) -> dict:
    from meromixis.column import solve_column
    from meromixis.inhibition import metallicity_crit
    from meromixis.regimes import transition_width

    base = params(row)
    nz = len(METALLICITIES)
    out = {k: np.full(nz, np.nan) for k in
           ("t_lo", "t_mid", "t_hi", "width", "lid_at_mid", "t_coldtrap_mid",
            "p_coldtrap_mid", "z_crit_at_mid", "t_deep_mid", "n_inhibited_zones",
            "nonfinite")}
    for i, z in enumerate(METALLICITIES):
        lo, mid, hi = transition_width(float(z), base)
        out["t_lo"][i], out["t_mid"][i], out["t_hi"][i] = lo, mid, hi
        out["width"][i] = hi - lo
        if not np.isfinite(mid):
            continue
        # column just on the cold-trapped side of the transition
        col = solve_column(base, teq=float(mid) - 5.0, metallicity=float(z))
        ct, layer = col.cold_trap, col.stable_layer
        out["lid_at_mid"][i] = float(layer.exists)
        out["t_coldtrap_mid"][i] = ct.T if ct.exists else np.nan
        out["p_coldtrap_mid"][i] = ct.p if ct.exists else np.nan
        out["t_deep_mid"][i] = float(col.T[-1])
        out["nonfinite"][i] = float(nonfinite_profile(col))
        # thermodynamic lid threshold at the cold-trap temperature (no RT involved)
        if ct.exists:
            out["z_crit_at_mid"][i] = metallicity_crit(ct.T, base.saturation)
        # VERIFY E2: how many separate inhibited zones does this column have?
        from meromixis.column import INHIBITED

        m = (col.regime == INHIBITED).astype(int)
        out["n_inhibited_zones"][i] = int(np.sum(np.diff(np.r_[0, m]) == 1))
    return out


def main() -> int:
    task, ntask, outdir = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows = draws()
    idx = np.array_split(np.arange(len(all_rows)), ntask)[task]
    recs = {k: [] for k in
            ("t_lo", "t_mid", "t_hi", "width", "lid_at_mid", "t_coldtrap_mid",
             "p_coldtrap_mid", "z_crit_at_mid", "t_deep_mid", "n_inhibited_zones",
             "nonfinite")}
    for n, i in enumerate(idx):
        r = one_draw(all_rows[i])
        for k in recs:
            recs[k].append(r[k])
        print(f"[{task}] draw {i} ({n + 1}/{len(idx)}) done", flush=True)
    ds = xr.Dataset(
        {k: (("draw", "metallicity"), np.array(v)) for k, v in recs.items()},
        coords={"draw": idx, "metallicity": METALLICITIES},
    )
    for j, nm in enumerate(NUISANCE_NAMES):
        ds[nm] = ("draw", all_rows[idx, j])
    ds.attrs["description"] = "cold-trap transition under grey-RT nuisance uncertainty"
    ds.to_netcdf(outdir / f"job_a_{task:04d}.nc")
    print(f"[{task}] wrote {outdir / f'job_a_{task:04d}.nc'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
