"""Job B: the ammonia-barrier kill test over its real uncertainty space.

The README quotes K_crit at one point in a four-dimensional space: one grey column,
one hand-estimated Phi_UV, one eta, one assumed K_layer.  This job maps
K_crit(Z) over the ensemble of grey columns AND over (Phi_UV, eta), so the
"barrier survives in roughly 12-100x solar" claim becomes a measured region with a
stated confidence instead of a single curve.

Also records, per (draw, Z):
  * thermal_factor and co_enhancement  -- the Leconte+2024 route, which is
    independent of Phi_UV and eta, so its spread is pure RT uncertainty;
  * in_domain for the NH3 quench level  -- VERIFY C5: when False, f_deep is just
    the model-bottom equilibrium value and K_crit inherits that arbitrariness;
  * the number of separate inhibited zones -- VERIFY E2.

Columns are solved once per (draw, Z) and reused across (Phi_UV, eta), which is
what makes the full grid affordable.

Usage: python job_b_barrier.py <task_id> <n_tasks> <outdir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nuisance import (  # noqa: E402
    METALLICITIES, NUISANCE_NAMES, draws, params, nonfinite_profile)
from uvflux import photospheric_photon_flux  # noqa: E402

N_DRAWS_B = 256
ETA_GRID = np.array([0.03, 0.1, 0.3, 1.0])
# Phi_UV spanning the blackbody photospheric floor up 4 decades (brackets both the
# README hand estimate, ~1.1 dex above the floor, and an active M dwarf).
PHI_FLOOR = photospheric_photon_flux(300.0)
PHI_GRID = PHI_FLOOR * 10.0 ** np.array([0.0, 1.0, 2.0, 3.0, 4.0])
# Teq held at a value inside the cold-trapped+lidded regime for all Z, so that a lid
# exists to test.  Chosen 15 K below the ensemble-median transition.
TEQ_TEST = 265.0


def one_draw(row: np.ndarray) -> dict:
    from meromixis.barrier import Mixing, Photolysis, ammonia_budget
    from meromixis.column import INHIBITED, solve_column

    base = params(row)
    shp = (len(METALLICITIES), len(PHI_GRID), len(ETA_GRID))
    out = {k: np.full(shp, np.nan) for k in
           ("k_crit", "k_crit_simple", "hidden_inh", "hidden_ref", "hidden_by_barrier")}
    out2 = {k: np.full(len(METALLICITIES), np.nan) for k in
            ("thermal_factor", "co_enhancement", "nh3_quench_x", "nh3_quench_p",
             "nh3_in_domain", "co_in_domain", "co_quench_p", "lid_dz", "lid_p_top",
             "lid_p_base", "n_zones", "supercriticality", "t_deep",
             "nonfinite")}
    mix = Mixing()
    for i, z in enumerate(METALLICITIES):
        col_i = solve_column(base, teq=TEQ_TEST, metallicity=float(z), inhibition=True)
        col_r = solve_column(base, teq=TEQ_TEST, metallicity=float(z), inhibition=False)
        layer = col_i.stable_layer
        m = (col_i.regime == INHIBITED).astype(int)
        out2["n_zones"][i] = int(np.sum(np.diff(np.r_[0, m]) == 1))
        out2["t_deep"][i] = float(col_i.T[-1])
        out2["nonfinite"][i] = float(nonfinite_profile(col_i) or nonfinite_profile(col_r))
        if out2["nonfinite"][i]:
            # Numerical failure in the deepest level(s): skip the kill test for
            # this (draw, Z) and report the count rather than averaging in a
            # meaningless number.  ~0.2% of the ensemble, all at Z = 300.
            continue
        if layer.exists:
            out2["lid_dz"][i] = layer.thickness
            out2["lid_p_top"][i] = layer.p_top
            out2["lid_p_base"][i] = layer.p_base
            out2["supercriticality"][i] = layer.supercriticality
        first = True
        for j, phi in enumerate(PHI_GRID):
            for k, eta in enumerate(ETA_GRID):
                b = ammonia_budget(base, mix, Photolysis(photon_flux=float(phi), eta=float(eta)),
                                   col_inh=col_i, col_ref=col_r)
                out["k_crit"][i, j, k] = b.k_crit
                out["k_crit_simple"][i, j, k] = b.k_crit_simple
                out["hidden_inh"][i, j, k] = float(b.sol_inh.hidden())
                out["hidden_ref"][i, j, k] = float(b.sol_ref.hidden())
                out["hidden_by_barrier"][i, j, k] = float(b.hidden_by_barrier)
                if first:  # Phi/eta-independent quantities
                    out2["thermal_factor"][i] = b.thermal_factor
                    out2["co_enhancement"][i] = b.co_enhancement
                    out2["nh3_quench_x"][i] = b.quench_inh.x
                    out2["nh3_quench_p"][i] = b.quench_inh.p
                    out2["nh3_in_domain"][i] = float(b.quench_inh.in_domain)
                    out2["co_in_domain"][i] = float(b.co_inh.in_domain)
                    out2["co_quench_p"][i] = b.co_inh.p
                    first = False
    return {**out, **out2}


def main() -> int:
    task, ntask, outdir = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows = draws(N_DRAWS_B)
    idx = np.array_split(np.arange(len(all_rows)), ntask)[task]
    k3 = ("k_crit", "k_crit_simple", "hidden_inh", "hidden_ref", "hidden_by_barrier")
    k1 = ("thermal_factor", "co_enhancement", "nh3_quench_x", "nh3_quench_p",
          "nh3_in_domain", "co_in_domain", "co_quench_p", "lid_dz", "lid_p_top",
          "lid_p_base", "n_zones", "supercriticality", "t_deep", "nonfinite")
    acc = {k: [] for k in k3 + k1}
    for n, i in enumerate(idx):
        r = one_draw(all_rows[i])
        for k in acc:
            acc[k].append(r[k])
        print(f"[{task}] draw {i} ({n + 1}/{len(idx)}) done", flush=True)
    ds = xr.Dataset(
        {**{k: (("draw", "metallicity", "phi_uv", "eta"), np.array(acc[k])) for k in k3},
         **{k: (("draw", "metallicity"), np.array(acc[k])) for k in k1}},
        coords={"draw": idx, "metallicity": METALLICITIES, "phi_uv": PHI_GRID, "eta": ETA_GRID},
    )
    for j, nm in enumerate(NUISANCE_NAMES):
        ds[nm] = ("draw", all_rows[idx, j])
    ds.attrs.update(teq_test=TEQ_TEST, phi_floor_blackbody=PHI_FLOOR,
                    description="ammonia-barrier kill test over grey-RT and photochemical uncertainty")
    ds.to_netcdf(outdir / f"job_b_{task:04d}.nc")
    print(f"[{task}] wrote {outdir / f'job_b_{task:04d}.nc'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
