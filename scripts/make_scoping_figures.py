#!/usr/bin/env python
"""Regenerate the emulator scoping figures in docs/figures/.

These are design aids: they show where the 1-D grey emulator puts the two boundaries so
the GCM grid can be placed sensibly.  They are not results and must not be quoted as such.

    python scripts/make_scoping_figures.py [--fast]
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import yaml

from meromixis import ColumnParams, solve_column
from meromixis.barrier import Mixing, Photolysis, ammonia_budget
from meromixis.plotting import plot_column, plot_killtest, plot_regime_map
from meromixis.regimes import regime_map

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"


def main(fast: bool = False):
    OUT.mkdir(parents=True, exist_ok=True)
    base = ColumnParams()

    # 1. column anatomy at the fiducial point
    col, ref = solve_column(base), solve_column(base, inhibition=False)
    plot_column(col, ref).savefig(OUT / "emulator_column.png")

    # 2. regime map
    teqs = np.arange(200.0, 420.1, 20.0 if fast else 5.0)
    mets = np.geomspace(1.0, 300.0, 8 if fast else 30)
    ds = regime_map(teqs, mets, base)
    ds.to_netcdf(OUT / "emulator_regime_map.nc")
    targets = yaml.safe_load((ROOT / "experiments" / "targets.yaml").read_text())["targets"]
    plot_regime_map(ds, targets).savefig(OUT / "emulator_regime_map.png")

    # 3. kill test
    zs = np.geomspace(12.0, 300.0, 6 if fast else 16)
    k_crit = {}
    for eta in (1.0, 0.1):
        k_crit[eta] = [ammonia_budget(replace(base, metallicity=float(z)), Mixing(), Photolysis(eta=eta)).k_crit for z in zs]
    ks = np.geomspace(1e-3, 1e2, 12 if fast else 40)
    fronts = {}
    for z in (15.0, 30.0, 100.0):
        p = replace(base, metallicity=z)
        ci, cr = solve_column(p), solve_column(p, inhibition=False)
        fronts[z] = (ks, [ammonia_budget(p, Mixing(k_layer=float(k)), Photolysis(eta=1.0), ci, cr).sol_inh.p_front for k in ks])
    plot_killtest(zs, k_crit, fronts).savefig(OUT / "emulator_killtest.png")
    print(f"wrote figures to {OUT}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="coarse grids (seconds instead of a minute)")
    main(ap.parse_args().fast)
