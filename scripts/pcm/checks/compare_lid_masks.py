"""Section-8 implementation: do the three masks behave as Leconte et al. 2024 describe?
import os

Leconte et al. (2024) find the stable layer "can be largely subsaturated and remain
stable", because "the medium does not need to be saturated for the vapor gradient to
have a stabilizing effect".  So as the vapour field is dried, the moist mask (q>q_cri)
should collapse while the Ledoux mask, which tests the actual mu gradient, should hold.
"""
import numpy as np

from meromixis import diagnostics3d as d3
from meromixis.column import ColumnParams, solve_column
from meromixis.pcm.io import standardize
from meromixis.synthetic import make_synthetic_run

for Z in (15.0, 30.0, 100.0):
    raw = make_synthetic_run(ColumnParams(teq=255.0, metallicity=Z), contrast=0.3)
    ds = standardize(raw, gravity=12.43, metallicity=Z)
    col = solve_column(ColumnParams(teq=255.0, metallicity=Z))
    sc = col.stable_layer.supercriticality if col.stable_layer.exists else float("nan")
    print(f"\n== Z = {Z:g}x solar   1-D q_deep/q_cri = {sc:.2f} ==")
    print(f"{'RH':>5} {'moist':>8} {'ledoux':>8} {'union':>8} {'agree':>7} "
          f"{'K(moist)':>10} {'K(union)':>10}")
    for rh in (1.0, 0.75, 0.5, 0.3, 0.2):
        d = ds.copy(deep=True)
        d["q_vap"] = ds["q_vap"] * rh
        inh = d3.inhibition_fields(d, criterion="union")
        am = float(inh["lid_area_fraction_moist"].mean())
        al = float(inh["lid_area_fraction_ledoux"].mean())
        au = float(inh["lid_area_fraction_union"].mean())
        ag = float(inh["lid_agreement"].mean())
        km = d3.layer_permeability(d, criterion="moist").get("kzz_layer", np.nan)
        ku = d3.layer_permeability(d, criterion="union").get("kzz_layer", np.nan)
        f = lambda v: f"{v:.3g}" if v is not None and np.isfinite(v) else "nan"
        print(f"{rh:>5.2f} {am:>8.3f} {al:>8.3f} {au:>8.3f} {ag:>7.3f} {f(km):>10} {f(ku):>10}")
