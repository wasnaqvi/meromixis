"""Reproduce the README's quoted emulator numbers at default parameters."""
import os
import numpy as np
from meromixis.column import ColumnParams, solve_column
from meromixis.regimes import transition_width
from meromixis.barrier import Mixing, Photolysis, ammonia_budget

print("== README claim: transition centred ~245 K at 3x solar, ~282 K at 100x; width ~20-25 K ==")
for z in (3.0, 100.0):
    lo, mid, hi = transition_width(z, ColumnParams())
    print(f"  Z={z:5.0f}x  T*(0.1,0.5,0.9) = {lo:.1f} {mid:.1f} {hi:.1f}   width={hi-lo:.1f} K")

print("\n== README claim: stable layer exists above ~11-17x solar ==")
for z in (8.0, 11.0, 14.0, 17.0, 20.0):
    col = solve_column(ColumnParams(metallicity=z))
    print(f"  Z={z:5.0f}x  lid={col.stable_layer.exists}  cold_trap={col.cold_trap.exists}")

print("\n== README claim: K_crit ~3 at 15x, ~1.5 at 30x, ~0.2 at 100x, <1e-3 at 300x ==")
print("   (README: hand UV flux 1e13 cm-2 s-1 = 1e17 m-2 s-1; 'fully irreversible' => eta=1)")
for z in (15.0, 30.0, 100.0, 300.0):
    b = ammonia_budget(ColumnParams(metallicity=z), Mixing(),
                       Photolysis(photon_flux=1.0e17, eta=1.0))
    print(f"  Z={z:5.0f}x  k_crit={b.k_crit:11.4g}  k_crit_simple={b.k_crit_simple:11.4g}"
          f"  thermal={b.thermal_factor:.3f}  coX={b.co_enhancement:.3g}"
          f"  in_domain={b.quench_inh.in_domain}")

print("\n== README claim: two inhibited zones at 0.13-0.26 bar and 4-7.5 bar (fiducial) ==")
from meromixis.column import INHIBITED
col = solve_column(ColumnParams())
m = (col.regime == INHIBITED).astype(int)
st = np.nonzero(np.diff(np.r_[0, m]) == 1)[0]; en = np.nonzero(np.diff(np.r_[m, 0]) == -1)[0]
print(f"  default column (teq=255, Z=100): {len(st)} zone(s)")
for a, b_ in zip(st, en):
    print(f"    {col.p[a]/1e5:.4g} - {col.p[b_]/1e5:.4g} bar")
