"""Is NH3 optically thick in its own photolysis band? The sharp-front model assumes so."""
import os
import numpy as np
from meromixis.barrier import NH3_PHOTOLYSIS as P, Mixing
from meromixis.chemistry import quench
from meromixis.column import ColumnParams, solve_column

for Z in (15.0, 30.0, 100.0):
    col = solve_column(ColumnParams(teq=255.0, metallicity=Z))
    q = quench(col, "NH3")
    # vertical NH3 column above each level, assuming the quenched mixing ratio below
    # the photolysis front (which is what the sharp-front solution asserts)
    o = np.argsort(col.z)
    z, n, p = col.z[o], col.n[o], col.p[o]
    good = np.isfinite(z) & np.isfinite(n)
    z, n, p = z[good], n[good], p[good]
    N_above = np.concatenate([np.cumsum((0.5*(n[:-1]+n[1:])*np.diff(z))[::-1])[::-1], [0.0]])
    tau = P.sigma * q.x * N_above
    for ptarget in (P.p_floor, 1.0e3, 1.0e2):
        i = int(np.argmin(np.abs(p - ptarget)))
        print(f"Z={Z:5.0f}x  p={p[i]/100:9.3g} mbar   f_NH3={q.x:.3g}   tau_NH3 = {tau[i]:.3g}")
    print()
