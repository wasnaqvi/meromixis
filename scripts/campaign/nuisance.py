"""Shared definitions for the gate-0 ensemble.

The grey emulator has four radiative/interior nuisance parameters that `rcm1d` has
not yet calibrated (docs/VERIFY.md E1).  Every headline emulator number in the
README is quoted at their default values, so the honest question is not "what is
T*?" but "how much of T* is the radiation scheme?".  This module defines the prior
we sample, so that spread becomes a measured error bar instead of a caveat.

Ranges (all uniform, log-uniform where marked) and why:

kappa_th    3e-3 .. 3e-2 m2/kg, log  Thermal opacity at p_ref = 1 bar.  H2-He CIA
                                     plus the water continuum spans roughly a
                                     decade over 250-400 K; the default 1e-2 sits
                                     in the middle.
kappa_pexp  0.7 .. 1.3               kappa ~ p**kappa_pexp.  Exactly 1 for pure
                                     collision-induced absorption; lower when line
                                     opacity (which is closer to p**0) contributes.
gamma_v     0.1 .. 1.0, log          Visible/thermal opacity ratio.  < 1 is a
                                     "transparent in the optical" atmosphere, which
                                     an M-dwarf-irradiated H2 envelope should be,
                                     but the value is not pinned without a real
                                     shortwave calculation.
tint        30 .. 80 K               Internal temperature.  Sets how convective the
                                     cloud region is, hence the lid's base.

Everything else (composition closure, saturation curve, q_cri) is thermodynamic and
is NOT sampled: those are the parts of the answer that do not depend on the RT.
"""
from __future__ import annotations

import numpy as np

NUISANCE = {
    "kappa_th": (3.0e-3, 3.0e-2, "log"),
    "kappa_pexp": (0.7, 1.3, "lin"),
    "gamma_v": (0.1, 1.0, "log"),
    "tint": (30.0, 80.0, "lin"),
}
NUISANCE_NAMES = list(NUISANCE)

# Science grid.  Z <= 300 because the CH4/H2O/NH3 closure runs out of H2 near 420x.
METALLICITIES = np.array([3.0, 6.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 300.0])
TEQ_GRID = np.arange(200.0, 421.0, 5.0)

N_DRAWS = 512


def draws(n: int = N_DRAWS, seed: int = 20261122) -> np.ndarray:
    """(n, 4) Sobol sample of the nuisance prior, mapped to physical units.

    Sobol rather than uniform random: we want the marginal spread of T* to converge
    with a few hundred draws, not a few thousand.
    """
    from scipy.stats import qmc

    m = int(np.ceil(np.log2(max(n, 2))))
    u = qmc.Sobol(d=len(NUISANCE), scramble=True, seed=seed).random_base2(m)[:n]
    out = np.empty_like(u)
    for j, name in enumerate(NUISANCE_NAMES):
        lo, hi, kind = NUISANCE[name]
        out[:, j] = 10.0 ** (np.log10(lo) + u[:, j] * np.log10(hi / lo)) if kind == "log" else lo + u[:, j] * (hi - lo)
    return out


def params(row: np.ndarray, **kw):
    """ColumnParams for one nuisance draw."""
    from meromixis.column import ColumnParams

    return ColumnParams(**dict(zip(NUISANCE_NAMES, (float(v) for v in row))), **kw)


def default_row() -> np.ndarray:
    """The README's default values, for reference inside the ensemble."""
    from meromixis.column import ColumnParams

    d = ColumnParams()
    return np.array([getattr(d, n) for n in NUISANCE_NAMES])


def nonfinite_profile(col) -> bool:
    """True if the column integration produced any non-finite temperature.

    This is a genuine numerical failure, not a physical regime.  Deep temperatures
    routinely exceed the water critical point (647.096 K) -- the ensemble median
    deep temperature is ~1400 K -- and the integrator handles that correctly by
    taking the dry branch, because psat is infinite there so the gas is never
    saturated.  The failure happens in the narrow band just BELOW the critical
    point, where at high enrichment the water partial pressure still exceeds psat
    so the moist branch is selected, and the moist adiabat is singular as the
    latent heat vanishes.  It affects about 0.2% of the ensemble, all at Z = 300,
    and only the deepest one or two levels; the cold trap and the lid sit four
    orders of magnitude higher in pressure and are unaffected.

    Columns flagged here are reported and excluded rather than averaged in.
    """
    import numpy as np

    return bool((~np.isfinite(np.asarray(col.T, float))).any())
