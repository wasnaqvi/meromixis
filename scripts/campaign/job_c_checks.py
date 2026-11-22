"""Job C: four numerical checks that settle open items in docs/VERIFY.md.

C1  Saturation-curve consistency (VERIFY P3).
    The emulator defaults to Wagner-Pruss/Murphy-Koop; the Generic PCM evaluates a
    constant-latent-heat Clausius-Clapeyron curve for generic condensables.  If T*
    differs between them, the emulator is not a valid 1-D baseline for the GCM no
    matter how well its opacities are calibrated.  Measures dT* for both curves.

C2  Model-bottom sensitivity of the NH3 reservoir (VERIFY C5).
    The repo notes the NH3 quench level falls below the 1000 bar bottom, so f_deep
    is really "the equilibrium value at whatever p_bot I chose".  Since K_crit is
    proportional to 1/f_deep at fixed geometry, that choice propagates straight into
    the headline number.  Measures how both move as p_bot goes 100 -> 3000 bar.

C3  Sharp-front model against a numerical solve of the SAME equation.
    barrier.steady_state is an analytic sharp-front construction;
    barrier.solve_selfshielded integrates the self-shielded diffusion-loss equation
    numerically by fixed-point iteration on the photolysis rate.  Because NH3 is
    optically thick in its own band, this (not an optically thin reference) is the
    like-for-like test of the machinery the kill test rests on.

C4  Two-zone census (VERIFY E2).
    The fiducial column has two inhibited zones; the upper one dominates the series
    resistance, so if it is a grey-RT artefact the headline K_crit is wrong.  Counts
    zones and the resistance each contributes across the nuisance ensemble.

Usage: python job_c_checks.py <outdir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nuisance import METALLICITIES, draws, params  # noqa: E402


def check_saturation() -> dict:
    from meromixis.column import ColumnParams
    from meromixis.regimes import transition_width
    from meromixis.thermo import Saturation

    wp, cc = Saturation("wp"), Saturation("cc")
    out = {"metallicity": METALLICITIES.tolist(), "t_mid_wp": [], "t_mid_cc": [],
           "width_wp": [], "width_cc": [], "psat_ratio": {}}
    for t in (250.0, 300.0, 400.0, 500.0):
        out["psat_ratio"][str(int(t))] = float(cc.psat(t) / wp.psat(t))
    for z in METALLICITIES:
        lo_w, mid_w, hi_w = transition_width(float(z), ColumnParams(saturation=wp))
        lo_c, mid_c, hi_c = transition_width(float(z), ColumnParams(saturation=cc))
        out["t_mid_wp"].append(mid_w); out["t_mid_cc"].append(mid_c)
        out["width_wp"].append(hi_w - lo_w); out["width_cc"].append(hi_c - lo_c)
    d = np.array(out["t_mid_cc"]) - np.array(out["t_mid_wp"])
    out["dT_star_cc_minus_wp"] = d.tolist()
    out["dT_star_median"] = float(np.nanmedian(d))
    out["dT_star_max_abs"] = float(np.nanmax(np.abs(d)))
    return out


def check_pbot() -> dict:
    from meromixis.barrier import Mixing, Photolysis, ammonia_budget
    from meromixis.column import ColumnParams

    out = {"p_bot_bar": [], "results": {}}
    for z in (30.0, 100.0):
        rows = []
        for pb in (1.0e7, 3.0e7, 1.0e8, 3.0e8):
            b = ammonia_budget(ColumnParams(teq=265.0, metallicity=z, p_bot=pb),
                               Mixing(), Photolysis())
            rows.append({"p_bot_bar": pb / 1e5, "nh3_x": b.quench_inh.x,
                         "nh3_p_bar": b.quench_inh.p / 1e5,
                         "in_domain": bool(b.quench_inh.in_domain),
                         "k_crit": b.k_crit, "thermal_factor": b.thermal_factor,
                         "t_deep": None})
        out["results"][f"Z{int(z)}"] = rows
    return out


def check_sharp_front() -> dict:
    """Sharp-front analytic solution vs a numerical solve of the SAME equation.

    The first version of this check compared steady_state against a sink fixed at the
    optically thin rate.  That was not a like-for-like test: NH3 is optically thick in
    its own band (tau = 6e2 to 5e5 between 0.3 bar and 1 mbar at the quenched
    abundance), so the two were solving different problems and the comparison said
    nothing about the analytic model.  barrier.solve_selfshielded now solves the
    self-shielded equation numerically by fixed-point iteration on J(z), which IS the
    equation steady_state approximates.
    """
    from meromixis.barrier import (
        Mixing, Photolysis, band_mean, solve_selfshielded, steady_state)
    from meromixis.chemistry import quench
    from meromixis.column import ColumnParams, solve_column

    mix = Mixing()
    out = {}
    for Z in (15.0, 30.0, 100.0):
        for eta, phi in ((0.3, 1.0e17), (1.0, 1.0e18)):
            col = solve_column(ColumnParams(teq=255.0, metallicity=Z))
            q = quench(col, "NH3")
            phot = Photolysis(photon_flux=phi, eta=eta)
            sol = steady_state(col, q.x, q.p, mix, phot)
            p_num, f_num, info = solve_selfshielded(col, q.x, q.p, mix, phot)
            an_band = sol.mean()
            nu_band = band_mean(p_num, f_num)
            out[f"Z{int(Z)}_eta{eta}_phi{phi:.0e}"] = {
                "f_deep": float(q.x),
                "analytic_p_front_mbar": float(sol.p_front / 100.0),
                "numerical_p_front_mbar": float(info["p_front"] / 100.0),
                "p_front_ratio": float(sol.p_front / info["p_front"]) if info["p_front"] else None,
                "analytic_band_mean": float(an_band),
                "numerical_band_mean": float(nu_band),
                "band_mean_ratio": float(nu_band / an_band) if an_band > 0 else None,
                "analytic_column_loss": float(sol.column_loss),
                "numerical_column_loss": float(info["column_loss"]),
                "loss_ratio": float(info["column_loss"] / sol.column_loss) if sol.column_loss else None,
                "analytic_hidden_1mbar": bool(sol.hidden(1.0e3)),
                "numerical_hidden_1mbar": bool(info["p_front"] > 1.0e3),
                "tau_at_numerical_front": float(info["tau_at_front"]),
                "iterations": int(info["iterations"]),
                "converged": bool(info["converged"]),
            }
    vals = [v for v in out.values()]
    agree = [v["analytic_hidden_1mbar"] == v["numerical_hidden_1mbar"] for v in vals]
    ratios = [v["p_front_ratio"] for v in vals if v["p_front_ratio"]]
    out["_summary"] = {
        "cases": len(agree),
        "hidden_verdict_agrees": int(sum(agree)),
        "p_front_ratio_min": min(ratios) if ratios else None,
        "p_front_ratio_max": max(ratios) if ratios else None,
    }
    return out


def check_two_zones(n_draws: int = 96) -> dict:
    from meromixis.barrier import Mixing, kzz_profile, layer_mask, resistance
    from meromixis.column import INHIBITED, solve_column

    rows = draws(512)[:n_draws]
    mix = Mixing()
    res = {"metallicity": METALLICITIES.tolist(), "n_draws": int(n_draws),
           "zone_count_hist": {}, "upper_zone_resistance_fraction": {}}
    for z in METALLICITIES:
        counts, fracs = [], []
        for row in rows:
            col = solve_column(params(row), teq=265.0, metallicity=float(z))
            m = (col.regime == INHIBITED).astype(int)
            starts = np.nonzero(np.diff(np.r_[0, m]) == 1)[0]
            ends = np.nonzero(np.diff(np.r_[m, 0]) == -1)[0]
            counts.append(len(starts))
            if len(starts) >= 2:
                # p is top->bottom, so zone 0 is the UPPER (lowest pressure) zone
                k = kzz_profile(col, mix)
                tot = resistance(col, k, col.p[ends[-1]], col.p[starts[0]])
                up = resistance(col, k, col.p[ends[0]], col.p[starts[0]])
                fracs.append(up / tot if tot > 0 else np.nan)
        vals, cnt = np.unique(counts, return_counts=True)
        res["zone_count_hist"][f"Z{int(z)}"] = {int(v): int(c) for v, c in zip(vals, cnt)}
        res["upper_zone_resistance_fraction"][f"Z{int(z)}"] = (
            {"median": float(np.nanmedian(fracs)), "n": len(fracs)} if fracs else {"median": None, "n": 0}
        )
    return res


def main() -> int:
    outdir = Path(sys.argv[1]); outdir.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, fn in (("C1_saturation", check_saturation), ("C2_pbot", check_pbot),
                     ("C3_sharp_front", check_sharp_front), ("C4_two_zones", check_two_zones)):
        print(f"=== {name} ===", flush=True)
        out[name] = fn()
        print(json.dumps(out[name], indent=2, default=str)[:1500], flush=True)
    (outdir / "job_c_checks.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"wrote {outdir / 'job_c_checks.json'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
