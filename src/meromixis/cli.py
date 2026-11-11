"""Command-line interface: ``meromixis <command> --help``."""
from __future__ import annotations

import argparse
import sys

import numpy as np


def _params(a):
    from .column import ColumnParams

    return ColumnParams(teq=a.teq, metallicity=a.metallicity, tint=a.tint, inhibition=not getattr(a, "no_inhibition", False))


def _add_column_args(p):
    p.add_argument("--teq", type=float, default=255.0, help="effective equilibrium temperature, K (albedo included)")
    p.add_argument("--metallicity", type=float, default=100.0, help="x solar")
    p.add_argument("--tint", type=float, default=50.0, help="internal temperature, K")


def cmd_column(a):
    from .column import REGIME_NAMES, solve_column

    col = solve_column(_params(a))
    ct, layer = col.cold_trap, col.stable_layer
    print(f"deep water: x = {col.comp.x_h2o_deep:.4g} mol/mol, q = {col.comp.q_h2o_deep:.4g} kg/kg, mu_dry = {col.comp.mu_dry * 1e3:.3f} g/mol")
    if ct.exists:
        print(f"cold trap : p = {ct.p / 100:.3g} mbar, T = {ct.T:.1f} K, x_above = {ct.x:.3g}, transmission = {ct.transmission:.3g}")
    else:
        print("cold trap : none (water never saturates: steam atmosphere)")
    if layer.exists:
        print(f"lid       : {layer.p_top / 1e5:.3g}-{layer.p_base / 1e5:.3g} bar, dT = {layer.delta_t:.1f} K, "
              f"dz = {layer.thickness / 1e3:.1f} km, q_deep/q_cri = {layer.supercriticality:.2f}")
    else:
        print("lid       : none")
    change = np.nonzero(np.diff(col.regime))[0]
    edges = [col.p[0], *col.p[change + 1], col.p[-1]]
    for i, r in enumerate([col.regime[0], *col.regime[change + 1]]):
        print(f"  {edges[i] / 1e5:10.3g} - {edges[i + 1] / 1e5:<10.3g} bar  {REGIME_NAMES[int(r)]}")
    print(f"T(1 bar) = {col.at(1e5):.1f} K, T(100 bar) = {col.at(1e7):.1f} K")
    if a.save:
        col.to_dataset().to_netcdf(a.save)
        print(f"saved {a.save}")


def cmd_transition(a):
    from .column import ColumnParams
    from .regimes import transition_width

    base = ColumnParams(tint=a.tint)
    print("metallicity   T(0.1)   T(0.5)   T(0.9)   width   [K, effective Teq; grey emulator]")
    for z in a.metallicity:
        lo, mid, hi = transition_width(z, base)
        print(f"{z:11g} {lo:8.1f} {mid:8.1f} {hi:8.1f} {hi - lo:7.1f}")


def cmd_killtest(a):
    from .barrier import Mixing, Photolysis, ammonia_budget

    b = ammonia_budget(_params(a), Mixing(k_layer=a.k_layer, k_rad=a.k_rad),
                       Photolysis(photon_flux=a.photon_flux, eta=a.eta))
    yr = 3.156e7
    print(f"quenched NH3      : {b.quench_ref.x:.3g} (no lid) -> {b.quench_inh.x:.3g} (lid); thermal factor {b.thermal_factor:.2f}"
          + ("" if b.quench_inh.in_domain else "   [quench level below model bottom]"))
    print(f"quenched CO       : {b.co_ref.x:.3g} -> {b.co_inh.x:.3g} (x{b.co_enhancement:.3g})")
    print(f"K_crit            : {b.k_crit:.3g} m2/s (series)   {b.k_crit_simple:.3g} m2/s (L dz / n f)")
    print(f"assumed K_layer   : {b.k_layer:.3g} m2/s   -> margin K_crit/K_layer = {b.margin:.3g}")
    print(f"photolysis front  : {b.sol_ref.p_front / 100:.3g} mbar (no lid)   {b.sol_inh.p_front / 100:.3g} mbar (lid)")
    print(f"timescales        : NH3 above the lid empties in {b.nh3_depletion_time / yr:.3g} yr; CO refills in {b.co_fill_time / yr:.3g} yr")
    print(f"verdict           : {b.verdict()}")
    print("(grey emulator + hand-estimated UV flux: a scoping number, not a result)")


def cmd_nodes(a):
    from .pcm.grid import teq_nodes

    print(" ".join(f"{t:g}" for t in teq_nodes(a.t_star, a.lo, a.hi, a.n, a.half_width)))


def cmd_pcm_init(a):
    from .pcm.grid import write_all

    dirs = write_all(a.config, a.out, with_1d=not a.no_1d)
    print(f"wrote {len(dirs)} run directories under {a.out} (table: {a.out}/runs.tsv)")


def _open(a):
    from .pcm.io import open_run

    varmap = dict(kv.split("=", 1) for kv in (a.varmap or []))
    return open_run(a.files if len(a.files) > 1 else a.files[0], a.gravity, a.metallicity, a.lon_substellar,
                    varmap, a.w_is_omega, a.last_n)


def _add_run_args(p):
    p.add_argument("files", nargs="+", help="diagfi segment(s)")
    p.add_argument("--gravity", type=float, required=True)
    p.add_argument("--metallicity", type=float, required=True)
    p.add_argument("--lon-substellar", type=float, default=0.0)
    p.add_argument("--varmap", nargs="*", help="standard=pcm name pairs, e.g. q_vap=H2O_vap")
    p.add_argument("--w-is-omega", action="store_true", help="vertical velocity is in Pa/s")
    p.add_argument("--last-n", type=int, default=None, help="use only the last n time samples")


def cmd_diag(a):
    from . import diagnostics3d as d3

    ds = _open(a)
    summ = d3.transmission_summary(ds)
    print("== Result I: cold trap ==")
    print(summ.to_dataframe().round(3).to_string())
    print(f"limb asymmetry (east - west): {summ.attrs['limb_asymmetry_dex']:+.3f} dex")
    print("\n== Result II: stable layer ==")
    if "w" in ds:
        for k, v in d3.layer_permeability(ds).items():
            print(f"{k:18s} {v}")
    else:
        inh = d3.inhibition_fields(ds)
        print(f"lid area fraction  {float(inh['lid_area_fraction'].mean()):.3f}   (no 'w' in file: K_eff not computed)")
    if a.out:
        summ.to_netcdf(a.out)


def cmd_monitor(a):
    from . import diagnostics3d as d3

    c = d3.convergence(_open(a), a.p_deep, a.intheat)
    print(" ".join(f"{k}={v:.5g}" for k, v in c.items()))


def cmd_synth(a):
    from .synthetic import make_synthetic_run

    make_synthetic_run(_params(a), contrast=a.contrast).to_netcdf(a.out)
    print(f"wrote {a.out} (synthetic, emulator columns; gravity 12.43, metallicity {a.metallicity:g})")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="meromixis", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("column", help="solve and summarise one emulator column")
    _add_column_args(p)
    p.add_argument("--no-inhibition", action="store_true")
    p.add_argument("--save", help="write the profile to NetCDF")
    p.set_defaults(func=cmd_column)

    p = sub.add_parser("transition", help="location and width of the cold-trap transition (emulator)")
    p.add_argument("--metallicity", type=float, nargs="+", default=[3, 10, 30, 100, 300])
    p.add_argument("--tint", type=float, default=50.0)
    p.set_defaults(func=cmd_transition)

    p = sub.add_parser("killtest", help="ammonia-barrier kill test for one planet (emulator)")
    _add_column_args(p)
    p.add_argument("--k-layer", type=float, default=1.0, help="m2/s in the stable layer")
    p.add_argument("--k-rad", type=float, default=100.0, help="m2/s in radiative zones")
    p.add_argument("--eta", type=float, default=0.3, help="irreversible fraction of NH3 photolysis")
    p.add_argument("--photon-flux", type=float, default=1e17, help="global-mean UV photons m-2 s-1")
    p.set_defaults(func=cmd_killtest)

    p = sub.add_parser("nodes", help="forcing labels clustered around a predicted transition")
    p.add_argument("--t-star", type=float, required=True)
    p.add_argument("--lo", type=float, default=250.0)
    p.add_argument("--hi", type=float, default=600.0)
    p.add_argument("-n", type=int, default=10)
    p.add_argument("--half-width", type=float, default=30.0)
    p.set_defaults(func=cmd_nodes)

    p = sub.add_parser("pcm-init", help="write Generic PCM run directories from an experiment YAML")
    p.add_argument("config")
    p.add_argument("out")
    p.add_argument("--no-1d", action="store_true", help="skip the rcm1d twins")
    p.set_defaults(func=cmd_pcm_init)

    p = sub.add_parser("diag", help="cold-trap and stable-layer diagnostics of a 3-D run")
    _add_run_args(p)
    p.add_argument("--out", help="write the summary to NetCDF")
    p.set_defaults(func=cmd_diag)

    p = sub.add_parser("monitor", help="one-line spin-up monitor (energy balance, deep-T drift)")
    _add_run_args(p)
    p.add_argument("--p-deep", type=float, default=1.0e6)
    p.add_argument("--intheat", type=float, default=0.0)
    p.set_defaults(func=cmd_monitor)

    p = sub.add_parser("synth", help="write a synthetic diagfi file to exercise the pipeline")
    _add_column_args(p)
    p.add_argument("--contrast", type=float, default=0.3)
    p.add_argument("--out", default="diagfi_synthetic.nc")
    p.set_defaults(func=cmd_synth)

    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
