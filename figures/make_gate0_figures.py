"""Figures for the gate-0 scoping campaign (A&A style).

    python figures/make_gate0_figures.py <gate0_run_dir> <out_dir>

Reads the Slurm array output of scripts/campaign/job_{a,b,c} and writes
PDF + PNG + SVG for each figure.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exoplot as xp  # noqa: E402

J = "aanda"
PCT = (2.5, 16, 50, 84, 97.5)
C_ENS = "#2c7bb6"
C_DEF = "#111111"
C_THERMO = "#b2182b"
KRES_BAND = (0.1, 10.0)  # m2/s, convection-resolving range quoted in the primer (VERIFY L1)

# JWST comparison sample (experiments/targets.yaml; several entries flagged verify: true)
TARGETS = [("LHS 1140 b", 215), ("K2-18 b", 258), ("TOI-732 c", 305),
           ("TOI-1231 b", 330), ("TOI-270 d", 352), ("LP 791-18 c", 355)]


def load_a(run: Path) -> xr.Dataset:
    fs = sorted((run / "job_a").glob("job_a_*.nc"))
    if not fs:
        raise SystemExit(f"no job_a output in {run / 'job_a'}")
    return xr.concat([xr.open_dataset(f) for f in fs], dim="draw").sortby("draw")


def load_b(run: Path) -> xr.Dataset:
    fs = sorted((run / "job_b").glob("job_b_*.nc"))
    if not fs:
        raise SystemExit(f"no job_b output in {run / 'job_b'}")
    return xr.concat([xr.open_dataset(f) for f in fs], dim="draw").sortby("draw")


def bands(da: xr.DataArray, dim: str = "draw") -> tuple:
    """(p2.5, p16, p50, p84, p97.5) ignoring NaN."""
    v = np.asarray(da.values, float)
    ax = da.dims.index(dim)
    return tuple(np.nanpercentile(v, q, axis=ax) for q in PCT)


# ---------------------------------------------------------------- figure 1
def fig_transition(a: xr.Dataset, out: Path):
    z = a.metallicity.values
    fig, axes = plt.subplots(1, 2, figsize=xp.figsize(J, "double", 0.42))

    ax = axes[0]
    lo2, lo1, med, hi1, hi2 = bands(a.t_mid)
    xp.model_band(ax, z, med, lo1, hi1, lo2, hi2, color=C_ENS,
                  label="ensemble (512 draws)")
    # default-parameter curve: the draw nearest the README defaults
    from meromixis.column import ColumnParams
    d = ColumnParams()
    dist = sum(((a[n].values - getattr(d, n)) / getattr(d, n)) ** 2
               for n in ("kappa_th", "kappa_pexp", "gamma_v", "tint"))
    ax.plot(z, a.t_mid.isel(draw=int(np.nanargmin(dist))).values, color=C_DEF, lw=1.1,
            ls="--", label="nearest to defaults")
    nudge = {"TOI-270 d": -7.0}  # sits 3 K below LP 791-18 c; separate the labels
    for name, teq in TARGETS:
        ax.axhline(teq, color="0.55", lw=0.5, zorder=0)
        ax.text(z[-1] * 0.92, teq + 1.5 + nudge.get(name, 0.0), name,
                fontsize=5.0, va="bottom", ha="right", color="0.35")
    ax.set_xscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(xp.label(r"Transition $T_\star$", "K", J))
    ax.set_xlim(z[0], z[-1])
    ax.legend(loc="upper left", frameon=False, fontsize=6)
    xp.panel_label(ax, "a", style="paren")

    ax = axes[1]
    lo2, lo1, med, hi1, hi2 = bands(a.width)
    xp.model_band(ax, z, med, lo1, hi1, lo2, hi2, color=C_ENS, label=None)
    ax.plot(z, a.width.isel(draw=int(np.nanargmin(dist))).values, color=C_DEF, lw=1.1, ls="--")
    ax.axhspan(20, 25, color=C_THERMO, alpha=0.12, zorder=0)
    ax.text(z[0] * 1.15, 22.5, "README quote", fontsize=5.5, color=C_THERMO, va="center")
    ax.set_xscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(xp.label(r"Transition width", "K", J))
    ax.set_xlim(z[0], z[-1])
    xp.panel_label(ax, "b", style="paren")

    fig.subplots_adjust(wspace=0.30)
    return xp.finalize(fig, str(out / "fig1_transition_ensemble"))


# ---------------------------------------------------------------- figure 2
def fig_sensitivity(a: xr.Dataset, out: Path):
    """Rank correlation of T* with each nuisance parameter, per metallicity."""
    from scipy.stats import spearmanr

    z = a.metallicity.values
    names = [("kappa_th", r"$\kappa_{\rm th}$"), ("gamma_v", r"$\gamma_{\rm v}$"),
             ("kappa_pexp", r"$\kappa$ exponent"), ("tint", r"$T_{\rm int}$")]
    cols = ["#2c7bb6", "#d7191c", "#fdae61", "#1a9641"]
    fig, ax = plt.subplots(figsize=xp.figsize(J, "single", 0.72))
    for (n, tex), c in zip(names, cols):
        pv = a[n].values
        rho = []
        for i in range(len(z)):
            t = a.t_mid.values[:, i]
            m = np.isfinite(t) & np.isfinite(pv)
            rho.append(spearmanr(pv[m], t[m]).statistic if m.sum() > 10 else np.nan)
        ax.plot(z, rho, color=c, lw=1.3, marker="o", ms=2.6, label=tex)
    ax.set_xscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(r"$\rho_{\rm Spearman}$ with $T_\star$")
    xp.zero_line(ax, 0)
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlim(z[0], z[-1])
    ax.legend(frameon=False, fontsize=6.5, ncol=4, loc="lower center")
    return xp.finalize(fig, str(out / "fig2_nuisance_sensitivity"))


# ---------------------------------------------------------------- figure 3
def fig_killtest(b: xr.Dataset, out: Path):
    z = b.metallicity.values
    phi = b.phi_uv.values
    fig, axes = plt.subplots(1, 2, figsize=xp.figsize(J, "double", 0.44))

    ax = axes[0]
    eta_i = int(np.argmin(np.abs(b.eta.values - 1.0)))
    shades = ["#bdd7e7", "#6baed6", "#2171b5"]
    pick = [0, 2, 4] if len(phi) >= 5 else list(range(len(phi)))
    for k, pi in enumerate(pick):
        kc = b.k_crit.isel(phi_uv=pi, eta=eta_i)
        v = np.asarray(kc.values, float)
        v = np.where(np.isfinite(v), v, np.nan)  # inf = hidden without a barrier
        lo1, med, hi1 = (np.nanpercentile(v, q, axis=0) for q in (16, 50, 84))
        dex = np.log10(phi[pi] / phi[0])
        ax.plot(z, med, color=shades[k], lw=1.4,
                label=rf"$\Phi_{{\rm UV}} = 10^{{{dex:.0f}}}\times$ floor")
        ax.fill_between(z, lo1, hi1, color=shades[k], alpha=0.25, lw=0)
    ax.axhspan(*KRES_BAND, color=C_THERMO, alpha=0.13, zorder=0)
    ax.text(z[0] * 1.1, 1.0, "convection-resolving $K$", fontsize=5.5, color=C_THERMO)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(xp.label(r"$K_{\rm crit}$", r"m$^2$ s$^{-1}$", J))
    ax.set_xlim(z[0], z[-1])
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    xp.panel_label(ax, "a", style="paren")

    ax = axes[1]
    frac = b.hidden_by_barrier.isel(eta=eta_i).mean("draw", skipna=True).values  # (Z, phi)
    dex = np.log10(phi / phi[0])
    pc = ax.pcolormesh(z, dex, frac.T, cmap="magma", vmin=0, vmax=1, shading="nearest")
    ax.set_xscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(r"$\log_{10}(\Phi_{\rm UV}/\Phi_{\rm bb})$")
    ax.yaxis.set_label_position("left")
    ax.text(0.04, 0.96, "(b)", transform=ax.transAxes, color="white",
            fontsize=8, fontweight="bold", va="top", ha="left")
    fig.subplots_adjust(wspace=0.38, right=0.88)
    cax = fig.add_axes([0.90, ax.get_position().y0, 0.015, ax.get_position().height])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label("fraction of ensemble hidden by the lid", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    return xp.finalize(fig, str(out / "fig3_killtest"))


# ---------------------------------------------------------------- figure 4
def fig_column(out: Path):
    from meromixis.barrier import Mixing, kzz_profile
    from meromixis.column import INHIBITED, ColumnParams, solve_column

    p0 = ColumnParams(teq=255.0, metallicity=100.0)
    ci = solve_column(p0, inhibition=True)
    cr = solve_column(p0, inhibition=False)
    fig, axes = plt.subplots(1, 2, figsize=xp.figsize(J, "double", 0.5))

    ax = axes[0]
    ax.plot(cr.T, cr.p / 1e5, color="0.55", lw=1.1, ls="--", label="no inhibition")
    ax.plot(ci.T, ci.p / 1e5, color=C_DEF, lw=1.4, label="with inhibition")
    m = (ci.regime == INHIBITED).astype(int)
    st = np.nonzero(np.diff(np.r_[0, m]) == 1)[0]
    en = np.nonzero(np.diff(np.r_[m, 0]) == -1)[0]
    for s, e in zip(st, en):
        ax.axhspan(ci.p[s] / 1e5, ci.p[e] / 1e5, color=C_THERMO, alpha=0.18, lw=0, zorder=0)
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlabel(xp.label("Temperature", "K", J))
    ax.set_ylabel(xp.label("Pressure", "bar", J))
    ax.set_ylim(1e3, 1e-4)
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    xp.panel_label(ax, "a", style="paren")

    ax = axes[1]
    k = kzz_profile(ci, Mixing())
    ax.plot(k, ci.p / 1e5, color=C_ENS, lw=1.4, label=r"$K(p)$, inhibited column")
    ax.plot(kzz_profile(cr, Mixing()), cr.p / 1e5, color="0.55", lw=1.1, ls="--",
            label=r"$K(p)$, no inhibition")
    for s, e in zip(st, en):
        ax.axhspan(ci.p[s] / 1e5, ci.p[e] / 1e5, color=C_THERMO, alpha=0.18, lw=0, zorder=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlabel(xp.label(r"Eddy diffusivity $K$", r"m$^2$ s$^{-1}$", J))
    ax.set_ylabel(xp.label("Pressure", "bar", J))
    ax.set_ylim(1e3, 1e-4)
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    xp.panel_label(ax, "b", style="paren")

    fig.subplots_adjust(wspace=0.32)
    return xp.finalize(fig, str(out / "fig4_column_anatomy"))


# ---------------------------------------------------------------- figure 5
def fig_checks(run: Path, out: Path):
    j = json.loads((run / "job_c" / "job_c_checks.json").read_text())
    c1, c2 = j["C1_saturation"], j["C2_pbot"]
    fig, axes = plt.subplots(1, 2, figsize=xp.figsize(J, "double", 0.44))

    ax = axes[0]
    z = np.array(c1["metallicity"], float)
    d = np.array(c1["dT_star_cc_minus_wp"], float)
    ax.plot(z, d, color=C_ENS, lw=1.4, marker="o", ms=3)
    xp.zero_line(ax, 0)
    ax.set_xscale("log")
    ax.set_xlabel(xp.label("Metallicity", r"$\times$ solar", J))
    ax.set_ylabel(xp.label(r"$T_\star$(PCM C--C) $-$ $T_\star$(Wagner-Pruss)", "K", J))
    ax.set_xlim(z[0], z[-1])
    xp.panel_label(ax, "a", style="paren")

    ax = axes[1]
    for key, c in (("Z30", "#2c7bb6"), ("Z100", "#d7191c")):
        rows = c2["results"][key]
        pb = np.array([r["p_bot_bar"] for r in rows], float)
        kc = np.array([r["k_crit"] for r in rows], float)
        ind = np.array([r["in_domain"] for r in rows], bool)
        lbl = key.replace("Z", r"$Z=") + r"\times$ solar"
        ax.plot(pb, kc, color=c, lw=1.4, marker="o", ms=3.4, mfc="white", label=lbl)
        if (~ind).any():
            ax.plot(pb[~ind], kc[~ind], color=c, marker="x", ms=5, ls="none")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xp.label(r"Model bottom $p_{\rm bot}$", "bar", J))
    ax.set_ylabel(xp.label(r"$K_{\rm crit}$", r"m$^2$ s$^{-1}$", J))
    ax.legend(frameon=False, fontsize=6)
    xp.panel_label(ax, "b", style="paren")

    fig.subplots_adjust(wspace=0.38)
    return xp.finalize(fig, str(out / "fig5_verify_checks"))


def main() -> int:
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    xp.use_style(fontsize=8, family="serif")
    a, b = load_a(run), load_b(run)
    print(f"job A: {a.sizes['draw']} draws   job B: {b.sizes['draw']} draws")
    for fn, args in ((fig_transition, (a, out)), (fig_sensitivity, (a, out)),
                     (fig_killtest, (b, out)), (fig_column, (out,)),
                     (fig_checks, (run, out))):
        try:
            paths = fn(*args)
            print("wrote", ", ".join(Path(p).name for p in paths))
        except Exception as e:  # keep going so one bad panel does not lose the rest
            print(f"FAILED {fn.__name__}: {type(e).__name__}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
