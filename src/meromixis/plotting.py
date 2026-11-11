"""Figures.  Every plotting function takes data and returns a matplotlib Figure."""
from __future__ import annotations

import numpy as np

from .column import DRY, INHIBITED, MOIST, RADIATIVE, Column
from .inhibition import q_crit

REGIME_COLORS = {RADIATIVE: "#8a8f98", MOIST: "#2a7fbf", DRY: "#c8642a", INHIBITED: "#7b2d8e"}
REGIME_TEXT = {RADIATIVE: "radiative", MOIST: "moist convection", DRY: "dry convection", INHIBITED: "inhibited (lid)"}


def _style():
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False, "axes.titlesize": 9,
    })


def _segments(ax, x, p, regime, lw=2.0, ls="-", alpha=1.0):
    change = np.nonzero(np.diff(regime))[0] + 1
    for a, b in zip(np.r_[0, change], np.r_[change, len(p)]):
        sl = slice(max(a - 1, 0), b)
        ax.plot(x[sl], p[sl] / 1e5, color=REGIME_COLORS[int(regime[a])], lw=lw, ls=ls, alpha=alpha)


def plot_column(col: Column, ref: Column | None = None):
    """Thermal structure and water profile of a column, coloured by regime."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    _style()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    if ref is not None:
        _segments(a, ref.T, ref.p, ref.regime, lw=1.2, ls="--", alpha=0.8)
        b.plot(ref.q, ref.p / 1e5, color="0.55", lw=1.0, ls="--")
    _segments(a, col.T, col.p, col.regime)
    b.plot(col.q, col.p / 1e5, color="k", lw=1.6, label="water vapour $q$")
    b.plot(np.minimum(q_crit(col.T, col.comp, col.params.saturation), 1.0), col.p / 1e5, color=REGIME_COLORS[INHIBITED],
           lw=1.0, ls=":", label=r"$q_{\rm cri}(T)$")
    ct, layer = col.cold_trap, col.stable_layer
    if ct.exists:
        for ax in (a, b):
            ax.axhline(ct.p / 1e5, color="#2a7fbf", lw=0.6, ls="-.")
        a.text(a.get_xlim()[1], ct.p / 1e5, "cold trap ", ha="right", va="bottom", fontsize=7, color="#2a7fbf")
    if layer.exists:
        for ax in (a, b):
            ax.axhspan(layer.p_top / 1e5, layer.p_base / 1e5, color=REGIME_COLORS[INHIBITED], alpha=0.10, lw=0)
    a.set(yscale="log", ylim=(col.p.max() / 1e5, 1e-4), xlabel="temperature (K)", ylabel="pressure (bar)")
    a.set_xlim(150, min(col.T.max(), 1600))
    b.set(xscale="log", xlim=(1e-5, 1.0), xlabel="water mass fraction")
    handles = [Line2D([], [], color=c, lw=2, label=REGIME_TEXT[k]) for k, c in REGIME_COLORS.items()]
    if ref is not None:
        handles.append(Line2D([], [], color="0.4", lw=1.2, ls="--", label="no inhibition"))
    a.legend(handles=handles, fontsize=7, loc="upper right")
    b.legend(fontsize=7, loc="upper left")
    p = col.params
    fig.suptitle(f"$T_{{\\rm eq}}$ = {p.teq:g} K, {p.metallicity:g}$\\times$ solar, $T_{{\\rm int}}$ = {p.tint:g} K"
                 "  (grey emulator)", fontsize=9)
    return fig


def plot_regime_map(ds, targets: list[dict] | None = None):
    """The two boundaries on one plane: water above the cold trap, and where the lid exists."""
    import matplotlib.pyplot as plt

    _style()
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    t, z = ds["teq"].values, ds["metallicity"].values
    logx = np.log10(ds["x_above"].values)
    pc = ax.pcolormesh(t, z, logx, cmap="YlGnBu", shading="nearest", vmin=-6, vmax=-0.5)
    fig.colorbar(pc, ax=ax, label=r"$\log_{10}$ H$_2$O mole fraction above the cold trap", pad=0.02)
    ax.contour(t, z, ds["transmission"].values, levels=[0.1, 0.5, 0.9], colors="k", linewidths=[0.6, 1.6, 0.6],
               linestyles=["--", "-", "--"])
    lid = ds["layer"].values
    if 0 < np.nanmean(lid) < 1:
        ax.contourf(t, z, lid, levels=[0.5, 1.5], colors="none", hatches=["///"])
        ax.contour(t, z, lid, levels=[0.5], colors="#7b2d8e", linewidths=1.4)
    ax.set(yscale="log", xlabel=r"effective $T_{\rm eq}$ (K), albedo included", ylabel=r"metallicity ($\times$ solar)")
    box = dict(facecolor="w", edgecolor="none", alpha=0.85, pad=1.5)
    ax.text(0.03, 0.30, "cold-trapped,\nconvecting", transform=ax.transAxes, va="top", fontsize=8, bbox=box)
    ax.text(0.97, 0.95, "steam", transform=ax.transAxes, va="top", ha="right", fontsize=8, color="w")
    ax.text(0.03, 0.95, "cold-trapped + lidded\n($q_{\\rm deep} > q_{\\rm cri}$)", transform=ax.transAxes, va="top",
            fontsize=8, color="#7b2d8e", bbox=box)
    for tg in targets or []:
        ax.axvline(tg["teq_a03"], ymin=0.0, ymax=0.04, color="k", lw=1.2)
        ax.text(tg["teq_a03"], z.min() * 1.25, tg["name"], rotation=90, fontsize=6, ha="right", va="bottom")
    ax.set_title("Emulator scoping map: solid = cold-trap transition (transmission 0.5), hatched = stable layer",
                 fontsize=8)
    return fig


def plot_killtest(metallicities, k_crit_by_eta: dict, fronts: dict, k_band=(0.1, 10.0)):
    """Left: critical layer diffusivity vs metallicity.  Right: photolysis front vs layer diffusivity."""
    import matplotlib.pyplot as plt

    _style()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.5))
    a.axhspan(*k_band, color="0.85", lw=0)
    a.text(0.98, np.sqrt(k_band[0] * k_band[1]), "convection-resolving\nrange (Leconte+2024)",
           transform=a.get_yaxis_transform(), fontsize=7, va="center", ha="right")
    for (eta, kc), ls in zip(k_crit_by_eta.items(), ["-", "--", ":"]):
        a.plot(metallicities, kc, ls, color="#7b2d8e", lw=1.8, label=rf"$\eta$ = {eta:g}")
    a.set(xscale="log", yscale="log", xlabel=r"metallicity ($\times$ solar)",
          ylabel=r"$K_{\rm crit}$ of the stable layer (m$^2$ s$^{-1}$)")
    from matplotlib.ticker import NullFormatter, ScalarFormatter

    a.set_xticks([15, 30, 60, 100, 300])
    a.xaxis.set_major_formatter(ScalarFormatter())
    a.xaxis.set_minor_formatter(NullFormatter())
    a.legend(title="irreversible fraction", fontsize=7, title_fontsize=7, loc="lower left")
    a.set_title("NH$_3$ hidden where the layer's $K$ lies below the curve")
    colors = ["#2a7fbf", "#7b2d8e", "#c8642a", "#3a8f4d"]
    for (z, (k, pf)), c in zip(fronts.items(), colors):
        b.plot(k, np.asarray(pf) / 100.0, color=c, lw=1.8, label=rf"{z:g}$\times$ solar")
    b.axhspan(1.0, 100.0, color="0.92", lw=0)
    b.text(0.02, 10.0, "probed in\ntransmission", transform=b.get_yaxis_transform(), fontsize=7, va="center")
    b.set(xscale="log", yscale="log", xlabel=r"$K$ in the stable layer (m$^2$ s$^{-1}$)",
          ylabel="NH$_3$ photolysis front (mbar)")
    b.invert_yaxis()
    b.legend(fontsize=7, loc="lower right")
    b.set_title("Front below the band = ammonia missing")
    fig.suptitle("Kill test (grey emulator, hand-estimated UV flux; not a result)", fontsize=9)
    fig.tight_layout()
    return fig
