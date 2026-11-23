"""Gate G1: does the PCM's own column show the three-layer structure, and does the
deep temperature rise when inhibition is on?"""
import os, sys
import numpy as np, xarray as xr
sys.path.insert(0, os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")) + "/src")
from meromixis.composition import Composition, q_to_x
from meromixis.inhibition import q_crit, find_stable_layer
from meromixis.thermo import Saturation

# Run directories; override to analyse a local copy of the outputs.
R = os.environ.get("MEROMIXIS_RUNS", f"/scratch/{os.environ.get('USER', 'wasi')}/pcm/runs")
sat, comp = Saturation("cc"), Composition(100.0)
out = {}

for tag in ("g1inh", "g1noinh"):
    f = f"{R}/T265_Z100_{tag}_1d/diagfi.nc"
    ds = xr.open_dataset(f, decode_times=False)
    nt = ds.sizes["Time"]
    g = lambda n: ds[n].isel(Time=-1).squeeze().values if n in ds else None
    T, p, qv = g("temperature"), g("p"), g("H2O_vap")
    rh, vt = g("RH_generic"), g("vteta")
    o = np.argsort(-p)                                   # bottom first
    T_, p_, q_ = T[o], p[o], qv[o]
    rh_ = rh[o] if rh is not None else None
    out[tag] = dict(T=T_, p=p_, q=q_, rh=rh_, vt=(vt[o] if vt is not None else None))

    print(f"\n{'='*64}\n== {tag} ==  ({nt} time samples)")
    ip = lambda P: float(np.interp(np.log(P), np.log(p_[::-1]), T_[::-1]))
    print(f"   T(0.1 bar)={ip(1e4):7.1f}   T(1 bar)={ip(1e5):7.1f}   "
          f"T(10 bar)={ip(1e6):7.1f}   T(80 bar)={T_[0]:7.1f} K at {p_[0]/1e5:.3g} bar")
    x_ = q_to_x(q_, comp.mu_dry)
    qc = q_crit(T_, comp, sat)
    sup = q_ / np.where(np.isfinite(qc) & (qc > 0), qc, np.inf)
    m = sup > 1.0
    print(f"   q_vap deep = {q_[0]:.4g} kg/kg   (emulator deep value {comp.q_h2o_deep:.4g})")
    print(f"   levels with q > q_cri : {int(m.sum())} of {len(q_)}", end="")
    if m.any():
        print(f"   spanning {p_[m].max()/1e5:.4g} - {p_[m].min()/1e5:.4g} bar"
              f"   peak q/q_cri = {np.nanmax(sup):.2f}")
    else:
        print()
    lay = find_stable_layer(p_, T_, x_, comp, 12.43, sat, float(x_[0]))
    print(f"   find_stable_layer: exists={lay.exists}", end="")
    print(f"  {lay.p_base/1e5:.4g}->{lay.p_top/1e5:.4g} bar, dT={lay.delta_t:.1f} K, "
          f"dz={lay.thickness/1e3:.1f} km, q_deep/q_cri={lay.supercriticality:.2f}"
          if lay.exists else "")
    if rh_ is not None and m.any():
        print(f"   RH inside the q>q_cri region: min={np.nanmin(rh_[m]):.3f} "
              f"median={np.nanmedian(rh_[m]):.3f} max={np.nanmax(rh_[m]):.3f}")
    # energy balance and drift
    for v in ("ASR", "OLR", "ISR"):
        if v in ds:
            print(f"   {v} = {float(ds[v].isel(Time=-1).squeeze().values):8.3f} W/m2", end="")
    print()
    Tl = ds["temperature"].values.squeeze()
    k = max(1, nt//4)
    print(f"   max |dT| over last {k} of {nt} samples: {np.abs(Tl[-1]-Tl[-k]).max():.3f} K")

if len(out) == 2:
    a, b = out["g1inh"], out["g1noinh"]
    print(f"\n{'='*64}\n== G1 VERDICT: inhibition on minus off ==")
    for P in (1e4, 1e5, 1e6, 1e7):
        ta = float(np.interp(np.log(P), np.log(a['p'][::-1]), a['T'][::-1]))
        tb = float(np.interp(np.log(P), np.log(b['p'][::-1]), b['T'][::-1]))
        print(f"   dT({P/1e5:7.3g} bar) = {ta-tb:+8.2f} K    (on {ta:7.1f}, off {tb:7.1f})")
    print(f"   dT(model bottom)   = {a['T'][0]-b['T'][0]:+8.2f} K"
          f"    (on {a['T'][0]:.1f}, off {b['T'][0]:.1f})")
