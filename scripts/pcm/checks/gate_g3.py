"""Gate G2/G3: the first 3-D runs.  Lid vs no-lid, with every number labelled by how
much of it the run can actually support (nday=50, two output samples)."""
import os, sys, json
import numpy as np, xarray as xr
sys.path.insert(0, os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")) + "/src")
from meromixis.pcm.io import open_run
from meromixis import diagnostics3d as d3
from meromixis.thermo import Saturation
from meromixis.composition import Composition, q_to_x
from meromixis.inhibition import q_crit

# Run directories; override to analyse a local copy of the outputs.
R = os.environ.get("MEROMIXIS_RUNS", f"/scratch/{os.environ.get('USER', 'wasi')}/pcm/runs")
G, Z = 12.4, 100.0
SAT = Saturation("cc", latent_heat=45038.2 / 0.01801528, t_ref=273.16, p_ref=611.657)
COMP = Composition(Z)
RUNS = {"g2dry": "T265_Z100_g2dry", "g3noinh": "T265_Z100_g3noinh", "g3lid": "T265_Z100_g3lid"}
PLEV = np.geomspace(7.5e6, 5e2, 24)          # 75 bar -> 5 mbar
res = {}

def openit(tag):
    """g2dry carries no water tracer, so standardize() would refuse it.  Inject a zero
    vapour field for that run only: it is the dry control and no water diagnostic is
    applied to it, but the energy-balance comparison needs it in the same form."""
    f = f"{R}/{RUNS[tag]}/diagfi.nc"
    raw = xr.open_dataset(f, decode_times=False)
    if "H2O_vap" not in raw:
        raw["H2O_vap"] = xr.zeros_like(raw["temperature"])
    from meromixis.pcm.io import standardize
    return standardize(raw, gravity=G, metallicity=Z, w_is_omega=False)

print("#" * 74)
print("# G2/G3: first 3-D Generic PCM runs.  nday=50, Time=[25,50] d, 64x48x40.")
print("# Rotation period 32.94 d, so the run covers 1.5 rotations: NOT spun up.")
print("#" * 74)

# ---------------------------------------------------------------- energy balance
print("\n=== 1. Energy balance and spin-up state ===")
print(f"{'run':10s} {'ASR':>8s} {'OLR':>8s} {'imbal':>8s} {'imbal/ASR':>10s} {'T(10bar)':>9s} {'dT/dt':>10s}")
for tag in RUNS:
    ds = openit(tag)
    asr = float(d3.wmean(ds["asr"], ds).mean("time"))
    olr = float(d3.wmean(ds["olr"], ds).mean("time"))
    td = d3.wmean(d3.to_pressure(ds, "T", [1e6]).isel(plev=0), ds)
    dt = float(td.isel(time=-1) - td.isel(time=0)) / float(ds.time[-1] - ds.time[0])
    print(f"{tag:10s} {asr:8.2f} {olr:8.2f} {asr-olr:8.2f} {(asr-olr)/asr:10.3f} "
          f"{float(td.isel(time=-1)):9.2f} {dt:+9.3f} K/d")
    res.setdefault(tag, {}).update(asr=asr, olr=olr, imbalance=asr-olr, t_10bar=float(td.isel(time=-1)),
                                   t_deep_drift_K_per_day=dt)
    ds.close()

# ---------------------------------------------------------------- the twin test
print("\n=== 2. Lid minus no-lid: does inhibition change the atmosphere at all? ===")
a, b = openit("g3lid"), openit("g3noinh")
Ta = d3.wmean(d3.to_pressure(a, "T", PLEV), a).mean("time").values
Tb = d3.wmean(d3.to_pressure(b, "T", PLEV), b).mean("time").values
qa = d3.wmean(d3.to_pressure(a, "q_vap", PLEV), a).mean("time").values
qb = d3.wmean(d3.to_pressure(b, "q_vap", PLEV), b).mean("time").values
print(f"{'p (bar)':>10s} {'T lid':>8s} {'T noinh':>8s} {'dT':>8s} {'q lid':>10s} {'q noinh':>10s}")
for p, t1, t2, q1, q2 in zip(PLEV, Ta, Tb, qa, qb):
    print(f"{p/1e5:10.4g} {t1:8.2f} {t2:8.2f} {t1-t2:+8.2f} {q1:10.3e} {q2:10.3e}")
dT = Ta - Tb
print(f"\n  max |dT| = {np.nanmax(np.abs(dT)):.2f} K at p = {PLEV[np.nanargmax(np.abs(dT))]/1e5:.4g} bar")
print(f"  dT at the bottom level  = {dT[0]:+.2f} K")
print(f"  identical?  {np.allclose(Ta, Tb, rtol=0, atol=1e-10)}")
res["twin"] = {"max_abs_dT": float(np.nanmax(np.abs(dT))), "dT_deep": float(dT[0]),
               "identical": bool(np.allclose(Ta, Tb, rtol=0, atol=1e-10))}

# ---------------------------------------------------------------- the lid
print("\n=== 3. The inhibited layer (three masks) ===")
for tag, ds in (("g3lid", a), ("g3noinh", b)):
    print(f"\n-- {tag}")
    for crit in ("moist", "ledoux", "union"):
        inh = d3.inhibition_fields(ds, sat=SAT, criterion=crit)
        f = float(inh["lid_area_fraction"].mean())
        if f > 0:
            pt = float(np.exp(d3.wmean(np.log(inh["p_lid_top"]), ds).mean("time")))
            pb = float(np.exp(d3.wmean(np.log(inh["p_lid_base"]), ds).mean("time")))
            print(f"   {crit:7s}: area fraction {f:6.3f}   top {pt/1e5:8.3g} bar   base {pb/1e5:8.3g} bar")
        else:
            print(f"   {crit:7s}: area fraction {f:6.3f}   (no lid anywhere)")
        res.setdefault(tag, {})[f"lid_area_{crit}"] = f
    print(f"   moist/union agreement: {float(inh['lid_agreement'].mean()):.3f}")
    sc = d3.inhibition_fields(ds, sat=SAT, criterion='union')['supercriticality']
    print(f"   supercriticality q_deep/q_cri: median {float(sc.median()):.3f}  max {float(sc.max()):.3f}")

# ---------------------------------------------------------------- sub-saturation
print("\n=== 4. Is the layer sub-saturated? (RH_generic inside the union mask) ===")
for tag, path in (("g3lid", RUNS["g3lid"]), ("g3noinh", RUNS["g3noinh"])):
    raw = xr.open_dataset(f"{R}/{path}/diagfi.nc", decode_times=False)
    ds = openit(tag)
    rh = raw["RH_generic"].rename({"Time": "time", "altitude": "lev",
                                   "latitude": "lat", "longitude": "lon"})
    if float(ds["p"].isel(lev=0).mean()) > float(ds["p"].isel(lev=-1).mean()):
        pass
    else:
        rh = rh.isel(lev=slice(None, None, -1))
    rh = rh.isel(lev=slice(None, None, -1)) if False else rh
    m = d3.inhibition_fields(ds, sat=SAT, criterion="union")["lid_union"]
    rh = rh.transpose(*m.dims)
    v = rh.values[m.values]
    if v.size:
        print(f"   {tag:9s} n={v.size:7d}  RH  median {np.nanmedian(v):.3f}  "
              f"p10 {np.nanpercentile(v,10):.3f}  p90 {np.nanpercentile(v,90):.3f}  "
              f"frac<0.95 {np.mean(v<0.95):.3f}")
        res.setdefault(tag, {})["rh_in_lid_median"] = float(np.nanmedian(v))
    else:
        print(f"   {tag:9s} no lid cells")
    raw.close(); ds.close()

# ---------------------------------------------------------------- cold trap
print("\n=== 5. Cold trap by region (geometric means) ===")
for tag, ds in (("g3lid", a), ("g3noinh", b)):
    ts = d3.transmission_summary(ds, sat=SAT)
    print(f"\n-- {tag}   limb asymmetry {ts.attrs['limb_asymmetry_dex']:+.3f} dex")
    print(f"   {'region':12s} {'log10 x_obs':>12s} {'log10 T_sat':>12s} {'log10 T_act':>12s} "
          f"{'T_trap (K)':>11s} {'p_trap (bar)':>13s}")
    for r in ts.region.values:
        s = ts.sel(region=r)
        print(f"   {str(r):12s} {float(s.log10_x_obs):12.3f} {float(s.log10_transmission_sat):12.3f} "
              f"{float(s.log10_transmission_actual):12.3f} {float(s.t_coldtrap):11.1f} "
              f"{float(s.p_coldtrap)/1e5:13.4g}")
    res.setdefault(tag, {})["coldtrap"] = {str(r): {k: float(ts[k].sel(region=r)) for k in ts.data_vars}
                                           for r in ts.region.values}

# ---------------------------------------------------------------- permeability
print("\n=== 6. Permeability of the layer (resolved flow only, lower bound) ===")
for tag, ds in (("g3lid", a), ("g3noinh", b)):
    for crit in ("moist", "union"):
        try:
            k = d3.layer_permeability(ds, criterion=crit)
        except Exception as e:
            print(f"   {tag:9s} {crit:6s}: FAILED {type(e).__name__}: {e}")
            continue
        if not k["has_lid"]:
            print(f"   {tag:9s} {crit:6s}: no lid")
            continue
        print(f"   {tag:9s} {crit:6s}: K_layer = {k['kzz_layer']:.4g} m2/s   "
              f"{k['p_base']/1e5:.3g} -> {k['p_top']/1e5:.3g} bar   "
              f"area {k['lid_area_fraction']:.3f}  agree {k['lid_agreement']:.3f}")
        res.setdefault(tag, {})[f"kzz_layer_{crit}"] = k["kzz_layer"]

OUT = os.path.join(os.path.dirname(R.rstrip("/")), "g3_results.json")
with open(OUT, "w") as f:
    json.dump(res, f, indent=2, default=str)
print(f"\n[wrote {OUT}]")
