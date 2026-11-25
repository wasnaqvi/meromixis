"""Stage A of the convecting inhibition test: rcm1d twins with grey kappa ~ p.

With constant grey kappa the PCM column never convected (dry or moist), so the
identical-start twins were bit-identical.  kappa ~ p lifts the deep radiative gradient
above the adiabat.  The emulator then puts a dry-convective deep below ~24 bar and a
moist-convective zone at 4-23 bar without inhibition, and a stable layer at 3.9-6.8 bar
with a deep 156 K hotter with it.

Writes the corrk table that carries kappa (pcm.greytable) and four run directories:
switch on/off x start column lidstart/mixstart.  lidstart is the decisive pair: with
inhibition off the layer should be mixed away within days.

NLEV (environment, default 40) sets the level count.  At 40 the layer is ~1.5 levels
thick and the adjustment cannot act on it; other counts need rcm1d built with that
many levels (jobs/rcm1d_build.sh), and their run directories get an _L<NLEV> suffix.
"""
import os
import sys
from dataclasses import replace
from pathlib import Path

import yaml

REPO = Path(os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")))
sys.path.insert(0, str(REPO / "src"))
from meromixis.pcm.deffiles import Planet, RunSpec, Star, one_d_twin, write_run_dir  # noqa: E402
from meromixis.pcm.greytable import write_power_law_dir  # noqa: E402

NDAY, PEXP, TABLE = 20_000, 1.0, "greyp1_H2He_20x25"
NLEV = int(os.environ.get("NLEV", 40))
SUF = "" if NLEV == 40 else f"_L{NLEV}"

cfg = yaml.safe_load((REPO / "experiments/gate_grey.yaml").read_text())
d = {k: (type(getattr(RunSpec(), k))(v) if type(getattr(RunSpec(), k)) in (int, float)
         and not isinstance(v, bool) else v)
     for k, v in cfg["defaults"].items()}
d.update(kappa_pexp=PEXP, corrkdir=TABLE, nlev=NLEV)
base = RunSpec(planet=Planet(**cfg["planet"]), star=Star(**cfg["star"]),
               teq=265.0, metallicity=100.0, **d)

pcm = Path("/scratch/" + os.environ.get("USER", "wasi") + "/pcm")
corrk = pcm / "datagcm/corrk_data"
t = write_power_law_dir(corrk / "null_H2He_20x25", corrk / TABLE, base.bands,
                        base.kappa_ir, base.kappa_vi, PEXP, base.mugaz_ref)
print(f"wrote {t}  kappa_IR = {base.kappa_ir:g}, kappa_VI = {base.kappa_vi:g} m2/kg, x (p/1 bar)^{PEXP:g};"
      f" mugaz_ref = {base.mugaz_ref:.4f}")

for start, init_inh in (("lidstart", True), ("mixstart", False)):
    for sw, inh in (("on", True), ("off", False)):
        spec = replace(base, inhibition=inh, virtual_theta=inh, init_inhibition=init_inh,
                       notes=f"kappa~p^{PEXP:g} twin: switch {sw}, start column {start}")
        one = replace(one_d_twin(replace(spec, run_id=f"T265_Z100_g1p{sw}_{start}")),
                      run_id=f"T265_Z100_g1p{sw}_{start}{SUF}_1d", tier="gate1d", nday=NDAY)
        p = write_run_dir(one, pcm / "runs", account=cfg["cluster"]["account"])
        print(f"wrote {p}  inhibition={inh}  init_inhibition={init_inh}  nday={NDAY}  nlev={NLEV}")
