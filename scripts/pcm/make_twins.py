"""Write identical-start twins for gates G1 (rcm1d) and G3 (3-D).

The first G1/G3 twins each started from the emulator column for their own switch, so
the on/off difference they reported was the 8 K they were given at the start.  Here
every on/off pair shares one start column, and the pair is run from both columns:

  lidstart  the column with inhibition (stable layer, warm deep).  Without inhibition
            the model should convect it away within days: the fast test.
  mixstart  the column without inhibition (convective).  With inhibition the deep can
            only warm at the internal-heat rate, ~500 yr per kelvin: the slow test.

3-D start files are not made here.  newstart's input does not depend on the switch,
so the 3-D runs reuse the start.nc/startfi.nc that make_starts.sh wrote for the old
g3lid (lid column) and g3noinh (mixed column); g3e_job.sh copies them and checks that
the profile.in they were made from matches.
"""
import os
import sys
from dataclasses import replace
from pathlib import Path

import yaml

REPO = Path(os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")))
sys.path.insert(0, str(REPO / "src"))
from meromixis.pcm.deffiles import Planet, RunSpec, Star, one_d_twin, write_run_dir  # noqa: E402

NDAY_1D, NDAY_3D = 100_000, 500

cfg = yaml.safe_load((REPO / "experiments/gate_grey.yaml").read_text())
d = {k: (type(getattr(RunSpec(), k))(v) if type(getattr(RunSpec(), k)) in (int, float)
         and not isinstance(v, bool) else v)
     for k, v in cfg["defaults"].items()}
base = RunSpec(planet=Planet(**cfg["planet"]), star=Star(**cfg["star"]),
               teq=265.0, metallicity=100.0, **d)
out = Path("/scratch/" + os.environ.get("USER", "wasi") + "/pcm/runs")

for start, init_inh in (("lidstart", True), ("mixstart", False)):
    for sw, inh in (("on", True), ("off", False)):
        spec = replace(base, inhibition=inh, virtual_theta=inh, init_inhibition=init_inh,
                       notes=f"identical-start twin: switch {sw}, start column {start}")
        one = replace(one_d_twin(replace(spec, run_id=f"T265_Z100_g1{sw}_{start}")),
                      run_id=f"T265_Z100_g1{sw}_{start}_1d", tier="gate1d", nday=NDAY_1D)
        three = replace(spec, run_id=f"T265_Z100_g3{sw}_{start}", tier="gate", nday=NDAY_3D)
        for s in (one, three):
            p = write_run_dir(s, out, account=cfg["cluster"]["account"])
            print(f"wrote {p}  inhibition={inh}  init_inhibition={init_inh}  nday={s.nday}")
