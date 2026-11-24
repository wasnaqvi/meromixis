"""Write the rcm1d twins for gate G1 (inhibition on and off)."""
import os
import sys
sys.path.insert(0, os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")) + "/src")
from pathlib import Path
import yaml
from dataclasses import replace
from meromixis.pcm.deffiles import Planet, RunSpec, Star, one_d_twin, write_run_dir

cfg = yaml.safe_load(Path("/home/wasi/meromixis/experiments/gate_grey.yaml").read_text())
d = {k: (type(getattr(RunSpec(), k))(v) if type(getattr(RunSpec(), k)) in (int, float)
         and not isinstance(v, bool) else v)
     for k, v in cfg["defaults"].items()}
base = RunSpec(planet=Planet(**cfg["planet"]), star=Star(**cfg["star"]),
               teq=265.0, metallicity=100.0, **d)
out = Path("/scratch/" + os.environ.get("USER", "wasi") + "/pcm/runs")
for tag, inhib in (("g1inh", True), ("g1noinh", False)):
    spec = replace(base, run_id=f"T265_Z100_{tag}", inhibition=inhib,
                   virtual_theta=inhib, tier="gate1d")
    twin = replace(one_d_twin(spec), run_id=f"T265_Z100_{tag}_1d", nday=19200)
    p = write_run_dir(twin, out, account=cfg["cluster"]["account"])
    print("wrote", p, " inhibition =", inhib, " ndt =", twin.nday)
