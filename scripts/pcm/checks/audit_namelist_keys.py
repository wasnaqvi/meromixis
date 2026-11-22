"""Definitive key audit against Generic PCM trunk r4609.
import os

Extracts every getin/getin_p key (either quote style) from the exported physics and
lon-lat dynamics, then checks every key meromixis writes.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, os.environ.get("MEROMIXIS_REPO", os.path.expanduser("~/meromixis")) + "/src")
from meromixis.pcm.deffiles import (  # noqa: E402
    RunSpec, render_callphys, render_rcm1d, render_run_def)

ROOT = Path("/scratch/" + os.environ.get("USER", "wasi") + "/lmd_src")
PAT = re.compile(r"""getin(?:_p)?\s*\(\s*['"]([A-Za-z0-9_]+)['"]""")
model, files = set(), 0
for p in ROOT.rglob("*"):
    if p.is_file() and p.suffix in (".F", ".F90", ".f90", ".h"):
        model |= set(PAT.findall(p.read_text(errors="replace")))
        files += 1
print(f"scanned {files} source files; model reads {len(model)} distinct keys")
print(f"sanity: 'nday' in model = {'nday' in model}; 'cpp_mugaz_mode' = {'cpp_mugaz_mode' in model}\n")

def keys_of(text):
    out = []
    for line in text.splitlines():
        line = line.split("#")[0].split("!")[0].strip()
        if not line or "=" not in line or line.startswith("INCLUDEDEF"):
            continue
        out.append(line.split("=")[0].strip())
    return out

variants = [
    ("callphys generic 3-D", render_callphys, RunSpec()),
    ("callphys generic 1-D", render_callphys, RunSpec(one_d=True)),
    ("callphys legacy", render_callphys, RunSpec(water_scheme="legacy")),
    ("callphys dry", render_callphys, RunSpec(water_scheme="dry")),
    ("callphys greybody", render_callphys, RunSpec(greybody=True)),
    ("callphys clouds", render_callphys, RunSpec(clouds_radiative=True)),
    ("run.def", render_run_def, RunSpec()),
    ("rcm1d.def", render_rcm1d, RunSpec(one_d=True)),
]
bad = {}
for label, fn, spec in variants:
    ks = keys_of(fn(spec))
    unknown = [k for k in ks if k not in model]
    print(f"{label:24s} {len(ks):3d} keys   "
          + ("all read by the model" if not unknown else f"{len(unknown)} NOT read"))
    for k in unknown:
        bad.setdefault(k, []).append(label)

print("\n=== keys the model never reads at r4609 ===")
for k, where in sorted(bad.items()) or [(None, None)]:
    if k is None:
        print("  (none)"); break
    near = sorted(m for m in model if k.lower() in m.lower() or m.lower() in k.lower())
    print(f"  {k:20s} ({len(where)} file variant(s))"
          + (f"   closest real keys: {near}" if near else "   no similar key exists"))
