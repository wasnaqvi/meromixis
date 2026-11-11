"""Expand an experiment YAML into run specifications and a branching (spin-up) plan."""
from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path

import numpy as np
import yaml

from .deffiles import Planet, RunSpec, Star, one_d_twin, write_run_dir

_SPEC_FIELDS = {f.name for f in fields(RunSpec)}
_SPEC_TYPES = {f.name: type(getattr(RunSpec(), f.name)) for f in fields(RunSpec)}


def _coerce(d: dict) -> dict:
    """Cast YAML scalars to the RunSpec field types (PyYAML reads ``8.0e6`` as a string)."""
    out = {}
    for k, v in d.items():
        typ = _SPEC_TYPES.get(k)
        out[k] = typ(v) if typ in (int, float) and not isinstance(v, bool) else v
    return out


def run_id(teq: float, metallicity: float, tag: str = "fid") -> str:
    return f"T{teq:03.0f}_Z{metallicity:03.0f}_{tag}"


def teq_nodes(t_star: float, t_lo: float, t_hi: float, n: int, half_width: float = 30.0) -> np.ndarray:
    """``n`` forcing labels between t_lo and t_hi, clustered around a predicted transition t_star.

    Half of the nodes fall within +-half_width of t_star; the rest thin out towards the
    ends, so that the transition is resolved without spending runs on the plateaus.
    """
    u = np.linspace(-1.0, 1.0, n)
    s = np.sinh(2.0 * u) / np.sinh(2.0)  # dense near 0
    span_lo, span_hi = t_star - t_lo, t_hi - t_star
    t = np.where(s < 0, t_star + s * span_lo, t_star + s * span_hi)
    core = np.abs(u) <= 0.5
    t = np.where(core, t_star + (u / 0.5) * half_width, t)
    return np.unique(np.round(np.clip(t, t_lo, t_hi) / 5.0) * 5.0)


def expand(config: dict) -> list[RunSpec]:
    """All 3-D runs described by a config dict (see ``experiments/steam_rising.yaml``)."""
    planet = Planet(**config.get("planet", {}))
    star = Star(**config.get("star", {}))
    defaults = _coerce(config.get("defaults", {}))
    unknown = set(defaults) - _SPEC_FIELDS
    if unknown:
        raise KeyError(f"unknown keys in defaults: {sorted(unknown)}")
    base = RunSpec(planet=planet, star=star, **defaults)
    runs: dict[str, RunSpec] = {}
    for tier, block in config.get("tiers", {}).items():
        if "runs" in block:  # explicit variants
            for item in block["runs"]:
                item = _coerce(dict(item))
                tag = item.pop("tag")
                spec = replace(base, tier=tier, **item)
                spec = replace(spec, run_id=run_id(spec.teq, spec.metallicity, tag))
                runs[spec.run_id] = spec
        else:  # cartesian product
            extra = _coerce({k: v for k, v in block.items() if k not in ("teq", "metallicity", "tag")})
            for z in block["metallicity"]:
                for t in block["teq"]:
                    spec = replace(base, tier=tier, teq=float(t), metallicity=float(z), **extra)
                    spec = replace(spec, run_id=run_id(t, z, block.get("tag", "fid")))
                    runs.setdefault(spec.run_id, spec)
    return _plan_branches(list(runs.values()), config.get("fiducial", {}))


def _plan_branches(runs: list[RunSpec], fiducial: dict) -> list[RunSpec]:
    """Each run starts from the restart of its nearest already-planned neighbour.

    The fiducial run (cold start from rcm1d) comes first; other fid-tagged runs are added
    in order of distance from it in (teq, log metallicity); variants branch from the
    fid run with the same forcing and composition.
    """
    fid_id = run_id(fiducial.get("teq", 300.0), fiducial.get("metallicity", 100.0), "fid")
    main = [r for r in runs if r.run_id.endswith("_fid")]
    variants = [r for r in runs if not r.run_id.endswith("_fid")]
    if not main:
        return runs
    root = next((r for r in main if r.run_id == fid_id), main[0])

    def dist(a: RunSpec, b: RunSpec) -> float:
        return abs(a.teq - b.teq) / 50.0 + abs(np.log10(a.metallicity / b.metallicity))

    planned = [root]
    todo = sorted((r for r in main if r is not root), key=lambda r: dist(r, root))
    for r in todo:
        parent = min(planned, key=lambda q: dist(r, q))
        planned.append(replace(r, branch_from=parent.run_id))
    by_key = {(r.teq, r.metallicity): r.run_id for r in planned}
    for v in variants:
        planned.append(replace(v, branch_from=by_key.get((v.teq, v.metallicity), root.run_id)))
    return planned


def load(path) -> list[RunSpec]:
    return expand(yaml.safe_load(Path(path).read_text()))


def write_all(path, out_root, with_1d: bool = True) -> list[Path]:
    """Write every run directory (and its rcm1d twin) plus a sweep table ``runs.tsv``."""
    cfg = yaml.safe_load(Path(path).read_text())
    account = cfg.get("cluster", {}).get("account", "def-ncowan")
    specs = expand(cfg)
    out_root = Path(out_root)
    dirs = []
    rows = ["run_id\ttier\tteq_K\tflux_W_m2\tFat1AU\tmetallicity\tqvap_deep\tbranch_from\tbuild"]
    for s in specs:
        dirs.append(write_run_dir(s, out_root, account))
        if with_1d and s.run_id.endswith("_fid"):
            dirs.append(write_run_dir(one_d_twin(s), out_root, account))
        rows.append(
            f"{s.run_id}\t{s.tier}\t{s.teq:g}\t{s.stellar_flux:.1f}\t{s.fat1au:.3f}\t{s.metallicity:g}\t"
            f"{s.qvap_deep:.5f}\t{s.branch_from or '(rcm1d)'}\t{s.build_command()}"
        )
    (out_root / "runs.tsv").write_text("\n".join(rows) + "\n")
    return dirs
