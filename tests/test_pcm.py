"""Run-directory generation for the Generic PCM."""
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

from meromixis.cli import main
from meromixis.pcm.deffiles import RunSpec, one_d_twin, pseudo_altitudes, render_callphys, render_traceur, write_run_dir
from meromixis.pcm.grid import expand, load, teq_nodes, write_all

CONFIG = Path(__file__).resolve().parents[1] / "experiments" / "steam_rising.yaml"


def test_forcing_arithmetic_matches_the_primer_table():
    for teq, flux in [(300, 1837), (350, 3403), (400, 5806), (500, 14175), (600, 29393)]:
        assert RunSpec(teq=teq).stellar_flux == pytest.approx(flux, rel=2e-3)
    s = RunSpec(teq=300.0)
    assert s.fat1au == pytest.approx(s.stellar_flux * 0.1429**2)
    assert s.planet.gravity == pytest.approx(12.4, abs=0.1)
    assert s.intheat == pytest.approx(5.670374419e-8 * 50.0**4)


def test_vertical_grid_is_log_pressure():
    s = RunSpec()
    z = pseudo_altitudes(s)
    assert len(z) == s.nlev and np.all(np.diff(z) > 0)
    assert np.allclose(np.diff(z), np.diff(z)[0])  # even in log p (Charnay+2021)
    assert 0 < z[0] < z[-1] < np.log10(s.p_bot / s.p_top)


def test_callphys_carries_the_post_r4146_keys():
    text = render_callphys(RunSpec())
    for key in ("cppd_ref", "mugaz_ref", "cppv_ref", "thermo_phy", "moist_convection_inhibition = .true.",
                "generic_condensation", "qvap_deep", "nosurf     = .true.", "tlocked    = .true.", "Fat1AU"):
        assert key in text, key
    assert "\ncpp " not in text and "\nmugaz " not in text  # the pre-r4146 names must not come back
    assert "moist_convection_inhibition = .false." in render_callphys(RunSpec(inhibition=False))
    dry = render_callphys(RunSpec(water_scheme="dry"))
    assert "generic_condensation        = .true." not in dry and "tracer     = .false." in dry
    assert "watercond  = .true." in render_callphys(RunSpec(water_scheme="legacy"))
    assert "global1d" in render_callphys(one_d_twin(RunSpec()))


def test_generic_tracers_avoid_the_legacy_name():
    """generic_tracer_index skips names containing 'h2o'; the GCS route must not use them."""
    names = RunSpec().tracer_names
    assert all("h2o" not in n for n in names) and names[0].endswith("_vap") and names[1].endswith("_ice")
    assert render_traceur(RunSpec()).startswith("#ModernTrac-v1\n2\n")
    assert render_traceur(RunSpec(water_scheme="legacy")).split() == ["2", "h2o_vap", "h2o_ice"]


def test_sweep_expansion_and_branching():
    specs = load(CONFIG)
    ids = [s.run_id for s in specs]
    assert len(ids) == len(set(ids)) and specs[0].run_id == "T300_Z100_fid" and specs[0].branch_from == ""
    by_id = {s.run_id: s for s in specs}
    seen = {specs[0].run_id}
    for s in specs[1:]:  # every parent is planned before its child
        assert s.branch_from in seen, s.run_id
        seen.add(s.run_id)
    assert by_id["T300_Z100_noinh"].branch_from == "T300_Z100_fid" and not by_id["T300_Z100_noinh"].inhibition
    assert by_id["T300_Z100_hires"].nlon == 128
    assert {s.metallicity for s in specs} == {3.0, 30.0, 100.0, 300.0}
    assert isinstance(by_id["T300_Z100_fid"].p_bot, float)
    with pytest.raises(KeyError):
        expand({"defaults": {"not_a_field": 1}})


def test_teq_nodes_cluster_on_the_transition():
    t = teq_nodes(285.0, 250.0, 600.0, 11)
    assert t.min() == 250 and t.max() == 600 and np.all(np.diff(t) > 0)
    assert np.sum(np.abs(t - 285) <= 30) >= np.sum(t > 400)


def test_write_run_directories(tmp_path):
    d = write_run_dir(RunSpec(), tmp_path)
    for name in ("callphys.def", "run.def", "gases.def", "traceur.def", "z2sig.def", "diagfi.def", "job.slurm",
                 "planet_constants.txt", "manifest.json"):
        assert (d / name).exists(), name
    man = json.loads((d / "manifest.json").read_text())
    assert man["derived"]["qvap_deep"] == pytest.approx(0.3445, abs=1e-3)
    job = (d / "job.slurm").read_text()
    assert "--account=def-ncowan" in job and "Everything is cool" in job and "sbatch job.slurm" in job
    assert (d / "z2sig.def").read_text().splitlines()[0] == "0.43429448"
    one = write_run_dir(one_d_twin(RunSpec()), tmp_path)
    assert (one / "rcm1d.def").exists() and not (one / "run.def").exists()
    dirs = write_all(CONFIG, tmp_path / "sweep")
    assert (tmp_path / "sweep" / "runs.tsv").read_text().count("\n") == len(load(CONFIG)) + 1
    assert len(dirs) > len(load(CONFIG))  # 1-D twins included


def test_init_inhibition_gives_twins_one_start(tmp_path):
    """On/off twins with a shared init_inhibition start from the same column."""
    base = replace(one_d_twin(RunSpec(teq=265.0, greybody=True)), nday=10)
    on = write_run_dir(replace(base, run_id="on", inhibition=True, init_inhibition=False), tmp_path)
    off = write_run_dir(replace(base, run_id="off", inhibition=False, init_inhibition=False), tmp_path)
    old = write_run_dir(replace(base, run_id="old", inhibition=True), tmp_path)
    assert (on / "profile").read_text() == (off / "profile").read_text()
    assert (on / "profile").read_text() != (old / "profile").read_text()
    assert "moist_convection_inhibition = .true." in (on / "callphys.def").read_text()
    assert "moist_convection_inhibition = .false." in (off / "callphys.def").read_text()


def test_cli_end_to_end(tmp_path, capsys):
    assert main(["column", "--teq", "255"]) == 0
    assert "cold trap" in capsys.readouterr().out
    assert main(["killtest", "--metallicity", "30", "--eta", "1"]) == 0
    assert "verdict" in capsys.readouterr().out
    syn = tmp_path / "syn.nc"
    assert main(["synth", "--teq", "270", "--out", str(syn)]) == 0
    assert main(["diag", str(syn), "--gravity", "12.43", "--metallicity", "100"]) == 0
    out = capsys.readouterr().out
    assert "Result I" in out and "Result II" in out
    assert main(["monitor", str(syn), "--gravity", "12.43", "--metallicity", "100", "--intheat", "0.35"]) == 0
    assert main(["pcm-init", str(CONFIG), str(tmp_path / "runs"), "--no-1d"]) == 0
    assert yaml.safe_load(CONFIG.read_text())["fiducial"]["teq"] == 300


def test_power_law_table(tmp_path):
    """kappa ~ p in Fortran order (T fastest, then p), constant elsewhere, in the model's units."""
    from meromixis.pcm.greytable import kappa_to_cross_section, power_law_table, write_power_law_dir

    k = power_law_table([-3.0, 0.0, 3.0, 6.0], n_t=3, n_q=2, n_band=4, n_g=5,
                        kappa_th=1e-2, kappa_pexp=1.0, mugaz_ref=2.5)
    assert k.shape == (3, 4, 2, 4, 5)
    assert np.allclose(k[:, 2], kappa_to_cross_section(1e-2, 2.5))  # 1e3 mbar = 1 bar = p_ref
    assert np.allclose(k[:, 3] / k[:, 2], 1e3) and np.ptp(k[:, 1]) == 0
    flat = k.ravel(order="F")
    assert flat[0] == flat[2] and flat[3] == pytest.approx(flat[0] * 1e3)  # T varies fastest

    ax = tmp_path / "null_axes"
    (ax / "20x25").mkdir(parents=True)
    (ax / "T.dat").write_text("3\n100.\n200.\n300.\n")
    (ax / "p.dat").write_text("4\n-3.0\n0.0\n3.0\n6.0\n")
    (ax / "Q.dat").write_text("2\nH2_\nH2O\n2\n1e-4\n1e-2\n")
    (ax / "g.dat").write_text("5\n.1\n.2\n.3\n.2\n.2\n")
    (ax / "20x25" / "narrowbands_IR.in").write_text("4\n10. 500.\n500. 1000.\n1000. 2000.\n2000. 2800.\n")
    # centres 2160, 2525, 2960, 3285, ...: the model treats the first three as infrared
    (ax / "20x25" / "narrowbands_VI.in").write_text(
        "6\n2020. 2300.\n2300. 2750.\n2750. 3170.\n3170. 3400.\n3400. 3700.\n3700. 4150.\n")
    with pytest.raises(ValueError):
        write_power_law_dir(ax, tmp_path / "null_bad", "20x25", 1e-2, 4e-3, 1.0, 2.5)
    d = write_power_law_dir(ax, tmp_path / "greyp1", "20x25", 1e-2, 4e-3, 1.0, 2.5, t_axis=(100., 500., 1500.))
    assert np.loadtxt(d / "20x25" / "corrk_gcm_IR.dat").size == 3 * 4 * 2 * 4 * 5
    vi = np.loadtxt(d / "20x25" / "corrk_gcm_VI.dat").reshape((3, 4, 2, 6, 5), order="F")
    # long-wave stellar bands carry kappa_IR's law, the rest kappa_VI's: ratio 0.4 at every p
    assert np.allclose(vi[:, :, :, 3:] / vi[:, :, :, :1], 0.4)
    assert (d / "p.dat").is_symlink() and (d / "T.dat").read_text().split()[-1] == "1500.0"


def test_kappa_pexp_changes_start_and_floor(tmp_path):
    base = replace(one_d_twin(RunSpec(teq=265.0, greybody=True)), nday=10)
    a = write_run_dir(replace(base, run_id="k0"), tmp_path)
    b = write_run_dir(replace(base, run_id="k1", kappa_pexp=1.0, corrkdir="greyp1_H2He_20x25"), tmp_path)
    assert (a / "profile").read_text() != (b / "profile").read_text()
    assert "kappa_IR   = 0.01" in (a / "callphys.def").read_text()
    assert "kappa_IR   = 1.0e-12" in (b / "callphys.def").read_text()
    assert "kappa_VI   = 1.0e-12" in (b / "callphys.def").read_text()
    assert "dtdryadj" in (b / "diagfi.def").read_text() and "dtmoistadj" in (b / "diagfi.def").read_text()


def test_reals_keep_a_decimal_point():
    """IOIPSL reads 1e-12 as 1e-13; every exponent must carry a decimal point."""
    from meromixis.pcm.deffiles import _real, check_reals

    assert _real(1e-12) == "1.0e-12" and _real(0.004) == "0.004" and _real(2.5e-7) == "2.5e-07"
    check_reals("ok.def", "kappa_IR   = 1.0e-12\n# 1e-12 in a comment is fine\n")
    with pytest.raises(ValueError):
        check_reals("bad.def", "kappa_IR   = 1e-12\n")


def test_rcm1d_ndt_is_days(tmp_path):
    """rcm1d multiplies ndt by day_step itself, so ndt must be nday, not steps."""
    d = write_run_dir(replace(one_d_twin(RunSpec()), nday=19200), tmp_path)
    assert "ndt        = 19200" in (d / "rcm1d.def").read_text()


def test_pcm_saturation_matches_the_condensable_table():
    """Start columns use the PCM's own water curve; keep it in step with the table row."""
    from meromixis.pcm.deffiles import PCM_H2O, pcm_saturation

    row = (Path(__file__).resolve().parents[1] / "experiments/table_tracers_condensable.H2O").read_text()
    line = next(ln for ln in row.splitlines() if ln.startswith("H2O "))
    kv = dict(tok.split("=") for tok in line.split()[1:])
    for k, v in PCM_H2O.items():
        assert float(kv[k]) == v, k
    assert pcm_saturation().psat_s(273.16) == pytest.approx(611.657)
