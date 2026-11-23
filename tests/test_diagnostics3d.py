"""3-D diagnostics on synthetic PCM-format output."""
import numpy as np
import pytest
import xarray as xr

from meromixis import ColumnParams, solve_column
from meromixis import diagnostics3d as d3
from meromixis.pcm.io import open_run, standardize
from meromixis.synthetic import make_mixing_test, make_synthetic_run

G, Z = 12.43, 100.0


@pytest.fixture(scope="module")
def run():
    raw = make_synthetic_run(ColumnParams(teq=270.0, metallicity=Z, n_lev=200, n_sub=2), contrast=0.4)
    return raw, standardize(raw, G, Z)


def test_standardize_layout(run, tmp_path):
    raw, ds = run
    assert ds["T"].dims == ("time", "lev", "lat", "lon")
    assert float(ds["p"].isel(lev=0).mean()) > float(ds["p"].isel(lev=-1).mean())  # bottom first
    assert float(ds["p"].max()) == pytest.approx(8.0e6, rel=1e-6)  # rebuilt from aps + bps * ps
    path = tmp_path / "diagfi1.nc"
    raw.to_netcdf(path)
    again = open_run(str(path), G, Z, last_n=2)
    assert again.sizes["time"] == 2 and np.allclose(again["T"], ds["T"].isel(time=slice(-2, None)))
    with pytest.raises(KeyError):
        standardize(raw.drop_vars("H2O_vap"), G, Z)


def test_region_masks_partition_the_sphere(run):
    _, ds = run
    m = d3.region_masks(ds)
    frac = {k: float(d3.wmean(v.astype(float), ds)) for k, v in m.items()}
    assert frac["day"] + frac["night"] == pytest.approx(1.0)
    assert frac["day"] == pytest.approx(0.5, abs=0.04)
    assert frac["limb_east"] == pytest.approx(frac["limb_west"], abs=0.01)
    assert frac["limb_east"] + frac["limb_west"] == pytest.approx(frac["terminator"])
    assert float(d3.area_weights(ds).sum()) == pytest.approx(1.0)


def test_interp_logp_is_exact_for_power_laws():
    p = np.geomspace(1e6, 10.0, 30)[None, :, None, None] * np.ones((2, 1, 3, 4))
    f = np.log(p) * 2.0 + 1.0
    out = d3.interp_logp(p, f, [3e5, 1e3, 1e8])
    assert np.allclose(out[:, 0], np.log(3e5) * 2 + 1) and np.allclose(out[:, 1], np.log(1e3) * 2 + 1)
    assert np.isnan(out[:, 2]).all()  # outside the column


def test_cold_trap_fields_match_the_column_model(run):
    _, ds = run
    ct = d3.coldtrap_fields(ds)
    j, i = int(np.argmin(np.abs(ds.lat.values))), int(np.argmin(np.abs(ds.lon.values - 180.0)))
    night = solve_column(ColumnParams(teq=270.0, metallicity=Z, flux_factor=0.6, n_lev=200, n_sub=2))
    got = float(ct["transmission_sat"].isel(time=0, lat=j, lon=i))
    assert got == pytest.approx(night.cold_trap.transmission, rel=0.35)  # 40 GCM levels vs 200
    # no transport in the synthetic run, so the two transmissions agree
    assert np.allclose(np.log10(ct["transmission_actual"]), np.log10(ct["transmission_sat"]), atol=0.15)


def test_transmission_summary_orders_regions(run):
    _, ds = run
    s = d3.transmission_summary(ds)
    x = s["log10_x_obs"]
    assert float(x.sel(region="substellar")) > float(x.sel(region="terminator")) >= float(x.sel(region="antistellar")) - 1e-6
    assert float(x.sel(region="day")) > float(x.sel(region="global")) > float(x.sel(region="night"))
    assert abs(s.attrs["limb_asymmetry_dex"]) < 0.02  # symmetric by construction


def test_uniform_planet_has_no_contrast():
    raw = make_synthetic_run(ColumnParams(teq=260.0, n_lev=200, n_sub=2), contrast=0.0, nlon=16, nlat=12)
    s = d3.transmission_summary(standardize(raw, G, Z))
    assert float(s["log10_x_obs"].max() - s["log10_x_obs"].min()) < 0.02


def test_inhibition_fields_find_the_lid(run):
    _, ds = run
    inh = d3.inhibition_fields(ds)
    assert 0.3 < float(inh["lid_area_fraction"].mean()) <= 1.0
    assert bool((inh["p_lid_base"] >= inh["p_lid_top"]).where(inh["has_lid"], True).all())
    n2 = d3.ledoux_n2(ds)
    assert float(n2.where(inh["lid"]).mean()) > 0  # the lid is Ledoux-stable
    # Above cloud base the lid must sit where the vapour still falls with height: in
    # the well-mixed region above the cold trap nothing is compositionally stabilised.
    x = d3.water_mole_fraction(ds)
    above_base = x < x.isel(lev=0) * (1.0 - 0.02)
    flat = d3.water_gradient(ds) <= inh.attrs["grad_tol"]
    for name in ("lid_moist", "lid_ledoux", "lid_union"):
        assert not bool((inh[name] & above_base & flat).any()), name
    sub = standardize(make_synthetic_run(ColumnParams(teq=230.0, metallicity=3.0, n_lev=200, n_sub=2), nlon=16, nlat=12), G, 3.0)
    assert float(d3.inhibition_fields(sub)["lid_area_fraction"].mean()) == 0.0


def test_flux_gradient_recovers_a_known_diffusivity():
    k_true = 5.0
    ds = make_mixing_test(k_true=k_true)
    out = d3.kzz_flux_gradient(ds, "chi", p_levels=np.geomspace(5e5, 3e2, 10))
    assert np.allclose(out["kzz"].values, k_true, rtol=0.03)
    assert np.all(out["flux"].values > 0)  # up-gradient tracer decays upward: flux is upward
    with pytest.raises(KeyError):
        d3.kzz_flux_gradient(ds.drop_vars("w"), "chi")


def test_convergence_monitor(run):
    _, ds = run
    c = d3.convergence(ds, p_deep=1e6, intheat=0.35)
    assert abs(c["toa_imbalance"]) < 5.0 and abs(c["t_deep_drift"]) < 1e-3
    assert isinstance(d3.global_profile(ds, "T", [1e5, 1e3]), xr.DataArray)
