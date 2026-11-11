"""The 1-D emulator: radiative limit, inhibition, cold trap and stable-layer diagnostics."""
import numpy as np
import pytest

from meromixis import ColumnParams, cold_trap, cold_trap_nd, find_stable_layer, solve_column, water_profile
from meromixis.coldtrap import skin_temperature
from meromixis.column import DRY, INHIBITED, MOIST, calibrate_grey, grey_temperature, optical_depth
from meromixis.regimes import LIDDED, STEAM, emulate_point, regime_map, transition_from_samples, transition_teq


@pytest.fixture(scope="module")
def fiducial():
    return solve_column(ColumnParams())


def test_radiative_zone_reproduces_guillot_profile(fiducial):
    """Above the first convective/inhibited level the integrator must return the analytic solution."""
    col = fiducial
    first = np.nonzero(col.regime != 0)[0][0]
    analytic = grey_temperature(optical_depth(col.p[:first], col.params), col.params)
    assert np.allclose(col.T[:first], analytic, rtol=2e-4)
    assert col.T[0] == pytest.approx(skin_temperature(col.params.teq, col.params.gamma_v, col.params.tint), rel=1e-3)


def test_profile_is_monotonic_and_hot_at_depth(fiducial):
    assert np.all(np.diff(fiducial.T) >= -1e-9)
    assert fiducial.regime[-1] == DRY
    assert 700 < fiducial.T[-1] < 3000


def test_inhibition_makes_a_lid_and_a_hotter_interior(fiducial):
    ref = solve_column(fiducial.params, inhibition=False)
    assert (fiducial.regime == INHIBITED).any() and not (ref.regime == INHIBITED).any()
    assert (ref.regime == MOIST).any()
    assert fiducial.T[-1] > ref.T[-1] + 100.0  # Leconte+2017/2024
    layer = fiducial.stable_layer
    assert layer.exists and layer.p_base > layer.p_top and layer.delta_t > 0 and layer.supercriticality > 1
    assert not ref.stable_layer.exists
    # the same answer from the profile-agnostic finder (as used on GCM columns)
    again = find_stable_layer(fiducial.p[::-1], fiducial.T[::-1], fiducial.x[::-1], fiducial.comp, fiducial.params.gravity)
    assert again.p_base == pytest.approx(layer.p_base) and again.thickness == pytest.approx(layer.thickness)


def test_no_lid_below_critical_metallicity():
    col = solve_column(ColumnParams(metallicity=3.0, teq=230.0))
    assert col.cold_trap.exists and not col.stable_layer.exists


def test_cold_trap_definition(fiducial):
    ct = fiducial.cold_trap
    xs = fiducial.params.saturation.x_sat(fiducial.T, fiducial.p)
    assert ct.exists and ct.x == pytest.approx(xs[fiducial.p > 3].min())
    assert 0 < ct.transmission < 1
    # rain-out profile: never increases upward, never exceeds saturation or the deep value
    x = fiducial.x
    assert np.all(np.diff(x) >= -1e-15) and np.all(x <= xs + 1e-15) and x[-1] == pytest.approx(fiducial.comp.x_h2o_deep)
    # input order must not matter
    assert np.allclose(water_profile(fiducial.p[::-1], fiducial.T[::-1], fiducial.comp.x_h2o_deep)[::-1], x)


def test_cold_trap_isothermal_analytic():
    """Isothermal column: x_sat = psat/p is smallest at the bottom, so the trap sits there."""
    p = np.geomspace(1e6, 1.0, 50)
    t = np.full_like(p, 250.0)
    ct = cold_trap(p, t, x_deep=0.5)
    from meromixis import Saturation

    assert ct.p == pytest.approx(1e6) and ct.x == pytest.approx(float(Saturation().psat(250.0)) / 1e6)


def test_cold_trap_nd_matches_column_loop(fiducial):
    cols = [solve_column(fiducial.params, teq=t, n_lev=120, n_sub=2) for t in (240.0, 270.0, 330.0)]
    p = np.stack([c.p for c in cols], axis=1)
    t = np.stack([c.T for c in cols], axis=1)
    out = cold_trap_nd(p, t, cols[0].comp.x_h2o_deep, axis=0, p_min=3.0)
    for i, c in enumerate(cols):
        assert out["transmission"][i] == pytest.approx(c.cold_trap.transmission)
    assert out["exists"].tolist() == [True, True, False]  # the hottest column is steam


def test_transition_is_monotonic_in_forcing_and_metallicity():
    base = ColumnParams(n_lev=200, n_sub=2)
    tr = [solve_column(base, teq=t).cold_trap.transmission for t in (230, 260, 290, 320)]
    assert all(b >= a for a, b in zip(tr, tr[1:])) and tr[0] < 0.1 and tr[-1] == 1.0
    t3, t100 = transition_teq(3.0, base, tol=1.0), transition_teq(100.0, base, tol=1.0)
    assert 200 < t3 < t100 < 350  # more water needs a warmer trap to get through


def test_regime_map_classification():
    ds = regime_map([240.0, 400.0], [3.0, 100.0], ColumnParams(n_lev=200, n_sub=2))
    assert ds["regime"].sel(teq=400.0).values.tolist() == [STEAM, STEAM]
    assert ds["regime"].sel(teq=240.0, metallicity=100.0) == LIDDED
    assert emulate_point(ColumnParams(teq=230.0, metallicity=3.0))["layer"] == 0.0


def test_calibration_recovers_known_parameters():
    truth = ColumnParams(kappa_th=2.0e-2, gamma_v=0.25, n_lev=200, n_sub=2)
    ref = solve_column(truth)
    best, rms = calibrate_grey(ref.p, ref.T, ColumnParams(n_lev=200, n_sub=2), p_range=(1e2, 1e5))
    assert rms < 1.0
    assert best.kappa_th == pytest.approx(truth.kappa_th, rel=0.15)


def test_transition_from_samples_agrees_with_bisection():
    base = ColumnParams(n_lev=200, n_sub=2)
    teqs = np.arange(250.0, 311.0, 5.0)
    tr = [solve_column(base, teq=t).cold_trap.transmission for t in teqs]
    lo, mid, hi = transition_from_samples(teqs, tr)
    assert lo < mid < hi
    assert mid == pytest.approx(transition_teq(100.0, base, tol=0.5), abs=2.5)
    assert np.isnan(transition_from_samples([300.0, 400.0], [1.0, 1.0])[1])  # never crosses
