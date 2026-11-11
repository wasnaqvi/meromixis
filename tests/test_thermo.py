"""Thermodynamics: saturation curve, composition closure, adiabats, inhibition threshold."""
import numpy as np
import pytest

from meromixis import Composition, Saturation, q_crit, q_to_x, x_to_q
from meromixis.constants import MOLAR_MASS, R_GAS
from meromixis.inhibition import metallicity_crit
from meromixis.thermo import dry_gradient, moist_gradient


def test_saturation_reference_points():
    sat = Saturation()
    assert sat.psat(273.16) == pytest.approx(611.657, rel=2e-3)  # triple point
    assert sat.psat(373.124) == pytest.approx(101325.0, rel=2e-3)  # normal boiling point
    assert sat.psat(230.0) == pytest.approx(8.95, rel=0.02)  # ice (Murphy & Koop)
    assert np.isinf(sat.psat(700.0))  # supercritical: no saturation
    assert sat.gamma(700.0) == 0.0


def test_scalar_and_array_paths_agree():
    sat = Saturation()
    for t in (180.0, 250.0, 273.16, 300.0, 450.0, 640.0):
        assert sat.psat_s(t) == pytest.approx(float(sat.psat(t)), rel=1e-12)
        assert sat.gamma_s(t) == pytest.approx(float(sat.gamma(t)), rel=1e-6)


def test_clausius_clapeyron_mode_is_the_pcm_form():
    sat = Saturation("cc", latent_heat=2.5e6)
    t = np.array([250.0, 300.0, 400.0])
    expect = sat.p_ref * np.exp(-(2.5e6 * MOLAR_MASS["H2O"] / R_GAS) * (1 / t - 1 / sat.t_ref))
    assert np.allclose(sat.psat(t), expect)
    assert np.allclose(sat.gamma(t), 2.5e6 * MOLAR_MASS["H2O"] / (R_GAS * t))


def test_latent_heat_from_curve_is_sane():
    assert Saturation().latent_heat_cc(300.0) == pytest.approx(2.44e6, rel=0.01)


def test_composition_closure():
    c = Composition(100.0)
    assert sum(c.x_deep.values()) == pytest.approx(1.0)
    assert sum(c.x_dry.values()) == pytest.approx(1.0)
    assert c.x_h2o_deep == pytest.approx(0.089, abs=0.002)
    assert Composition(1.0).mu_dry * 1e3 == pytest.approx(2.31, abs=0.02)
    assert c.mu_dry > Composition(1.0).mu_dry  # CH4 and NH3 load the dry gas
    with pytest.raises(ValueError):
        Composition(500.0)


def test_mole_mass_fraction_round_trip():
    mu_d = Composition(10.0).mu_dry
    x = np.array([1e-6, 1e-2, 0.3, 0.9])
    assert np.allclose(q_to_x(x_to_q(x, mu_d), mu_d), x)


def test_q_crit_matches_leconte_2017():
    # Leconte+2017 quote q_cri ~ 0.07 for water in H2 near 300 K, i.e. ~10x solar.
    assert float(q_crit(300.0, Composition(1.0))) == pytest.approx(0.065, abs=0.005)
    assert 8.0 < metallicity_crit(300.0) < 14.0
    # constant-L closed form R T / ((M_v - mu_d) L)
    c, sat = Composition(1.0), Saturation("cc")
    closed = R_GAS * 350.0 / ((MOLAR_MASS["H2O"] - c.mu_dry) * sat.latent_heat)
    assert float(q_crit(350.0, c, sat)) == pytest.approx(closed, rel=1e-12)


def test_moist_gradient_limits():
    c, sat = Composition(1.0), Saturation()
    # dilute limit -> dry adiabat of the dry gas
    assert moist_gradient(150.0, 1e5, c, sat) == pytest.approx(dry_gradient(150.0, 0.0, c), rel=1e-3)
    # steam limit -> the vapour-pressure curve
    t = 400.0
    assert moist_gradient(t, sat.psat_s(t) * 0.999999, c, sat) == pytest.approx(1.0 / sat.gamma_s(t), rel=1e-3)
    # in between: condensation flattens the profile
    assert moist_gradient(300.0, 1e5, c, sat) < dry_gradient(300.0, 0.0, c)
