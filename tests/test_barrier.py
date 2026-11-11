"""Chemistry and the ammonia-barrier machinery."""
import numpy as np
import pytest

from meromixis import ColumnParams, solve_column
from meromixis.barrier import (
    Mixing,
    Photolysis,
    ammonia_budget,
    k_crit,
    k_layer_crit,
    kzz_profile,
    layer_mask,
    max_supply,
    photon_flux_from_spectrum,
    resistance,
    solve_diffusion_loss,
    steady_state,
)
from meromixis.chemistry import co_equilibrium, k_eq, nh3_equilibrium, quench


def test_equilibrium_constants_against_reference_values():
    assert k_eq("CH4+H2O=CO+3H2", 1000.0) == pytest.approx(26.5, rel=0.03)  # JANAF: ~26
    assert k_eq("N2+3H2=2NH3", 700.0) == pytest.approx(8.5e-5, rel=0.1)


def test_equilibrium_limits_and_mass_balance():
    # cold, dense: all N in NH3; hot, thin: almost none
    assert nh3_equilibrium(400.0, 1e7, 1e-3, 0.8) == pytest.approx(1e-3, rel=1e-3)
    assert nh3_equilibrium(1500.0, 1e5, 1e-3, 0.8) < 1e-5
    x_co = co_equilibrium(2500.0, 1e4, 1e-3, 5e-3, 0.8)
    assert x_co == pytest.approx(1e-3, rel=1e-2)  # C is the limiting element
    assert co_equilibrium(500.0, 1e7, 1e-3, 5e-3, 0.8) < 1e-12
    # the returned root satisfies the mass-action law
    t, p, xc, xo, xh = 1100.0, 3e7, 2e-2, 8e-2, 0.7
    x = co_equilibrium(t, p, xc, xo, xh)
    assert x * xh**3 * (p / 1e5) ** 2 / ((xc - x) * (xo - x)) == pytest.approx(k_eq("CH4+H2O=CO+3H2", t), rel=1e-8)


def test_diffusion_loss_solver_against_analytic_solution():
    """Uniform K n and loss: f = cosh((L - z)/lam) / cosh(L/lam)."""
    z = np.linspace(0.0, 1.0e5, 801)
    n = np.full_like(z, 1e25)
    kk, loss = 10.0, 1e-8
    lam = np.sqrt(kk / loss)
    f = solve_diffusion_loss(z, kk * n, n, np.full_like(z, loss), 1.0)
    assert np.allclose(f, np.cosh((z[-1] - z) / lam) / np.cosh(z[-1] / lam), rtol=2e-3)


def test_solver_conserves_tracer_without_loss():
    z = np.linspace(0.0, 5e4, 60)
    n = 1e25 * np.exp(-z / 2e4)
    kn = n * np.where(z < 2e4, 1e3, 1e-2)  # sharp K contrast
    assert np.allclose(solve_diffusion_loss(z, kn, n, np.zeros_like(z), 3e-4), 3e-4)


@pytest.fixture(scope="module")
def cols():
    p = ColumnParams(metallicity=30.0)
    return p, solve_column(p), solve_column(p, inhibition=False)


def test_resistances_add_and_the_lid_dominates(cols):
    _, ci, _ = cols
    layer = ci.stable_layer
    k = kzz_profile(ci, Mixing(k_layer=0.1))
    assert set(np.unique(k[layer_mask(ci)])) == {0.1}
    p_up = 0.5 * layer.p_top
    r_all = resistance(ci, k, 1e7, p_up)
    r_lid = resistance(ci, k, layer.p_base, layer.p_top)
    assert r_lid <= r_all and r_lid / r_all > 0.9
    assert max_supply(ci, k, 1e-3, 1e7, p_up) == pytest.approx(1e-3 / r_all)
    assert resistance(ci, k, 1e7, 1e5) + resistance(ci, k, 1e5, p_up) == pytest.approx(r_all, rel=0.02)


def test_k_crit_simple_scaling():
    assert k_crit(1e17, 5e4, 2.4e24, 1e-4) == pytest.approx(2.08e1, rel=0.01)  # the explainer's numbers, in SI
    assert k_crit(2e17, 5e4, 2.4e24, 1e-4) == pytest.approx(2 * k_crit(1e17, 5e4, 2.4e24, 1e-4))


def test_front_regimes(cols):
    p, ci, _ = cols
    q = quench(ci, "NH3")
    leaky = steady_state(ci, q.x, q.p, Mixing(k_layer=50.0), Photolysis(eta=1.0))
    tight = steady_state(ci, q.x, q.p, Mixing(k_layer=1e-3), Photolysis(eta=1.0))
    assert not leaky.supply_limited and leaky.column_loss == pytest.approx(1e17)  # every photon used
    assert tight.supply_limited and tight.column_loss == pytest.approx(tight.supply_max, rel=1e-3)
    assert tight.p_front > leaky.p_front and tight.hidden() and not leaky.hidden()
    assert np.all(np.diff(leaky.f) >= -1e-18)  # mixing ratio never increases upward (p is top-down)
    assert tight.depletion < 1e-2 < leaky.depletion


def test_k_layer_crit_is_the_visibility_threshold(cols):
    p, ci, _ = cols
    q = quench(ci, "NH3")
    phot = Photolysis(eta=1.0)
    kc = k_layer_crit(ci, q.x, q.p, Mixing(), phot)
    assert np.isfinite(kc) and kc > 0
    assert kc == pytest.approx(k_layer_crit(ci, q.x, q.p, Mixing(k_layer=37.0), phot))  # independent of the guess
    assert steady_state(ci, q.x, q.p, Mixing(k_layer=kc / 3), phot).hidden()
    assert not steady_state(ci, q.x, q.p, Mixing(k_layer=kc * 3), phot).hidden()


def test_budget_separates_thermal_and_transport(cols):
    p, ci, cr = cols
    b = ammonia_budget(p, Mixing(k_layer=0.05), Photolysis(eta=1.0), ci, cr)
    assert 0 < b.thermal_factor <= 1.0 + 1e-9  # hotter interior never makes more NH3
    assert b.co_enhancement >= 1.0 - 1e-9  # ... and never less CO
    assert b.hidden_by_barrier and "hides" in b.verdict()
    assert b.co_fill_time > 0 and b.nh3_depletion_time > 0
    hot = ammonia_budget(ColumnParams(metallicity=100.0), Mixing(k_layer=0.05))
    assert hot.thermal_factor < 0.8 and hot.co_enhancement > 2.0
    none = ammonia_budget(ColumnParams(metallicity=3.0))
    assert np.isnan(none.k_crit) and "no stable layer" in none.verdict()


def test_photon_flux_integration():
    wl = np.linspace(100.0, 300.0, 2001)
    flux = np.full_like(wl, 1e-3)  # W m-2 nm-1
    got = photon_flux_from_spectrum(wl, flux, band=(160.0, 230.0), geometry=1.0)
    h, c = 6.62607015e-34, 2.99792458e8
    expect = 1e-3 * 1e-9 / (h * c) * 0.5 * (230.0**2 - 160.0**2)
    assert got == pytest.approx(expect, rel=1e-3)
