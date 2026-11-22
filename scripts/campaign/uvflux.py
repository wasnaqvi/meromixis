"""A self-consistent floor for the NH3-dissociating photon flux.

docs/VERIFY.md C4 flags the kill test's Phi_UV = 1e13 photons cm-2 s-1 as a hand
estimate.  We can do better than a guess without a measured spectrum, by noting that
for a blackbody the band-to-bolometric ratio is fixed by Teff alone, and the
bolometric flux is already set by the forcing label (F_sub = 4 sigma teq^4).  That
gives the *photospheric* flux exactly, with no stellar radius needed.

The result is a hard LOWER bound, not an estimate: M dwarfs have chromospheric and
transition-region emission at 160-230 nm that exceeds the photospheric blackbody by
two to four orders of magnitude.  Because K_crit scales with eta * Phi_UV, a lower
Phi_UV means a lower K_crit, i.e. it makes the ammonia barrier HARDER to maintain.
So the blackbody floor is conservative against the barrier, and the ensemble spans
upward from it to a MUSCLES-like active M dwarf.
"""
from __future__ import annotations

import numpy as np

H_PLANCK, C_LIGHT, K_BOLTZ = 6.62607015e-34, 2.99792458e8, 1.380649e-23
SIGMA_SB = 5.670374419e-8
NH3_BAND = (160.0, 230.0)  # nm: NH3 photodissociation continuum


def planck_lambda(wl_nm, teff: float):
    """B_lambda in W m-2 sr-1 nm-1."""
    wl = np.asarray(wl_nm, float) * 1e-9
    a = 2.0 * H_PLANCK * C_LIGHT**2 / wl**5
    return a / np.expm1(H_PLANCK * C_LIGHT / (wl * K_BOLTZ * teff)) * 1e-9


def photospheric_photon_flux(teq: float, teff: float = 3457.0, band=NH3_BAND,
                             geometry: float = 0.25) -> float:
    """Global-mean photon flux (m-2 s-1) in ``band`` for a blackbody star.

    Normalises a Planck spectrum to the substellar bolometric flux 4 sigma teq^4,
    so no stellar radius or orbital distance is required.
    """
    from meromixis.barrier import photon_flux_from_spectrum

    wl = np.geomspace(50.0, 50000.0, 20000)
    b = planck_lambda(wl, teff)
    f_sub = 4.0 * SIGMA_SB * teq**4
    spec = b * (f_sub / np.trapezoid(b, wl))  # W m-2 nm-1 at the planet
    return photon_flux_from_spectrum(wl, spec, band=band, geometry=geometry)


def flux_grid(teq: float = 300.0, teff: float = 3457.0, n: int = 6,
              decades_above: float = 4.0) -> np.ndarray:
    """Phi_UV grid from the blackbody floor up by ``decades_above`` decades."""
    floor = photospheric_photon_flux(teq, teff)
    return floor * 10.0 ** np.linspace(0.0, decades_above, n)


if __name__ == "__main__":
    for t in (250.0, 300.0, 350.0):
        f = photospheric_photon_flux(t)
        print(f"teq={t:5.0f} K  photospheric NH3-band flux = {f:.3e} m-2 s-1 "
              f"= {f * 1e-4:.3e} cm-2 s-1")
    print(f"\nREADME hand estimate: 1e13 cm-2 s-1 = 1e17 m-2 s-1")
