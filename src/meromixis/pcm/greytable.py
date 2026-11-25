"""Grey opacity that varies with pressure, carried by a correlated-k table.

Generic PCM's grey mode (``graybody = .true.``) adds a constant ``kappa_IR`` to the
gaseous k-coefficients.  A constant kappa caps the deep radiative gradient at
dlnT/dlnp = 1/4, below an H2 adiabat (R/cp ~ 2/7), so such a column never convects
at depth and the inhibition switch has nothing to act on.  kappa ~ p**n lifts the cap
to (1 + n)/4.

The model has no key for that, but it does not need one.
``rad_correlatedk_read_opacity_tables.F90`` reads ``corrk_gcm_IR.dat`` for any
corrkdir not starting with "null", adds kappa_IR, and interpolates log10(k) in
log10(p) with 4-point Lagrange polynomials.  A table holding kappa_th (p/p_ref)**n,
the same at every T, water amount, band and g-point, is therefore grey and is
reproduced exactly between its pressure nodes.  kappa_IR is then only a floor.

The stellar opacity scales the same way (kappa_vi (p/p_ref)**n, kappa_VI a floor), as
in the emulator, whose visible optical depth is gamma_v times the thermal one.  With
kappa_VI constant the top of the model absorbs starlight but can hardly emit, and in
the first try it heated past 700 K within days.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

M_P = 1.672621e-27  # kg, the factor the model itself uses to convert kappa_IR
KAPPA_FLOOR = 1.0e-12  # m2 kg-1, kappa_IR when the table carries the opacity


def kappa_to_cross_section(kappa, mugaz_ref):
    """m2 kg-1 -> cm2 per molecule, as the model converts kappa_IR (mugaz_ref in g/mol)."""
    return np.asarray(kappa, float) * 1.0e4 * mugaz_ref * M_P


def power_law_table(log10p_mbar, n_t, n_q, n_band, n_g, kappa_th, kappa_pexp, mugaz_ref, p_ref=1.0e5):
    """k(T, p, q, band, g) in cm2/molecule, constant except along p.  Write it with order="F"."""
    p = 10.0 ** np.asarray(log10p_mbar, float) * 100.0  # mbar -> Pa
    k = kappa_to_cross_section(kappa_th * (p / p_ref) ** kappa_pexp, mugaz_ref)
    return np.broadcast_to(k[None, :, None, None, None], (n_t, len(p), n_q, n_band, n_g)).copy()


def _counted(path: Path) -> list[str]:
    """Tokens of an axis file whose first token is the count."""
    tok = path.read_text().split()
    return tok[1: 1 + int(tok[0])]


T_AXIS = (40.0, 55.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 200.0, 240.0, 280.0, 320.0,
          360.0, 400.0, 450.0, 500.0, 600.0, 700.0, 850.0, 1000.0, 1200.0, 1500.0)


def write_power_law_dir(axes_dir, out_dir, bands, kappa_th, kappa_vi, kappa_pexp, mugaz_ref,
                        p_ref=1.0e5, t_axis=T_AXIS) -> Path:
    """A corrk directory with the p, q, g and band files of ``axes_dir`` (linked) and power-law k data.

    ``out_dir``'s name must not start with "null", or the model zeroes the table.  It
    gets its own T axis: the model stops if any level is hotter than the last node
    (700 K in the borrowed Jupiter axes, below the 714 K bottom of a lidded column),
    and the k-values do not depend on T, so a wider axis costs nothing.
    """
    axes_dir, out_dir = Path(axes_dir), Path(out_dir)
    if out_dir.name.startswith("null"):
        raise ValueError("a corrkdir starting with 'null' is zeroed by the model")
    (out_dir / bands).mkdir(parents=True, exist_ok=True)
    for f in ("p.dat", "Q.dat", "g.dat"):
        _link(axes_dir / f, out_dir / f)
    if (out_dir / "T.dat").is_symlink():
        (out_dir / "T.dat").unlink()
    (out_dir / "T.dat").write_text(f"{len(t_axis)}\n" + "".join(f"{t:.1f}\n" for t in t_axis))
    for f in ("narrowbands_IR.in", "narrowbands_VI.in"):
        _link(axes_dir / bands / f, out_dir / bands / f)

    n_t = len(t_axis)
    logp = [float(v) for v in _counted(axes_dir / "p.dat")]
    q = (axes_dir / "Q.dat").read_text().split()
    n_q = int(q[1 + int(q[0])])  # ngas, gas names, then the number of water values
    n_g = len(_counted(axes_dir / "g.dat"))

    # The model splits bands at IR_VI_WNLIMIT by band centre and adds kappa_IR to the
    # long-wave side of *both* band sets (3 of the 25 stellar bands at 20x25) and
    # kappa_VI to the short-wave side.  The table follows the same split.
    for name, kind in (("corrk_gcm_IR.dat", "IR"), ("corrk_gcm_VI.dat", "VI")):
        centres = band_centres(axes_dir / bands / f"narrowbands_{kind}.in")
        nb, nl = len(centres), n_longwave(centres)
        k = power_law_table(logp, n_t, n_q, nb, n_g, kappa_th, kappa_pexp, mugaz_ref, p_ref)
        k[:, :, :, nl:, :] = power_law_table(logp, n_t, n_q, nb - nl, n_g, kappa_vi, kappa_pexp, mugaz_ref, p_ref)
        np.savetxt(out_dir / bands / name, k.ravel(order="F"), fmt="%.6e")
    return out_dir


IR_VI_WNLIMIT = 3000.0  # cm-1, rad_correlatedk_read_opacity_tables.F90 (graybody)


def band_centres(path: Path) -> np.ndarray:
    """Band centres (cm-1) of a narrowbands_*.in file: the count, then one edge pair per band."""
    tok = path.read_text().split()
    edges = np.asarray(tok[1: 1 + 2 * int(tok[0])], float).reshape(-1, 2)
    return edges.mean(axis=1)


def n_longwave(centres) -> int:
    """Bands the model treats as infrared: those before the first centre above IR_VI_WNLIMIT."""
    above = np.nonzero(np.asarray(centres) > IR_VI_WNLIMIT)[0]
    return int(above[0]) if len(centres) > 1 and len(above) else len(centres)


def _link(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    os.symlink(src, dst)
