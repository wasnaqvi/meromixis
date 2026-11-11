"""Bulk composition bookkeeping for an H2-dominated envelope.

The envelope is described by a single enrichment factor ``metallicity`` (x solar,
applied equally to C, N and O).  Below the clouds all O is assumed to be in H2O, all
C in CH4 and all N in NH3 unless the ``f_*`` fractions say otherwise; H2 takes
whatever hydrogen is left.  This is the closure behind the hand estimates in the
"Ammonia Barrier" explainer, made explicit so that the dry-gas mean molecular weight
(which enters q_cri) grows with metallicity as it should.

"Dry gas" always means *everything except water*.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .constants import CP_MOLAR, CP_T_GRID, MOLAR_MASS, SOLAR_X_H

M_V = MOLAR_MASS["H2O"]


def cp_molar(species: str, T):
    """Molar heat capacity (J mol-1 K-1), linear interpolation on the NASA-polynomial table."""
    return np.interp(T, CP_T_GRID, CP_MOLAR[species])


@dataclass(frozen=True)
class Composition:
    """Deep (below-cloud) composition of the envelope.

    Parameters
    ----------
    metallicity : enrichment of C, N, O relative to solar (linear, not dex).
    f_o_in_h2o, f_c_in_ch4, f_n_in_nh3 : fraction of each element carried by the
        hydride below the clouds.  The remainder is left out of the mass budget (it
        would sit in CO/CO2/N2), which is adequate at the few-percent level below
        ~800 K.
    """

    metallicity: float = 100.0
    f_o_in_h2o: float = 1.0
    f_c_in_ch4: float = 1.0
    f_n_in_nh3: float = 1.0
    he_per_h: float = SOLAR_X_H["He"]

    def __post_init__(self):
        if self.metallicity <= 0:
            raise ValueError("metallicity must be > 0 (linear, x solar)")
        if self._n["H2"] <= 0:
            raise ValueError(
                f"metallicity={self.metallicity:g}x solar leaves no H2 in the CH4/H2O/NH3 closure "
                "(limit ~420x). Use an explicit composition for metal-dominated envelopes."
            )

    @cached_property
    def _n(self) -> dict:
        """Molecules per H nucleus."""
        z = self.metallicity
        n = {
            "H2O": z * SOLAR_X_H["O"] * self.f_o_in_h2o,
            "CH4": z * SOLAR_X_H["C"] * self.f_c_in_ch4,
            "NH3": z * SOLAR_X_H["N"] * self.f_n_in_nh3,
            "He": self.he_per_h,
        }
        n["H2"] = 0.5 * (1.0 - 2 * n["H2O"] - 4 * n["CH4"] - 3 * n["NH3"])
        return n

    # ---- deep mole fractions (water included) ---------------------------------
    @cached_property
    def x_deep(self) -> dict:
        tot = sum(self._n.values())
        return {k: v / tot for k, v in self._n.items()}

    @property
    def x_h2o_deep(self) -> float:
        return self.x_deep["H2O"]

    # ---- dry gas (everything but water) ---------------------------------------
    @cached_property
    def x_dry(self) -> dict:
        n = {k: v for k, v in self._n.items() if k != "H2O"}
        tot = sum(n.values())
        return {k: v / tot for k, v in n.items()}

    @cached_property
    def mu_dry(self) -> float:
        """Mean molar mass of the dry gas, kg mol-1."""
        return sum(x * MOLAR_MASS[k] for k, x in self.x_dry.items())

    @cached_property
    def _cp_dry_table(self) -> np.ndarray:
        return sum(x * np.asarray(CP_MOLAR[k]) for k, x in self.x_dry.items())

    def cp_dry_molar(self, T):
        """Molar cp of the dry gas, J mol-1 K-1."""
        return np.interp(T, CP_T_GRID, self._cp_dry_table)

    @property
    def varpi(self) -> float:
        """Reduced mean-molecular-weight contrast, (M_v - mu_d)/M_v (Leconte+2017)."""
        return 1.0 - self.mu_dry / M_V

    @property
    def q_h2o_deep(self) -> float:
        """Deep water *mass* fraction of the gas (specific humidity)."""
        return float(x_to_q(self.x_h2o_deep, self.mu_dry))

    # ---- helpers for a gas with water mole fraction x --------------------------
    def mu(self, x_h2o):
        """Mean molar mass (kg mol-1) of dry gas + water at mole fraction x."""
        return x_h2o * M_V + (1.0 - x_h2o) * self.mu_dry

    def cp_molar_mix(self, T, x_h2o):
        return x_h2o * cp_molar("H2O", T) + (1.0 - x_h2o) * self.cp_dry_molar(T)


def x_to_q(x, mu_dry):
    """Water mole fraction -> mass fraction of the gas."""
    x = np.asarray(x, dtype=float)
    return x * M_V / (x * M_V + (1.0 - x) * mu_dry)


def q_to_x(q, mu_dry):
    """Water mass fraction of the gas -> mole fraction."""
    q = np.asarray(q, dtype=float)
    return (q / M_V) / (q / M_V + (1.0 - q) / mu_dry)
