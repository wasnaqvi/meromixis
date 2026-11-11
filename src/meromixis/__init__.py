"""meromixis: cold traps and chemoclines in temperate sub-Neptune atmospheres.

A meromictic lake never overturns: a light upper layer rides on a dense, solute-rich
deep layer, separated by a chemocline.  A water-rich hydrogen atmosphere can do the
same thing under its cloud deck.  This package holds the shared machinery for two
questions about that structure:

I.  the cold-trap transition (when does water get past the trap?), and
II. the ammonia barrier (how leaky is the chemocline to deep gas?).
"""
from .coldtrap import ColdTrap, cold_trap, cold_trap_nd, water_profile
from .column import Column, ColumnParams, solve_column
from .composition import Composition, q_to_x, x_to_q
from .inhibition import StableLayer, find_stable_layer, q_crit, x_crit
from .thermo import Saturation

__version__ = "0.1.0"
__all__ = [
    "ColdTrap", "Column", "ColumnParams", "Composition", "Saturation", "StableLayer",
    "cold_trap", "cold_trap_nd", "find_stable_layer", "q_crit", "q_to_x", "solve_column",
    "water_profile", "x_crit", "x_to_q",
]
