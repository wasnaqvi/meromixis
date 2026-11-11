"""Physical constants and reference data (SI throughout unless stated)."""
from __future__ import annotations

R_GAS = 8.314462618  # J mol-1 K-1
K_B = 1.380649e-23  # J K-1
N_A = 6.02214076e23  # mol-1
SIGMA_SB = 5.670374419e-8  # W m-2 K-4
G_NEWTON = 6.67430e-11  # m3 kg-1 s-2
R_EARTH = 6.3781e6  # m
M_EARTH = 5.9722e24  # kg
AU = 1.495978707e11  # m
S_EARTH = 1361.0  # W m-2
BAR = 1.0e5  # Pa

# Molar masses, kg mol-1
MOLAR_MASS = {
    "H2": 2.01588e-3,
    "He": 4.002602e-3,
    "H2O": 18.01528e-3,
    "CH4": 16.0425e-3,
    "NH3": 17.03052e-3,
    "CO": 28.0101e-3,
    "N2": 28.0134e-3,
}

# Solar photospheric abundances, log10(N_X / N_H) + 12  (Asplund, Amarsi & Grevesse 2021)
SOLAR_LOG_EPS = {"He": 10.914, "C": 8.46, "N": 7.83, "O": 8.69}
SOLAR_X_H = {k: 10.0 ** (v - 12.0) for k, v in SOLAR_LOG_EPS.items()}

# Molar heat capacities cp(T), J mol-1 K-1, tabulated from NASA polynomials
# (GRI-Mech 3.0 thermo via Cantera 3.2; regenerate with scripts/dev/make_thermo_tables.py).
CP_T_GRID = (200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 3000.0)
CP_MOLAR = {
    "H2": (27.527, 28.851, 29.277, 29.298, 29.255, 29.616, 30.163, 32.279, 34.257, 37.066),
    "CH4": (33.439, 35.761, 40.530, 46.494, 52.730, 63.999, 73.617, 90.414, 100.436, 111.613),
    "NH3": (33.702, 35.699, 38.650, 41.994, 45.336, 51.252, 56.529, 66.416, 72.865, 79.051),
    "H2O": (33.350, 33.596, 34.268, 35.214, 36.321, 38.733, 41.295, 47.291, 51.752, 56.791),
    "N2": (28.793, 29.075, 29.320, 29.636, 30.087, 31.394, 32.762, 34.805, 35.988, 37.028),
    "CO": (29.133, 29.143, 29.375, 29.811, 30.415, 31.894, 33.163, 35.211, 36.247, 37.209),
    "He": (20.786,) * 10,
}

# Water: triple and critical points
T_TRIPLE_H2O = 273.16  # K
P_TRIPLE_H2O = 611.657  # Pa
T_CRIT_H2O = 647.096  # K
P_CRIT_H2O = 22.064e6  # Pa
