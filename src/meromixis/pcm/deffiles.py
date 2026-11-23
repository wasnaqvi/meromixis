"""Render Generic PCM run directories for the sweep.

Every key written here carries a provenance tag in the generated file:

``[r4146]``  seen in the LMDZ.GENERIC trunk at changeset 4146 (19 Mar 2026: the
             "thermodynamics update" that renamed ``cpp``/``mugaz`` to
             ``cppd_ref``/``mugaz_ref`` and added ``thermo_phy``) or in ``deftank``;
``[wiki]``   documented on the LMD Planets wiki / the project primer;
``[verify]`` written from memory of older versions.  Diff against the ``deftank``
             example closest to your checkout before the first run (docs/VERIFY.md).

The generator never invents physics: it only translates a :class:`RunSpec` into the
text files the model reads, and writes a JSON manifest so each run is reproducible.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from ..composition import Composition, cp_molar
from ..constants import AU, G_NEWTON, M_EARTH, MOLAR_MASS, R_EARTH, SIGMA_SB


def _b(flag: bool) -> str:
    return ".true." if flag else ".false."


@dataclass(frozen=True)
class Planet:
    name: str = "K2-18b-like"
    radius_rearth: float = 2.61
    mass_mearth: float = 8.63
    rotation_period_days: float = 32.94  # synchronous: equals the orbital period
    semi_major_axis_au: float = 0.1429

    @property
    def radius(self) -> float:
        return self.radius_rearth * R_EARTH

    @property
    def gravity(self) -> float:
        return G_NEWTON * self.mass_mearth * M_EARTH / self.radius**2

    @property
    def omega(self) -> float:
        return 2.0 * np.pi / (self.rotation_period_days * 86400.0)


@dataclass(frozen=True)
class Star:
    teff: float = 3457.0
    model: str = "blackbody"  # "blackbody" | "startype" | "file"
    startype: int = 1
    spectrum_file: str = ""


@dataclass(frozen=True)
class RunSpec:
    """One simulation.  ``teq`` is a *label for the forcing*: F = 4 sigma teq^4 (zero albedo,
    full redistribution).  The model computes its own albedo."""

    run_id: str = "T300_Z100_fid"
    tier: str = "spine"
    teq: float = 300.0
    metallicity: float = 100.0
    tint: float = 50.0
    planet: Planet = field(default_factory=Planet)
    star: Star = field(default_factory=Star)
    # grid
    one_d: bool = False
    nlon: int = 64
    nlat: int = 48
    nlev: int = 40
    p_bot: float = 8.0e6  # Pa (Charnay+2021: 80 bar)
    p_top: float = 20.0  # Pa (0.2 mbar)
    # physics
    water_scheme: str = "generic"  # "generic" (GCS route, has inhibition) | "legacy" | "dry"
    condensable_name: str = "H2O"  # tracer stem for the GCS route; see VERIFY P2
    inhibition: bool = True
    virtual_theta: bool = True
    mass_redistrib: bool = False
    clouds_radiative: bool = False
    greybody: bool = False
    kappa_ir: float = 1.0e-2
    kappa_vi: float = 4.0e-3
    corrkdir: str = "SET_ME_H2-He-H2Ovar"  # no public table has x_H2O > 1e-2 in H2/He
    bands: str = "38x36"
    datadir: str = "./datagcm"
    # time stepping / output
    nday: int = 500
    day_step: int = 480
    iphysiq: int = 10
    iradia: int = 4
    output_every_days: float = 25.0
    segments: int = 40
    branch_from: str = ""  # run_id whose restart initialises this run
    notes: str = ""

    # ---- derived quantities --------------------------------------------------
    @property
    def composition(self) -> Composition:
        return Composition(self.metallicity)

    @property
    def stellar_flux(self) -> float:
        """Substellar flux at the planet, W m-2."""
        return 4.0 * SIGMA_SB * self.teq**4

    @property
    def fat1au(self) -> float:
        """Stellar flux scaled to 1 AU.  The sweep varies this at fixed orbit and rotation."""
        return self.stellar_flux * self.planet.semi_major_axis_au**2

    @property
    def intheat(self) -> float:
        return SIGMA_SB * self.tint**4

    @property
    def mugaz_ref(self) -> float:
        """Dry-gas molar mass, g mol-1."""
        return self.composition.mu_dry * 1.0e3

    @property
    def cppd_ref(self) -> float:
        """Dry-gas specific heat at 300 K, J kg-1 K-1."""
        return float(self.composition.cp_dry_molar(300.0)) / self.composition.mu_dry

    @property
    def cppv_ref(self) -> float:
        return float(cp_molar("H2O", 300.0)) / MOLAR_MASS["H2O"]

    @property
    def qvap_deep(self) -> float:
        return 0.0 if self.water_scheme == "dry" else self.composition.q_h2o_deep

    @property
    def physics_steps_per_day(self) -> float:
        return self.day_step / self.iphysiq

    @property
    def diagfi_output_rate(self) -> int:
        return max(1, int(round(self.output_every_days * self.physics_steps_per_day)))

    @property
    def dims(self) -> str:
        return f"{self.nlev}" if self.one_d else f"{self.nlon}x{self.nlat}x{self.nlev}"

    @property
    def tracer_names(self) -> tuple[str, ...]:
        if self.water_scheme == "generic":
            return (f"{self.condensable_name}_vap", f"{self.condensable_name}_ice")
        if self.water_scheme == "legacy":
            return ("h2o_vap", "h2o_ice")
        return ()

    def build_command(self, arch: str = "FIR") -> str:
        target = "rcm1d" if self.one_d else "gcm"
        par = "" if self.one_d else " -parallel mpi_omp"
        return f"./makelmdz_fcm -arch {arch} -p generic -d {self.dims} -b {self.bands}{par} {target}"


# ---- individual files --------------------------------------------------------------------

def comments_on_own_lines(text: str) -> str:
    """Move every trailing ``#`` comment onto its own line above its key.

    IOIPSL's ``getin`` silently drops some keys that carry a trailing comment, and
    whether it does depends on the length of the key: at r4609
    ``virtual_theta_correction = .true.  # ...`` is read, while
    ``moist_convection_inhibition = .true.  # ...`` on the very next line is not, and
    the model falls back to the default with only a ``USING DEFAULTS`` line in the log
    to say so.  That silently disabled the central switch of this project -- the
    inhibited and non-inhibited runs came out bit-identical -- so no generated .def
    file may carry a trailing comment on a key line.

    Comments on their own lines are parsed correctly, so this keeps the provenance
    annotations without risking the key.
    """
    out = []
    for raw in text.splitlines():
        line = raw.rstrip()
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#") or "=" not in line:
            out.append(line)
            continue
        i = line.find("#")
        if i > 0:
            out.append(line[i:].rstrip())
            out.append(line[:i].rstrip())
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def render_callphys(s: RunSpec) -> str:
    gen, leg = s.water_scheme == "generic", s.water_scheme == "legacy"
    star = {
        "blackbody": f"stelbbody = .true.                 # [wiki]\nstelTbb   = {s.star.teff:.1f}              # [wiki]",
        "startype": f"stelbbody = .false.\nstartype  = {s.star.startype}                       # [r4146]",
        "file": f"stelbbody = .false.\nstartype  = 0                       # [verify]\nstelspec_file = {s.star.spectrum_file}   # [wiki]",
    }[s.star.model]
    lines = f"""# callphys.def generated by meromixis for run {s.run_id}
# forcing label Teq = {s.teq:g} K  ->  F = {s.stellar_flux:.1f} W m-2 ;  [M/H] = {s.metallicity:g} x solar ;  Tint = {s.tint:g} K
# provenance tags: [r4146] trunk changeset 4146 / deftank, [wiki] LMD Planets wiki, [verify] check before use

## Orbit / general options
tracer     = {_b(s.water_scheme != "dry")}         # [r4146]
diurnal    = .false.                # [wiki] tidally locked: fixed substellar point
tlocked    = .true.                 # [wiki]
season     = .false.                # [verify]
lwrite     = .false.                # [r4146]
callstats  = .false.                # [r4146]
enertest   = .false.                # [r4146] on, it prints the energy budget every step: 1.6 GB per 1-D run
diagfi_output_rate = {s.diagfi_output_rate}          # [r4146] in physics steps (= every {s.output_every_days:g} days)

## External data
datadir    = {s.datadir}            # [r4146]

## Radiative transfer
callrad    = .true.                 # [r4146]
iradia     = {s.iradia}                      # [r4146] radiation every iradia physics steps
corrk      = .true.                 # [r4146]
corrkdir   = {s.corrkdir}   # [r4146] must match -b {s.bands} and gases.def
continuum_CIA      = {_b(not s.greybody)}         # [r4609] H2-H2, H2-He collision-induced absorption
continuum_far_wing = {_b(not s.greybody)}         # [r4609] far-wing (H2O) continuum
# NB: a bare `continuum` key is NOT read by the model; it was silently ignored.
# In graybody mode the gaseous k-coefficients are zeroed and kappa_IR/kappa_VI carry
# the opacity, so the continuum must be off or it would be counted on top of it
# (and it needs datagcm/continuum/CIA_data, which a null corrk run does not ship).
# Aerosols.  aeroco2 and aeroh2o default to .true. in the model and each adds one
# to the aerosol count, which then demands Mars/CO2-ice optical property files that
# have no business in an H2 envelope.  Turn them off explicitly; naerkind must equal
# the number actually switched on (aerosol_global_variables.F90 aborts otherwise).
aeroco2    = .false.                # [r4609]
aeroh2o    = .false.                # [r4609]
# aerosol_global_variables.F90 forces ia=1 with noaero=.true. when nothing is
# switched on, so naerkind is 1 even for a run with no radiatively active aerosol.
# The CO2-ice optical property files are then read but never used.
naerkind   = 1                      # [r4609] must match the model's own count
callgasvis = .true.                 # [r4146]
rayleigh   = .true.                 # [r4146]
graybody   = {_b(s.greybody)}                # [deftank r596] double-grey gate test
kappa_IR   = {s.kappa_ir:g}                  # [deftank r596] m2 kg-1, only if graybody
kappa_VI   = {s.kappa_vi:g}                  # [deftank r596] m2 kg-1, only if graybody
specOLR    = .false.                # [wiki]
meanOLR    = .false.                # [r4146]
kastprof   = .false.                # [r4146]
varactive  = {_b(s.water_scheme != "dry")}         # [r4146] variable gas (water) radiatively active
varfixed   = .false.                # [r4146]

## Star
{star}
Fat1AU     = {s.fat1au:.4f}             # [wiki] W m-2; the ONLY quantity the Teq sweep changes

## Reference thermodynamics (renamed at r4146: cpp -> cppd_ref, mugaz -> mugaz_ref)
thermo_phy     = thermo_uni_ideal   # [r4146] the only mode available with the lon-lat core
cpp_mugaz_mode = 1                  # [r4146] 1: read the two values below
mugaz_ref      = {s.mugaz_ref:.4f}             # [r4146] g mol-1, dry gas (H2, He, CH4, NH3)
cppd_ref       = {s.cppd_ref:.1f}            # [r4146] J kg-1 K-1, dry gas at 300 K
cppv_ref       = {s.cppv_ref:.1f}             # [r4146] J kg-1 K-1, required if generic_condensation

## Lower boundary: no surface, internal heat flux
nosurf     = .true.                 # [r4146]
callsoil   = .false.                # [r4146] required by nosurf
intheat    = {s.intheat:.4f}               # [r4146] W m-2 (= sigma Tint^4)

## Mixing
calladj    = .true.                 # [r4146] dry convective adjustment
calldifv   = .true.                 # [r4146]
UseTurbDiff = .true.                # [r4146]
"""
    if gen:
        lines += f"""
## Water as a generic condensable species (the route that carries convective inhibition)
water      = .false.
generic_condensation        = .true.     # [r4146]
generic_rain                = .true.     # [r4146]
evap_prec_generic           = .true.     # [r4146]
precip_scheme_generic       = 1          # [r4146]
moistadjustment_generic     = .true.     # [r4146]
moist_convection_inhibition = {_b(s.inhibition)}    # [r4146] q > q_cri switches moist adjustment off
virtual_theta_correction    = {_b(s.virtual_theta)}    # [r4146] mu effect in the dry adjustment
metallicity                 = 0.0        # [r4146] log10; only shifts psat of refractory species
qvap_deep                   = {s.qvap_deep:.6f}   # [r4146] kg/kg, vapour held at this value at the bottom
mass_redistrib              = {_b(s.mass_redistrib)}    # [r4146] pressure change from condensation (variant run)
"""
    elif leg:
        lines += f"""
## Legacy water cycle (Charnay+2021-style control; no moist_convection_inhibition here)
water      = .true.                 # [r4146]
watercond  = .true.                 # [r4146]
waterrain  = .true.                 # [r4146]
moistadjustment = .true.            # [r4146]
qvap_deep  = {s.qvap_deep:.6f}            # [r4146] read by largescale
mass_redistrib = {_b(s.mass_redistrib)}           # [r4146]
"""
    else:
        lines += "\n## Dry control (Innes & Pierrehumbert 2022 benchmark): no water\nwater = .false.\n"
    if s.one_d:
        lines += """
## 1-D global-mean column
global1d   = .true.                 # [r4146]
szangle    = 60.                    # [r4146] mean stellar zenith angle
"""
    lines += f"""
## Clouds
# Radiatively active condensate clouds need is_rgcs=1 in traceur.def, an aerogeneric count
# here and a build with -t <n>.  First pass keeps them passive (clouds_radiative = {str(s.clouds_radiative).lower()}).
aerogeneric = {1 if s.clouds_radiative else 0}                     # [wiki]
"""
    return lines


def render_run_def(s: RunSpec) -> str:
    return f"""# run.def generated by meromixis for run {s.run_id}   (keys from deftank/run.def.* unless tagged)
planet_type = generic
nday        = {s.nday}
day_step    = {s.day_step}
iperiod     = 5
iconser     = {s.day_step}
idissip     = 5
lstardis    = .true.
hybrid      = .true.
nitergdiv   = 1
nitergrot   = 2
niterh      = 2
tetagdiv    = 4000.
tetagrot    = 5000.
tetatemp    = 5000.
coefdis     = 0.
purmats     = .false.
iphysiq     = {s.iphysiq}
clon        = 0.
clat        = 0.
grossismx   = 1.
grossismy   = 1.
fxyhypb     = .false.
dzoomx      = 0.
dzoomy      = 0.
taux        = 2.
tauy        = 2.
# upper sponge: treat as a parameter to test (jets are sensitive to it)        [verify]
iflag_top_bound = 1
mode_top_bound  = 1
tau_top_bound   = 1.e-5
INCLUDEDEF=callphys.def
"""


def pseudo_altitudes(s: RunSpec) -> np.ndarray:
    """Layer pseudo-altitudes for z2sig.def with H = 1/ln(10): z = log10(p_bot / p)."""
    edges = np.linspace(0.0, np.log10(s.p_bot / s.p_top), s.nlev + 1)
    return 0.5 * (edges[1:] + edges[:-1])


def render_z2sig(s: RunSpec) -> str:
    z = pseudo_altitudes(s)
    return "0.43429448\n" + "\n".join(f"{v:.6f}" for v in z) + "\n"


def render_gases(s: RunSpec) -> str:
    xd = s.composition.x_dry
    tot = xd["H2"] + xd["He"]
    names, fracs = ["H2_", "He_"], [xd["H2"] / tot, xd["He"] / tot]
    if s.water_scheme != "dry":
        names.append("H2O" if s.water_scheme == "legacy" else s.condensable_name)
        fracs.append(-1.0)
    body = "\n".join(names) + "\n" + "\n".join("-1" if f < 0 else f"{f:.6f}" for f in fracs)
    # su_gases.F90 does a bare `read(90,*)` and then reads ngasmx from the NEXT
    # line, so gases.def may carry exactly ONE header line.  Two comment lines make
    # the model abort with "error reading number of gases".
    return (
        "# meromixis: H2/He renormalised, water variable (-1); order must match the corrk Q.dat\n"
        f"{len(names)}\n{body}\n"
    )


def render_traceur(s: RunSpec) -> str:
    names = s.tracer_names
    if s.water_scheme == "generic":
        opt = " is_condensable=1" + (" is_rgcs=1" if s.clouds_radiative else "")
        return "#ModernTrac-v1\n" + f"{len(names)}\n" + "\n".join(n + opt for n in names) + "\n"
    return f"{len(names)}\n" + "".join(n + "\n" for n in names)


def render_diagfi(s: RunSpec) -> str:
    """Variables to write to diagfi.

    "temperature", not "temp": the generic physics writes
    ``write_output("temperature",...)`` (physiq_mod.F90) and has no "temp", so the
    old name silently produced output files with no temperature in them.

    ``vteta`` and ``RH_generic`` are written by the model itself and are exactly the
    two fields the lid diagnosis needs: the virtual potential temperature is the
    Ledoux criterion Leconte et al. (2024) recast in their Eq. 21, and the relative
    humidity settles whether the layer is sub-saturated without having to infer it.
    """
    base = ["temperature", "teta", "vteta", "ps", "p", "u", "v", "w",
            "ASR", "OLR", "ISR", "tsurf"]
    if s.water_scheme == "generic":
        base.append("RH_generic")
    elif s.water_scheme == "legacy":
        base.append("RH")
    return "\n".join(base + list(s.tracer_names)) + "\n"


def initial_column(s: RunSpec) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Emulator equilibrium on the model's layers: (p_layer [Pa], T_layer [K], T_bottom [K], q_vap [kg/kg]).

    rcm1d's own initial state is isothermal at ``tref`` with water only in the bottom
    layer.  An isothermal column is stable, so convection never starts and water never
    leaves the bottom layer: after 13 000 days the G1 twins had x_H2O ~ 1e-15 above
    10 bar and no cloud, cold trap or lid at all, and inhibition had nothing to act
    on.  The deep column cannot spin up from there either -- its heat capacity is
    ~5e9 J m-2 K-1 against 0.35 W m-2 of internal heat, ~500 yr per kelvin.  So the
    column starts from the emulator's equilibrium for the same forcing, composition,
    inhibition switch and opacity law (constant kappa in graybody mode).
    """
    from ..column import ColumnParams, solve_column
    from ..composition import x_to_q

    col = solve_column(ColumnParams(
        teq=s.teq, metallicity=s.metallicity, tint=s.tint, gravity=s.planet.gravity,
        kappa_th=s.kappa_ir, gamma_v=s.kappa_vi / s.kappa_ir,
        kappa_pexp=0.0 if s.greybody else 1.0, inhibition=s.inhibition,
        p_bot=max(1.0e8, 1.5 * s.p_bot),
    ))
    p_lay = s.p_bot * 10.0 ** (-pseudo_altitudes(s))  # sigma = exp(-z/H), H = 1/ln 10
    lp, lpc = np.log(p_lay), np.log(col.p)
    t_lay = np.interp(lp, lpc, col.T)
    t_bot = float(np.interp(np.log(s.p_bot), lpc, col.T))
    q = np.asarray(x_to_q(col.x, col.comp.mu_dry), float)
    q_lay = np.minimum(np.interp(lp, lpc, q), s.qvap_deep) if s.qvap_deep > 0 else np.zeros_like(p_lay)
    return p_lay, t_lay, t_bot, q_lay


def render_profile(s: RunSpec) -> str:
    """rcm1d ``profile`` for ichoice = 8: the bottom (surface) temperature, then one per layer."""
    _, t, t_bot, _ = initial_column(s)
    return "\n".join(f"{v:.4f}" for v in (t_bot, *t)) + "\n"


def render_profile_vap(s: RunSpec) -> str:
    """rcm1d ``profile_<tracer>``: the surface reservoir (kg m-2; none here), then q per layer."""
    *_, q = initial_column(s)
    return "0.0\n" + "\n".join(f"{v:.6e}" for v in q) + "\n"


def render_newstart_input(s: RunSpec) -> str:
    """Answers for ``newstart`` (option 0, from a start_archive.nc), one per line.

    The lon-lat core cannot cold-start generic physics: with ``read_start=.false.``
    ``iniacademic`` only fills fields when ``iflag_phys >= 2``, so a 3-D run with no
    start.nc stops in ``dynetat0``.  Option 1 (from start/startfi) does no regridding,
    so any seed on another grid has to go through ``start2archive`` and option 0.
    The prompt order follows newstart.F and tabfi_mod.F90 at r4609.  The state is
    flat, at rest, horizontally uniform at the emulator column (``profile.in`` and
    ``profile_<vap>``, written alongside), with p_s = p_bot everywhere.
    """
    p = s.planet
    rcp = 8.314463 / (s.mugaz_ref / 1.0e3) / s.cppd_ref
    consts = [
        ("obliquit", "0."), ("peri_day", "0."),
        ("periastr", f"{p.semi_major_axis_au:.5f}"), ("apoastr", f"{p.semi_major_axis_au:.5f}"),
        ("year_day", f"{p.rotation_period_days:.4f}"), ("rad", f"{p.radius:.1f}"),
        ("omeg", f"{p.omega:.6e}"), ("g", f"{p.gravity:.4f}"),
        ("mugaz_ref", f"{s.mugaz_ref:.4f}"), ("rcp_ref", f"{rcp:.6f}"), ("daysec", "86400."),
    ]
    lines = ["0",                                  # from start_archive.nc
             "y", f"{s.p_bot:.1f}", f"{s.p_bot / 200.0:.1f}"]  # preff, pa
    for k, v in consts:
        lines += [k, v]
    lines += ["",                                  # end of tab_cntrl changes
              "none", "0.1", "1000.",              # no surface file; albedo, inertia (nosurf: unused)
              "0"]                                 # lect_start_archive: which stored date (the seed has day 0)
    lines += ["flat", "y", f"{s.p_bot:.1f}",       # no topography, uniform p_s
              "t=profile", "q=0"]
    if s.water_scheme != "dry":
        lines += ["q=profile", "1"]                # tracer 1 is the vapour (traceur.def order)
    lines += [""]                                  # end of changes: write restart.nc, restartfi.nc
    return "\n".join(lines) + "\n"


def render_rcm1d(s: RunSpec) -> str:
    return f"""# rcm1d.def generated by meromixis for run {s.run_id}  (keys from deftank/earth1d.rcm1d.def; [verify] the rest)
daysec     = 86400.
g          = {s.planet.gravity:.4f}
year_day   = {s.planet.rotation_period_days:.4f}
periastr   = {s.planet.semi_major_axis_au:.5f}
apoastr    = {s.planet.semi_major_axis_au:.5f}
peri_day   = 0.
obliquit   = 0.
psurf      = {s.p_bot:.1f}
day0       = 0
time       = 12
day_step   = 48
ndt        = {s.nday * 48}
diagfi_output_rate = {48 * 50}
# 8: read T from the file "profile" (emulator equilibrium; see initial_column)
ichoice    = 8
tref       = {max(s.teq, 200.0) * 1.5:.1f}
isin       = 0
rad        = {s.planet.radius:.1f}             # [verify]
"""


def render_planet_constants(s: RunSpec) -> str:
    p = s.planet
    return f"""Planet constants for run {s.run_id}: enter these in newstart (they live in startfi.nc / start.nc)

radius            {p.radius:.1f} m
gravity           {p.gravity:.4f} m s-2
rotation period   {p.rotation_period_days:.4f} d   (omega = {p.omega:.6e} rad s-1), synchronous
day length        86400 s as the model "day" unit            [verify: convention for tlocked runs]
year              {p.rotation_period_days:.4f} days
periastr=apoastr  {p.semi_major_axis_au:.5f} AU  ({p.semi_major_axis_au * AU / 1e9:.3f} million km)
obliquity         0
bottom pressure   {s.p_bot:.1f} Pa
mugaz             {s.mugaz_ref:.4f} g mol-1      cp {s.cppd_ref:.1f} J kg-1 K-1    (must agree with callphys.def)
deep water        q = {s.qvap_deep:.6f} kg/kg   (x = {s.composition.x_h2o_deep:.5f} mol/mol)
initial state     {'restart of ' + s.branch_from if s.branch_from else 'T(p), q(p) from the matching rcm1d run'}
build             {s.build_command()}
"""


def render_slurm(s: RunSpec, account: str = "def-ncowan", walltime: str = "24:00:00", ntasks: int = 24,
                 cpus_per_task: int = 1, exe: str = "gcm.e") -> str:
    # One OpenMP thread per rank: generic_cloud_common_h reads table_tracers_condensable
    # on the fixed unit 117 from every thread, and two threads race on it.
    if s.one_d:
        run = "./rcm1d.e > lrun${next} 2>&1"
        ntasks, cpus_per_task = 1, 1
    else:
        run = f"srun ./{exe} > lrun${{next}} 2>&1"
    return f"""#!/bin/bash
#SBATCH --account={account}
#SBATCH --job-name={s.run_id}
#SBATCH --time={walltime}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem-per-cpu=2G
#SBATCH --output=logs/%x-%j.out
# Self-chaining Generic PCM job for DRAC Fir, generated by meromixis.
# One submission = one segment of {s.nday} days; it resubmits itself until num_run == NUMMAX.
set -euo pipefail
NUMMAX={s.segments}

module load StdEnv/2023 gcc openmpi netcdf-fortran-mpi    # check with `module spider`
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=400M
ulimit -s unlimited

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs segments
num=$(cat num_run 2>/dev/null || echo 0)
if [ "$num" -ge "$NUMMAX" ]; then echo "{s.run_id}: all $NUMMAX segments done"; exit 0; fi
next=$((num + 1))

{run}
if ! grep -q "Everything is cool" "lrun${{next}}"; then
    echo "{s.run_id}: segment ${{next}} did not finish cleanly; not resubmitting" >&2
    exit 1
fi

cp restart.nc   "segments/restart${{next}}.nc"
cp restartfi.nc "segments/restartfi${{next}}.nc"
mv diagfi.nc    "segments/diagfi${{next}}.nc"
mv restart.nc start.nc
mv restartfi.nc startfi.nc
echo "$next" > num_run

# spin-up log: energy balance and deep temperature, one line per segment
if command -v meromixis >/dev/null 2>&1; then
    meromixis monitor "segments/diagfi${{next}}.nc" --gravity {s.planet.gravity:.4f} --metallicity {s.metallicity:g} \\
        --intheat {s.intheat:.4f} >> spinup.log || true
fi

if [ "$next" -lt "$NUMMAX" ]; then sbatch job.slurm; fi
"""


def manifest(s: RunSpec) -> dict:
    d = asdict(s)
    d["derived"] = {
        "stellar_flux_W_m2": s.stellar_flux,
        "Fat1AU": s.fat1au,
        "intheat_W_m2": s.intheat,
        "gravity": s.planet.gravity,
        "mugaz_ref": s.mugaz_ref,
        "cppd_ref": s.cppd_ref,
        "cppv_ref": s.cppv_ref,
        "qvap_deep": s.qvap_deep,
        "x_h2o_deep": s.composition.x_h2o_deep,
        "build": s.build_command(),
    }
    return d


def write_run_dir(s: RunSpec, root, account: str = "def-ncowan") -> Path:
    """Create ``root/<run_id>/`` with every input file.  Existing files are overwritten."""
    from .. import __version__

    d = Path(root) / s.run_id
    d.mkdir(parents=True, exist_ok=True)
    files = {
        "callphys.def": comments_on_own_lines(render_callphys(s)),
        "gases.def": render_gases(s),
        "traceur.def": render_traceur(s),
        "z2sig.def": render_z2sig(s),
        "diagfi.def": render_diagfi(s),
        "planet_constants.txt": render_planet_constants(s),
        "job.slurm": render_slurm(s, account=account),
    }
    if s.one_d:
        files["rcm1d.def"] = comments_on_own_lines(render_rcm1d(s))
        files["profile"] = render_profile(s)
        if s.water_scheme != "dry":
            # rcm1d reads profile_<name> for every *_vap tracer except legacy h2o_vap,
            # which has its own branch reading profile_h2o_vap: same file name either way.
            files[f"profile_{s.tracer_names[0]}"] = render_profile_vap(s)
        # rcm1d still reads run.def, and aborts unless run_1d is set there.
        files["run.def"] = (
            "# 1-D run: rcm1d reads this as well as rcm1d.def\n"
            "run_1d=.true.\n"
            "restart=.false.\n"
            "INCLUDEDEF=callphys.def\n"
            "INCLUDEDEF=rcm1d.def\n"
        )
    else:
        files["run.def"] = comments_on_own_lines(render_run_def(s))
        if not s.branch_from:
            # inputs for newstart (see render_newstart_input); same file formats as rcm1d's
            files["profile.in"] = render_profile(s)
            files["newstart.in"] = render_newstart_input(s)
            if s.water_scheme != "dry":
                files[f"profile_{s.tracer_names[0]}"] = render_profile_vap(s)
    for name, text in files.items():
        (d / name).write_text(text)
    m = manifest(s)
    m["meromixis_version"] = __version__
    (d / "manifest.json").write_text(json.dumps(m, indent=2, default=str) + "\n")
    return d


def one_d_twin(s: RunSpec) -> RunSpec:
    """The rcm1d baseline that accompanies a 3-D run (same forcing and composition)."""
    return replace(s, run_id=s.run_id + "_1d", one_d=True, tier="baseline1d", nday=3000, segments=1, branch_from="")
