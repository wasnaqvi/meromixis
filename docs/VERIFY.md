# Verification checklist

Everything below was written without access to a compiled model or to the papers' full
text. Tick each item against your checkout or the source before relying on it.

## Generic PCM interface

| ID | Item | Where | How to check |
|---|---|---|---|
| P1 | Which thermodynamics the lon-lat core supports (`thermo_uni_ideal` only at r4146) | `callphys.def` | `changelog.txt` at your revision; ask the LMD team (gate G0) |
| P2 | Name of water as a generic condensable (`condensable_name: H2O`) and its entry in `datadir/table_tracers_condensable` (delta_vapH, Tref, Pref, mass). `generic_tracer_index` skips tracer names containing lowercase `h2o` or `co2` | `traceur.def`, `gases.def` | read `generic_cloud_common_h.F90` and the table |
| P3 | Emulator saturation curve matches the PCM's: set `Saturation("cc", latent_heat, t_ref, p_ref)` from the P2 table values when comparing | `thermo.py` | compare psat at 250, 300, 400 K |
| P4 | Output variable names (`temp`, `ps`, `p`, `w`, `ASR`, `OLR`, tracer names) and whether `aps`/`bps` are written | `pcm/io.py` `VARMAP` | `ncdump -h diagfi.nc` |
| P5 | `z2sig.def`: values are mid-layer pseudo-altitudes, count equals `nlev` | `pcm/deffiles.py` | compare with a `deftank` z2sig file and the first lines of the run log |
| P6 | Sponge keys in `run.def` (`iflag_top_bound`, `mode_top_bound`, `tau_top_bound`) | `run.def` | `LMDZ.COMMON/libf/dyn3d*/conf_gcm.F90` |
| P7 | `rcm1d.def` keys beyond the Earth example (`rad`), and `global1d`/`szangle` placement | `rcm1d.def`, `callphys.def` | `deftank/jupiter1d/` |
| P8 | Build flag `-p generic` and executable name | `job.slurm`, `planet_constants.txt` | `ls LMDZ.COMMON/bin` after building |
| P9 | Substellar longitude for `tlocked` runs (assumed 0°) and day-length convention | `diagnostics3d.stellar_mu`, `planet_constants.txt` | map of `ISR` from the first run |
| P10 | Sign and units of `w` (assumed m s⁻¹, positive up) | `diagnostics3d.kzz_flux_gradient` | variable attributes; `--w-is-omega` otherwise |
| P11 | Radiatively active generic clouds: `aerogeneric`, `is_rgcs`, `-t` | `clouds` variants | wiki page "Radiative Generic Condensable Specie" |

## Physics inputs

| ID | Item | Where | How to check |
|---|---|---|---|
| C1 | Quench timescales `t_CO = 1.5e-6 p⁻¹ m⁻⁰·⁷ exp(42000/T)` s and `t_NH3 = 1.0e-7 p⁻¹ exp(52000/T)` s (p in bar) | `chemistry.py` | Zahnle & Marley (2014), the equations for CO and NH3 |
| C2 | Solar abundances (Asplund, Amarsi & Grevesse 2021: O 8.69, C 8.46, N 7.83, He 10.914) | `constants.py` | the paper's Table 2 |
| C3 | Gibbs-energy fits and cp tables | `constants.py`, `chemistry.py` | `python scripts/dev/make_thermo_tables.py`; the test suite pins Kp(1000 K) and Kp(700 K) |
| C4 | NH3 UV photon flux (10¹³ cm⁻² s⁻¹), cross-section, irreversible fraction η, UV floor pressure | `barrier.Photolysis` | integrate a MUSCLES spectrum with `photon_flux_from_spectrum`; η from a kinetics code |
| C5 | NH3 quench level falls below the 1000 bar model bottom in the fiducial column; the bottom equilibrium value is used | `chemistry.quench` (`in_domain`) | extend `p_bot` with a real-gas EOS, or take the deep adiabat from an interior model |

## Literature numbers

| ID | Item | Where |
|---|---|---|
| L1 | K_zz = 0.1-10 m² s⁻¹ in the stable layer of Leconte et al. (2024), taken from the project primer | `barrier.Mixing`, kill-test figure |
| L2 | Target parameters flagged `verify: true` | `experiments/targets.yaml` |
| L3 | Novelty: repeat the ADS queries | `docs/science_case.md` section 6 |

## Emulator

| ID | Item |
|---|---|
| E1 | Grey opacities (`kappa_th`, `kappa_pexp`, `gamma_v`) are uncalibrated defaults. Calibrate against `rcm1d` before quoting any emulator number |
| E2 | Two inhibited zones appear in the fiducial column (0.13-0.26 bar and 4-7.5 bar). Check with `rcm1d` whether the upper one is real: it dominates the series resistance at 100× solar |
| E3 | Ideal gas and constant-composition dry gas to 1000 bar; no cloud radiative effect; no CO2 or N2 in the mass budget |
