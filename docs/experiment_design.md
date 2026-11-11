# Experiment design

The sweep is defined in `experiments/steam_rising.yaml` and expanded by
`meromixis pcm-init`. This note explains the choices and fixes the acceptance criteria in
advance, so that "converged" and "the lid exists" are measurements.

## 1. Fixed choices

| Choice | Value | Reason |
|---|---|---|
| Planet | K2-18 b-like: 2.61 R⊕, 8.63 M⊕, g = 12.4 m s⁻² | continuity with Charnay+2021, Leconte+2024 |
| Rotation | synchronous, 32.94 d, **fixed across the sweep** | resolves primer open item 2: only `Fat1AU` changes, so thermodynamics is not confounded with rotation; the Kepler-consistent case is a variant |
| Forcing label | `teq`, with F = 4σ teq⁴ | the model computes its own albedo; realised Teq is an output |
| Grid | 64×48×40, 80 bar to 0.2 mbar, even in log p | Charnay+2021 |
| Deep water | `qvap_deep` fixed at the bottom level | the "deep composition" knob |
| Internal heat | Tint = 50 K (`intheat` = 0.35 W m⁻²) | varied in one variant; it controls how convective the cloud region is |
| Water scheme | generic condensable species | the only route with `moist_convection_inhibition` and `virtual_theta_correction` |

## 2. Gates (do these in order; each can stop the campaign cheaply)

**G0, desk check.** Ask the Generic PCM team what the lon-lat core supports at the head of
trunk. At r4146 the changelog states that `thermo_binary_ideal` works only with DYNAMICO or
in 1-D. Decide the core for Result II on that answer. Result I does not depend on it: the
cold trap sits where water is dilute.

**G1, column (`*_1d` twins, `rcm1d`).** With the generic route at 100× solar and
teq = 250-300 K: (a) a three-layer structure appears with inhibition on and not with it
off; (b) `inhibition.find_stable_layer` on the output finds it; (c) the deep temperature
rises with inhibition on. Then calibrate the emulator (`column.calibrate_grey`) so later
1-D against 3-D comparisons use the PCM's own column as the baseline. This also closes
primer open item 3 in practice: it fails immediately if the corrk table's water axis stops
short of `x_h2o_deep` (0.089 at 100×, 0.30 at 300×).

**G2, dry 3-D (`T300_Z100_dry`).** Reproduce Innes & Pierrehumbert (2022): weak horizontal
temperature gradients, day-night overturning, superrotation aloft, kinetic energy
equilibrating in ~10⁴ days.

**G3, lid survival (`T300_Z100_grey`, then `T275_Z100_fid` for 2000 days).** Initialise
from the G1 column. Pass if after 2000 days `diagnostics3d.inhibition_fields` still reports
a lid over most of the area and `ledoux_n2` is positive inside it. Fail if resolved
overturning has erased it: the dynamics then do not feel the composition gradient, and
Result II needs the other core. Run the emulator twin (`meromixis column`) alongside.

## 3. Tiers

| Tier | Runs | Purpose |
|---|---|---|
| spine | Z = 100×, ten forcings 250-600 K, dense at 275-325 K | Result I backbone |
| ribs | Z = 3, 30, 300× at 250-350 K | metallicity dependence; 3× is the no-lid control for Result II |
| variants | no inhibition, legacy water, radiative clouds, mass redistribution, 128×96, 80 levels, Tint = 80 K | attribution and robustness |

Re-centre the forcing nodes once G1 has located the 1-D transition:
`meromixis nodes --t-star <K> --lo 250 --hi 600 -n 10`.

Each run branches from the restart of its nearest neighbour (`branch_from` in
`runs.tsv`); only the fiducial run starts from the column profile.

## 4. Convergence criteria

A run enters the analysis when, over the last 2000 days:

* global TOA imbalance `|ASR + intheat − OLR|` < 1% of ASR;
* global-mean temperature at 10 bar drifts by < 0.1 K per 100 days;
* terminator `log10_x_obs` changes by < 0.05 dex between the two halves of the window;
* (Result II) lid top and base pressures change by < 5%.

`job.slurm` appends `meromixis monitor` output to `spinup.log` after every segment, so
these are plots, not guesses. Deep layers at 80 bar may not meet the temperature criterion
in affordable time; report the residual drift rather than hide it, and show that the
diagnostics above 10 bar are insensitive to it (the `tint80` variant brackets this).

## 5. Output needed

`diagfi.def` requests `temp ps p u v w ASR OLR ISR` and the two water tracers. `w` is
required for the flux-gradient K_zz. For that estimate sample at least every 5 days over
the analysis window (`output_every_days`), since it is a covariance.

## 6. Analysis: figure-to-function map

| Quantity | Function |
|---|---|
| column cold trap, saturation-limited and actual transmission | `diagnostics3d.coldtrap_fields` |
| region table (substellar, limbs, night) and limb asymmetry | `diagnostics3d.transmission_summary` |
| transition location and width from a set of runs | `regimes.transition_from_samples` on terminator `transmission_actual` against realised Teq (same definition as the emulator's `transition_teq`) |
| lid mask, top and base, area fraction | `diagnostics3d.inhibition_fields` |
| stratification seen by the dynamics | `diagnostics3d.ledoux_n2` |
| resolved-flow K_zz(p), harmonic mean across the lid | `diagnostics3d.kzz_flux_gradient`, `layer_permeability` |
| K_crit for the same column | `barrier.k_layer_crit` on the global-mean profile |
| thermal against transport split | `barrier.ammonia_budget` |
| spectra | Pytmosph3R on `diagfi`, then a 1-D retrieval; the gap between input and retrieved abundance is itself a result |

## 7. A better permeability measurement (optional patch)

The flux-gradient K of water vapour is available from standard output but is contaminated
by condensation near the cloud base. The cleaner measurement is a pair of passive tracers
with idealised chemistry: one relaxed to 1 below the lid and destroyed with a fixed
lifetime above it, one the reverse. This needs a few lines in the physics (a relaxation
tendency on two extra tracers) and is not written here because it must be coded against
your checkout. `barrier.solve_diffusion_loss` is the 1-D analogue to validate it against,
and `kzz_flux_gradient(ds, tracer="<name>")` already accepts any tracer.

## 8. Cost

Not estimated: no timing exists for this configuration on Fir. Day 4 of the primer's plan
(scaling test on 1, 2 and 4 nodes) gives simulated days per core-hour; multiply by
`segments × nday` per run in `runs.tsv`. Budget the spine first.
