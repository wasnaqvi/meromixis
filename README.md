# meromixis

**Cold traps and chemoclines in temperate sub-Neptune atmospheres: one Generic PCM sweep, two results.**

A *meromictic* lake never overturns. A light upper layer (the mixolimnion) rides on a dense,
solute-rich deep layer (the monimolimnion), and the thin *chemocline* between them lets
almost nothing through. A water-rich hydrogen atmosphere can arrange itself the same way
under its cloud deck, because in H2 the humid gas is the heavy gas. This repository holds
the shared machinery for two questions about that structure:

| | Question | Deliverable |
|---|---|---|
| **I. The cold-trap transition** | At what forcing does water get past the cold trap, how sharp is the switch, and does 3-D circulation move it? | T★(Z) and its width in 3-D against 1-D, with limb-resolved water |
| **II. The ammonia barrier** | Is the stable layer under the clouds tight enough to starve the upper atmosphere of NH3, so that a deep envelope imitates a shallow surface? | The layer's resolved-flow permeability against the critical value K_crit(Z), with CO as the discriminator |

Both are read off the same (forcing × metallicity) grid of Generic PCM runs with the same
column diagnostics: I lives on the boundary of the cold-trapped region, II lives inside the
part of it where `q_deep > q_cri`.

![Emulator scoping map](docs/figures/emulator_regime_map.png)

*Scoping map from the built-in 1-D grey emulator. It is a design aid for placing GCM runs,
not a result.*

## Status

| Component | State |
|---|---|
| Thermodynamics, cold-trap and stable-layer diagnostics, 1-D emulator | implemented, 46 tests passing |
| Barrier kill test (supply against photolysis, thermal against transport split) | implemented and tested; UV flux and quench timescales are placeholders to verify |
| 3-D diagnostics (limb-resolved cold trap, lid mask, flux-gradient K_zz) | implemented, tested on synthetic PCM-format output only |
| Generic PCM run-directory generator and Fir job chain | implemented; keys checked against trunk r4146 where tagged, **no run has been executed** |
| GCM results | none yet |

Read [`docs/VERIFY.md`](docs/VERIFY.md) before the first run. The most important open item:
as of trunk r4146 (March 2026) the variable-composition thermodynamics
(`thermo_binary_ideal`) is available only with DYNAMICO or in 1-D, so the lon-lat core
carries the lid only through `virtual_theta_correction` and `moist_convection_inhibition`
in the physics. Gate G1 in [`docs/experiment_design.md`](docs/experiment_design.md) tests
whether that is enough for Result II.

## Quick start

```bash
git clone <your-remote>/meromixis.git && cd meromixis
pip install -e ".[dev]"
pytest                                   # 46 tests, a few seconds

meromixis column --teq 255 --metallicity 100       # one emulator column
meromixis transition --metallicity 3 30 100        # where the cold trap opens (emulator)
meromixis killtest --metallicity 30 --eta 1        # ammonia-barrier kill test (emulator)
python scripts/make_scoping_figures.py             # regenerate docs/figures/

meromixis pcm-init experiments/steam_rising.yaml runs/   # 36 run directories + rcm1d twins
meromixis synth --teq 270 --out diagfi_synthetic.nc      # fake PCM output ...
meromixis diag diagfi_synthetic.nc --gravity 12.43 --metallicity 100   # ... through the real pipeline
```

On real output, `meromixis diag segments/diagfi*.nc --gravity ... --metallicity ... --last-n 40`
prints the region-by-region cold-trap table (Result I) and the lid's area fraction and
resolved-flow K_zz (Result II).

## Layout

```
src/meromixis/
  composition.py   metallicity -> mole/mass fractions, dry-gas mu and cp
  thermo.py        water saturation (Murphy-Koop + Wagner-Pruss, or the PCM's constant-L form), adiabats
  inhibition.py    q_cri, Ledoux N^2, stable-layer finder
  coldtrap.py      rain-out profile, cold-trap level and transmission (1-D and n-D)
  column.py        1-D grey radiative-convective emulator with inhibited layers; calibrate_grey()
  regimes.py       steam / cold-trapped / lidded classification, transition finder
  chemistry.py     NH3 and CO equilibrium and quenching
  barrier.py       K(p), series resistance, K_crit, sharp-front profile, thermal/transport budget
  diagnostics3d.py the same diagnostics on (time, lev, lat, lon) output; flux-gradient K_zz
  synthetic.py     PCM-format synthetic output for testing the pipeline
  pcm/             .def renderers, sweep expansion and branching plan, Slurm chain, output I/O
experiments/       steam_rising.yaml (the sweep), targets.yaml (JWST sample for over-plotting)
scripts/           scoping figures; dev/ regenerates the thermochemical tables
docs/              science case, experiment design, Fir runbook, verification checklist
tests/             unit and end-to-end tests
```

## What the emulator already says (and what it cannot)

Run on a K2-18 b-like planet with default grey opacities:

* The cold trap opens over a ~20-25 K interval of effective equilibrium temperature, centred
  near 245 K at 3× solar and 282 K at 100× solar. That sits between K2-18 b and TOI-270 d,
  where the population papers place an observational "critical temperature", and **below
  the 300-600 K range of the original sweep plan**: with zero albedo every one of those
  points is steam in 1-D. The sweep in `experiments/steam_rising.yaml` therefore starts at
  250 K and clusters around 275-325 K.
* The stable layer exists above ~11-17× solar and only on the cold-trapped side.
* The barrier is **marginal**. With the hand-estimated UV flux of the explainer, K_crit is
  ~3 m² s⁻¹ at 15× solar, ~1.5 at 30×, ~0.2 at 100× and <10⁻³ at 300×, against the
  0.1-10 m² s⁻¹ quoted for convection-resolving simulations. It is lower than the explainer's
  back-of-envelope 2×10⁵ cm² s⁻¹ (20 m² s⁻¹) because the self-consistent layer is deeper
  (larger n) and quenched NH3 is more abundant than 10⁻⁴. If only 10% of photolysis is
  irreversible, divide by ten. The idea survives in a window of roughly 12-100× solar for
  fully irreversible photolysis, narrowing to below ~40× at 10%.
* CO behaves as hoped: it is enhanced by the hotter interior and refills the upper
  atmosphere in 10²-10⁵ yr for layer diffusivities between 1 and 10⁻³ m² s⁻¹.

None of this is quotable. Radiation is double-grey and uncalibrated, clouds do nothing,
the gas is ideal to 1000 bar, and two chemical inputs are unverified. The emulator exists
to place runs, test diagnostics and provide a 1-D baseline once it is calibrated against
`rcm1d` (`meromixis.column.calibrate_grey`).

## Documents

* [`docs/science_case.md`](docs/science_case.md): literature map, the two gaps, hypotheses, outcome matrix, paper plan
* [`docs/experiment_design.md`](docs/experiment_design.md): gates, tiers, convergence criteria, figure-to-function map
* [`docs/pcm_on_fir.md`](docs/pcm_on_fir.md): build and run on DRAC Fir
* [`docs/VERIFY.md`](docs/VERIFY.md): everything written from memory or not yet checked against a run

## Citing

See `CITATION.cff`. The Generic PCM itself is developed at LMD and collaborators and is not
distributed here; cite its papers and Charnay et al. (2021), Leconte et al. (2017, 2024)
for the physics this project leans on.

## License

MIT. See `LICENSE`.
