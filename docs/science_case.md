# Science case

Two results from one grid. This note records where each sits in the literature, the
hypotheses that can fail, and what gets written up in either outcome. The literature
survey was done on 12 November 2026; repeat the ADS queries in section 6 before submission.

## 1. Shared physics

In an H2-dominated gas water vapour is about eight times heavier than the background.
Three consequences follow from one water profile q(p):

1. **Cold trap.** The water that reaches the upper atmosphere is the minimum of the
   saturation mixing ratio along the column. `coldtrap.cold_trap` returns the level, its
   temperature and the transmission `x_above / x_deep`.
2. **Convective inhibition.** Where q exceeds `q_cri = 1/(ϖ γ)` moist convection shuts off
   (Guillot 1995; Leconte et al. 2017, arXiv:1610.05506). A stable layer forms between the
   cloud base and the level where q falls to q_cri. `inhibition.find_stable_layer`.
3. **Hotter interior and throttled exchange.** Heat must cross that layer radiatively, so
   the deep envelope is hotter; gas must cross it by whatever weak mixing survives.

The control parameters are the forcing (instellation and albedo) and the deep water
abundance (metallicity). In that plane there are three regimes and two boundaries:

| Regime | Cold trap | Lid | Emulator location (K2-18 b-like) |
|---|---|---|---|
| steam | no | no | effective Teq above ~245-290 K |
| cold-trapped, convecting | yes | no | cooler, Z below ~11-17× solar |
| cold-trapped and lidded | yes | yes | cooler, Z above ~11-17× solar |

Result I is the boundary between steam and the rest. Result II is a property of the lidded
regime. A lid without a cold trap cannot occur, which already gives a population-level
prediction: planets whose water is visible should not be affected by the barrier.

## 2. Literature map

**Convective inhibition.** Leconte et al. (2017) derive the criterion and the three-layer
structure. Leconte et al. (2024, A&A 686; arXiv:2401.06608) run the Generic PCM physics in a
convection-resolving box for K2-18 b: the layer forms, can be sub-saturated, has low eddy
diffusivity, and makes the interior hot enough that CO, CO2 and CH4 become detectable while
NH3 falls. arXiv:2409.18217 (ApJ 2025) finds the same suppression with CM1: resolved and
sub-grid K_zz drop by orders of magnitude in the layer. Cavalié et al. (2017,
arXiv:1703.04358) show for Uranus and Neptune that such a layer changes what deep abundances
one infers from quenched CO.

**3-D models of temperate sub-Neptunes.** Charnay et al. (2021, arXiv:2011.11553) run
K2-18 b in the LMD Generic GCM at several metallicities and find the coldest point of the
0.5-5 mbar region at the terminator. Innes & Pierrehumbert (2022, arXiv:2112.11108) give the
dry dynamics. arXiv:2502.12234 presents a new convection scheme for such GCMs and notes
that reduced vertical mixing from inhibition has been put forward to explain K2-18 b's
chemistry. arXiv:2511.07546 runs GCMs of a Hycean K2-18 b with different surface pressures
and albedos. Each of these treats one planet.

**The cold trap in the JWST sample.** K2-18 b shows CH4 and CO2 without H2O or NH3. Hu et
al. (2025, arXiv:2507.12622) turn the water upper limit into a cold-trap temperature below
175 K (1σ) or 215 K (2σ), which needs a substantial Bond albedo. 1-D models reach the same
tension: arXiv:2501.18477 needs an extra albedo source to cold-trap the water in a
mini-Neptune K2-18 b; arXiv:2603.19803 invokes apoastron clouds on the eccentric orbit;
arXiv:2407.09009 maps how condensation alters surface identification. The PNAS review
(arXiv:2509.19247) reads the rising H2O abundance from K2-18 b to TOI-270 d as evidence for
a critical temperature separating cold-trapped from steamy atmospheres. TOI-732 c
(arXiv:2512.15844), LP 791-18 c (arXiv:2512.10876) and TOI-1231 b (arXiv:2607.00973) now
populate the 300-360 K interval.

**Ammonia as a surface diagnostic.** Hu et al. (2021, arXiv:2108.04419), Yu et al. (2021)
and Tsai et al. (2021) propose missing NH3 as the mark of a shallow surface. Competing
explanations for K2-18 b's missing NH3 are an ocean, dissolution in a magma ocean
(discussed in arXiv:2501.18477), and the hotter interior under a stable layer (Leconte et
al. 2024). Hu et al. (2025) find that photochemistry alone, with standard mixing, cannot
remove NH3 from a massive atmosphere.

## 3. Gap I: the cold-trap transition in 3-D

*What exists:* an observational claim of a critical temperature; 1-D models that place the
trap but need ad hoc albedo or reduced irradiation for K2-18 b; 3-D models of single
planets.

*What is missing:* the location **and width** of the transition as a function of
metallicity from a 3-D model with non-dilute moist physics, and what the two limbs look
like while the planet crosses it.

*Hypotheses.*

* **H1a (offset).** The 3-D transition sits at higher forcing than the 1-D one computed
  from the same model's column mode, because the terminator and nightside trap water that
  the global-mean column does not. If true, 3-D relaxes the albedo that 1-D models need for
  K2-18 b. Falsified if `T★_3D − T★_1D` is zero within the node spacing.
* **H1b (width).** The transition is broader in 3-D: the substellar column opens first, the
  limbs last. Measured as the forcing interval between limb transmission 0.1 and 0.9.
* **H1c (asymmetry).** East-west limb contrast in water peaks inside the transition band
  and vanishes on both plateaus, giving an observable (ingress/egress asymmetry) that
  identifies planets in transition.
* **H1d (transport).** `transmission_actual / transmission_sat` departs from 1: advection
  and sub-saturation, not local temperature alone, set the water aloft.

*Deliverables:* T★(Z), width(Z), limb asymmetry(Teq, Z), each for 3-D and for `rcm1d`, with
the JWST targets over-plotted. Functions: `diagnostics3d.transmission_summary`,
`regimes.transition_teq`.

## 4. Gap II: is the lid a barrier?

*What exists:* the thermal route to low NH3 (hot interior, Leconte et al. 2024); evidence
that mixing collapses in the layer (convection-resolving boxes); a passing mention that
reduced mixing could matter for chemistry.

*What is missing:* (a) a quantitative supply-against-loss criterion for a photochemically
fragile gas; (b) the leak through the layer by the global circulation, which a local box
cannot contain; (c) a clean separation of the thermal and transport contributions;
(d) the implication for shallow-surface inference, with a discriminator.

*The kill test.* NH3 is hidden at the levels spectra probe when the series resistance from
the quench level up to that level satisfies `R ≥ f_deep / (η Φ_UV)`; equivalently when the
layer diffusivity is below `K_crit` (`barrier.k_layer_crit`). The GCM supplies the layer
geometry and a **lower bound** on its permeability: `diagnostics3d.layer_permeability` gives
the flux-gradient K of the resolved flow, which excludes plumes, waves and sub-grid mixing.

* If `K_resolved > K_crit` over the plausible range of η and Φ_UV, the circulation alone
  defeats the barrier. That is a clean negative result: the shallow-surface diagnostic is
  robust against convective inhibition, and only the thermal effect remains.
* If `K_resolved < K_crit`, the barrier survives the large-scale leak and the question
  passes to the small scales, where the convection-resolving numbers (0.1-10 m² s⁻¹ as
  quoted in the primer) are compared with `K_crit(Z)`.

*Scoping numbers* (emulator; hand-estimated UV flux of 10¹³ photons cm⁻² s⁻¹; see README):
`K_crit` falls from ~3 m² s⁻¹ at 15× solar to ~0.2 at 100× and below 10⁻³ at 300×. The
barrier is therefore a low-to-intermediate metallicity phenomenon: plausible if K2-18 b's
percent-level CH4 reflects an enrichment of a few tens times solar, implausible at several
hundred times solar. At 100× solar the thermal effect alone lowers quenched NH3 by a factor ~2 and
raises quenched CO by a factor ~20 in the emulator.

*Discriminator.* CO has no efficient sink above the layer and refills on 10²-10⁵ yr, so
"NH3 missing, CO present" separates a lidded deep envelope from a shallow surface. It does
not by itself separate the barrier from a magma ocean; nitrogen carriers that bypass NH3
(N2, HCN) are the next thing to compute and need a real photochemical model.

*Limits of this repository for Result II.* The sharp-front model replaces a photochemical
network with two numbers (η, Φ_UV). It decides whether the idea is alive; the paper needs
a kinetics code (for example VULCAN or EPACRIS) run on the GCM-derived K(p) and T(p).
`barrier.kzz_profile` and `Column.to_dataset` produce those inputs.

## 5. Outcome matrix and paper plan

| | Result I | Result II |
|---|---|---|
| Positive | 3-D shifts and broadens the transition; limb asymmetry flags planets in transition | A lidded window exists where NH3 is starved; NH3-without-CO-loss is a false positive for a surface |
| Negative | 3-D agrees with 1-D: population studies may keep using columns, with a measured error bar | The resolved circulation leaks faster than K_crit: missing NH3 cannot be blamed on the lid |

Either row is publishable because each cell answers a question the literature currently
answers by assumption.

**Paper I** (A&A or ApJ, full length): *the cold-trap transition of temperate sub-Neptunes
in three dimensions*. Figures: (1) regime plane with JWST targets; (2) T(p) and q(p) at
substellar, limbs, antistellar across the transition; (3) T★ and width against metallicity,
3-D and 1-D; (4) limb asymmetry against forcing; (5) transmission_actual / transmission_sat
maps; (6) sensitivity: clouds, resolution, rotation, internal heat; (7) synthetic spectra
through Pytmosph3R for three forcings.

**Paper II** (A&A Letter or ApJL): *can a stable layer hide ammonia?* Figures: (1) column
anatomy with and without inhibition; (2) lid area fraction and thickness over the grid;
(3) resolved K_zz(p) through the lid against K_crit; (4) kill-test plane (Z against K) with
the convection-resolving band; (5) NH3 and CO profiles from a kinetics code for the three
scenarios (deep, shallow, deep and lidded).

## 6. ADS queries to repeat before submission

```
abs:("cold trap" AND "sub-Neptune") year:2025-2027
abs:("convective inhibition" OR "convection inhibition" OR "moist convective inhibition") abs:(ammonia OR NH3) year:2024-2027
abs:("stable layer" AND "sub-Neptune" AND (quench OR "eddy diffusion" OR Kzz)) year:2024-2027
abs:("critical temperature" AND "sub-Neptune" AND water) year:2025-2027
citations(arXiv:2401.06608) abs:(GCM OR "general circulation")
citations(arXiv:2011.11553) year:2025-2027
```
