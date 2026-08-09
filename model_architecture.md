# Model architecture (current version)

This describes the model as implemented on `main` — the C++ engine in
`src/micropoint_new2.cpp` (types in `src/micropointheaders3.h`) and
its R interface (`R/NewMicroPoint.R`, `R/Ectotherm.R`,
`R/PlantGeometry.R`).

## Overview

The core is a multilayer canopy–soil surface energy balance model.
Given hourly above-canopy weather forcing, vegetation structure, and
soil properties, it solves for temperature, humidity, wind speed and
radiation at any height above or below ground, and at any time step,
by iterating the coupled energy and water balance equations to
convergence at each step. Soil heat and moisture are tracked as a
column of layered state that persists between time steps, so a run
effectively integrates forward through the supplied time series.

An ectotherm sub-model consumes the microclimate output and solves a
separate steady-state energy balance for an animal's body temperature
at each requested point.

R functions are thin wrappers: they build/validate the C++ struct
inputs (`vegpstruct`, `soilpstruct`, `climstruct`, ... in
`micropointheaders3.h`), call an `// [[Rcpp::export]]` entry point,
and reshape the returned `Rcpp::List`/`DataFrame` into the shape
documented for the R function.

## C++ engine sections (in file order, `micropoint_new2.cpp`)

1. **Solar model** — Julian day, solar position (zenith/azimuth) for a
   given lat/long/time.
2. **Radiation model** — Two-stream shortwave (direct + diffuse, split
   by sunlit/shaded leaf fraction) and longwave radiation exchange
   through canopy layers; ground and canopy absorption.
3. **Wind model** — Zero-plane displacement, roughness length,
   Monin–Obukhov stability corrections (`dpsim`/`dpsih`/`dphih`, smooth
   tanh-tapered stable branch as of 2026-08-09, not a hard clamp), and
   the within/above-canopy wind profile. Exported as
   `zeroplanedisCpp2`, `roughlengthCpp2`, `dpsimCpp2`, `dpsihCpp2`,
   `dphihCpp2` (the last two are internal-only in the current build;
   only the three `Ectotherm.R` needs are exported). `windmodelCpp`'s
   own internal `uf`/`L` iteration is Aitken-damped.
4. **Plant model** — Leaf energy balance: Penman-Monteith latent/
   sensible heat exchange, stomatal conductance from the Eller et al.
   (2020) hydraulic-optimisation scheme (uses `vegpstruct` fields like
   `Vcmax25`, `Kxmx`, `psi50`, `Dcrit`), leaf temperature. `leafgs`
   lazily caches `vegp.rpmin` (whole-plant hydraulic resistance) on
   first use, same pattern as `vegp.apsi`.
5. **Soil model** — Layered soil heat and Campbell-model water
   transport (`soilwatermod`/`soilmod`); the Campbell hydraulic
   conductivity exponent is always derived as `2*b + 3`, never stored.
   Solves a Thomas (tridiagonal) system per step, with Aitken-
   accelerated iteration (`WAitkenState`, and `SoilHeatCpp`'s own
   surface-temperature iterate) to convergence.
6. **Below-canopy Lagrangian model** — Propagates temperature/humidity/
   wind through canopy layers using a Lagrangian (near-field/far-field)
   dispersion approach, given source/sink strengths from the plant and
   soil models.
7. **Run model for one step** (`onestep`/`onestepbare` structs) —
   Couples radiation, wind, plant, soil, and canopy-dispersion pieces
   and iterates them to a converged single time step, for vegetated
   and bare-ground cases respectively. `OneStepBelow` Aitken-damps `H`
   (sensible heat) between outer passes before feeding it to the wind
   model — see "Convergence & performance" below.
8. **Above canopy** (`Tabove`, `RHabove`, `Uabove`) — Extrapolates
   temperature/humidity/wind from canopy top up to the reference
   height `zref`, assuming `zref` is above canopy height (see
   `CLAUDE.md`'s parked-issues note).
9. **Bigleaf model** — A simplified single-layer ("big leaf") version
   used only to spin up/initialize the soil heat and water state
   before a full multilayer run (`BigLeafCpp2`, `BigLeafBareCpp`,
   called from `weatherhgt_adjust`/`InitailizeWater`).
10. **Ectotherm model** — Separate steady-state body-heat-balance
    solver (`Ectotherm`/`EctothermM`), given ambient conditions and an
    animal geometry/physiology parameter list (`create_ectop()`
    output).
11. **R wrappers** — The step-8/9/10 machinery is looped over full
    time series / vertical profiles here: `profilebareR`, `profileR`,
    `RunBareR`, `RunModelR`, `RunBelowFullBare`, `RunBelowFull`,
    `WeatherhgtCpp2`. These are what the R-level `return_profile`,
    `RunMicro`, `RunModelFull`, `weatherhgt_adjust` functions actually
    call.
12. **Other useful functions** — `expand_outputCpp` (spline
    upsampling for plotting), `geometricCpp` (soil layer depth
    geometry), clear-sky radiation and solar altitude helpers.

## R-level API map

| Purpose | Function(s) | Defined in |
|---|---|---|
| Build vegetation params | `createvegp()` | `NewMicroPoint.R` |
| Build soil params | `createsoilc()` | `NewMicroPoint.R` |
| Foliage vertical profile | `PAIgeometry()` (forest), `PAIgrass()` (grass), `vegpforgrass()` | `PlantGeometry.R` |
| Adjust weather to reference height | `weatherhgt_adjust()` | `NewMicroPoint.R` |
| Spin up soil water state | `InitailizeWater()` | `NewMicroPoint.R` |
| CO2 ppm from year | `Cafromyear()` | `NewMicroPoint.R` |
| Single-hour vertical profile | `return_profile()` | `NewMicroPoint.R` |
| Time series at one height | `RunMicro()` | `NewMicroPoint.R` |
| Full time × height matrices | `RunModelFull()` | `NewMicroPoint.R` |
| Upsample output matrices | `expand_output()` | `NewMicroPoint.R` |
| Ectotherm parameters | `create_ectop()` | `Ectotherm.R` |
| Ectotherm profile / time series / full | `profile_ecto()`, `timeseries_ecto()`, `full_ecto()` | `Ectotherm.R` |

Every one of these — run against a full year of hourly `climdata`, for
both forest and grassland, canopy and bare ground, and the ectotherm
paths — was exercised end-to-end on 2026-08-09 with no errors or
warnings (see git log for the fix that was needed: exporting the
wind-profile Cpp2 functions `Ectotherm.R` depends on).

## Data flow (single time step, vegetated case)

```
weather (climdata row) + vegp + soilc
        |
        v
  solar position, radiation model  ---->  Rdirdown/Rdifdown/Rswup/Rlwdown/Rlwup
        |
        v
  wind model (d, zm, stability)   ---->  uz profile
        |
        v
  plant model (per canopy layer)  <---->  Lagrangian canopy dispersion
        |                                  (iterated to convergence)
        v
  soil model (heat + water, layered, persists across time steps)
        |
        v
  onestep result: tair/tleaf/rh/uz profiles, H, L, Et, Ev, soil state
```

The ectotherm model is a separate pass that takes an `onestep`-style
microclimate result (from `return_profile`/`RunMicro`/`RunModelFull`)
plus `ectop` (animal parameters) and solves the animal's body energy
balance independently — it does not feed back into the microclimate
solve.

## Convergence & performance

Full detail, methodology, and measurements: `performance_assessment.md`.
Short version:

- Benchmarked at **20 canopy layers**; runtime is superlinear in layer
  count (`LangrangianOne`'s pairwise-layer loop is O(n²)) — going from
  20→160 layers costs ~46x, not ~8x.
- Fixed: undamped oscillations in `windmodelCpp`'s `uf`/`L` and
  `OneStepBelow`'s `H`, `SoilHeatCpp`'s surface-temperature iterate —
  all now Aitken-damped, all verified non-regressive.
- Found but not fixed: some hours need a large, genuinely monotonic
  jump in `H` between consecutive hours (not an oscillation — damping
  can't accelerate it; needs real extrapolation or a smarter warm
  start).
- `plantmodelCpp`/`leafgs` micro-optimised (redundant `satvapCpp2` call,
  cached `rpmin`): ~29% faster end-to-end at n=20. Also fixed a real
  bug found in passing: the woody-vegetation evaporation formula had a
  misplaced bracket (changes model output for vegetation with a woody
  fraction).
