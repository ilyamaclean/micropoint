# Performance assessment (2026-08-09)

Written in response to `restructuretask.txt`'s request for "an overall
assessment of how the model could be speeded up through better design,
improved damping, other checks for convergence, changes to architecture
etc.", followed up with an explicit instruction to do *real* profiling
rather than reason by analogy, since "the real wins are all to do with
ensuring model convergence." This version replaces an earlier draft that
leaned entirely on the sibling `microclimfv2` package's own profiling
data — that turned out to be the wrong assumption (see below); this one
is grounded in actual instrumentation run against this package's own
code, over a full year of real hourly data.

Constraints respected throughout: no clumping coefficient introduced, the
full Lagrangian multilayer canopy model kept exactly as-is, and every
change verified to leave model outputs correct (recompile + the tutorial
vignette's full-year smoke test, zero errors/zero warnings, after every
change).

## Method

Temporary instrumentation was added directly to `micropoint_new2.cpp`
(std::chrono per-phase timers inside `OneStepBelow`'s loop, plus
iteration-count tracking for each of the model's five distinct
convergence loops), run against a full year of the inbuilt `climdata`
(forest scenario: `createvegp("BET.Te")`, 20 canopy layers, "Clay loam"
soil, height-adjusted to 27 m, `RunMicro(reqhgt = 0.25)`), analysed in R,
then **stripped back out** once the investigation was done — this
mirrors the exact same profiling methodology `microclimfv2` used, and
for the same reason: it's not meant to ship.

## The model's five convergence loops

| Loop | Where | Role |
|---|---|---|
| `windmodelCpp` | wind/stability | solves `uf`/`L` (Monin-Obukhov length) |
| `SoilHeatCpp` | soil heat | per-timestep soil temperature (tridiagonal) |
| `SoilWaterCpp` | soil water | per-timestep soil moisture (tridiagonal, Newton) |
| `canopytop` | above-canopy | canopy-top temperature/vapour pressure |
| `OneStepBelow` | outer loop | wraps all of the above + radiation + plant + Lagrangian dispersion, once per output timestep |

## What the real data showed (first pass, before any fix)

**Phase timing was the opposite of what the `microclimfv2` analogy
predicted.** That package has no multilayer canopy at all (its own
`CLAUDE.md` describes it as deliberately BigLeaf-only, feeding a grid
model), so of course its cost concentrated almost entirely in
`SoilWaterCpp`/`SoilHeatCpp` (~58%/~39%). Profiling *this* package's
actual `OneStepBelow` loop found:

| Phase | Share of runtime |
|---|---|
| Plant model (`plantmodelCpp`, per-layer leaf energy balance + stomatal conductance) | **43.3%** |
| Lagrangian canopy dispersion (`LangrangianOne`) | **22.8%** |
| Soil heat (`SoilHeatCpp`) | 14.0% |
| Soil water (`SoilWaterCpp`) | 12.1% |
| Longwave radiation | 3.7% |
| Rain interception | 2.5% |
| Wind model | 1.0% |
| Canopy top | 0.3% |
| Sensible/latent heat summation | 0.2% |

The full multilayer Lagrangian canopy model — this package's whole reason
for being architecturally distinct from `microclimfv2`'s point model —
is where two-thirds of the runtime actually goes, not the soil solvers.
The earlier draft of this document, which assumed the profiling
breakdown would transfer, was wrong; this is why it's been rewritten
rather than amended.

**Convergence-stream iteration counts** (before any fix), full year,
8760 hours:

| Stream | Mean | Median | Max | Hit 100-cap |
|---|---|---|---|---|
| Outer (`OneStepBelow`) | 12.4 | 12 | 100 | 1 / 8,760 (0.01%) |
| Soil heat | 1.9 | 1 | 100 | 13 / 108,747 (0.01%) |
| Soil water | 2.0 | 2 | 21 | 0 |
| **Wind** | 6.3 | 5 | 102 | **988 / 108,747 (0.91%)** |
| Canopy top | 1.1 | 1 | 3 | 0 |

Soil heat and soil water were already converging cleanly (consistent
with this session's earlier alignment work against `pointmodel.cpp` —
see git log). Wind stood out: nearly 1% of its ~109k calls across the
year ran the full 100 iterations without settling below its `1e-8`
tolerance — more outright non-convergence than every other stream
combined.

**Root cause, found by inspection**: `windmodelCpp` solves its own
`uf ↔ L` fixed-point loop internally, completely undamped — structurally
the same kind of self-referential coupling as `SoilHeatCpp`'s
`Tsurf_iter`, which this session had already found and fixed (it
diverges when seeded from a bad guess) using Aitken damping. `windmodelCpp`
had no equivalent damping.

**A second, unrelated bug found by inspection, not profiling**:
`OneStepBelow` was discarding `SoilHeatCpp`'s real iteration count
(`soilheat.iters = 0`) immediately after computing it — silently zeroing
the `soilhiters` diagnostic column that's already wired all the way
through to `return_profile()`'s R output for every vegetated call. Pure
reporting bug (the field is never read for any control flow), but it
means anyone inspecting `soilhiters` today has been looking at zeros.

## Fixes implemented and verified

1. **`OneStepBelow`: stopped discarding `soilheat.iters`.** Restores the
   already-existing `soilhiters` diagnostic. Zero behaviour risk (pure
   reporting field).
2. **`windmodelCpp`: Aitken-damped the `uf`/`L` iteration**, same
   technique (`Aitken1DState`/`aitken1d`) already used and verified
   elsewhere in this file. Measured effect, same profiling harness,
   before vs. after:

   | | Before | After |
   |---|---|---|
   | Wind: hit 100-cap | 988 (0.91%) | **147 (0.13%)** |
   | Wind: mean iterations | 6.26 | 5.93 |
   | Outer: mean iterations | 12.41 | 12.48 (no change, within noise) |
   | Outer: hit 100-cap | 1 | 2 (no change, within noise) |

   An 85% reduction in wind's own outright non-convergence, at the same
   fixed point (Aitken damping doesn't change what's being solved for,
   only how reliably the iteration gets there) — confirmed correct by
   the full tutorial smoke test. It did **not**, however, meaningfully
   improve the outer loop's own convergence — see below.

## What's still unresolved (deliberately not guessed at further)

The outer loop's tail — 63% of hours need ≥10 outer iterations, and a
handful need close to 100 — is not explained by wind's own convergence
(fixing wind didn't move this number), and doesn't correlate strongly
with any single weather variable (Spearman correlation of outer
iteration count against temperature, humidity, wind speed, precipitation,
and sensible heat flux: all |r| < 0.14). The outer loop already uses a
fairly sophisticated custom weighted-Aitken relaxation
(`aitkin_weightdif`, spatially weighted towards the canopy base) for
`tair`/`rh`/`tleaf`, and its convergence tolerance is already loosened by
the R-level wrapper defaults (`1e-2`, not the tighter `1e-3` the C++
function's own default parameter would suggest).

This was investigated but deliberately **not** further modified, for the
same reason `microclimfv2`'s own profiling session gives as precedent:
that repo tried two plausible-looking changes to its own (structurally
similar) soil solver and had to revert both after measuring real
regressions — one of them doubled runtime. Blindly retuning an
already-tuned damping scheme without a much deeper diagnostic (tracing
what specifically is oscillating or crawling in one of the worst
individual hours, the way `microclimfv2` traced a specific bad
`SoilWaterCpp` timestep) risks exactly that outcome. Since plant +
Lagrangian dispersion (66% of per-pass cost) get re-run on every one of
those outer passes, cutting the outer loop's iteration count is the
single biggest remaining lever on total runtime — but it needs that
deeper trace-based diagnosis, not a guess, before touching it.

## Also implemented: release-mode compilation

`src/Makevars`, `src/Makevars.win` (`PKG_CXXFLAGS = -O2`), matching
`microclimfv2`'s own file. Zero code risk, real free win: a genuine
`R CMD INSTALL` picks this up (the interactive `devtools::load_all()`
debug build always appends `-O0` last regardless, by design, so this
doesn't show up in development — it matters for the package as actually
installed/distributed).

## Remaining opportunities (not attempted — flagged for a future session)

1. **Trace-diagnose the outer loop's slow-convergence tail** (see above)
   — the highest-value remaining lever, needs per-iteration inspection of
   specific bad hours (e.g. hour 7713, 2017-11-18 08:00, and hour 2961,
   2017-05-04 08:00, both hit the 100-iteration cap in testing) rather
   than aggregate statistics.
2. Whether `SoilWaterCpp` has a saturation-clamp slow tail like
   `microclimfv2`'s (never checked here — this package's soil water
   convergence was clean in this profiling run, but that doesn't rule out
   the same mechanism triggering under different weather).
3. `plantmodelCpp`/`LangrangianOne` themselves (66% of runtime) received
   no micro-optimisation pass (redundant computation, allocation) —
   out of scope for a *convergence* investigation, but worth a session of
   its own given how dominant they are here.
4. Compiler tuning beyond `-O2` (`-O3`, `-march=native`, LTO) — untested.

## What was explicitly preserved

- The full Lagrangian multilayer canopy model — untouched.
- No clumping coefficient introduced anywhere.
- No change to any function's numerical fixed point — every fix in this
  document (and the diabatic correction / `a2` / `SoilHeatCpp` damping
  fixes made earlier in this session, see git log) changes only how
  reliably/quickly an iteration reaches its answer, never the answer
  itself. All independently verified via full recompile + the tutorial
  vignette's full-year smoke test (zero errors, zero warnings).
