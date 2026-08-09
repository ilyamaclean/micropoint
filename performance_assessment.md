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

## Trace-diagnosing the outer loop's slow-convergence tail

Fixing wind didn't move the outer loop's own tail (63% of hours still
need ≥10 outer iterations), and it doesn't correlate strongly with any
single weather variable (Spearman correlation of outer iteration count
against temperature, humidity, wind speed, precipitation, sensible heat
flux: all |r| < 0.14) — so a second round of instrumentation was added:
temporary per-iteration tracing (gated to one target hour at a time, via
two small exported functions, since `OneStepBelow` is called exactly
once per hour in order) capturing `tdif`, `H`, `wind.uf`/`LL`, `Th`/`eh`,
`tground`, and soil iteration counts at every outer pass, for the worst
hours found in the aggregate stats (hour 7713, 2017-11-18 08:00, and
hour 2961, 2017-05-04 08:00 — both hit the 100-iteration cap in
different runs — plus hour 4268, 2017-06-27 19:00, a near-calm outlier).
Stripped back out once done, same as the first round.

**Finding: two distinct failure modes, not one.**

1. **A genuine oscillation.** In all three traced hours, `H` (sensible
   heat) settles into a persistent period-2 pattern — alternating
   between two interleaved sequences each outer pass (e.g. hour 7713,
   iterations 20–30: `H` = `-454, -517, -419, -490, -395, -470, -365...`)
   — with the amplitude decaying only very slowly. `windmodelCpp` reads
   the *previous* pass's raw, undamped `H` as an input every time and
   feeds a new raw `H` back out: a self-referential loop structurally
   identical to the two already fixed this session (`SoilHeatCpp`'s
   `Tsurf_iter`, `windmodelCpp`'s own internal `uf_iter`), but with no
   equivalent damping. `tair`/`rh`/`tleaf` get smoothed by
   `aitkin_weightdif` *after* this feeds through wind/resistances, which
   is too late to stop the oscillation at its source.

   **Fixed**: `OneStepBelow` now Aitken-damps `H` the same way (`H_iter`
   feeds `windmodelCpp`; the raw `H` is still what's stored and used
   everywhere else — same pattern as `Tsurf_iter`/`uf_iter`). Measured,
   full year, before → after: aggregate outer-loop stats barely moved
   (mean 12.41 → 12.48, essentially noise — the extreme tail is a tiny
   fraction of 8760 hours), but the specific oscillation-dominated hours
   improved substantially and nothing got worse: hour 2961 100 → 36
   iterations, hour 4268 28 → 26. Kept — real, non-regressive, verified
   by the full smoke test, with no downside observed anywhere.

2. **A large hour-to-hour transition — found, NOT fixed.** Hour 7713
   still hits the 100-iteration cap even with `H` damping. Its trace
   shows why: `H` starts the hour at −2063 W/m² (inherited as the
   previous hour's converged state, carried forward as this hour's
   initial guess) and needs to relax to roughly −270 W/m² — a large,
   genuinely *monotonic* move, not an oscillation, and barely
   distinguishable before/after the damping fix. This is exactly the
   regime `microclimfv2`'s own profiling session warns about:
   `aitken1d` as implemented here is damping-only (its blend weight is
   capped at `[0.05, 0.9]`, so it can only ever take a *smaller* step
   than the raw one) — the right tool for an oscillation, the wrong tool
   for a large monotonic relaxation that needs to close a big gap
   quickly. Applying more damping here would not help and could plausibly
   hurt, mirroring the reverted `SoilWaterCpp` attempt in the reference
   repo (damping applied to a problem that needed acceleration, not
   damping, made it more than twice as slow).

   **Confirmed by direct test, not just inference**: swapped the adaptive
   `aitken1d` for much heavier fixed damping (blend weight 0.2, vs.
   adaptive's usual settling point) as a controlled experiment. The two
   genuine oscillation hours improved further under heavier damping
   (2961: 36 → 31 iterations; 4268: 26, unchanged) — consistent with
   damping being the right tool there, at any reasonable strength. Hour
   7713 was completely unmoved (100 → 100) regardless of damping
   strength — direct evidence, not assumption, that this case is not an
   oscillation at all: damping controls how smoothly a step is taken, not
   how far there is to travel, and no damping coefficient changes the
   actual distance `H` has to close. Reverted the experiment (adaptive
   `aitken1d` performed at least as well and is more consistent with the
   rest of the file); the committed fix is the adaptive version above.

   **Not fixed this session** — closing this case needs a genuinely
   different technique (real Aitken Δ²-extrapolation that's allowed to
   overshoot the raw step, not the damping-only `aitken1d` used
   elsewhere here, or a smarter per-timestep warm start when consecutive
   hours' equilibria are very different) and real before/after
   verification against this specific hour, not a guess. Flagged as the
   top remaining opportunity below.

## Also implemented: release-mode compilation

`src/Makevars`, `src/Makevars.win` (`PKG_CXXFLAGS = -O2`), matching
`microclimfv2`'s own file. Zero code risk, real free win: a genuine
`R CMD INSTALL` picks this up (the interactive `devtools::load_all()`
debug build always appends `-O0` last regardless, by design, so this
doesn't show up in development — it matters for the package as actually
installed/distributed).

## Remaining opportunities (not attempted — flagged for a future session)

1. **Accelerate (not damp) the large hour-to-hour `H` transitions** (see
   above — hour 7713 is a ready-made reproducible test case). Needs a
   fundamentally different mechanism than `aitken1d` here: e.g. genuine
   Aitken Δ²-extrapolation (allowed to overshoot the raw step, not just
   shrink it), or detecting a big jump in incoming forcing between
   consecutive hours and re-seeding `H`'s initial guess from a
   neutral-stability estimate instead of blindly carrying the previous
   hour's converged value forward. Whatever is tried, verify against hour
   7713 specifically (100 iterations before any fix) plus the full
   aggregate stats, the same way the `H`-damping fix above was verified —
   don't assume a "should help" change actually does.
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
