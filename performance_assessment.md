# Performance assessment (2026-08-09)

Written in response to `restructuretask.txt`'s request for "an overall
assessment of how the model could be speeded up through better design,
improved damping, other checks for convergence, changes to architecture
etc.", with the explicit constraints that overall functionality must stay
the same and no clumping coefficient should be introduced. Both are
respected below — nothing in this document changes model outputs; the one
change actually implemented is a build-configuration flag.

## Method

Rather than re-run a fresh profiling investigation from scratch, this
assessment leans on one already done, days earlier, on the near-identical
`SoilHeatCpp`/`SoilWaterCpp` pair in the sibling `microclimfv2` package
(`src/pointmodel.cpp`, commit `a428791`, documented in that repo's
`CLAUDE.md` under "PERFORMANCE PROFILING (2026-08-04)"). The two
functions are close to line-for-line identical between the packages after
this session's soil-model alignment work (see git log) — the profiling
picture almost certainly transfers.

## Implemented this session

**Release-mode compilation** (`src/Makevars`, `src/Makevars.win`, both
new, `PKG_CXXFLAGS = -O2`) — matches the reference package's own file
exactly. Zero code changed, zero behaviour risk. The reference repo
measured this as a real, ~33% runtime reduction versus the `-O0` debug
build that `devtools::load_all()`/`pkgload` use by default for
interactive development — a real end user's `R CMD INSTALL` /
`install_github()` already picks up R's own default `-O2` on most
platforms, but pinning it explicitly guards against a personal
`~/.R/Makevars` overriding that to something slower, and makes the
intent explicit for anyone building this package from source.

## The main cost driver (not changed this session — see below)

The reference profiling found `SoilWaterCpp` (~58% of runtime) and
`SoilHeatCpp` (~39%) account for essentially all of it (everything else
combined: <4%, confirmed not worth optimising). In this package, both are
called **inside `OneStepBelow`'s own outer convergence loop**
(`R/../src/micropoint_new2.cpp`, "Run model for one step" section): each
outer pass fully re-solves both the soil heat and soil water nonlinear
systems from scratch (each with its own up-to-`maxIter` inner Newton
loop), even though only `OneStepBelow`'s boundary conditions (`Rabs`,
`rHa`, `Th`, `eh`) shift modestly between outer passes as the canopy
profile converges. This is very likely what `restructuretask.txt` meant
by "several of them might be more effectively embedded in one loop."

**This was investigated but not touched.** Two reasons:

1. It is a materially bigger, riskier change than anything else in this
   session's work — it touches the coupling between the canopy
   convergence loop and the two most expensive functions in the model,
   not just a self-contained formula. Verifying it doesn't change
   outputs needs a proper before/after run against real data with
   `witers`/`iters`/`G`/`Tground` comparisons, not just "it still
   compiles and the tutorial runs" (this session's verification bar for
   the smaller, self-contained fixes).
2. The reference repo already tried two related ideas on exactly this
   pair of functions and reverted both after measuring them — directly
   relevant precedent against assuming an "obviously wasteful"-looking
   pattern here is actually free to fix:
   - **Per-layer Aitken relaxation** on `SoilWaterCpp`'s `psiw[i]`
     update, replacing the current dryness-based `lambda` damping:
     runtime **more than doubled** (2.16s → 4.82s), 100-iteration cap
     hits went 34 → 449. `aitken1d` is damping-only (can only shrink a
     step); the slow cases it was applied to needed larger extrapolated
     steps, not smaller ones. **Do not try this on this package's
     `SoilWaterCpp` either, for the same reason.**
   - **Excluding a saturation-clamped layer's residual from the
     `massBalance` convergence check**: reduced cap-hits (34 → 26) but
     gave no measurable net runtime improvement, while changing the
     actual convergence stopping criterion. Reverted.

Given that evidence, restructuring the outer/inner loop coupling without
the ability to profile-and-verify against real data in this session would
be exactly the kind of speculative change the reference repo's own
experience warns against. Recommendation: revisit with a proper
before/after harness (timing + output-diff against a known-good run),
not as a blind refactor.

## Remaining opportunities (untested here, flagged for a future session)

Carried over from the reference repo's own "not attempted" list, equally
applicable to this package's near-identical soil functions:

1. A *complete* fix for `SoilWaterCpp`'s saturation-clamp slow-convergence
   tail — making a clamped layer a proper Dirichlet row in the
   tridiagonal system, rather than the damping-only or
   stopping-criterion tweaks already tried and reverted. Bigger change to
   the linear system assembly; not attempted anywhere yet.
2. Whether `SoilHeatCpp` has a similar slow-tail pattern to
   `SoilWaterCpp`'s clamp issue — never actually investigated (only its
   damping was fixed).
3. `SoilHeatCpp`/`SoilWaterCpp` allocate several fresh `std::vector`s on
   every call (tens of thousands of calls per full-year run) — a
   plausible but never-quantified allocation cost. Would need real
   profiling to know if it matters.
4. Compiler tuning beyond `-O2` (`-O3`, `-march=native`, LTO) — untested,
   and platform-dependent (this package ships as source, built on the
   end user's machine, so `-march=native` in particular needs care).
5. The `OneStepBelow` outer/inner loop coupling described above.

## What was explicitly preserved

- The full Lagrangian multilayer canopy model (`LangrangianOne`,
  `plantmodelCpp`, the below-canopy dispersion machinery) — untouched,
  per `restructuretask.txt`'s explicit instruction not to drop it in
  favour of a bulk BigLeaf approach.
- No clumping coefficient was introduced anywhere.
- No change to any function's numerical behaviour beyond the diabatic
  correction, `a2`, and `SoilHeatCpp` damping fixes documented separately
  in git log (each independently verified via full recompile + the
  tutorial vignette's full-year smoke test).
