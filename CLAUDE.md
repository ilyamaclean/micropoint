# CLAUDE.md

Guidance for working in this repository with Claude Code.

## What this is

`micropoint` is an R package (Rcpp-backed) implementing a mechanistic,
multilayer canopy–soil point microclimate model, plus an ectotherm
body-temperature model that consumes its output. See
`model_architecture.md` for the technical design.

## Branch layout

This repository was split into two branches on 2026-08-09:

- **`main`** — the *current* model. Everything here is "the new
  version": `R/NewMicroPoint.R`, `R/Ectotherm.R`, `R/PlantGeometry.R`,
  `src/micropoint_new2.cpp` / `src/micropointheaders3.h`. This is what
  gets developed going forward.
- **`legacy`** — the old model, frozen. `R/BigLeaf.R`, `R/Ground.R`,
  `R/Radiation.R`, `R/Runfunctions.R`, `src/micropoint.cpp` /
  `src/micropointheaders.h`. Kept for backward compatibility /
  reference only; do not add new features here.

The two versions do not share any R or C++ source files. The inbuilt
datasets (`data/*.rda`, documented in `R/data.R`) are duplicated
identically across both branches since both versions use a subset of
them.

## Build workflow

After changing any `.cpp`/`.h` file under `src/`, or adding/removing
`// [[Rcpp::export]]` functions, or changing roxygen `#'` comments in
`R/`, regenerate the generated files in this order:

```r
Rcpp::compileAttributes(".")                          # src/RcppExports.cpp, R/RcppExports.R
roxygen2::roxygenise(".", roclets = c("collate", "namespace", "rd"))  # NAMESPACE, man/*.Rd
```

`roxygenise()` also triggers a recompile via `pkgload`/`devtools`, so
it doubles as a compile check. To load and smoke-test interactively:

```r
devtools::load_all(".")
```

Compiled build artifacts (`*.o`, `*.dll`) are checked into this repo
(unusual for an R package, but that's the existing convention here —
don't gitignore them without discussing it first).

## Known parked issues

- A round of numerical fixes to the soil water / heat solver
  (`micropointheaders3.h`, `micropoint_new2.cpp`) and to
  `newsoilparamstable` was committed as a checkpoint
  (`45a7758`) but **has not been independently verified for
  correctness** — only that the package compiles and the tutorial
  vignette runs without error. Validating the numerics against
  known-good output is a separate, not-yet-started task.
- `roughlengthCpp2()` (used by `Ectotherm.R`'s `profile_ecto()`) takes
  no diabatic-correction argument, unlike the legacy `roughlengthCpp()`
  it replaced — it always uses a fixed roughness-sublayer constant
  (see the comment above `roughlengthCpp2` in
  `src/micropoint_new2.cpp`). This was a required signature match, not
  a stylistic choice.
- The above-canopy wind/temperature profile functions (`Tabove`,
  `RHabove`, `Uabove`) assume `zref` is above the canopy top. Calling
  them with `zref` below canopy height (e.g. a tall forest with
  `zref = 2`) produces `NaN` — seen during testing when the tutorial's
  own example parameters were combined incorrectly, not a repo bug,
  but worth being aware of if you see the same warning.

## Vignettes

`vignettes/new_model_tutorial.Rmd` and `vignettes/modelequations.Rmd`
document the current model. Both were manually verified to run
end-to-end against this branch on 2026-08-09 (see git log for the
commit that fixed the one structural gap found:
`Ectotherm.R` calling unexported wind-profile functions).

Note `.Rbuildignore` excludes all `vignettes/*.Rmd` from `R CMD build`
— vignettes here are maintained/rendered externally (there's an
`rsconnect`/RPubs publish history in `vignettes/rsconnect/`, itself
gitignored) rather than built as part of the package.
