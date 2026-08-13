# Handover: stable-side Obukhov-length safeguard missing in `micropoint`'s `clipMOlength`

**STATUS: not started.** The identical bug was found and fixed in `microclimfv2` on 2026-08-13 (see that
repo's `CLAUDE.md`, dated PRIORITY entry same date, and commit `b887102`, "clipMOlength: add the missing
stable-side Obukhov-length safeguard"). This doc is the handover for doing the equivalent fix here, per
direct instruction from Ilya: fix cleanly in both packages, but have a separate Claude Code session action
the `micropoint` side.

**For a fresh Claude Code session to action** -- read this whole doc before touching code, then implement,
verify, and commit in this repo (`git commit`, not just edit-and-leave).

**Target file**: `C:\micro\micropoint\src\micropoint_new2.cpp`.

## The mechanism, traced directly in the real source (both repos)

`clipMOlength` (`micropoint_new2.cpp:512-538`) is meant to keep the Monin-Obukhov length `L` within a
physically sensible bound of the neutral log-law profile. Every one of its 4 call sites here wraps the
result in what looks like a proper two-sided clamp:

```cpp
if (H > 0) { if (LL < Lsafe) LL = Lsafe; }
else       { if (LL > Lsafe) LL = Lsafe; }
```

But `clipMOlength` itself only ever modifies `L` when `L < 0` (the `if (L < 0.0 && psim < psim_min)` guard
at line 519) -- for `L >= 0` (stable) it falls through to `return L;` at line 537, unchanged. So
`Lsafe = clipMOlength(LL, ...)` comes back numerically identical to whatever `LL` was fed in, and
`if (LL > Lsafe)` becomes `if (LL > LL)` -- structurally impossible to fire. Looks symmetric at the call
site; isn't, underneath. Confirmed by reading the code directly (this file, this session), not assumed from
`microclimfv2`'s own equivalent finding -- `micropoint`'s `clipMOlength` is byte-for-byte the same
unstable-only logic, just with `dpsimCpp2` in place of `dpsimCpp`.

**Practical effect**: under strong stable stratification and near-calm wind, resistance can diverge the same
way it can on the unstable/hot side -- producing an implausibly *cold* surface/canopy temperature, with no
safeguard catching it. Ilya's own framing when this was raised: "Doesn't seem to be causing major problems
but good to fix cleanly" -- i.e. not an emergency, but a real, worth-closing gap.

## The fix that was applied in `microclimfv2` -- port this, don't reinvent it

`microclimfv2`'s `src/utils.cpp` extended `clipMOlength` to handle `L >= 0` symmetrically, reusing an
already-existing function, `lStableFinalCpp` (`src/utils.cpp:128-177` in that repo), for the inversion.
**`micropoint` does not have an equivalent of `lStableFinalCpp` at all** (confirmed: no match for
`lStableFinalCpp`/`recoverLCpp` anywhere in this repo's `src/`) -- it needs porting in first, not just
`clipMOlength` itself. `lStableFinalCpp` has no dependencies beyond `dpsimCpp`/`dpsimCpp2` (this package's own
stable/unstable psi_m formula), so the port should be close to verbatim: copy the function body, rename
`dpsimCpp` calls to `dpsimCpp2`, keep everything else identical (including its own internal comments, which
explain the fold-back below).

The fixed `clipMOlength` (`microclimfv2` `src/utils.cpp`, current state):

```cpp
double clipMOlength(double L, double zref, double d, double zm, double beta)
{
    const double ln_z = std::log((zref - d) / zm);
    const double psim_min = -beta * ln_z;
    const double tol = 1e-4;
    double psim = dpsimCpp(zm / L) - dpsimCpp((zref - d) / L);
    if (L < 0.0 && psim < psim_min) {
        // ... unchanged unstable-branch bisection, not reproduced here ...
        return L_high;
    }
    if (L >= 0.0) {
        const double psim_max = beta * ln_z;
        const double L_bound = lStableFinalCpp(psim_max, zref, d, zm);
        if (L < L_bound) {
            return L_bound;
        }
    }
    return L;
}
```

Port this shape directly (`dpsimCpp`->`dpsimCpp2`, `lStableFinalCpp` ported in as above).

## The gotcha that will bite a naive port -- read this before writing any code

The OBVIOUS first implementation -- literally mirroring the unstable branch's own `if (psim < psim_min)`
check as `if (psim > psim_max)` -- **is wrong on the stable branch specifically**, and was caught in
`microclimfv2` by a standalone test *before* it shipped, not assumed correct because it compiled.

`psim(L)` for `L >= 0` has a genuine interior fold-back (this is exactly what `lStableFinalCpp`'s own doc
comment in `microclimfv2` describes, for the *inverse* direction): as `L` decreases from `+infinity` toward
`0+`, `psim(L)` first RISES from 0 up to a peak at some `Lpeak`, then FALLS back down toward 0 again as `L`
continues toward `0+`. A deliberately extreme test case, `L=1e-6`, gives `psim(L)` essentially **0** --
indistinguishable from a perfectly neutral length -- despite being nonsensically over-stable. A naive
`psim(L) > psim_max` check silently treats this as "already safe" and lets the single worst case straight
through uncapped.

**The fix**: compare `L` itself against `L_bound = lStableFinalCpp(psim_max, ...)`, not `psim(L)` against
`psim_max`. On the safe monotonic branch (`L > Lpeak`), `psim(L)` decreases as `L` increases, so
`L < L_bound` and `psim(L) > psim_max` are equivalent there. But for any `L <= Lpeak` (including the
wrapped-past-the-peak case), `L` is always `< L_bound` too (since `L_bound` is itself `> Lpeak` whenever
`psim_max` is reachable at all), so the direct-`L` comparison catches it correctly either way. This is
exactly the code shape above -- don't "simplify" it back to a `psim`-based comparison.

**A separate thing to expect, not a bug if you see it**: for small roughness lengths relative to `zref`
(large `ln_z`), `beta*ln_z` can exceed the stable-branch formula's own achievable peak `psim` value
entirely. In that regime `lStableFinalCpp` correctly saturates to `Lpeak` (its own documented fallback --
"a target psi_m at or beyond what this equation can produce at all... saturates to Lpeak") rather than a
value that reproduces `psim_max` exactly. Confirmed genuine saturation, not a bug, in `microclimfv2` by
checking `Lsafe` stays byte-identical across wildly different `beta` values (0.9/2.0/5.0) for the affected
`(zref,zm)` combinations -- if it were a bug, a bigger target would have moved the answer. Worth re-running
the same sanity check here after porting, not assuming it transfers automatically.

## A second, real issue found while tracing the call sites here -- needed for the fix to actually bind, not just exist

All 4 of `micropoint`'s own call sites compute `Lsafe` **once, before** the convergence loop that actually
uses it, from whatever `LL` happened to be at that point -- not fresh, from the current iteration's own
`LL`, every pass. Confirmed by reading each site directly:

- **`windmodelCpp`, line ~611-614**: `double LL = 1e99; if (H > 0.0) LL = (...); double Lsafe =
  clipMOlength(LL, zref, d, zm);` -- for `H <= 0` (stable/neutral), the seed `LL` fed into `clipMOlength`
  is the neutral **sentinel `1e99`**, never a real physical value, because this line uses `H > 0.0` (not
  `H != 0.0`). Even with `clipMOlength` itself fixed, `clipMOlength(1e99, ...)` returns `1e99` unchanged
  (since `1e99` is never `< L_bound`), so `Lsafe` stays at the sentinel for the *entire* convergence loop
  that follows (lines ~631-650) -- the stable-side fix would be silently inert at this specific call site
  without also fixing this. Worth noting: `microclimfv2`'s own `windmodelCpp` had this exact `H > 0.0` ->
  `H != 0.0` bug fixed back on 2026-08-04 (see that repo's `CLAUDE.md`), and that fix's own commit message
  claims it was "mirrored in `micropoint`" -- confirmed by direct inspection (this session, not left as a
  TODO) that the mirroring was incomplete, not lost: this specific line (612) still has the old `H > 0.0`
  form, while sibling guards in the very same file (lines 631, 652, 2370) already correctly use
  `H != 0.0` -- the fix landed almost everywhere it needed to except this one seed line, most likely just
  missed during that pass rather than deliberately reverted.
- **Bare-soil surface iteration, line ~2353-2356**: `LL` computed unconditionally (no `H>0`/`H!=0` guard, so
  no sentinel problem) but still only once, before the `while` loop at line ~2368, using `H = onestepin.H`
  (last timestep's converged value, not this iteration's live one).
- **`solveonestep` (canopy), line ~2715-2718**: inside its own loop, `Lsafe` is computed from the *previous*
  iteration's `LL`, then `LL` itself is recomputed fresh immediately after -- i.e. `Lsafe` is always one
  iteration behind, not from a sentinel but not current either.
- **`solveonestepbare`, line ~2806-2809**: same one-iteration-behind pattern as `solveonestep`.

**Recommended fix, uniform across all 4**: move `Lsafe = clipMOlength(...)` to *immediately after* each
loop's own fresh per-iteration `LL` recompute, replacing the pre-loop `Lsafe` entirely (delete the pre-loop
`LL`/`Lsafe` pair once confirmed unused elsewhere -- check each site individually, since the exact
surrounding code differs slightly between the four; at `windmodelCpp` specifically, the pre-loop `LL`
computation at line 611-612 appears to be used for nothing except seeding the now-to-be-deleted pre-loop
`Lsafe`, so removing it should also retire that `H>0.0`-only bug as a side effect rather than needing a
separate fix -- but verify this before deleting, don't assume it from this description alone). This is a
mechanical, low-risk change (relocate one line, delete its now-dead predecessor) -- `clipMOlength`/
`lStableFinalCpp`'s own cost is cheap (confirmed negligible in `microclimfv2`'s real-data benchmark), so
recomputing every iteration instead of once isn't a performance concern.

## Suggested validation

1. **Self-consistency, mirroring `microclimfv2`'s own verification**: a standalone test (this file's
   `clipMOlength`/`lStableFinalCpp` have no Rcpp dependency once isolated, so a plain g++ harness works --
   `microclimfv2` used exactly this approach) confirming `psim(Lsafe)` reproduces `psim_max` on the safe
   monotonic branch, and that a deliberately extreme `L` (e.g. `1e-6`) still gets clamped correctly (not
   silently passed through) on the wrapped-past-the-peak branch.
2. **Real-data before/after**: run a full real year (or whatever this repo's own standard regression driver
   is) and confirm no new non-convergence / NaN introduced.
3. **A deliberately extreme stable/cold/near-calm synthetic scenario** (mirroring `microclimfv2`'s own check:
   cold air, near-zero incoming shortwave, low incoming longwave, ~0.15 m/s wind, sustained several hours) --
   confirm it converges cleanly with a physically ordinary cold-night profile, not a divergence, and ideally
   confirm (e.g. via a temporary counter, removed before committing) that the stable-side clamp actually
   fires somewhere in that scenario -- i.e. that the fix isn't just present but inert.
4. Full existing test-suite / vignette pass if this repo has one, per its own normal commit standards.

Once implemented and verified, commit in this repo and update this repo's own `CLAUDE.md` with a dated entry
(this repo already has that convention -- see e.g. `multilayer_TL_blowup_handover.md`'s own STATUS-line
pattern for how a completed handover gets marked done here).
