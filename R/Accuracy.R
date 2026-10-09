# What each accuracy level sets:
#   tolerance          convergence tolerance on canopy air temperature (K)
#   minpass, maxiter   fewest and most passes an hour takes
#   canopylayers       most canopy layers the run is made with
#   soillayers         most soil layers the run is made with
#   leafsteps          refinements of each leaf temperature
#   soilheatloosening  factor by which the soil heat tolerance is relaxed
#   toptolerance       tolerance on canopy-top air temperature (K)
#   groundrow          ground adjusted with the canopy air (2) or held while it is solved (0)
#   nearlocal          near field from a layer's own foliage and its two neighbours only
#   leafclosed         leaf response to air temperature from a formula, not a second solve
ACCURACY_LEVELS <- list(
  "very high" = list(tolerance = 0.003, minpass = 3, maxiter = 200, canopylayers = Inf, soillayers = Inf,
                     leafsteps = 4, soilheatloosening = 1, toptolerance = 1e-8,
                     groundrow = 0, nearlocal = FALSE, leafclosed = FALSE),
  "high"      = list(tolerance = 0.1, minpass = 3, maxiter = 75, canopylayers = Inf, soillayers = Inf,
                     leafsteps = 4, soilheatloosening = 1, toptolerance = 1e-8,
                     groundrow = 0, nearlocal = FALSE, leafclosed = FALSE),
  "moderate"  = list(tolerance = 0.3, minpass = 1, maxiter = 50, canopylayers = 10, soillayers = Inf,
                     leafsteps = 2, soilheatloosening = 100, toptolerance = 1e-4,
                     groundrow = 2, nearlocal = FALSE, leafclosed = FALSE),
  "fast"      = list(tolerance = 0.5, minpass = 1, maxiter = 50, canopylayers = 5, soillayers = Inf,
                     leafsteps = 2, soilheatloosening = 100, toptolerance = 1e-4,
                     groundrow = 2, nearlocal = TRUE, leafclosed = TRUE))
# Settings for a run from its accuracy level. tolerance and maxiter override the level if
# given. A list given as accuracy is used over the settings of "high".
.accuracy_settings <- function(accuracy = "high", tolerance = NULL, maxiter = NULL) {
  if (is.list(accuracy)) {
    st <- utils::modifyList(ACCURACY_LEVELS[["high"]], accuracy)
    if (is.null(st$level)) st$level <- "custom"
  } else {
    accuracy <- match.arg(accuracy, names(ACCURACY_LEVELS))
    st <- ACCURACY_LEVELS[[accuracy]]
    st$level <- accuracy
  }
  if (!is.null(tolerance)) st$tolerance <- tolerance
  if (!is.null(maxiter)) st$maxiter <- maxiter
  st
}
# Passes the solver settings to the C++ code.
.set_solver <- function(st) {
  o <- function(k, d) if (is.null(st[[k]])) d else st[[k]]
  SolverSettingsCpp(st$minpass, st$leafsteps, st$soilheatloosening, st$toptolerance,
                    as.integer(o("groundrow", 0L)), o("nearlocal", FALSE), o("leafclosed", FALSE))
  invisible(NULL)
}
# One warning for the run if some hours did not converge.
.warn_unsettled <- function(error, st) {
  k <- sum(error > st$tolerance, na.rm = TRUE)
  if (k > 0) warning(sprintf("The hour did not converge within %d passes in %d of %d hours; results there are those of the steadiest pass",
                             as.integer(st$maxiter), k, length(error)), call. = FALSE)
  invisible(k)
}
# A canopy re-binned onto nmax layers of equal depth, keeping total plant area and its
# distribution with height. NULL if it already has nmax layers or fewer.
.coarsen_canopy <- function(paii, Lfrac, nmax) {
  n <- length(paii)
  if (!is.finite(nmax) || n <= nmax) return(NULL)
  below <- function(x) stats::approx((0:n) / n, c(0, cumsum(x)), xout = (0:nmax) / nmax)$y
  p <- pmax(diff(below(paii)), 0)
  live <- diff(below(paii * Lfrac))
  centre <- stats::approx((seq_len(n) - 0.5) / n, Lfrac, xout = (seq_len(nmax) - 0.5) / nmax, rule = 2)$y
  L <- ifelse(p > 0, live / p, centre)
  list(paii = p, Lfrac = pmin(pmax(L, 0), 1))
}
# A canopy given as one layer is split into two halves: the model needs at least two layers.
.at_least_two_layers <- function(paii, Lfrac) {
  if (length(paii) > 1) return(list(paii = paii, Lfrac = Lfrac))
  list(paii = rep(paii / 2, 2), Lfrac = rep(Lfrac, 2))
}
# Interpolates canopy profiles from the layers the run used to the n layers the user gave.
# Above the highest layer centre, values run to their canopy-top values; below the lowest,
# air temperature runs to ground temperature with the logarithm of height.
.canopy_to_user <- function(mout, n, h, zm) {
  nc <- ncol(mout$tair)
  if (nc == n) return(mout)
  zc <- c((seq_len(nc) - 0.5) / nc, 1)
  zu <- (seq_len(n) - 0.5) / n
  lo <- findInterval(zu, zc, all.inside = TRUE)
  w <- pmax((zu - zc[lo]) / (zc[lo + 1] - zc[lo]), 0)
  onto <- function(M, top = M[, nc]) {
    Mx <- cbind(M, top)
    out <- sweep(Mx[, lo, drop = FALSE], 2, 1 - w, `*`) + sweep(Mx[, lo + 1, drop = FALSE], 2, w, `*`)
    dimnames(out) <- NULL
    out
  }
  lowest <- mout$tair[, 1]
  for (v in c("Rdirdown", "Rdifdown", "Rswup", "Rlwdown", "Rlwup", "tleaf", "leafwet")) mout[[v]] <- onto(mout[[v]])
  mout$windspeed <- onto(mout$windspeed, mout$ucanopy)
  mout$tair <- onto(mout$tair, mout$tcanopy)
  mout$relhum <- onto(mout$relhum, mout$rhcanopy)
  zh <- zm / 5
  low <- which(zu < zc[1] & zu * h > zh)
  if (length(low) > 0) {
    tground <- mout$tsoil[, 1]
    share <- log(zu[low] * h / zh) / log(zc[1] * h / zh)
    mout$tair[, low] <- tground + outer(lowest - tground, share)
  }
  mout
}
# Depths of the soil nodes (m), surface first.
.soil_depths <- function(soilc) geometricCpp(soilc$nLayers, soilc$totalDepth)[1:(soilc$nLayers + 1)]
# A soil re-described on nmax layers. NULL if it already has nmax layers or fewer.
.coarsen_soil <- function(soilc, nmax) {
  n <- soilc$nLayers
  if (!is.finite(nmax) || n <= nmax) return(NULL)
  zu <- .soil_depths(soilc)
  out <- soilc
  out$nLayers <- nmax
  zc <- .soil_depths(out)
  for (nm in names(soilc)) {
    if (is.numeric(soilc[[nm]]) && length(soilc[[nm]]) == n + 1) out[[nm]] <- stats::approx(zu, soilc[[nm]], xout = zc, rule = 2)$y
  }
  list(soilc = out, zu = zu, zc = zc)
}
# Interpolates soil profiles from the nodes the run used to the user's nodes.
.soil_to_user <- function(M, cs) {
  lo <- findInterval(cs$zu, cs$zc, all.inside = TRUE)
  w <- (cs$zu - cs$zc[lo]) / (cs$zc[lo + 1] - cs$zc[lo])
  out <- sweep(M[, lo, drop = FALSE], 2, 1 - w, `*`) + sweep(M[, lo + 1, drop = FALSE], 2, w, `*`)
  dimnames(out) <- NULL
  out
}
