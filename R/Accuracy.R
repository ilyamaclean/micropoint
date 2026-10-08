#' What each accuracy level sets
#'
#' One setting, `accuracy`, governs how closely each hour is solved and how much
#' the run costs. Each level fixes:
#' * `tolerance`: how far a pass may still move the canopy air (K) when the hour
#'   is accepted;
#' * `minpass`: the least number of passes an hour takes;
#' * `maxiter`: the most, set at about five times what one hour in a hundred
#'   needs at the level's tolerance, so that an hour reaching it is one that is
#'   not settling;
#' * `canopylayers`: the most canopy layers the run itself is made with. A
#'   canopy given on more is run on that many and its results returned on the
#'   layers the user gave; one given on that many or fewer is run as given,
#'   however few;
#' * `soillayers`: the most soil layers the run is made with;
#' * `leafsteps`: successive approximations in each solution of a leaf's energy
#'   balance;
#' * `soilheatloosening`, `toptolerance`: how closely the soil heat and
#'   canopy-top balances inside each pass are closed.
#'
#' The levels were placed by comparing each with the model solved very tightly,
#' for forest and grass in five kinds of year, against two targets each: the
#' error of air, leaf and ground temperature over a year (0.1, 0.5 and 1 K for
#' very high, high and moderate) and in the worst hour (five times that).
#' @noRd
ACCURACY_LEVELS <- list(
  "very high" = list(tolerance = 0.003, minpass = 3, maxiter = 200, canopylayers = Inf, soillayers = Inf,
                     leafsteps = 4, soilheatloosening = 1, toptolerance = 1e-8),
  "high"      = list(tolerance = 0.1, minpass = 3, maxiter = 75, canopylayers = Inf, soillayers = Inf,
                     leafsteps = 4, soilheatloosening = 1, toptolerance = 1e-8),
  "moderate"  = list(tolerance = 0.3, minpass = 1, maxiter = 50, canopylayers = 10, soillayers = Inf,
                     leafsteps = 2, soilheatloosening = 100, toptolerance = 1e-4))
#' Settings for a run from its accuracy level
#'
#' @param accuracy one of `"high"`, `"very high"` or `"moderate"`. A list of
#'   settings is taken as given, over those of `"high"`.
#' @param tolerance,maxiter where supplied, used in place of the level's own.
#' @return a list of settings; see `ACCURACY_LEVELS`.
#' @noRd
accuracy_settings <- function(accuracy = "high", tolerance = NULL, maxiter = NULL) {
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
#' Hands the engine the parts of a level that are not arguments of a run
#' @noRd
set_solver <- function(st) {
  SolverSettingsCpp(st$minpass, st$leafsteps, st$soilheatloosening, st$toptolerance)
  invisible(NULL)
}
#' One warning for a run in which some hours did not settle
#'
#' @param error the move of canopy air in the last pass of each hour (K).
#' @param st the run's settings.
#' @noRd
warn_unsettled <- function(error, st) {
  k <- sum(error > st$tolerance, na.rm = TRUE)
  if (k > 0) warning(sprintf("The hour did not converge within %d passes in %d of %d hours; results there are those of the steadiest pass",
                             as.integer(st$maxiter), k, length(error)), call. = FALSE)
  invisible(k)
}
#' A canopy described on fewer layers
#'
#' Layers are of equal depth, so the plant area of each new layer is the plant
#' area of the old ones lying within it, counting a split layer by the share of
#' its depth on either side; total plant area and its distribution with height
#' are kept. The live fraction of a new layer is that of the plant area in it.
#' @param paii,Lfrac plant area and live fraction of each layer, from the ground up.
#' @param nmax the most layers the level runs a canopy on.
#' @return `NULL` if the canopy is to be run as given, otherwise a list of
#'   `paii` and `Lfrac` on `nmax` layers.
#' @noRd
coarsen_canopy <- function(paii, Lfrac, nmax) {
  n <- length(paii)
  if (!is.finite(nmax) || n <= nmax) return(NULL)
  below <- function(x) stats::approx((0:n) / n, c(0, cumsum(x)), xout = (0:nmax) / nmax)$y
  p <- pmax(diff(below(paii)), 0)
  live <- diff(below(paii * Lfrac))
  centre <- stats::approx((seq_len(n) - 0.5) / n, Lfrac, xout = (seq_len(nmax) - 0.5) / nmax, rule = 2)$y
  L <- ifelse(p > 0, live / p, centre)
  list(paii = p, Lfrac = pmin(pmax(L, 0), 1))
}
#' A canopy given as a single layer, as the two halves of that layer
#'
#' Each layer exchanges longwave radiation, heat and vapour with the others, so
#' the canopy equations need at least two. One layer is taken as two of half
#' its depth and half its plant area, which describes the same canopy.
#' @param paii,Lfrac plant area and live fraction of each layer.
#' @return a list of `paii` and `Lfrac`, unchanged if there are two layers or more.
#' @noRd
at_least_two_layers <- function(paii, Lfrac) {
  if (length(paii) > 1) return(list(paii = paii, Lfrac = Lfrac))
  list(paii = rep(paii / 2, 2), Lfrac = rep(Lfrac, 2))
}
#' Canopy profiles through time, from the layers a run was made with to those the user gave
#'
#' Every profile returned is a value at the centre of each layer, so each is
#' interpolated in height between the centres the run used. Above the highest of
#' those, air temperature, humidity and wind run on to their canopy-top values,
#' which the model carries. Below the lowest centre, air temperature runs on
#' towards the ground's as the logarithm of height, reaching it at the soil's
#' roughness height for heat, since most of the difference between the ground
#' and the air lies within a short distance of the soil. Everything else beyond
#' the outermost centres takes the value of the nearest.
#' @param mout output of the full run; matrices of hours by layers.
#' @param n the number of layers wanted.
#' @param h canopy height (m).
#' @param zm roughness length of the soil (m).
#' @noRd
canopy_to_user <- function(mout, n, h, zm) {
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
  for (v in c("Rdirdown", "Rdifdown", "Rswup", "Rlwdown", "Rlwup", "tleaf")) mout[[v]] <- onto(mout[[v]])
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
#' Depths of a soil's nodes (m), the surface first
#' @noRd
soil_depths <- function(soilc) geometricCpp(soilc$nLayers, soilc$totalDepth)[1:(soilc$nLayers + 1)]
#' A soil described on fewer layers
#'
#' The same depth of soil with nodes spaced as [createsoilc()] would space
#' them, each property taken from the user's soil at the depth of the new node.
#' @param soilc soil parameter list.
#' @param nmax the most layers wanted.
#' @return `NULL` if the soil has no more than `nmax` layers, otherwise a list
#'   of the soil on `nmax` layers (`soilc`) and the node depths of the user's
#'   soil (`zu`) and of the new one (`zc`).
#' @noRd
coarsen_soil <- function(soilc, nmax) {
  n <- soilc$nLayers
  if (!is.finite(nmax) || n <= nmax) return(NULL)
  zu <- soil_depths(soilc)
  out <- soilc
  out$nLayers <- nmax
  zc <- soil_depths(out)
  for (nm in names(soilc)) {
    if (is.numeric(soilc[[nm]]) && length(soilc[[nm]]) == n + 1) out[[nm]] <- stats::approx(zu, soilc[[nm]], xout = zc, rule = 2)$y
  }
  list(soilc = out, zu = zu, zc = zc)
}
#' Soil profiles through time, from the nodes a run was made with to the user's
#' @param M matrix of hours by nodes.
#' @param cs output of `coarsen_soil()`.
#' @noRd
soil_to_user <- function(M, cs) {
  lo <- findInterval(cs$zu, cs$zc, all.inside = TRUE)
  w <- (cs$zu - cs$zc[lo]) / (cs$zc[lo + 1] - cs$zc[lo])
  out <- sweep(M[, lo, drop = FALSE], 2, 1 - w, `*`) + sweep(M[, lo + 1, drop = FALSE], 2, w, `*`)
  dimnames(out) <- NULL
  out
}
