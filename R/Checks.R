#' Check and tidy model inputs
#'
#' Tests the weather, vegetation and soil inputs against the ranges they can
#' physically take and the ranges they plausibly take, before a model run.
#' Impossible values with a single sensible correction are corrected, with a
#' warning; impossible values without one stop the run; possible but unusual
#' values give a warning. Every model run function calls this first.
#'
#' @param climdata data.frame of hourly weather, named as in the inbuilt
#' dataset `climdata`.
#' @param vegp vegetation parameter list as returned by [createvegp()], or `NA`
#' for bare ground.
#' @param soilc soil parameter list as returned by [createsoilc()].
#' @param lat latitude (decimal degrees), used to place the sun.
#' @param long longitude (decimal degrees), used to place the sun.
#' @param paii optional vector of plant area index for each canopy layer.
#' @param Lfrac optional vector of the live leaf fraction of each canopy layer.
#'
#' @details
#' **Direct and diffuse shortwave.** The direct beam is recorded as the part of
#' `swdown` not in `difrad`, as it arrives on a horizontal surface. Dividing by
#' the cosine of the solar zenith angle gives the beam as it arrives normal to
#' the sun, which is what leaves and slopes intercept. That cannot exceed what
#' the clearest sky lets through: the solar constant reduced by the depth of
#' atmosphere the beam has crossed, which is many times greater with the sun
#' near the horizon than overhead. The limit used is
#' `1352 * 0.76^((m * pres / 101.3)^0.678)` W m^-2, where `m` is the relative
#' air mass of Kasten and Young (1989). The form is that of Meinel and Meinel
#' (1976), in which each further unit of air mass removes a little less than
#' the one before, with the transmittance of a very clear atmosphere (0.76 for
#' an overhead sun at sea level) in place of their typical one. It is a
#' ceiling, not an estimate of the beam: no account is taken of season, water
#' vapour or aerosol.
#' Hourly records sum an hour of radiation, but the sun's position is computed
#' for a single instant, and the two disagree most near the horizon: the implied
#' beam can be impossibly strong, or present while the sun is below the horizon.
#' In those hours the excess is treated as diffuse, so total shortwave is always
#' the recorded value and only its partition changes.
#'
#' @return A list with elements `climdata`, `vegp` and `soilc`, with any
#' corrections applied.
#' @references Kasten F, Young AT (1989) Revised optical air mass tables and
#'   approximation formula. Applied Optics 28: 4735-4738.
#'
#'   Meinel AB, Meinel MP (1976) Applied Solar Energy. Addison Wesley.
#' @export
runchecks <- function(climdata, vegp, soilc, lat, long, paii = NA, Lfrac = NA) {
  # ~~~~ Weather ~~~~
  need <- c("obs_time", "temp", "relhum", "pres", "swdown", "difrad", "lwdown",
            "windspeed", "winddir", "precip")
  miss <- setdiff(need, names(climdata))
  if (length(miss) > 0) stop("climdata is missing: ", paste(miss, collapse = ", "))
  if (anyNA(climdata[, need])) stop("climdata contains NAs")
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  if (length(tme) > 1) {
    step <- as.numeric(diff(as.POSIXct(tme)), units = "secs")
    if (any(step != 3600)) stop("climdata must be hourly and without gaps")
  }
  if (min(climdata$temp) < -65) warning("Minimum temperature seems low")
  if (max(climdata$temp) > 65) warning("Maximum temperature seems high")
  if (min(climdata$relhum) < 0) stop("Relative humidity cannot be negative")
  if (max(climdata$relhum) > 100) {
    warning("Relative humidity above 100% set to 100%")
    climdata$relhum <- pmin(climdata$relhum, 100)
  }
  if (median(climdata$pres) > 200) stop("pres must be in kPa")
  if (min(climdata$pres) < 50) warning("Minimum pressure seems low. OK if site at high altitude")
  if (max(climdata$pres) > 108.5) warning("Maximum pressure seems high")
  if (min(climdata$swdown) < 0 || min(climdata$difrad) < 0) {
    warning("Negative shortwave radiation set to zero")
    climdata$swdown <- pmax(climdata$swdown, 0)
    climdata$difrad <- pmax(climdata$difrad, 0)
  }
  s <- which(climdata$difrad > climdata$swdown)
  if (length(s) > 0) {
    warning("Diffuse radiation above total shortwave in ", length(s), " hours. Set to total shortwave")
    climdata$difrad[s] <- climdata$swdown[s]
  }
  # Beam normal to the sun cannot exceed what the clearest sky passes (see details)
  obstime <- data.frame(year = tme$year + 1900, month = tme$mon + 1, day = tme$mday,
                        hour = tme$hour)
  alt <- solaltCpp(obstime, lat, long)
  cosz <- sin(alt * pi / 180)
  beamh <- climdata$swdown - climdata$difrad
  beammax <- clearestbeam(alt, climdata$pres) * pmax(cosz, 0)
  s <- which(beamh > beammax)
  if (length(s) > 0) {
    warning("Direct beam stronger than the sun can deliver, or with the sun below the horizon, in ",
            length(s), " hours. Excess treated as diffuse")
    climdata$difrad[s] <- climdata$swdown[s] - beammax[s]
  }
  lwmax <- 1.1 * 5.67e-8 * (climdata$temp + 273.15)^4
  if (min(climdata$lwdown) < 0) stop("Longwave radiation cannot be negative")
  if (any(climdata$lwdown > lwmax)) warning("Downward longwave radiation seems high given air temperature")
  if (min(climdata$windspeed) < 0) stop("Wind speed cannot be negative")
  if (max(climdata$windspeed) > 55) warning("Maximum wind speed seems high")
  climdata$winddir <- climdata$winddir %% 360
  if (min(climdata$precip) < 0) stop("Precipitation cannot be negative")
  if (max(climdata$precip) > 150) warning("Maximum precipitation rate seems high")
  # ~~~~ Vegetation ~~~~
  if (!(length(vegp) == 1 && is.na(vegp[1]))) {
    if (vegp$h < 0) stop("vegp$h cannot be negative")
    if (vegp$pai < 0) stop("vegp$pai cannot be negative")
    if (vegp$pai > 20) warning("vegp$pai seems high")
    if (vegp$x <= 0) stop("vegp$x must be positive")
    for (nm in c("lref", "ltra", "lrefp", "ltrap", "vegem", "pTAW")) {
      if (vegp[[nm]] < 0 || vegp[[nm]] > 1) stop("vegp$", nm, " must be between 0 and 1")
    }
    if (vegp$lref + vegp$ltra > 1) stop("vegp$lref + vegp$ltra cannot exceed 1")
    if (vegp$lrefp + vegp$ltrap > 1) stop("vegp$lrefp + vegp$ltrap cannot exceed 1")
    if (vegp$len <= 0 || vegp$wid <= 0) stop("vegp$len and vegp$wid must be positive")
    if (vegp$mwft < 0) stop("vegp$mwft cannot be negative")
    if (vegp$Vcmx25 < 0) stop("vegp$Vcmx25 cannot be negative")
    if (vegp$Tup <= vegp$Tlow) stop("vegp$Tup must exceed vegp$Tlow")
    if (vegp$psi50 > 0) stop("vegp$psi50 must be negative")
    if (is.null(vegp$root50) || is.null(vegp$root95))
      stop("vegp needs root50 and root95, the depths (m) above which 50% and 95% of roots lie; rootskew is no longer used")
    if (vegp$root50 <= 0) stop("vegp$root50 must be positive")
    if (vegp$root95 <= vegp$root50) stop("vegp$root95 must exceed vegp$root50")
    if (!(length(paii) == 1 && is.na(paii[1]))) {
      if (any(paii < 0)) stop("paii cannot be negative")
      # Plant area includes the stems and branches that hold foliage up, so every
      # layer of a real canopy holds some. A layer with none above foliage means
      # the canopy height is above the vegetation; one below foliage means a
      # crown with nothing beneath it.
      if (any(paii == 0))
        stop("paii must be greater than zero in every canopy layer: plant area includes the stems, ",
             "so an empty layer means either a canopy height above the vegetation or foliage with nothing holding it up")
      if (abs(sum(paii) - vegp$pai) > 0.01 * max(vegp$pai, 1e-6))
        warning("sum(paii) differs from vegp$pai")
    }
    if (!(length(Lfrac) == 1 && is.na(Lfrac[1]))) {
      if (any(Lfrac < 0 | Lfrac > 1)) stop("Lfrac must be between 0 and 1")
      if (!(length(paii) == 1 && is.na(paii[1])) && length(Lfrac) != length(paii))
        stop("paii and Lfrac must have one value per canopy layer")
    }
  }
  # ~~~~ Soil ~~~~
  for (nm in c("gref", "grefPAR", "groundem")) {
    if (soilc[[nm]] < 0 || soilc[[nm]] > 1) stop("soilc$", nm, " must be between 0 and 1")
  }
  if (soilc$slope < 0 || soilc$slope > 90) stop("soilc$slope must be between 0 and 90 degrees")
  soilc$aspect <- soilc$aspect %% 360
  nn <- soilc$nLayers + 1
  for (nm in c("Smax", "Smin", "Ksat", "b", "psi_e", "rho")) {
    if (length(soilc[[nm]]) != nn) stop("soilc$", nm, " must have nLayers + 1 values")
  }
  if (any(soilc$Smin < 0 | soilc$Smax > 1)) stop("soilc$Smin and soilc$Smax must be between 0 and 1")
  if (any(soilc$Smin >= soilc$Smax)) stop("soilc$Smin must be below soilc$Smax")
  if (any(soilc$Ksat <= 0)) stop("soilc$Ksat must be positive")
  if (any(soilc$b <= 0)) stop("soilc$b must be positive")
  return(list(climdata = climdata, vegp = vegp, soilc = soilc))
}
#' Strongest direct beam the clearest sky can deliver
#'
#' @param alt solar altitude (degrees above the horizon).
#' @param pres atmospheric pressure (kPa).
#' @return Beam irradiance normal to the sun (W m^-2); zero with the sun at or
#'   below the horizon.
#' @details
#' The solar constant reduced by a very clear atmosphere's transmittance. The
#' depth of air the beam crosses is the relative air mass of Kasten and Young
#' (1989), which stays finite at the horizon, scaled by pressure for the site's
#' elevation. Transmittance falls with that depth as Meinel and Meinel (1976)
#' have it: each further unit of air mass removes less than the one before,
#' because the wavelengths most easily scattered are already gone. A constant
#' transmittance per unit air mass would cut the beam of a low sun to half its
#' real strength. The base, 0.76 for an overhead sun at sea level (about 1030
#' W m^-2 at the ground), is that of an exceptionally clear day where they give
#' 0.7 for a typical one; it is chosen to bound the beam from above, not to
#' estimate it.
#' @noRd
clearestbeam <- function(alt, pres) {
  up <- alt > 0
  m <- rep(Inf, length(alt))
  m[up] <- 1 / (sin(alt[up] * pi / 180) + 0.50572 * (alt[up] + 6.07995)^-1.6364)
  ifelse(up, 1352 * 0.76^((m * pres / 101.325)^0.678), 0)
}
