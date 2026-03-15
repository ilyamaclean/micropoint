#' @title creates ectotherm parameters
#'
#' @description The function `create_ectop` creates a list of geometric,
#' thermal and physiological parameters for an ectotherm from body dimensions
#' and organism type. Body volume is calculated assuming an ellipsoid and
#' surface area is estimated using the Knud Thomsen approximation. Default
#' values are assigned for evaporative skin fraction, thermal conductivity,
#' conductive contact fraction, body tilt and metabolic parameters.
#'
#' @param height numeric value giving body height (cm)
#' @param width numeric value giving body width (cm)
#' @param len numeric value giving body length (cm)
#' @param refl numeric value giving shortwave reflectance of the body surface (unitless)
#' @param adir numeric value giving body orientation relative to north (decimal degrees - see details)
#' @param position character string describing assumptions made about body position (see details)
#' @param type character string giving ectotherm type; one of `"arthropod"`,
#' `"reptile"`, `"amphibian"` or `"other"`.
#'
#' @details
#' `position` can be one of the following:
#' - `"fixed"`: The body position is constant, determined by the specified direction (`adir`) and tilt (`atilt`).
#' - `"randomdir"`: The body tilt is fixed, but the direction the organism faces is random.
#' - `"random"`: Both the direction and tilt of the organism are random.
#' - `"max"`: The organism orients itself to maximize radiation absorption.
#' - `"min"`: The organism orients itself to minimize radiation absorption.
#'
#' @return a list of ectotherm parameters:
#' \describe{
#'   \item{height}{Body height (cm)}
#'   \item{width}{Body width (cm)}
#'   \item{len}{Body length (cm)}
#'   \item{refl}{Shortwave reflectance of the body surface (unitless)}
#'   \item{confrac}{Fraction of body surface in conductive contact with the substrate (unitless)}
#'   \item{skinwetfrac}{Fraction of body surface capable of evaporative water loss (unitless)}
#'   \item{em}{Longwave emissivity of the body surface (unitless)}
#'   \item{volume}{Estimated body volume assuming an ellipsoid (m^3)}
#'   \item{area}{Estimated body surface area using the Knud Thomsen ellipsoid approximation (m^2)}
#'   \item{rho}{Body density used for mass estimation (kg m^-3)}
#'   \item{a0}{Metabolic normalization constant in the relationship \eqn{M = a_0 m^b Q_{10}^{(T-T_{ref})/10}} (W kg^-b)}
#'   \item{b}{Mass scaling exponent for metabolic rate (unitless)}
#'   \item{Q10}{Temperature sensitivity coefficient for metabolic rate (unitless)}
#'   \item{Tref}{Reference temperature for metabolic scaling (°C)}
#'   \item{k}{Thermal conductivity of body tissue (W m^-1 K^-1)}
#'   \item{adir}{Body orientation relative to the direct solar beam (degrees)}
#'   \item{atilt}{Body tilt angle relative to the horizontal (degrees)}
#'   \item{position}{Descriptor of body position in the environment (character)}
#' }
#' @export
create_ectop <- function(height, width, len, refl = 0.1, adir = 0, position = "randomdir",
                         type = c("arthropod", "reptile", "amphibian", "other")) {
  type <- match.arg(type)
  # Estimate typical wet skin fraction
  wetfrac <- switch(type, arthropod = 0.01, reptile = 0.05, amphibian = 0.9, other = 0.03)
  # Volume and surface area
  AA <- len / 200
  BB <- width / 200
  CC <- height / 200
  V <- (4 / 3) * pi * AA * BB * CC
  P <- 1.6075
  A <- 4 * pi * (((AA^P * BB^P + AA^P * CC^P + BB^P * CC^P) / 3)^(1 / P))
  m <- V * 1000
  # variables for computing metabolic rate
  a0 <- switch(type, arthropod = 0.282, reptile = 0.124, amphibian = 0.203, other = 0.124)
  b <- switch(type, arthropod = 0.75, reptile = 0.768, amphibian = 0.884, other = 0.75)
  Q10 <- switch(type, arthropod = 2.0, reptile = 2.44, amphibian = 2.21, other = 2.3)
  Tref <- switch(type, arthropod = 25, reptile = 20, amphibian = 20, other = 20)
  # compute oxygen consumption at reference temperature (mL O2 s^-1)
  VO2ref <- a0 * m^b / 20.1
  # metabolic heat function (W)
  Mfun <- function(Tbody) a0 * m^b * Q10^((Tbody - Tref)/10)
  # oxygen consumption function (mL O2 s^-1)
  VO2fun <- function(Tbody) VO2ref * Q10^((Tbody - Tref)/10)
  # compute other variables
  k <- switch(type, arthropod = 0.2, reptile = 0.5, amphibian = 0.5, other = 0.4)
  confrac <- switch(type, arthropod = 0.01, reptile = 0.2, amphibian = 0.1, other = 0.02)
  atilt <- switch(type, arthropod = 0, reptile = 5, amphibian = 20, other = 5)
  return(list(height = height, width = width, len = len, refl = refl, confrac = confrac,
              skinwetfrac = wetfrac, em = 0.97, volume = V, area = A, rho = 1000,
              a0 = a0, b = b, Q10 = Q10, Tref = Tref, k = k,
              adir = adir, atilt = atilt, position = position))
}
#' Return a vertical ectotherm body-temperature profile for a specified hour
#'
#' Runs the ectotherm model across the vertical microclimate profile returned by
#' [return_profile()] and returns predicted body temperature through the soil,
#' canopy, ground surface, and air above the canopy for a single hour.
#'
#' @param moutprofile List of vertical microclimate variables for a single hour,
#'   as returned by [return_profile()].
#' @param climdata data.frame of meteorological forcing data using the same
#'   naming conventions as the inbuilt dataset `climdata`.
#' @param hr Integer. Hour to evaluate, indexed from 1 to `nrow(climdata)`.
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param ectop List of ectotherm parameters required by [Ectotherm()].
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) above ground of the meteorological forcing data.
#' @param maxiter Integer. Maximum number of iterations used to achieve
#'   convergence in the ectotherm model.
#' @param tolerance Numeric. Convergence tolerance for the ectotherm model.
#'
#' @details
#' The function evaluates ectotherm body temperature separately below ground, at
#' the ground surface, within the canopy, and above the canopy, then combines
#' these into a single vertical profile.
#'
#' Within the canopy, body temperature is calculated by passing the profile
#' microclimate variables in `moutprofile` to [Ectotherm()]. At the ground
#' surface, body temperature is calculated separately using surface soil
#' temperature and an estimate of surface relative humidity derived from soil
#' water potential. Above the canopy, body temperature is calculated using the
#' above-canopy air temperature and humidity profiles together with wind speeds
#' estimated from the canopy wind profile. Below ground, body temperature is
#' assumed to equal soil temperature.
#'
#' @return A data.frame with two columns:
#' \describe{
#'   \item{z}{Vertical position (m), with negative values below ground, `0` at
#'   the soil surface, and positive values above ground.}
#'   \item{Tbody}{Predicted ectotherm body temperature (degrees C).}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
profile_ecto <- function(moutprofile, climdata, hr, vegp, soilc, ectop, lat, long, zref, maxiter = 100, tolerance = 1e-2) {
  # ** Compute ectotherm body temperature at ground level
  n <- 2
  psie <- -abs(soilc$psi_e[1])
  psiw <- psie * (moutprofile$theta / soilc$Smax[1])^-soilc$b[1]
  Tk <- moutprofile$Soiltemp[1] + 273.15
  surfrh <- exp(0.018015 * psiw / (8.314 * Tk))
  tme <- as.POSIXlt(climdata$obs_time[hr], tz = "UTC")
  obstime <- data.frame(year = rep(tme$year +1900, n),
                        month = rep(tme$mon + 1, n),
                        day = rep(tme$mday, n),
                        hour = rep(tme$hour, n))
  climin <- data.frame(Rdirdown = rep(moutprofile$Rdirdown[1], n),
                       Rdifdown = rep(moutprofile$Rdifdown[1], n),
                       Rswup = rep(moutprofile$Rswup[1], n),
                       Rlwdown = rep(moutprofile$Rlwdown[1], n),
                       Rlwup = rep(moutprofile$Rlwup[1], n),
                       uz = rep(moutprofile$uz[1], n),
                       wdir = rep(climdata$winddir[hr], n),
                       Ta = rep(moutprofile$tair[1], n),
                       Ts = rep(moutprofile$Soiltemp[1], n),
                       rh = rep(moutprofile$rh[1], n),
                       pk = rep(climdata$pres[hr], n),
                       surfrh = rep(surfrh, n))
  Tbodyg <- Ectotherm(obstime, climin, ectop, lat, long, vegp$ltra, maxiter, tolerance)[1]
  # ** Compute ectotherm body temperature below ground
  Tbodyb <- rev(moutprofile$Soiltemp)
  n2 <- length(moutprofile$zb) - 1
  zb <- moutprofile$zb[1:n2]
  # ** Compute ectotherm body temperature above ground and within canopy
  # Create obstime data frame
  n <- length(moutprofile$Rswup)
  obstime <- data.frame(year = rep(tme$year +1900, n),
                        month = rep(tme$mon + 1, n),
                        day = rep(tme$mday, n),
                        hour = rep(tme$hour, n))
  # Create input climate data
  climin <- data.frame(Rdirdown = moutprofile$Rdirdown,
                       Rdifdown = moutprofile$Rdifdown,
                       Rswup = moutprofile$Rswup,
                       Rlwdown = moutprofile$Rlwdown,
                       Rlwup = moutprofile$Rlwup,
                       uz = moutprofile$uz,
                       wdir = rep(climdata$winddir[hr], n),
                       Ta = moutprofile$tair,
                       Ts = moutprofile$tair,
                       rh = moutprofile$rh,
                       pk = rep(climdata$pres[hr], n),
                       surfrh = rep(1, n))
  if (class(vegp) == "logical") {
    ectop$confrac <- 0
    Tbodyca <- Ectotherm(obstime, climin, ectop, lat, long, vegp$ltra, maxiter, tolerance)
    z <- c(rev(-zb), 0, moutprofile$z)
  } else {
    climin$Ts <- moutprofile$tleaf
    # Compute ectotherm temperature within canopy
    Tbodyc <- Ectotherm(obstime, climin, ectop, lat, long, vegp$ltra, maxiter, tolerance)
    # Compute bpdy temperature above canopy
    # wind speed
    za <- moutprofile$za
    n2 <- length(za)
    n <- length(moutprofile$Rswup)
    d <- zeroplanedisCpp(vegp$h, vegp$pai)
    zm <- roughlengthCpp(vegp$h, vegp$pai, d, moutprofile$psih)
    uf <- (0.41 * climdata$windspeed[hr]) / (log((zref - d)/zm) + moutprofile$psim)
    psim <- 0
    for (i in 1:length(za)) psim[i] <- dpsimCpp(zm / moutprofile$LL) - dpsimCpp((za[i] - d) / moutprofile$LL)
    uza <- (uf / 0.41) * (log((moutprofile$za - d) / zm) + psim)
    obstime <- data.frame(year = rep(tme$year +1900, n2),
                          month = rep(tme$mon + 1, n2),
                          day = rep(tme$mday, n2),
                          hour = rep(tme$hour, n2))
    climin <- data.frame(Rdirdown = rep(moutprofile$Rdirdown[n], n2),
                         Rdifdown = rep(moutprofile$Rdifdown[n], n2),
                         Rswup = rep(moutprofile$Rswup[n], n2),
                         Rlwdown = rep(moutprofile$Rlwdown[n], n2),
                         Rlwup = rep(moutprofile$Rlwup[n], n2),
                         uz = uza,
                         wdir = rep(climdata$winddir[hr], n2),
                         Ta = moutprofile$tair_above,
                         Ts = rep(moutprofile$tleaf[n], n2),
                         rh = moutprofile$rh_above,
                         pk = rep(climdata$pres[hr], n2),
                         surfrh = rep(1, n2))
    ectop$confrac <- 0
    Tbodya <- Ectotherm(obstime, climin, ectop, lat, long, vegp$ltra, maxiter, tolerance)
    Tbodyca <- c(Tbodyc, Tbodya)
    z <- c(rev(-zb), 0, moutprofile$z, moutprofile$za)
  }
  Tbody <- c(Tbodyb, Tbodyg, Tbodyca)
  out <- data.frame(z = z, Tbody = Tbody)
  return(out)
}
#' Return time-series of ectotherm body-temperature for a specified height
#'
#' Runs the ectotherm model through time using microclimate data returned by
#' [RunMicro()] and returns predicted body temperature for very hour.
#'
#' @param mout List of microclimate variables as returned by [RunMicro()]
#' @param climdata data.frame of meteorological forcing data using the same
#'   naming conventions as the inbuilt dataset `climdata`.
#' @param reqhgt height above or below ground for which body temperature is wanted
#' (negative is below ground). Must match value supplied to [RunMIcro()].
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param ectop List of ectotherm parameters required by [Ectotherm()].
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param maxiter Integer. Maximum number of iterations used to achieve
#'   convergence in the ectotherm model.
#' @param tolerance Numeric. Convergence tolerance for the ectotherm model.
#'
#' @details
#' If `confrac` in `ectop` is > 0, and `0 < reqhgt < vegp$h` the ectotherm is
#' assumed to be sitting on leaves. If `reqhgt = 0` the ectotherm is assumed
#' to be sitting on the ground. Below ground, body temperature is
#' assumed to equal soil temperature.
#'
#' @return A data.frame with two columns:
#' \describe{
#'   \item{obs_time}{POSIXlt object of dates and times associated with eahc body
#'   temperature estimate}
#'   \item{Tbody}{Predicted ectotherm body temperature (degrees C).}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
timeseries_ecto <- function(mout, climdata, reqhgt, vegp, soilc, ectop, lat, long, maxiter = 100, tolerance = 1e-2) {
  if (reqhgt < 0) {  # below ground
    Tbody <- mout$tsoil
  } else {
    tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
    obstime <- data.frame(year = tme$year +1900,
                          month = tme$mon + 1,
                          day = tme$mday,
                          hour = tme$hour)
    n <- length(tme)
    climin <- data.frame(Rdirdown = mout$Rdirdown,
                         Rdifdown = mout$Rdifdown,
                         Rswup = mout$Rswup,
                         Rlwdown = mout$Rlwdown,
                         Rlwup = mout$Rlwup,
                         uz = mout$windspeed,
                         wdir = climdata$winddir,
                         Ta = mout$tair,
                         Ts = mout$tair,
                         rh = mout$relhum,
                         pk = climdata$pres,
                         surfrh = rep(1, n))
    if (reqhgt == 0) {   # on ground surface
      psie <- -abs(soilc$psi_e[1])
      psiw <- psie * (mout$soilwater / soilc$Smax[1])^-soilc$b[1]
      Tk <- mout$tground + 273.15
      climin$surfrh <- exp(0.018015 * psiw / (8.314 * Tk))
    } else { # above ground
      if (class(vegp) != "logical") {
        if (reqhgt < vegp$h) {  # vegetation present and below canopy
          climin$Ts <- mout$tleaf
        } else { # vegetation present and above canopy
          ectop$confrac <- 0
        }
      } else { # no vegetation present
        ectop$confrac <- 0
      }
    }
    # Calculate transmission of surface
    if (class(vegp) != "logical") {
      ltra <- vegp$ltra
    } else ltra = 0
    Tbody <- Ectotherm(obstime, climin, ectop, lat, long, ltra, maxiter, tolerance)
  }
  out <- data.frame(obs_time = climdata$obs_time, Tbody = Tbody)
  return(out)
}
#' Run the full ectotherm model through time
#'
#' Runs the full ectotherm model for all time steps and heights using microclimate
#' data returned by [RunModelFull()]

#' @param mout List of microclimate variables as returned by [RunModelFull()]
#' @param climdata data.frame of meteorological forcing data using the same
#'   naming conventions as the inbuilt dataset `climdata`.
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param ectop List of ectotherm parameters required by [Ectotherm()].
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param maxiter Integer. Maximum number of iterations used to achieve
#'   convergence in the ectotherm model.
#' @param tolerance Numeric. Convergence tolerance for the ectotherm model.
#'
#' @details
#' Below ground, body temperature is assumed to equal soil temperature.
#'
#' @return A list of two Numeric Matrices:
#' \describe{
#'   \item{Tbodya}{Estimated body temperature above ground (degrees C).}
#'   \item{Tbody}{Estimated body temperature below ground (degrees C).}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
full_ecto <- function(mout, climdata, vegp, soilc, ectop, lat, long, maxiter = 100, tolerance = 1e-2) {
  # Calculate for canopy bit
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  # Create List of climate input variables
  climin <- list(Rdirdown = mout$Rdirdown,
                 Rdifdown = mout$Rdifdown,
                 Rswup = mout$Rswup,
                 Rlwdown = mout$Rlwdown,
                 Rlwup = mout$Rlwup,
                 uz = mout$windspeed,
                 Ta = mout$tair,
                 Ts = mout$tair,
                 rh = mout$relhum,
                 pk = climdata$pres,
                 wdir = climdata$winddir)
  if (class(vegp) != "logical") {
    climin$Ts <- mout$tleaf
    ltra <- vegp$ltra
  } else {
    ectop$confrac <- 0
    ltra <- 0
  }
  Tbodya <- EctothermM(obstime, climin, ectop, lat, long, ltra, maxiter, tolerance)
  Tbodyb <- mout$tsoil
  return(list(Tbodya = Tbodya, Tbodyb = Tbodyb))
}
