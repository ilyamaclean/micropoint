#' Create list of vegetation parameters for running the microclimate model
#'
#' @param vegtype Character string giving the plant function type for which
#' parameters are to be created (see details)
#' @details
#' Recognised plant functional types are `BET.Tr`	- Broadleaf Evergreen Tree (Tropical),
#' `BET.Te` - Broadleaf Evergreen Tree (Temperate), 	`BDT`	-  Broadleaf Deciduous Tree,
#' `NET` - Needleleaf Evergreen Tree,	`NDT`	- Needleleaf Deciduous Tree, `C3` - C3 grass,
#' 	`C4` - C4 grass,	`ESh` - Evergreen Shrub and	`DSh`- Deciduous Shrub.
#'
#' @return a list of vegetation parameters:
#' \describe{
#'   \item{h}{Vegetation height (m)}
#'   \item{pai}{Plant Area Index (m^2 m^-2)}
#'   \item{x}{Campbell leaf angle distribution parameter (unitless)}
#'   \item{lref}{Leaf reflectance (shortwave) (unitless)}
#'   \item{ltra}{Leaf transmittance (shortwave) (unitless)}
#'   \item{lrefp}{Leaf reflectance (PAR) (unitless)}
#'   \item{ltrap}{Leaf transmittance (PAR) (unitless)}
#'   \item{Vcmx25}{Maximum Rubisco carboxylation rate at 25 °C (µmol m^-2 s^-1)}
#'   \item{Tup}{Upper temperature limit for photosynthesis (°C)}
#'   \item{Tlow}{Lower temperature limit for photosynthesis (°C)}
#'   \item{alpha}{Photosynthetic quantum efficiency (mol CO2 mol^-1 PAR)}
#'   \item{Kxmx}{Maximum xylem hydraulic conductance (µmol m^-2 s^-1)}
#'   \item{hv}{Huber value (cm^2 m^-2)}
#'   \item{fd}{Fraction of Vcmax converted to leaf dark respiration (unitless)}
#'   \item{psi50}{Water potential causing 50\% loss of maximum hydraulic conductance (MPa)}
#'   \item{apsi}{Xylem hydraulic conductance parameter (calculated if \code{< 0}) (unitless)}
#'   \item{len}{Leaf length (m)}
#'   \item{wid}{Leaf width (m)}
#'   \item{vegem}{Vegetation emissivity (unitless)}
#'   \item{mwft}{Maximum water film thickness (mm)}
#'   \item{pTAW}{Fraction of total available soil water that can be depleted before stress (unitless)}
#'   \item{root50}{Depth above which half of all roots lie (m)}
#'   \item{root95}{Depth above which 95\% of all roots lie (m)}
#' }
#' Root depths are from Schenk & Jackson (2002), Table 4.
#' @references Schenk, H.J. & Jackson, R.B. (2002) The global biogeography of
#' roots. Ecological Monographs, 72, 311-328.
#' @export
createvegp <- function(vegtype = "BET.Te") {
  nms <- names(PFTparams)
  s <- which(nms == vegtype)
  v <- PFTparams[,s] * PFTparams$multiplier
  vnms <- PFTparams$varname
  vegp <- as.list(v)
  names(vegp) <- vnms
  return(vegp)
}
#' Create soil parameter list
#'
#' Creates a list of soil physical parameters for a selected soil type and
#' replicates them across soil layers, with optional adjustment for slope,
#' aspect, surface organic matter, and deep boundary saturation.
#'
#' @param soiltype Character string giving the soil type. See details. Default is \code{"Clay loam"}.
#' @param nlayers Integer. Number of soil layers.
#' @param totalDepth Numeric. Total soil depth (m).
#' @param slope Numeric. Surface slope (degrees).
#' @param aspect Numeric. Surface aspect (degrees).
#' @param surface_organicmu Numeric. Multiplication factor controlling the
#'   vertical profile of soil organic matter near the surface.
#' @param FreeDrain Logical. If `TRUE`, a free-drainage lower boundary is
#' assumed, allowing water to drain out of the bottom of the soil profile under
#' gravity. If `FALSE`, the lower boundary is treated as saturated, meaning the
#' soil is assumed to be sitting above a water table, which limits drainage from
#' the deepest layer.
#'
#' @details
#' The selected soil parameter values are replicated across \code{nlayers + 1}
#' soil nodes. Organic matter content is adjusted vertically with corresponding
#' rescaling of mineral and related soil constituents to preserve total volume.
#' Bulk density is also adjusted slightly with depth.
#'
#' Available soil types (argument \code{soiltype}) are:
#' \itemize{
#'   \item Sand
#'   \item Loamy sand
#'   \item Sandy loam
#'   \item Loam
#'   \item Silt loam
#'   \item Sandy clay loam
#'   \item Clay loam
#'   \item Silty clay loam
#'   \item Sandy clay
#'   \item Silty clay
#'   \item Clay
#' }
#' @return A list of soil parameters:
#' \describe{
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3), repeated for each soil layer}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3), repeated for each soil layer}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3), repeated for each soil layer}
#'   \item{Vq}{Quartz fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vm}{Mineral fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vo}{Organic fraction by volume (m^3 m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{Mc}{Mass fraction of clay (kg kg^-1), repeated for each soil layer}
#'   \item{b}{Campbell soil water retention parameter (unitless), repeated for each soil layer}
#'   \item{psi_e}{Air-entry water potential (J/kg, numerically equal to kPa), repeated for each soil layer}
#'   \item{gref}{Ground shortwave reflectance (unitless)}
#'   \item{groundem}{Ground longwave emissivity (unitless)}
#'   \item{grefPAR}{Ground reflectance for photosynthetically active radiation (unitless)}
#'   \item{slope}{Surface slope (degrees)}
#'   \item{aspect}{Surface aspect (degrees)}
#'   \item{nLayers}{Number of soil layers}
#'   \item{totalDepth}{Total soil depth (m)}
#'   \item{FreeDrain}{Logical indicating whether lower boundary layer is free draining}
#' }
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
createsoilc <- function(soiltype = "Clay loam", nlayers = 15, totalDepth = 3, slope = 0, aspect = 180,
                        surface_organicmu = 3, FreeDrain = TRUE) {
  s <- which(newsoilparamstable$Soil.type == soiltype)
  # Campbell parameters, selected by name
  campbellcols <- c("Smax", "Smin", "Ksat", "Vq", "Vm", "Vo", "Mc", "b", "psi_e")
  v <- as.numeric(newsoilparamstable[s, campbellcols])
  names(v) <- campbellcols
  soilc <- as.list(v)
  for (nm in campbellcols) {
    soilc[[nm]] <- rep(v[[nm]], nlayers + 1)
  }
  # Adjust organic
  mu <- rev(geometricCpp(nlayers, surface_organicmu))[2:(nlayers + 2)]  # multiplication factor for organic
  soilc$Vo <- soilc$Vo * mu
  dif <- soilc$Vo - v[["Vo"]]
  sm <- v[["Vq"]] + v[["Vm"]]
  mu <- (sm - dif) / sm # multiplication factor for other soil constituents
  soilc$Vm <- soilc$Vm * mu
  soilc$Vq <- soilc$Vq * mu
  soilc$Mc <- soilc$Mc * mu
  soilc$gref <- 0.15
  soilc$groundem <- 0.97
  soilc$grefPAR <- 0.75 * soilc$gref
  soilc$slope <- slope
  soilc$aspect <- aspect
  soilc$nLayers <- nlayers
  soilc$totalDepth <-  totalDepth
  soilc$FreeDrain <- FreeDrain
  return(soilc)
}
#' Estimate atmospheric CO2 concentration from year
#'
#' Estimates atmospheric CO2 concentration (ppm) from calendar year using
#' polynomial fits to historical observations and projections.
#'
#' @param year Numeric. Calendar year.
#'
#' @details
#' The function uses piecewise quadratic relationships fitted to atmospheric
#' CO2 data for different periods. Estimates are valid for years between
#' 1750 and 2100. Values outside this range produce an error.
#'
#' @return Numeric. Estimated atmospheric CO2 concentration (ppm).
#'
#' @export
Cafromyear <- function(year) {
  if (year > 2100) {
    stop("Cannot estimate CO2 for year > 2100\n")
  } else if (year > 2026) {
    Ca <- -0.011905 * year^2 + 51.571429 * year - 55191.666667
  } else if (year > 1962) {
    Ca <- 0.013401 * year^2 - 51.718345 * year + 50201.979441
  } else if (year >= 1750) {
    Ca <- 0.001040 * year^2 - 3.677483 * year + 3528.996720
  } else {
    stop("Cannot estimate CO2 for year < 1750\n")
  }
  return(Ca)
}
#' Convert standard meteorological observations to above-canopy forcing
#'
#' Converts meteorological observations measured at height `zin`
#' (typically ~2 m above ground) to height `zout` so they can be used as
#' meteorological forcing when vegetation height exceeds the measurement
#' height and the observations therefore lie within the canopy rather than
#' above it.
#'
#' @param climdata data.frame of weather data matching the format of the inbuilt
#' dataset `climdata`
#' @param zin Numeric. Height (m) at which the meteorological variables
#' were measured.
#' @param zout Numeric. Height (m) used as meteorological forcing.
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param SoilTempIni Optional numeric vector of initial soil temperatures.
#' If `NA`, a default profile is generated.
#' @param ThetaIni Optional numeric vector of initial volumetric soil water
#'   contents. If `NA`, a default profile is generated.
#' @param CO2ppm Optional numeric. Atmospheric CO₂ concentration (ppm).
#' If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#' (°C) used to initialise the soil temperature profile. If `NA`, it is
#' estimated from `climdata$temp`.
#' @param accuracy Character. Controls the trade-off between model run time and
#' numerical accuracy. One of "high" (the default), "very high", "moderate" or "fast".
#' @param maxiter Optional integer. Maximum number of iterations within each
#'   hour. Set by `accuracy` if not supplied.
#' @param tolerance Optional numeric. Convergence tolerance (K) for air
#'   temperature. Set by `accuracy` if not supplied.
#'
#' @details
#' Standard meteorological observations are typically made about 2 m above the
#' ground. If vegetation height exceeds the measurement height, these observations
#' lie within the canopy. The microclimate model requires meteorological forcing
#' above the canopy, so the observations are adjusted from height `zin` to height
#' `zout` before being used to force the model.
#'
#' @return A data.frame identical to `climdata` but with `temp`, `relhum`,
#' `windspeed`, and `pres` converted to height `zout`.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
weatherhgt_adjust <- function(climdata, zin, zout, lat, long, SoilTempIni = NA, ThetaIni = NA,
                              CO2ppm = NA, boundaryT = NA, accuracy = "high", maxiter = NULL, tolerance = NULL) {
  st <- .accuracy_settings(accuracy, tolerance, maxiter)
  n <- dim(climdata)[1]
  if (class(boundaryT) == "logical") {
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  if (class(SoilTempIni) == "logical") {
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(15, boundaryT - surfaceT) + surfaceT)[1:16]
  }
  mutheta <- seq(0.6, 1, length.out =  16)
  if (class(ThetaIni) == "logical")  ThetaIni <- 0.46 * mutheta
  vegp <-  createvegp(vegtype = "C3")
  vegp$h <- 0.12
  paii <- PAIgrass(0.12, 10)
  vegp$pai <- sum(paii)
  Lfrac <- seq(0.7, 0.9, length.out = length(paii))
  soilc <- createsoilc(soiltype = "Clay loam")
  climdata <- runchecks(climdata, NA, soilc, lat, long)$climdata
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  .set_solver(st)
  hgtvars <- WeatherhgtCpp2(obstime, climdata, soilc, vegp, paii, Lfrac, zin, zout, lat, long,
                            SoilTempIni, ThetaIni, CO2ppm, st$maxiter, st$tolerance)
  clim2 <- climdata
  clim2$temp <- hgtvars$new_temp
  clim2$relhum <- hgtvars$new_relhum
  clim2$windspeed <- hgtvars$new_windspeed
  clim2$pres <- hgtvars$new_pressure
  return(clim2)
}
# Lowest height (m) at which weather forcing can be measured for a canopy:
# d + c (h - d), where c is the factor by which heat diffusivity at canopy top
# exceeds its surface-layer value. Zero for bare ground.
.reference_height <- function(vegp, zm = 0.004) {
  if (class(vegp) == "logical") return(0)
  referenceheightCpp(vegp$h, vegp$pai, zm)
}
# Adjusts weather measured too close to a canopy to the height the model needs, rounded up
# to a whole metre. Returns the weather and the height it then refers to.
.forcing_at_reference <- function(climdata, vegp, zref, lat, long, CO2ppm = NA, boundaryT = NA, zm = 0.004,
                                 st = .accuracy_settings()) {
  zreq <- .reference_height(vegp, zm)
  if (zref >= zreq - 1e-6) return(list(climdata = climdata, zref = zref))
  zreq <- ceiling(zreq - 1e-6)   # adjusted to a whole number of metres
  message(sprintf("Weather at %.4g m is too close to a canopy of this height and density; translated to %.4g m.",
                  zref, zreq))
  # forcing already checked by the caller
  climdata <- suppressWarnings(weatherhgt_adjust(climdata, zref, zreq, lat, long,
                                                 CO2ppm = CO2ppm, boundaryT = boundaryT, accuracy = st))
  list(climdata = climdata, zref = zreq)
}
# A canopy shorter than this multiple of the soil roughness length is run as bare ground,
# as is one with no height or no plant area.
SHORT_CANOPY_MULTIPLE <- 2
.asbareifshort <- function(vegp, zm = 0.004) {
  if (class(vegp) == "logical") return(vegp)
  if (vegp$h <= 0 || vegp$pai <= 0) return(NA)   # no canopy
  if (vegp$h < SHORT_CANOPY_MULTIPLE * zm) {
    warning("Canopy height is below ", SHORT_CANOPY_MULTIPLE, " times the soil roughness length (",
            signif(SHORT_CANOPY_MULTIPLE * zm * 1000, 3), " mm), so the vegetation is no taller than the soil is rough. ",
            "Running as bare ground.
", call. = FALSE)
    return(NA)
  }
  vegp
}
# Plant area given per unit horizontal area, re-expressed per unit area of the sloping ground.
.pai_on_slope <- function(vegp, paii, soilc, paiFlat) {
  if (paiFlat && class(vegp) != "logical") {
    mu <- cos(soilc$slope * pi / 180)
    vegp$pai <- vegp$pai * mu
    paii <- paii * mu
  }
  list(vegp = vegp, paii = paii)
}
#' Runs model quickly for a period of time to get an initial soil water profile
#'
#' @param climdata data.frame of weather conditions that uses the same naming conventions
#'   as in the inbuilt dataset `climdata`.
#' @param vegp Vegetation parameter list as returned by [createvegp()]. If `NA`,
#'   bare ground is assumed. Canopies shorter than twice the soil roughness
#'   length are also run as bare ground.
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param Lfrac Fraction of plant area that is living vegetation: a single value
#'   or one for each canopy layer. Required unless `vegp` is `NA`.
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) of the meteorological forcing data.
#'   Default is 2 (see Details).
#' @param zmr Numeric. Roughness length (m) of the soil surface.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#'   If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#'   used to initialise the soil temperature profile. If `NA`, it is estimated
#'   from `climdata$temp`, which usually requires roughly a full year of data.
#' @param maxiter Integer. Maximum number of iterations used to achieve convergence
#'   within each time step.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#' @param accuracy Character. Controls the trade-off between model run time and
#'   numerical accuracy. One of "high" (the default), "very high",
#'   "moderate" or "fast". Used only where weather is adjusted to a greater
#'   height.
#' @param paiFlat Logical. If `TRUE`, plant area index is given per unit
#'   horizontal area. If `FALSE` (the default), it is given per unit area of the
#'   sloping ground surface. Makes no difference on flat ground.
#' @details
#' The soil water model is run over the supplied climate data, and the final
#' soil water profile is used as the initial condition. Where antecedent climate
#' data are not available, use data from the corresponding period in the current
#' year ending immediately before the time of interest.
#'
#' @return Numeric vector of initial volumetric soil water contents in each layer
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
InitailizeWater <- function(climdata, vegp, soilc, lat, long, zref = 2, zmr = 0.004,
                            CO2ppm = NA, boundaryT = NA, maxiter = 100, C3 = TRUE,
                            Lfrac = NA, accuracy = "high", paiFlat = FALSE) {
  st <- .accuracy_settings(accuracy)
  chk <- runchecks(climdata, vegp, soilc, lat, long, Lfrac = Lfrac)
  climdata <- chk$climdata; soilc <- chk$soilc
  vegp <- .pai_on_slope(vegp, NA, soilc, paiFlat)$vegp
  soilc$zm <- zmr   # soil roughness length
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  vegp <- .asbareifshort(vegp, zmr)
  if (class(vegp) == "logical") {
    blm <- BigLeafBareCpp(obstime, climdata, soilc, zref, zmr, lat, long,
                          boundaryT, maxiter)
  } else {
    if (length(Lfrac) == 1 && is.na(Lfrac)) {
      stop("Lfrac must be supplied when vegp describes a canopy\n")
    }
    if (length(Lfrac) == 1) Lfrac <- rep(Lfrac, 10)
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    fr <- .forcing_at_reference(climdata, vegp, zref, lat, long, CO2ppm, boundaryT, zmr, st)
    climdata <- fr$climdata; zref <- fr$zref
    blm <- BigLeafCpp2(obstime, climdata, soilc, vegp, Lfrac, zref, CO2ppm,
                       lat, long, boundaryT, maxiter, C3)
  }
  return(blm$SoilThetaIni)
}


#' Return a vertical microclimate profile for a specified hour
#'
#' Runs the microclimate model up to the specified hour and then returns the resulting
#' vertical profile for the hour, through the soil, canopy, and optionally the air above the
#' canopy. Profiles can be returned for bare ground or vegetated conditions.
#'
#' @param climdata data.frame of weather conditions that uses the same naming conventions
#' as in the inbuilt dataset `climdata`.
#' @param hr Integer. Hour to return, indexed from 1 to `nrow(climdata)`.
#' @param vegp Vegetation parameter list as returned by [createvegp()]. If `NA`,
#'   bare ground is assumed. Canopies shorter than twice the soil roughness
#'   length are also run as bare ground.
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param paii Numeric vector of plant area index values for each canopy layer
#' as returned by [PAIgeometry()] or [PAIgrass()] (see details).
#' @param Lfrac Numeric vector giving the fraction of plant area that is living vegetation
#' in each canopy layer (see details).
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) of the meteorological forcing data.
#'   Default is 2 (see details).
#' @param SoilTempIni Optional numeric vector of initial soil temperatures.
#'   If `NA`, a default profile is generated.
#' @param ThetaIni Optional numeric vector of initial volumetric soil water
#'   contents. If `NA`, a default profile is generated.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#'   If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#'   used to initialise the soil temperature profile. If `NA`, it is estimated
#'   from `climdata$temp`, which requires roughly a full year of data.
#' @param zm Numeric. Roughness length (m) of the soil surface.
#' @param maxiter Optional integer. Maximum number of iterations within each
#'   hour. Set by `accuracy` if not supplied.
#' @param tolerance Optional numeric. Convergence tolerance (K) for air
#'   temperature. Set by `accuracy` if not supplied.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#' @param plotout Logical. If `TRUE`, plot the requested profile.
#' @param accuracy Character. One of `"high"` (the default), `"very high"`,
#'   `"moderate"` or `"fast"`. See details of [RunMicro()].
#' @param varn Character. Variable to plot. One of `"temp"`, `"relhum"`,
#'   or `"radiation"`.
#' @param paiFlat Logical. If `TRUE`, plant area index is given per unit
#'   horizontal area. If `FALSE` (the default), it is given per unit area of the
#'   sloping ground surface. Makes no difference on flat ground.
#'
#' @details
#' The length of `paii` and `Lfrac` must be the same, and `sum(paii)` must equal
#' `vegp$pai`.
#'
#' The model is run up to hour `hr`, then the full profile for that hour is
#' returned. If `vegp` is supplied, a vegetated profile is calculated;
#' otherwise a bare-ground profile is returned.
#'
#' The model requires weather data from above the canopy, typically at least
#' 1.3 to 1.9 times canopy height for a dense canopy and somewhat lower for a
#' sparse one. If `zref` is too low, the weather data are adjusted automatically
#' to the required height. If using the function repeatedly for different
#' hours, it is quicker to adjust `climdata` once beforehand using
#' [weatherhgt_adjust()].
#'
#' If `plotout = TRUE`, the selected profile is plotted. `"temp"` plots
#' temperature, `"relhum"` relative humidity and soil moisture, and
#' `"radiation"` shortwave and longwave radiation.
#'
#' @return A list of model outputs for the hour, including:
#' \describe{
#'   \item{z}{Heights (m) of the canopy layers or, for bare ground, of the profile returned.}
#'   \item{tair, rh, uz}{Air temperature (°C), relative humidity (\%) and wind speed (m s^-1) at heights `z`.}
#'   \item{tleaf}{Leaf temperature (°C) of each canopy layer. Vegetated runs only.}
#'   \item{leafwet}{Fraction of the leaf surface that is wet in each canopy layer. Vegetated runs only.}
#'   \item{Rdirdown, Rdifdown, Rswup, Rlwdown, Rlwup}{Shortwave and longwave radiation (W m^-2).}
#'   \item{zb, Soiltemp, thetas}{Depths (m) of the soil nodes, and soil temperature (°C) and volumetric water content at them.}
#'   \item{H, L, G}{Sensible, latent and ground heat fluxes (W m^-2).}
#'   \item{za, tair_above, rh_above}{Heights (m) above the canopy, and air temperature and relative humidity at them. Vegetated runs only.}
#'   \item{tcanopy, rhcanopy, ucanopy}{Air temperature, relative humidity and wind speed at canopy top. Vegetated runs only.}
#'   \item{zref}{Height (m) of the weather data used. Vegetated runs only.}
#'   \item{iters, error}{Number of iterations and final convergence error (K).}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
return_profile <- function(climdata, hr, vegp, soilc, paii, Lfrac, lat, long, zref = 2,
                           SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA, zm = 0.004,
                           maxiter = NULL, tolerance = NULL, C3 = TRUE, plotout = TRUE, varn = "temp",
                           accuracy = "high", paiFlat = FALSE) {
  st <- .accuracy_settings(accuracy, tolerance, maxiter)
  chk <- runchecks(climdata, vegp, soilc, lat, long, paii, Lfrac)
  climdata <- chk$climdata; soilc <- chk$soilc
  ps <- .pai_on_slope(vegp, paii, soilc, paiFlat); vegp <- ps$vegp; paii <- ps$paii
  soilc$zm <- zm   # soil roughness length
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  nlay <- soilc$nLayers
  if (class(SoilTempIni) == "logical") {
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(nlay + 1)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  vegp <- .asbareifshort(vegp, zm)
  if (class(vegp) == "logical") {  # Bare ground assumed
    z <- seq(0, zref, length.out = 1000)
    mout <- profilebareR(hr - 1, obstime, climdata, soilc, z, zref, lat, long, SoilTempIni, ThetaIni,
                         zm, st$maxiter, st$tolerance)
    zref2 <- zref
  }  else {
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    # A canopy given as one layer is returned as its two halves
    two <- .at_least_two_layers(paii, Lfrac); paii <- two$paii; Lfrac <- two$Lfrac
    # The hours before the one returned are run on at most twenty layers;
    # a canopy given on very few layers is run as given.
    n <- length(paii)
    nlead <- min(20, st$canopylayers)
    if (n >= 4) {
      paii20 <- spline(x = seq_len(n), y = paii, xout = seq(1, n, length.out = nlead),
                       method = "natural")$y
      paii20 <- pmax(paii20, 0)
      mu <- sum(paii) / sum(paii20)
      paii20 <- paii20 * mu
      Lfrac20 <- spline(x = seq_len(n), y = Lfrac, xout = seq(1, n, length.out = nlead),
                        method = "natural")$y
      Lfrac20 <- pmin(pmax(Lfrac20, 0), 1)
    } else {
      paii20 <- paii
      Lfrac20 <- Lfrac
    }
    # Forcing measured too close to the canopy is translated up to the height the model needs
    fr <- .forcing_at_reference(climdata, vegp, zref, lat, long, CO2ppm, boundaryT, zm, st)
    climdata2 <- fr$climdata; zref2 <- fr$zref
    .set_solver(st)
    mout <- profileR(hr - 1, obstime, climdata2, soilc, vegp, paii20, paii,
                     Lfrac20, Lfrac, zref2, CO2ppm, lat, long, SoilTempIni, ThetaIni,
                     st$maxiter, st$tolerance, -1, 1.25, C3)
  }
  if (class(vegp) != "logical") {
    na <- max(2, round((2 / vegp$h) * length(paii), 0))
    za <- seq(vegp$h, zref2, length.out = na)
    # get canopy top values
    Th <- mout$tcanopy
    Rh <- mout$rhcanopy
    tabove <- 0
    rhabove <- 0
    for (i in 1:length(za)) {
      tabove[i] <- Tabove(za[i], zref2, Th, climdata2$temp[hr], vegp$h, vegp$pai, mout$LL, zm, mout$fade)
      rhabove[i] <- RHabove(za[i], zref2, Rh, Th, climdata2$temp[hr], tabove[i], climdata2$relhum[hr], vegp$h, vegp$pai, mout$LL, zm, mout$fade)
    }
    mout$za <- za
    mout$tair_above <- tabove
    mout$rh_above <- rhabove
    # the height the forcing was referred to, and the wind there, for the hour returned
    mout$zref <- zref2
    mout$uref <- climdata2$windspeed[hr]
    # plant area index as the model used it
    mout$pai <- vegp$pai
  }
  if (plotout) {
    tp <- paste0("Hour: ", hr)
    if (varn == "temp") {
      nb <- length(mout$Soiltemp)
      if (class(vegp) != "logical") {
        zall <- c(rev(-mout$zb[1:nb]), mout$z, za)
        tair <- c(rev(mout$Soiltemp), mout$tair, tabove)
        xmn <- floor(min(tair, mout$tleaf))
        xmx <- ceiling(max(tair, mout$tleaf))
      } else {
        zall <- c(rev(-mout$zb[1:nb]), mout$z)
        tair <- c(rev(mout$Soiltemp), mout$tair)
        xmn <- floor(min(tair))
        xmx <- ceiling(max(tair))
      }
      plot(zall ~ tair, type = "l", xlim = c(xmn, xmx),  ylim = c(min(zall), zref2),
           col = "blue", xlab = "Temperature", ylab = "Height(m)", lwd = 2, main = tp)
      if (class(vegp) != "logical") {
        par(new = T)
        plot(mout$z ~ mout$tleaf, type = "l", xlim = c(xmn, xmx),  ylim = c(min(zall), zref2),
             col = "darkgreen", xlab = "", ylab = "", lwd = 2)
        abline(h = vegp$h, lty = 2)
      }
      abline(h = 0, lty = 2)
    } else if (varn == "relhum") {
      nb <- length(mout$thetas)
      soilrh <- (mout$thetas / soilc$Smax) * 100
      if (class(vegp) != "logical") {
        zall <- c(rev(-mout$zb[1:nb]), mout$z, za)
        rh <- c(rev(soilrh), mout$rh, rhabove)
      } else {
        zall <- c(rev(-mout$zb[1:nb]), mout$z)
        rh <- c(rev(soilrh), mout$rh)
      }
      plot(zall ~ rh, type = "l", xlim = c(0, 100),  ylim = c(min(zall), zref2),
           col = "blue", xlab = "Relative Humidity (%)", ylab = "Height(m)", lwd = 2, , main = tp)
      abline(h = 0, lty = 2)
      if (class(vegp) != "logical") abline(h = vegp$h, lty = 2)
    } else if (varn == "radiation") {
      xmn1 <- floor(min(mout$Rswup, mout$Rdirdown, mout$Rdifdown))
      xmx1 <- ceiling(max(mout$Rswup, mout$Rdirdown, mout$Rdifdown))
      xmn2 <- floor(min(mout$Rlwup, mout$Rlwdown))
      xmx2 <- ceiling(max(mout$Rlwup, mout$Rlwdown))
      par(mfrow = c(1, 2))
      plot(mout$z ~ mout$Rdirdown, type = "l", xlim = c(xmn1, xmx1),  ylim = c(0, vegp$h),
           col = "red", xlab = "Shortwave radiation", ylab = "Height(m)", lwd = 2, , main = tp)
      par(new = T)
      plot(mout$z ~ mout$Rdifdown, type = "l", xlim = c(xmn1, xmx1),  ylim = c(0, vegp$h),
           col = "blue", xlab = "", ylab = "", lwd = 2)
      par(new = T)
      plot(mout$z ~ mout$Rswup, type = "l", xlim = c(xmn1, xmx1),  ylim = c(0, vegp$h),
           col = "black", xlab = "", ylab = "", lwd = 2)
      plot(mout$z ~ mout$Rlwdown, type = "l", xlim = c(xmn2, xmx2),  ylim = c(0, vegp$h),
           col = "red", xlab = "Longwave radiation", ylab = "Height(m)", lwd = 2)
      par(new = T)
      plot(mout$z ~ mout$Rlwup, type = "l", xlim = c(xmn2, xmx2),  ylim = c(0, vegp$h),
           col = "black", xlab = "", ylab = "", lwd = 2)
    } else (stop("varn not reconised\n"))
  }
  mout$thetas <- .reportedtheta(mout$thetas, soilc)
  return(mout)
}
#' Run the microclimate model at a specified height
#'
#' Runs the microclimate model for each time step in `climdata` and returns
#' the simulated microclimate at the requested height `reqhgt`.
#'
#' @param climdata data.frame of weather conditions that uses the same naming conventions
#' as in the inbuilt dataset `climdata`.
#' @param reqhgt Numeric. Height (m) at which output is required. Positive
#' values are above ground, `0` is the soil surface, and negative values are
#' below ground.
#' @param vegp Vegetation parameter list as returned by [createvegp()]. If `NA`,
#'   bare ground is assumed.
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param paii Numeric vector of plant area index values for each canopy layer
#' as returned by [PAIgeometry()] or [PAIgrass()] (see details).
#' @param Lfrac Numeric vector giving the fraction of plant area that is living vegetation
#' in each canopy layer (see details).
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) of the meteorological forcing data.
#' Default is 2 (see details).
#' @param SoilTempIni Optional numeric vector of initial soil temperatures.
#' If `NA`, a default profile is generated.
#' @param ThetaIni Optional numeric vector of initial volumetric soil water
#'   contents. If `NA`, a default profile is generated.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#' If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#' used to initialise the soil temperature profile. If `NA`, it is estimated
#' from `climdata$temp`, which requires roughly a full year of data.
#' @param zm Numeric. Roughness length (m) of the soil surface.
#' @param maxiter Optional integer. Maximum number of iterations within each
#'   hour. Set by `accuracy` if not supplied.
#' @param tolerance Optional numeric. Convergence tolerance (K) for air
#'   temperature. Set by `accuracy` if not supplied.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#' @param accuracy Character. Controls the trade-off between model run time and
#' numerical accuracy. One of `"high"` (the default), `"very high"`,
#' `"moderate"` or `"fast"`. `"very high"` is slower and more accurate, while
#' `"moderate"` and `"fast"` progressively favour speed over accuracy. `"fast"`
#' is intended primarily for applications where the model is run repeatedly
#' over many grid cells.
#' @param paiFlat Logical. If `TRUE`, plant area index is given per unit
#'   horizontal area. If `FALSE` (the default), it is given per unit area of the
#'   sloping ground surface. Makes no difference on flat ground.
#'
#' @details
#' If `vegp` is `NA`, bare-ground conditions are assumed and the bare-ground
#' model is used. If `reqhgt <= 0`, output is treated as below-ground and soil
#' temperature and soil water content are returned instead of air temperature
#' and relative humidity.
#'
#' @return A data.frame containing model outputs for each time step.
#'
#' For bare-ground runs, the returned columns are:
#' \describe{
#'   \item{year, month, day, hour}{Date-time components of the simulation.}
#'   \item{Rdirdown}{Direct shortwave radiation downward.}
#'   \item{Rdifdown}{Diffuse shortwave radiation downward.}
#'   \item{Rswup}{Upward shortwave radiation.}
#'   \item{Rlwdown}{Downward longwave radiation.}
#'   \item{Rlwup}{Upward longwave radiation.}
#'   \item{tair}{Air temperature at the requested height.}
#'   \item{tground}{Ground surface temperature.}
#'   \item{relhum}{Relative humidity at the requested height.}
#'   \item{windspeed}{Wind speed at the requested height.}
#'   \item{H}{Sensible heat flux.}
#'   \item{L}{Latent heat flux.}
#'   \item{G}{Ground heat flux.}
#'   \item{iters}{Number of iterations required for convergence.}
#'   \item{error}{Final convergence error.}
#' }
#'
#' For vegetated runs, the returned data.frame contains all of the columns
#' listed above, plus:
#' \itemize{
#'   \item `tleaf`: leaf temperature at the requested height within the canopy,
#'   or `-999.99` where not applicable.
#'   \item `leafwet`: fraction of the leaf surface that is wet, or `-999.99`
#'   where not applicable.
#'   \item `Evaporation`: soil evaporation.
#'   \item `Transpiration`: plant transpiration.
#' }
#'
#' Both also return `Evaporation` from the soil (mm), `surfwater`, the volumetric
#' water content of the soil surface, `runoff` from the soil
#' surface (mm), `drainage` from the base of the soil (mm) and `witers`, the
#' number of iterations of the soil water model.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
RunMicro <- function(climdata, reqhgt, vegp, soilc, paii, Lfrac, lat, long, zref = 2,
                     SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA, zm = 0.004,
                     maxiter = NULL, tolerance = NULL, C3 = TRUE, accuracy = "high", paiFlat = FALSE) {
  st <- .accuracy_settings(accuracy, tolerance, maxiter)
  chk <- runchecks(climdata, vegp, soilc, lat, long, paii, Lfrac)
  climdata <- chk$climdata; soilc <- chk$soilc
  ps <- .pai_on_slope(vegp, paii, soilc, paiFlat); vegp <- ps$vegp; paii <- ps$paii
  soilc$zm <- zm   # soil roughness length
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  # Soil layers used for the run, if the accuracy level allows fewer than were given.
  cs <- .coarsen_soil(soilc, st$soillayers)
  if (!is.null(cs)) {
    if (class(SoilTempIni) != "logical") SoilTempIni <- stats::approx(cs$zu, SoilTempIni, xout = cs$zc, rule = 2)$y
    if (class(ThetaIni) != "logical") ThetaIni <- stats::approx(cs$zu, ThetaIni, xout = cs$zc, rule = 2)$y
    soilc <- cs$soilc
  }
  nlay <- soilc$nLayers
  if (class(SoilTempIni) == "logical") {
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(nlay + 1)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  vegp <- .asbareifshort(vegp, zm)
  if (class(vegp) == "logical") {
    mout <- RunBareR(reqhgt, obstime, climdata, soilc, zref, lat, long, SoilTempIni,
                     ThetaIni, zm, st$maxiter, st$tolerance)
  } else {
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    fr <- .forcing_at_reference(climdata, vegp, zref, lat, long, CO2ppm, boundaryT, zm, st)
    climdata <- fr$climdata; zref <- fr$zref
    # Canopy layers used for the run, if the accuracy level allows fewer than were given.
    cc <- .coarsen_canopy(paii, Lfrac, st$canopylayers)
    if (!is.null(cc)) { paii <- cc$paii; Lfrac <- cc$Lfrac }
    two <- .at_least_two_layers(paii, Lfrac); paii <- two$paii; Lfrac <- two$Lfrac
    .set_solver(st)
    mout <- RunModelR(reqhgt, obstime, climdata, soilc, vegp, paii, Lfrac, zref, CO2ppm, lat, long, SoilTempIni,
                      ThetaIni, st$maxiter, st$tolerance, -1, 1.25, C3)
  }
  # One warning for the run if hours ended far from settled, or if the soil's water did not
  # balance in more than 1% of hours.
  .warn_unsettled(mout$error, st)
  stuck <- sum(mout$witers >= max(st$maxiter, 50))
  if (stuck > 0.01 * length(mout$witers)) warning(sprintf("The soil's water did not balance in %d of %d hours; soil water there is approximate", stuck, length(mout$witers)))
  mout$surfwater <- .reportedtheta(mout$surfwater, soilc)
  # Rename variables if below ground
  if (reqhgt <= 0) {
    mout$tsoil <- mout$tair; mout$tair <- NULL
    mout$soilwater <- .reportedtheta(mout$relhum, soilc); mout$relhum <- NULL
  }
  return(mout)
}
#' Run the full microclimate model through time
#'
#' Runs the full microclimate model for all time steps in `climdata` and returns
#' simulated microclimate profiles for all heights above ground and all soil
#' nodes below ground.
#'
#' @param climdata data.frame of weather conditions using the same naming
#' conventions as the inbuilt dataset `climdata`.
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param vegp Vegetation parameter list as returned by [createvegp()]. If `NA`,
#'   bare ground is assumed. Canopies shorter than twice the soil roughness
#'   length are also run as bare ground.
#' @param paii Numeric vector of plant area index values for each canopy layer,
#' ordered from bottom to top.
#' @param Lfrac Numeric vector giving the fraction of plant area in each canopy
#' layer that is living leaf tissue.
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) of the meteorological forcing data.
#' Default is 2.
#' @param SoilTempIni Optional numeric vector of initial soil temperatures
#' (degrees C), ordered from the soil surface downward and including the lower
#' boundary node. If `NA`, a default profile is generated.
#' @param ThetaIni Optional numeric vector of initial soil water contents,
#' ordered from the soil surface downward and including the lower boundary
#' node. If `NA`, a default profile is generated.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm). If `NA`,
#' it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual temperature used
#' to initialise the soil temperature profile. If `NA`, it is estimated from
#' `climdata$temp`, which requires approximately a full year of data.
#' @param zm Numeric. Roughness length (m) of the soil surface.
#' @param maxiter Optional integer. Maximum number of iterations within each
#'   hour. Set by `accuracy` if not supplied.
#' @param tolerance Optional numeric. Convergence tolerance (K) for air
#'   temperature. Set by `accuracy` if not supplied.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#' @param accuracy Character. One of `"high"` (the default), `"very high"`,
#'   `"moderate"` or `"fast"`. See details of [RunMicro()].
#' @param paiFlat Logical. If `TRUE`, plant area index is given per unit
#'   horizontal area. If `FALSE` (the default), it is given per unit area of the
#'   sloping ground surface. Makes no difference on flat ground.
#'
#' @details
#' The model requires weather data from above the canopy, typically at least
#' 1.3 to 1.9 times canopy height for a dense canopy and somewhat lower for a
#' sparse one. If `zref` is too low, the weather data are adjusted automatically
#' to the required height using [weatherhgt_adjust()] before the model is run.
#'
#' If `vegp` is `NA`, the model runs in bare-ground mode. Otherwise the full
#' canopy model is used with vegetation structure defined by `vegp`, `paii`,
#' and `Lfrac`.
#'
#' The model is run for every time step in `climdata` and returns profiles
#' through the canopy and soil rather than output at a single requested height.
#'
#' @return A list of numeric matrices giving simulated microclimate variables
#' for all time steps and all canopy heights (above ground) or soil nodes
#' (below ground). Rows correspond to time steps and columns correspond to
#' canopy heights or soil nodes.
#'
#' The heights (m) of the columns are returned as `z` and the depths (m) of the
#' soil nodes as `zsoil`. For vegetated runs `z` gives the middle of each canopy
#' layer; the air between canopy top and `zref` is returned separately at
#' heights `za`.
#'
#' For **bare-ground runs** (`vegp = NA`) the returned elements are:
#' \describe{
#'   \item{tair}{Air temperature profile (°C).}
#'   \item{relhum}{Relative humidity profile (\%).}
#'   \item{windspeed}{Wind speed profile (m s^-1).}
#'   \item{tsoil}{Soil temperature profile (°C).}
#'   \item{theta}{Volumetric soil water content profile (m^3 m^-3).}
#'   \item{runoff}{Runoff from the soil surface in each time step (mm). A vector.}
#'   \item{drainage}{Drainage from the base of the soil in each time step (mm). A vector.}
#'   \item{witers}{Number of iterations of the soil water model in each time step. A vector.}
#'   \item{iters}{Number of iterations in each time step. A vector.}
#'   \item{error}{Final convergence error in each time step (K). A vector.}
#' }
#'
#' For **vegetated runs** the returned elements are:
#' \describe{
#'   \item{Rdirdown}{Direct shortwave radiation downward (W m^-2).}
#'   \item{Rdifdown}{Diffuse shortwave radiation downward (W m^-2).}
#'   \item{Rswup}{Upward shortwave radiation (W m^-2).}
#'   \item{Rlwdown}{Downward longwave radiation (W m^-2).}
#'   \item{Rlwup}{Upward longwave radiation (W m^-2).}
#'   \item{tair}{Air temperature profile (°C).}
#'   \item{tleaf}{Leaf temperature profile (°C).}
#'   \item{leafwet}{Fraction of the leaf surface that is wet.}
#'   \item{relhum}{Relative humidity profile (\%).}
#'   \item{windspeed}{Wind speed profile (m s^-1).}
#'   \item{tsoil}{Soil temperature profile (°C).}
#'   \item{theta}{Volumetric soil water content profile (m^3 m^-3).}
#'   \item{runoff}{Runoff from the soil surface in each time step (mm). A vector.}
#'   \item{drainage}{Drainage from the base of the soil in each time step (mm). A vector.}
#'   \item{witers}{Number of iterations of the soil water model in each time step. A vector.}
#'   \item{tcanopy}{Air temperature at canopy top (°C). A vector.}
#'   \item{rhcanopy}{Relative humidity at canopy top (\%). A vector.}
#'   \item{ucanopy}{Wind speed at canopy top (m s^-1). A vector.}
#'   \item{tair_above}{Air temperature profile above the canopy (°C).}
#'   \item{relhum_above}{Relative humidity profile above the canopy (\%).}
#'   \item{windspeed_above}{Wind speed profile above the canopy (m s^-1).}
#'   \item{iters}{Number of iterations in each time step. A vector.}
#'   \item{error}{Final convergence error in each time step (K). A vector.}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
RunModelFull <- function(climdata, soilc, vegp, paii, Lfrac, lat, long, zref = 2,
                       SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA,
                       zm = 0.004, maxiter = NULL, tolerance = NULL, C3 = TRUE, accuracy = "high",
                       paiFlat = FALSE) {
  st <- .accuracy_settings(accuracy, tolerance, maxiter)
  chk <- runchecks(climdata, vegp, soilc, lat, long, paii, Lfrac)
  climdata <- chk$climdata; soilc <- chk$soilc
  ps <- .pai_on_slope(vegp, paii, soilc, paiFlat); vegp <- ps$vegp; paii <- ps$paii
  soilc$zm <- zm   # soil roughness length
  # Derive boundary temperature if not supplied
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  # Soil layers used for the run, if the accuracy level allows fewer than were given.
  soiluser <- soilc
  cs <- .coarsen_soil(soilc, st$soillayers)
  if (!is.null(cs)) {
    if (class(SoilTempIni) != "logical") SoilTempIni <- stats::approx(cs$zu, SoilTempIni, xout = cs$zc, rule = 2)$y
    if (class(ThetaIni) != "logical") ThetaIni <- stats::approx(cs$zu, ThetaIni, xout = cs$zc, rule = 2)$y
    soilc <- cs$soilc
  }
  nlay <- soilc$nLayers
  if (class(SoilTempIni) == "logical") {
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(1+nlay)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  vegp <- .asbareifshort(vegp, zm)
  if (class(vegp) == "logical") {  # Bare ground assumed
    z <- (c(1:20)/ 20) * zref
    mout <- RunBelowFullBare(obstime, climdata, soilc, z, zref, zm, lat, long,
                             SoilTempIni, ThetaIni, st$maxiter, st$tolerance)
    mout$z <- z
  } else {
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    fr <- .forcing_at_reference(climdata, vegp, zref, lat, long, CO2ppm, boundaryT, zm, st)
    climdata <- fr$climdata; zref <- fr$zref
    # Canopy layers used for the run, if the accuracy level allows fewer than were given.
    nuser <- length(paii)
    cc <- .coarsen_canopy(paii, Lfrac, st$canopylayers)
    if (!is.null(cc)) { paii <- cc$paii; Lfrac <- cc$Lfrac }
    two <- .at_least_two_layers(paii, Lfrac); paii <- two$paii; Lfrac <- two$Lfrac
    .set_solver(st)
    mout <- RunBelowFull(obstime, climdata, soilc, vegp, paii, Lfrac, zref, CO2ppm,
                         lat, long, SoilTempIni, ThetaIni, st$maxiter, st$tolerance,
                         -1, 1.25, C3)
    mout <- .canopy_to_user(mout, nuser, vegp$h, zm)
    # Air between canopy top and the height of the weather data, at about the spacing of the canopy layers
    mout$z <- (seq_len(nuser) - 0.5) * vegp$h / nuser
    mout$za <- seq(vegp$h, zref, length.out = min(max(2, ceiling((zref - vegp$h) * nuser / vegp$h) + 1), 50))
    above <- AboveCanopyCpp(mout$za, zref, mout$tcanopy, mout$rhcanopy, mout$ucanopy, climdata$temp,
                            climdata$relhum, climdata$windspeed, vegp$h, vegp$pai, mout$LL, zm, mout$fade)
    mout$tair_above <- above$tair; mout$relhum_above <- above$relhum; mout$windspeed_above <- above$windspeed
    mout$LL <- NULL; mout$fade <- NULL
  }
  mout$zsoil <- .soil_depths(soiluser)
  if (!is.null(cs)) {
    mout$tsoil <- .soil_to_user(mout$tsoil, cs)
    mout$theta <- .soil_to_user(mout$theta, cs)
  }
  # One warning for the run if hours ended far from settled, or if the soil's water did not
  # balance in more than 1% of hours.
  .warn_unsettled(mout$error, st)
  stuck <- sum(mout$witers >= max(st$maxiter, 50))
  if (stuck > 0.01 * length(mout$witers)) warning(sprintf("The soil's water did not balance in %d of %d hours; soil water there is approximate", stuck, length(mout$witers)))
  mout$theta <- .reportedtheta(mout$theta, soiluser)
  return(mout)
}
# Soil water content as returned to the user: no lower than the residual water content.
.reportedtheta <- function(theta, soilc) {
  if (is.null(theta)) return(theta)
  pmax(theta, min(soilc$Smin))
}
#' R wrapper for c++ expand_outputCpp function
#'
#' Linearly interpolates a Numeric Matrix for better visual representation
#'
#' @param mat matrix to expand.
#' @param nout number of rows of expanded matrix
#' @return an expanded matrix
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
expand_output <- function(mat, nout) {
  mout <- expand_outputCpp(mat, nout)
  return(mout)
}

