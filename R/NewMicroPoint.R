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
#'   \item{Dcrit}{Vapour pressure deficit at onset of strong stomatal limitation (kPa)}
#'   \item{alpha}{Photosynthetic quantum efficiency (mol CO2 mol^-1 PAR)}
#'   \item{Kxmx}{Maximum xylem hydraulic conductance (µmol m^-2 s^-1)}
#'   \item{hv}{Huber value (cm^2 m^-2)}
#'   \item{f0}{Ratio of ci/(ci − Γ*) under low vapour pressure deficit from Jacobs (1994) (unitless)}
#'   \item{fd}{Fraction of Vcmax converted to leaf dark respiration (unitless)}
#'   \item{psi50}{Water potential causing 50\% loss of maximum hydraulic conductance (MPa)}
#'   \item{apsi}{Xylem hydraulic conductance parameter (calculated if \code{< 0}) (unitless)}
#'   \item{len}{Leaf length (m)}
#'   \item{wid}{Leaf width (m)}
#'   \item{vegem}{Vegetation emissivity (unitless)}
#'   \item{mwft}{Maximum water film thickness (mm)}
#'   \item{pTAW}{Fraction of total available soil water that can be depleted before stress (unitless)}
#'   \item{rootskew}{Root distribution skew towards the top of the soil profile (unitless)}
#' }
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
#' @param totalDepth Numeric. Total soil depth (m). Internally this is multiplied
#'   by 1.5 to allow for a boundary layer.
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
#'   \item{n}{Pore size distribution parameter (unitless), repeated for each soil layer}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3), repeated for each soil layer}
#'   \item{Vq}{Quartz fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vm}{Mineral fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vo}{Organic fraction by volume (m^3 m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{Mc}{Coarse fragment or gravel fraction (m^3 m^-3), repeated for each soil layer}
#'   \item{rho}{Bulk density (kg m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{b}{Campbell soil water retention parameter (unitless), repeated for each soil layer}
#'   \item{psi_e}{Air-entry water potential (m), repeated for each soil layer}
#'   \item{gref}{Ground shortwave reflectance (unitless)}
#'   \item{groundem}{Ground longwave emissivity (unitless)}
#'   \item{grefPAR}{Ground reflectance for photosynthetically active radiation (unitless)}
#'   \item{slope}{Surface slope (degrees)}
#'   \item{aspect}{Surface aspect (degrees)}
#'   \item{nLayers}{Number of soil layers}
#'   \item{totalDepth}{Total soil depth used internally (m), after multiplying input depth by 1.5}
#'   \item{FreeDrain}{Logical indicating whether lower boundary layer is free draining}
#' }
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
createsoilc <- function(soiltype = "Clay loam", nlayers = 15, totalDepth = 2, slope = 0, aspect = 180,
                        surface_organicmu = 3, FreeDrain = TRUE) {
  # Adjust total depth to allow for boundary layer
  totalDepth <- 1.5 * totalDepth
  s <- which(newsoilparamstable$Soil.type == soiltype)
  v <- as.numeric(newsoilparamstable[s,2:13])
  nms <- names(newsoilparamstable)[2:13]
  soilc <- as.list(v)
  for (i in 1:12) {
    soilc[[i]] <- rep(v[i], nlayers + 1)
  }
  names(soilc) <- nms
  # Adjust organic
  mu <- rev(geometricCpp(nlayers, surface_organicmu))[2:(nlayers + 2)]  # multiplication factor for organic
  soilc$Vo <- soilc$Vo * mu
  dif <- soilc$Vo - v[8]
  sm <- v[6] + v[7]
  mu <- (sm - dif) / sm # multiplication factor for other soil constituents
  soilc$Vm <- soilc$Vm * mu
  soilc$Vq <- soilc$Vq * mu
  soilc$Mc <- soilc$Mc * mu
  mu <- seq(0.9, 1.1, length.out = nlayers + 1)  # multiplication factor for bulk density
  soilc$rho <- soilc$rho * mu
  soilc$gref <- 0.15
  soilc$groundem <- 0.97
  soilc$grefPAR <- 0.75 * soilc$gref
  soilc$slope <- slope
  soilc$aspect <- aspect
  soilc$nLayers <- nlayers
  soilc$totalDepth <-  totalDepth
  soilc$FreeDrain <- FreeDrain
  soilc$alpha <- NULL
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
#' contents. If `NA`, all layers are initialised to 0.419.
#' @param CO2ppm Optional numeric. Atmospheric CO₂ concentration (ppm).
#' If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#' (°C) used to initialise the soil temperature profile. If `NA`, it is
#' estimated from `climdata$temp`.
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
                              CO2ppm = NA, boundaryT = NA) {
  n <- dim(climdata)[1]
  if (class(boundaryT) == "logical") {
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  # Create SoilTemp and Theta vectors if not provided
  if (class(SoilTempIni) == "logical") {
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(15, boundaryT - surfaceT) + surfaceT)[1:16]
  }
  mutheta <- seq(0.6, 1, length.out =  16)
  if (class(ThetaIni) == "logical")  ThetaIni <- 0.46 * mutheta
  # Create vegp inputs
  vegp <-  createvegp(vegtype = "C3")
  vegp$h <- 0.12
  paii <- PAIgrass(0.12, 10)
  vegp$pai <- sum(paii)
  Lfrac <- seq(0.7, 0.9, length.out = length(paii))
  # Create soilc inputs
  soilc <- createsoilc(soiltype = "Clay loam")
  # Create obstime
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  # Estimate CO2
  if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  # Run model
  hgtvars <- WeatherhgtCpp2(obstime, climdata, soilc, vegp, paii, Lfrac, zin, zout, lat, long,
                            SoilTempIni, ThetaIni, CO2ppm)
  clim2 <- climdata
  clim2$temp <- hgtvars$new_temp
  clim2$relhum <- hgtvars$new_relhum
  clim2$windspeed <- hgtvars$new_windspeed
  clim2$pres <- hgtvars$new_pressure
  return(clim2)
}
#' Runs model quickly for a period of time to get an initial soil water profile
#'
#' @param climdata data.frame of weather conditions that uses the same naming conventions
#'   as in the inbuilt dataset `climdata`.
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#'   If `NA`, bare-ground conditions are assumed.
#' @param soilc Soil parameter list as returned by [createsoilc()].
#' @param paii Numeric vector of plant area index values for each canopy layer,
#'   as returned by [PAIgeometry()] or [PAIgrass()] (see Details).
#' @param Lfrac Single numeric value or numeric vector giving the fraction of plant area
#'   that is living vegetation in each canopy layer.
#' @param lat Numeric. Latitude in decimal degrees.
#' @param long Numeric. Longitude in decimal degrees.
#' @param zref Numeric. Height (m) of the meteorological forcing data.
#'   Default is 2 (see Details).
#' @param zmr Numeric. Roughness length for bare-ground calculations.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#'   If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#'   used to initialise the soil temperature profile. If `NA`, it is estimated
#'   from `climdata$temp`, which usually requires roughly a full year of data.
#' @param maxiter Integer. Maximum number of iterations used to achieve convergence
#'   within each time step.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#'
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
                            CO2ppm = NA, boundaryT = NA, maxiter = 100, C3 = TRUE) {
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
  if (class(vegp) == "logical") {
    blm <- BigLeafBareCpp(obstime, climdata, soilc, zref, zmr, lat, long,
                          boundaryT, maxiter)
  } else {
    if (length(Lfrac) == 1) Lfrac=rep(Lfrac, 10)
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
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
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' If `NA`, bare-ground conditions are assumed.
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
#'   contents in each layer. If `NA`, all plausible profile is computed
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#'   If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#'   used to initialise the soil temperature profile. If `NA`, it is estimated
#'   from `climdata$temp`, which requires roughly a full year of data.
#' @param zm Numeric. Roughness length for bare-ground calculations.
#' @param maxiter Integer. Maximum number of iterations used to achieve
#'   convergence.
#' @param tolerance Numeric. Convergence tolerance.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#' @param plotout Logical. If `TRUE`, plot the requested profile.
#' @param varn Character. Variable to plot. One of `"temp"`, `"relhum"`,
#'   or `"radiation"`.
#'
#' @details
#' The length of vectors `paii` and `Lfrac` must be the same and `sum(paii)` must
#' equal `vegp$pai`.
#'
#' The function first runs the model up to hour `hr`, then returns the full
#' profile for that hour. If `vegp` is provided, a vegetated profile is
#' calculated; otherwise a bare-ground profile is returned.
#'
#' Standard meteorological observations are often measured about 2 m above the
#' ground. If vegetation height exceeds this measurement height, the forcing
#' data lie within the canopy rather than above it. Because the microclimate
#' model requires above-canopy forcing, the weather data are first adjusted to
#' a height above the canopy before the profile is calculated. If using the function
#' repeatedly to generate height profiles for different hours, it is much quicker
#' to pre height-adjust `climdata` using [weatherhgt_adjust()].
#'
#' If `plotout = TRUE`, the selected profile is plotted. For `"temp"`, soil,
#' air, and leaf temperature profiles are plotted where available. For
#' `"relhum"`, soil moisture is converted to relative humidity below ground and
#' joined to the air humidity profile above ground. For `"radiation"`, vertical
#' profiles of shortwave and longwave radiation are plotted.
#'
#' @return A data.frame containing model outputs for each time step. For
#' bare-ground runs, the returned columns are:
#'
#' \describe{
#'   \item{year, month, day, hour}{Date-time components of the simulation.}
#'   \item{Rdirdown}{Direct shortwave radiation downward.}
#'   \item{Rdifdown}{Diffuse shortwave radiation downward.}
#'   \item{Rswup}{Upward shortwave radiation.}
#'   \item{Rlwdown}{Downward longwave radiation.}
#'   \item{Rlwup}{Upward longwave radiation.}
#'   \item{tair}{Air temperature at the requested height.}
#'   \item{tground}{Ground surface temperature.}
#'   \item{relhum}{Relative humidity.}
#'   \item{windspeed}{Wind speed.}
#'   \item{H}{Sensible heat flux.}
#'   \item{L}{Latent heat flux.}
#'   \item{G}{Ground heat flux.}
#'   \item{iters}{Number of iterations required for convergence.}
#'   \item{error}{Final convergence error.}
#' }
#'
#' For vegetated runs, the returned data.frame contains all of the columns
#' listed above, plus:
#'
#' \itemize{
#'   \item `tleaf`: leaf temperature.
#'   \item `Evaporation`: soil evaporation.
#'   \item `Transpiration`: plant transpiration.
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
return_profile <- function(climdata, hr, vegp, soilc, paii, Lfrac, lat, long, zref = 2,
                           SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA, zm = 0.004,
                           maxiter = 100, tolerance = 1e-2, C3 = TRUE, plotout = TRUE, varn = "temp") {
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  # Create SoilTemp and Theta vectors if not provided
  nlay <- soilc$nLayers
  if (class(SoilTempIni) == "logical") {
    boundaryT <- mean(climdata$temp)
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(nlay + 1)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  # Create obstime
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  if (class(vegp) == "logical") {  # Bare ground assumed
    z <- seq(0, zref, length.out = 1000)
    mout <- profilebareR(hr - 1, obstime, climdata, soilc, z, zref, lat, long, SoilTempIni, ThetaIni,
                         zm, maxiter, tolerance)
    zref2 <- zref
  }  else {
    # Estimate CO2
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    # Create vegp inputs
    n <- length(paii)
    paii20 <- spline(x = seq_len(n), y = paii, xout = seq(1, n, length.out = 20),
                     method = "natural")$y
    mu <- sum(paii) / sum(paii20)
    paii20 <- paii20 * mu
    Lfrac20 <- spline(x = seq_len(n), y = paii, xout = seq(1, n, length.out = 20),
                      method = "natural")$y
    # Check whether height adjustment is required
    if (vegp$h > (zref - 1)) {
      zref2 <- vegp$h + 2
      climdata2 <- weatherhgt_adjust(climdata, zref, zref2, lat, long, SoilTempIni, ThetaIni, CO2ppm, boundaryT)
    } else {
      zref2 <- zref
      climdata2 <- climdata
    }
    mout <- profileR(hr - 1, obstime, climdata2, soilc, vegp, paii20, paii,
                     Lfrac20, Lfrac, zref2, CO2ppm, lat, long, SoilTempIni, ThetaIni,
                     maxiter, tolerance, 0.25, 1.25, C3)
  }
  if (class(vegp) != "logical") {
    na <- round((2 / vegp$h) * length(paii), 0)
    za <- seq(vegp$h, zref2, length.out = na)
    # get canopy top values
    Th <- mout$tair[length(mout$tair)]
    Rh <- mout$rh[length(mout$rh)]
    tabove <- 0
    rhabove <- 0
    for (i in 1:length(za)) {
      tabove[i] <- Tabove(za[i], zref2, Th, climdata2$temp[hr], vegp$h, vegp$pai)
      rhabove[i] <- RHabove(za[i], zref2, Rh, Th, climdata2$temp[hr], tabove[i], climdata2$relhum[hr], vegp$h, vegp$pai)
    }
    mout$za <- za
    mout$tair_above <- tabove
    mout$rh_above <- rhabove
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
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' If `NA`, bare-ground conditions are assumed.
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
#' contents. If `NA`, all layers are initialised to 0.419.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm).
#' If `NA`, it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual air temperature
#' used to initialise the soil temperature profile. If `NA`, it is estimated
#' from `climdata$temp`, which requires roughly a full year of data.
#' @param zm Numeric. Roughness length for bare-ground calculations.
#' @param maxiter Integer. Maximum number of iterations used to achieve
#' convergence.
#' @param tolerance Numeric. Convergence tolerance.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#'
#' @details
#' Standard meteorological observations are typically made about 2 m above the
#' ground. If vegetation is taller than the measurement height, these data lie
#' within the canopy rather than above it. The microclimate model requires
#' above-canopy meteorological forcing, so when vegetation is present the model
#' uses the canopy structure defined by `vegp`, `paii` and `Lfrac` to simulate
#' conditions at the requested height `reqhgt`.
#'
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
#'   \item `Evaporation`: soil evaporation.
#'   \item `Transpiration`: plant transpiration.
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
RunMicro <- function(climdata, reqhgt, vegp, soilc, paii, Lfrac, lat, long, zref = 2,
                     SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA, zm = 0.004,
                     maxiter = 100, tolerance = 1e-2, C3 = TRUE) {
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  nlay <- soilc$nLayers
  # Create SoilTemp and Theta vectors if not provided
  if (class(SoilTempIni) == "logical") {
    boundaryT <- mean(climdata$temp)
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(nlay + 1)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  # Create obstime
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  if (class(vegp) == "logical") {
    mout <- RunBareR(reqhgt, obstime, climdata, soilc, zref, lat, long, SoilTempIni,
                     ThetaIni, zm, maxiter, tolerance)
  } else {
    # Estimate CO2
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    mout <- RunModelR(reqhgt, obstime, climdata, soilc, vegp, paii, Lfrac, zref, CO2ppm, lat, long, SoilTempIni,
                      ThetaIni, maxiter, tolerance, 0.25, 1.25, C3)
  }
  # Rename variables if below ground
  if (reqhgt <= 0) {
    mout$tsoil <- mout$tair; mout$tair <- NULL
    mout$soilwater <- mout$relhum; mout$relhum <- NULL
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
#' @param vegp Vegetation parameter list as returned by [createvegp()].
#' If `NA`, bare-ground conditions are assumed.
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
#' node. If `NA`, all layers are initialised to 0.419.
#' @param CO2ppm Optional numeric. Atmospheric CO2 concentration (ppm). If `NA`,
#' it is estimated from the year of the first observation using [Cafromyear()].
#' @param boundaryT Optional numeric. Approximate mean annual temperature used
#' to initialise the soil temperature profile. If `NA`, it is estimated from
#' `climdata$temp`, which requires approximately a full year of data.
#' @param zm Numeric. Roughness length used for bare-ground calculations.
#' @param maxiter Integer. Maximum number of iterations allowed for convergence.
#' @param tolerance Numeric. Convergence tolerance.
#' @param C3 Logical. If `TRUE`, use C3 photosynthesis; otherwise C4.
#'
#' @details
#' The meteorological forcing data must represent above-canopy conditions.
#' If the supplied forcing height `zref` is below the canopy, the model applies
#' the required height adjustment internally before running the simulation.
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
#' For **bare-ground runs** (`vegp = NA`) the returned elements are:
#' \describe{
#'   \item{tair}{Air temperature profile (°C).}
#'   \item{relhum}{Relative humidity profile (\%).}
#'   \item{windspeed}{Wind speed profile (m s^-1).}
#'   \item{tsoil}{Soil temperature profile (°C).}
#'   \item{theta}{Volumetric soil water content profile (m^3 m^-3).}
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
#'   \item{relhum}{Relative humidity profile (\%).}
#'   \item{windspeed}{Wind speed profile (m s^-1).}
#'   \item{tsoil}{Soil temperature profile (°C).}
#'   \item{theta}{Volumetric soil water content profile (m^3 m^-3).}
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
RunModelFull <- function(climdata, soilc, vegp, paii, Lfrac, lat, long, zref = 2,
                       SoilTempIni = NA, ThetaIni = NA, CO2ppm = NA, boundaryT = NA,
                       zm = 0.004, maxiter = 100, tolerance = 1e-2, C3 = TRUE) {
  # Derive boundary temperature if not supplied
  if (class(boundaryT) == "logical") {
    n <- length(climdata$temp)
    if (n < 8750) stop("Incomplete year. Need to provide boundaryT as ~ mean annual temperature\n")
    boundaryT <- mean(climdata$temp)
  }
  nlay <- soilc$nLayers
  # Create SoilTemp and Theta vectors if not provided
  if (class(SoilTempIni) == "logical") {
    boundaryT <- mean(climdata$temp)
    surfaceT <- climdata$temp[1]
    SoilTempIni = (geometricCpp(nlay, boundaryT - surfaceT) + surfaceT)[1:(1+nlay)]
  }
  mutheta <- seq(0.6, 1, length.out =  nlay + 1)
  if (class(ThetaIni) == "logical")  ThetaIni <- mutheta * soilc$Smax
  # Create obstime
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year +1900,
                        month = tme$mon + 1,
                        day = tme$mday,
                        hour = tme$hour)
  if (class(vegp) == "logical") {  # Bare ground assumed
    # Create z vector
    z <- (c(1:20)/ 20) * zref
    mout <- RunBelowFullBare(obstime, climdata, soilc, z, zref, zm, lat, long,
                             SoilTempIni, ThetaIni, maxiter, tolerance)
  } else {
    # Estimate CO2
    if (is.na(CO2ppm)) CO2ppm <- Cafromyear(tme$year[1] + 1900)
    mout <- RunBelowFull(obstime, climdata, soilc, vegp, paii, Lfrac, zref, CO2ppm,
                         lat, long, SoilTempIni, ThetaIni, maxiter, tolerance,
                         0.25, 1.25, C3)

  }
  return(mout)
}
#' R wrapper for c++ expand_outputCpp function
#'
#' Spline interpolates a Numeric Matrix for better visual representation
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

