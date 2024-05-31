#' Cap values
.lim<-function(x,l,up=FALSE) {
  if (length(l)==1) l<-x*0+l
  if (up) {
    s<-which(x>l)
    x[s]<-l[s]
  } else {
    s<-which(x<l)
    x[s]<-l[s]
  }
  x
}
.phair<- function(tc=15,pk=101.3) {
  tk<-tc+273.15
  ph<-44.6*(pk/101.3)*(273.15/tk)
  ph
}
.cpair<-function(tc) {
  cp<-2e-05*tc^2+0.0002*tc+29.119
  cp
}
#' Calculate zero-plane displacement
.zeroplanedis<-function(h,pai) {
  pai[pai<0.001]<-0.001
  d<-(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*h
  d
}
#' Calculate roughness length
.roughlength<-function(h,pai,d=NA,psi_h=0) {
  if (class(d)=="logical") d<-.zeroplanedis(h,pai)
  Be<-sqrt(0.003+(0.2*pai)/2)
  zm<-(h-d)*exp(-0.4/Be)*exp(psi_h)
  zm[zm<0.0005]<-0.0005
  zm
}
#' Calculate integrated diabatic correction coefficient for momentum
.dpsim<-function(ze) {
  # unstable
  x<-(1-15*ze)^0.25
  psi_m<-log(((1+x)/2)^2*(1+x^2)/2)-2*atan(x)+pi/2
  # stable
  s<-which(ze>=0)
  p2<- -4.7*ze
  psi_m[s]<-p2[s]
  psi_m[psi_m < -4] <- -4
  psi_m[psi_m > 3] <- 3
  psi_m
}
#' Calculate integrated diabatic correction coefficient for heat
.dpsih<-function(ze) {
  # unstable
  y<-(1-9*ze)^0.5
  psi_h<-log(((1+y)/2)^2)
  # stable
  s<-which(ze>=0)
  p2<- -(4.7*ze)/0.74
  psi_h[s]<-p2[s]
  psi_h[psi_h < -4] <- -4
  psi_h[psi_h > 3] <- 3
  psi_h
}
.dphih<-function(ze) {
  # unstable
  phim<-1/((1-16*ze)^0.25)
  phih<-phim^2
  # stable
  phis<-1+((6*ze)/(1+ze))
  s<-which(ze>0)
  phih[s]<-phis[s]
  phih
}
#' Calculate free convection
.gfree<-function(leafd,H) {
  d<-0.71*leafd
  dT<-0.7045388*(d*H^4)^(1/5)
  gha<-0.0375*(dT/d)^0.25
  gha[gha<0.1]<-0.1
  gha
}
#' Calculate molar conductance above canopy
.gturb<-function(uf,d,zm,z1,z0=NA,ph=43,psi_h=0,gmin=0.03) {
  if (is.na(z0)) z0<-0.2*zm+d
  ln<-log((z1-d)/(z0-d))
  g<-(0.4*ph*uf)/(ln+psi_h)
  g<-.lim(g,gmin)
  g<-.lim(g,0.0001)
  g
}
#' Stomatal conductance (set g50 to 300 and multiply by 3 to give bulk surface)
.stomcond <- function(Rsw, gsmax, q50 = 100) {
  rpar<-Rsw*4.6
  gs<-(gsmax * rpar)/(rpar + q50)
  gs
}
.satvap <- function(tc, method = "Tetens") {
  if (method == "Buck") {
    es<-0.61121*exp((18.678-tc/234.5)*(tc/(257.14+tc)))
    ei<-0.61115*exp((23.036-tc/333.7)*(tc/(279.82+tc)))
    s<-which(tc<0)
    es[s]<-ei[s]
  } else {
    es<-0.61078*exp(17.27*tc/(tc+237.3))
    ei<-0.61078*exp(21.875*tc/(tc+265.5))
    s<-which(tc<0)
    es[s]<-ei[s]
  }
  es
}
.dewpoint <- function(tc, ea) {
  # Dew point
  e0 <- 611.2/1000
  L <- (2.501*10^6) - (2340 * tc)
  T0 <- 273.15
  Rv <- 461.5
  it <- 1/T0 - (Rv/L) * log(ea/e0)
  Tdew <- 1/it - 273.15
  # Frost point
  e0 <- 610.78/1000
  L <- 2.834*10^6
  T0 <- 273.16
  it <- 1/T0 - (Rv/L) * log(ea/e0)
  Tfrost <- 1/it - 273.15
  Tdew[Tdew < 0] <- Tfrost[Tdew < 0]
  Tdew
}
#' @title Run big-leaf model
#' @description Runs big leaf model iteratively to derive ground and canopy
#' heat exchange surface temperatures
#' @param weather a data.frame of hourly weather as in the example dataset `climdata`
#' @param vegp an object of class vegparams formatted as for the inbuilt example dataset `vegparams`
#' @param groundp an object of class vegparams formatted as for the inbuilt example dataset `groundparams`
#' @param soilm a vector of volumetric soil moisture fractions in the top 10 cm of soil.
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long longitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @param dTmx maximum by which vegetation or soil surface temperature can exceed air temperature (deg C, set to ensure model convergence)
#' @param zref height above ground (m) of measurements in `weather` (see details)
#' @param uref height above gorund (m) of wind speed measurements in `climdata` (if different from zref)
#' @param merid Longitude of local time zone meridian (decimnal degrees), Default: 0 (Greenwich Mean Time)
#' @param dst Daylight saving hours. E.g. with `merid = 0` 0 for UTC or 1  for British Summer Time.
#' @param maxiter Maximum number of iterations over which to run the model
#' @param bwgt backward weighting to apply when iteratively running the model (default 0.5)
#' @param tol Error margin (deg C) deemed acceptable for model convergence (default 0.5).
#' @param gmn optional minimum convective conductance value (mol/m^2/s). Lower values are more
#' physically realistic, but can reduce likelihood of model convergence.
#' @param plotout optional logical indicating whether to plot a time-series of ground surface
#' temperatures on each iteration
#' @param swmethod method used to calculate fraction of ground surface acting as free water surface
#' (0 = based on rainfall, 1 computed from soil effective relative humidity, 2 computed from soil moisture fraction)
#' @param yearG optional logical indicating whether or not to calculate and account for annual ground heat flux cycle
#' @return an object of class pointmicro, namely a list of the following:
#' \describe{
#'  \item{tme}{POSIXlt object of times and dates corresponding to model outputs}
#'  \item{Tc}{a vector of canopy heat exchange surface temperatures (deg C)}
#'  \item{Tg}{a vector of ground surface temperatures (deg C)}
#'  \item{H}{a vector of sensible heat fluxes (W/m^2)}
#'  \item{G}{a vector of ground heat fluxes (W/m^2)}
#'  \item{psih}{diabatic correction factors for heat}
#'  \item{psim}{diabatic correction factors for momentum}
#'  \item{phih}{diabatic influencing factors for heat}
#'  \item{OL}{Obukhov length}
#'  \item{error.mar}{maximum temperature difference between ultimate and penultimate iteration of model}
#' }
#' @details When running the model `zref` must be greater than he height of vegetation.
#' If running for a location with tall vegetation, function [WeatherHeight()] can
#' be used to adjust e.g. temperature and wind speed to a user-specified height.
#'
#' @rdname RunBigLeaf
#' @export
#'
#' @examples
#' # Run Big leaf model with inbuilt parameters
#' bigleafp <- RunBigLeaf(climdata, vegparams, groundparams, lat = 49.96807,
#' long = -5.215668)
#' # Plot ground and canopy temperature
#' tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
#' par(mar = c(6, 6, 3, 3))
#' # Ground temperature
#' plot(bigleafp$Tg ~ tme, type = "l", cex.axis = 2, cex.lab = 2,
#'      xlab = "Month", ylab = "Temperature", ylim = c(-5, 40),
#'      col = rgb(1, 0, 0, 0.5))
#' par(new = TRUE)
#' # Canopy temperature
#' plot(bigleafp$Tc ~ tme, type = "l", cex.axis = 2, cex.lab = 2,
#'      xlab = "", ylab = "", ylim = c(-5, 40), col = rgb(0, 0.5, 0, 0.5))
RunBigLeaf<-function(climdata,  vegp, groundp, lat, long, zref = 2, uref = zref, soilm = NA, surfwet = NA, dTmx = 25,
                     maxiter = 20) {
  # Create date data.frame of obstime
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour+tme$min/60+tme$sec/3600)
  climdata$obs_time<-NULL
  # Create vectors of vegp and groundp
  vegpp<-as.vector(unlist(vegp))
  groundpp<-as.vector(unlist(groundp))
  if (vegp$h > zref) {
    climdata<-weatherhgtCpp(obstime, climdata, zref, uref, vegp$h, lat, long)
    zref<-vegp$h
  }
  climdata$windspeed[climdata$windspeed<0.5]<-0.5
  if (class(soilm)=="logical") {
    soiltype<-.soilmatch(groundp)
    spa<-micropoint::soilparams
    ii<-which(spa$Soil.type==soiltype)
    soilm<-soilmCpp(climdata,spa$rmu[ii],spa$mult[ii],spa$pwr[ii],spa$Smax[ii],spa$Smin[ii],spa$Ksat[ii],spa$a[ii])
    soilm<-stats::spline(soilm,n=length(climdata$temp))$y
  }
  # check whether time sequence is for whole year
  yearG<-TRUE
  nn<-dim(climdata)[1]
  if (nn < 8760) yearG<-FALSE
  # Run big Leaf model
  microp<-BigLeafCpp(obstime,climdata,vegpp,groundpp,soilm,lat,long,dTmx,zref,maxiter,0.5,0.5,0.1,yearG)
  class(bigleafp)<-"pointmicro"
  return(microp)
}


#'
#' @title Adjust weather data for height
#' @description Adjust weather data to user-specified height above ground
#' @param weather a data.frame of hourly weather as in the example dataset `climdata`
#' @param zin height above ground of temperature and relative humidity measurements (m)
#' @param uzin height above ground of wind speed measurements (m)
#' @param zout height above ground for which observations are needed
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long latitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @return a data.frame of hourly weather
#'
#' @details When running [RunBigLeaf()], `zref` the height at which weather
#' observations are assumed to have been taken must be greater than the
#' height of vegetation. This function derives estimates of temperature,
#' relative humidity and wind speed at a user specified height above ground,
#' by assuming that weather observations were taken in an open area with
#' with a vegetation height of 12 cm and a plant area index of 1, consistent
#' with WMO guidelines for locating weather stations.
#' @examples
#' weather<-weather_hgt(climdata, 2, 2, 10, 49.96807, -5.215668)
#' # Temperature comparison
#' plot(weather$temp ~ climdata$temp, xlab = "Temperature at 2 m",
#'      ylab = "Temperature at 10m", pch = 15)
#' abline(a = 0, b = 1, lwd = 2, col = "red")
#' # Relative humidity comparison
#' plot(weather$relhum ~ climdata$relhum, xlab = "Relative humidity at 2 m",
#'      ylab = "Relative humidity at 10m", pch = 15)
#' abline(a = 0, b = 1, lwd = 2, col = "red")
#' # Wind speed comparison
#' plot(weather$windspeed ~ climdata$windspeed, xlab = "Wind speed at 2 m",
#'      ylab = "Wind speed at 10m", pch = 15)
#' abline(a = 0, b = 1, lwd = 2, col = "red")
#' @rdname weather_hgt
#' @export
weather_hgt<-function(climdata, zin = 2, uzin = zin, zout = 10, lat, long) {
  # Create date data.frame of obstime
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour+tme$min/60+tme$sec/3600)
  climdata<-weatherhgtCpp(obstime, climdata, zin, uzin, zout, lat, long)
  return(climdata)
}
