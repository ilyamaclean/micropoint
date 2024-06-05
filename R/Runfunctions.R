#' @title runs checks on microclimate input data
#' @description The function `runchecksfun` runs checks on the microclimate model
#' inputs
#' @param climdata a data.frame of hourly weather as in the example dataset
#' `climdata` (see details)
#' @param vegp an object of class vegparams or forestparams formatted as for the inbuilt example
#' datasets `vegparams` or `forestparams` (see details).
#' @param groundp an object of class groundparams formatted as for the inbuilt
#' example dataset `groundparams`.
#' @param latitude (decimal degrees) - used for radiation checks
#' @param longitude (decimal degrees) - used for radiation checks
#' @details this function runs checks on provided data testing whether provided data lie within plausable ranges
#' and returning errors or warnings accordingly. The original data are returned as a list with
#' some values capped if checks fail.
runchecksfun <- function(climdata, vegp, groundp, lat, long) {
  # run checks in vegp
  if (vegp$h < 0) stop("vegp$h must be positive")
  if (vegp$pai < 0) stop("vegp$pai must be positive")
  if (vegp$pai > 50) warning("vegp$pai seems too high")
  if (vegp$x < 0) stop("vegp$x must be positive")
  if (vegp$clump < 0) stop("vegp$clump must be between 0 and 1")
  if (vegp$clump > 1) stop("vegp$clump must be between 0 and 1")
  if (vegp$lref < 0) stop("vegp$lref must be between 0 and 1")
  if (vegp$lref > 1) stop("vegp$lref must be between 0 and 1")
  if (vegp$ltra < 0) stop("vegp$ltra must be between 0 and 1")
  if (vegp$ltra > 1) stop("vegp$ltra must be between 0 and 1")
  if (vegp$ltra+vegp$lref > 1) stop("vegp$ltra + vegp$lref cannot be greater than 1")
  if (vegp$leafd > 5) warning("vegp$leafd seems too large")
  if (vegp$em < 0) stop("vegp$em must be between 0 and 1")
  if (vegp$em > 1) stop("vegp$em must be between 0 and 1")
  if (vegp$gsmax < 0) stop("vegp$gsmax cannot be negative")
  if (vegp$gsmax < 0.01) warning("vegp$gsmax seems low")
  if (vegp$gsmax > 2) warning("vegp$gsmax seems high")
  if (vegp$q50 < 0) stop("vegp$q50 cannot be negative")
  if (vegp$q50 > 500) warning("vegp$q5 seems high")
  if (length(vegp$skew) > 0) {
    if (vegp$skew<0) stop("vegp$skew cannot be negative")
    if (vegp$skew>10) stop("vegp$skew cannot exceed 10")
  }
  if (length(vegp$spread) > 0) {
    if (vegp$spread<0) stop("vegp$spread cannot be negative")
    if (vegp$spread>100) stop("vegp$spread cannot exceed 100")
  }
  # run checks in groundp
  if (groundp$gref < 0.0001) stop("groundp$gref cannot be less than 0.0001")
  if (groundp$gref > 1) stop("groundp$gref cannot be greater than 1")
  if (groundp$slope < 0) stop("groundp$slope cannot be negative")
  if (groundp$slope > 180) stop("groundp$slope cannot be grewater than 180")
  if (groundp$slope > 90) warning("model not tested with groundp$slope > 90")
  if (groundp$aspect < 0) {
    warning("aspect adjusted to range 0-360 degrees")
    groundp$aspect<-groundp$aspect%%360
  }
  if (groundp$aspect > 360) {
    warning("aspect adjusted to range 0-360 degrees")
    groundp$aspect<-groundp$aspect%%360
  }
  if (groundp$em < 0) stop("groundp$em must be between zero and 1")
  if (groundp$em > 1) stop("groundp$em must be between zero and 1")
  if (groundp$rho < 0) stop("groundp$rho cannot be negative")
  if (groundp$rho > 2.5) warning("groundp$rho seems high. Check units")
  if (groundp$Vm < 0) stop("groundp$Vm must be between 0 and 1")
  if (groundp$Vm > 1) stop("groundp$Vm must be between 0 and 1")
  if (groundp$Vq < 0) stop("groundp$Vq must be between 0 and 1")
  if (groundp$Vq > 1) stop("groundp$Vq must be between 0 and 1")
  if (groundp$Mc < 0) stop("groundp$Mc must be between 0 and 1")
  if (groundp$Mc > 1) stop("groundp$Mc must be between 0 and 1")
  if (groundp$b < 0) stop("groundp$b must be positive")
  if (groundp$b > 10) warning("groundp$b seems high")
  if (groundp$Psie > 0) stop("groundp$Psie must be negative")
  if (groundp$Smax < 0) stop("groundp$Smax must be between 0 and 1")
  if (groundp$Smax > 1) stop("groundp$Smax must be between 0 and 1")
  if (groundp$Smin < 0) stop("groundp$Smin must be between 0 and 1")
  if (groundp$Smin > 1) stop("groundp$Smin must be between 0 and 1")
  if (groundp$Smin > groundp$Smax) stop("groundp$Smin cannot be greater than groundp$Smax")
  if (groundp$alpha < 0) stop("groundp$alpha must be positive")
  if (groundp$alpha > 0.2) warning("groundp$alpha seems high")
  if (groundp$n < 0) stop("groundp$n must be positive")
  if (groundp$n > 5) warning("groundp$n seems high")
  if (groundp$Ksat < 0) stop("groundp$Ksat must be positive")
  if (groundp$Ksat > 2000) warning("groundp$Ksat seems high. Check units")
  # Check climate data
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  if (min(climdata$temp) < -65) warning("Minumum temperature seems low")
  if (max(climdata$temp) > 65) warning("Maximum temperature seems high")
  if (min(climdata$relhum) < 15) warning("Minumum relative humidity seems low")
  if (min(climdata$relhum) < 0) stop("Relative humidity cannot be negative")
  if (max(climdata$relhum) > 100) {
    warning("Maximum relative humidity cannot exceed 100. Capping at 100")
    s<-which(climdata$relhum>100)
    climdata$relhum[s]<-100
  }
  if (max(climdata$relhum) > 100) {
    warning("Maximum relative humidity cannot exceed 100. Capping at 100")
    s<-which(climdata$relhum>100)
    climdata$relhum[s]<-100
  }
  if (min(climdata$pres) < 80) warning("Minumum pressure seems low. OK if site at high altitude")
  if (max(climdata$pres) > 108.4) warning("Maximum pressure seems high")
  lt<-tme$hour+tme$min/60+tme$sec/3600
  csr<-clearskyradCpp(tme$year+1900,tme$mon+1,tme$mday,lt,lat,long,climdata$temp,climdata$relhum,climdata$pres)
  s<-which(climdata$swdown > csr)
  if (length(s) > 0) {
    warning("Downward shortwave radiation higher than expected clear-sky radiation. Capping at clear-sky value")
    climdata$swdown[s]<-csr[s]
  }
  s2<-which(climdata$difrad > climdata$swdown)
  if (length(s2) > 0) {
    warning("Diffuse radiation higher than total downward radiation. Capping at downward radiation")
    climdata$difrad[s2]<-climdata$swdown[s2]
  }
  lwup<-5.67*10^-8*(climdata$temp+273.15)^4
  s<-which(climdata$lwdown > (lwup*1.1))
  if (length(s) > 0) warning("Downward longwave radiation seems high given air temperature")
  if (min(climdata$windspeed) < 0) stop("Wind speed cannot be negative")
  if (max(climdata$windspeed) > 55) warning("Maximum wind speed seems high")
  if (min(climdata$precip) < 0) stop("Precipitation cannot be negative")
  if (max(climdata$precip) > 150) warning("Maximum precipitation rate seems high")
  return(list(climdata=climdata,vegp=vegp,gorundp=groundp))
}

#' @title runs point microclimate model
#' @description The function `runpointmodel` runs the point microclimate model
#' above or below ground.
#' @param climdata a data.frame of hourly weather as in the example dataset
#' `climdata` (see details)
#' @param reqhgt height (m) above (positive) or below (negative) ground for
#' which model output is wanted.
#' @param vegp an object of class vegparams or forestparams formatted as for the inbuilt example
#' datasets `vegparams` or `forestparams` (see details).
#' @param paii optional vector of plant area index values for each canopy layer
#' as returned by [PAIgeometry()] (see details).
#' @param groundp an object of class groundparams formatted as for the inbuilt
#' example dataset `groundparams`.
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long longitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @param zref height above ground (m) of e.g. temperature measurements in `climdata`
#' @param uref height above gorund (m) of wind speed measurements in `climdata` (if different from zref)
#' @param soilm optional vector of hourly volumetric soil moisture fractions in
#' the top 10 cm of soil. Calculated if not supplied.
#' @param surfwet optional vector of canopy surface acting like a saturated
#' water surface. Calculated if not supplied.
#' @param dTmx = maximum by which vegetation or soil surface temperature can
#' exceed air temperature (deg C, set to ensure model convergence).
#' @param maxiter = maximum number of iterations over which to run the model to
#' achieve model convergence.
#' @param n number of canopy layers used to run model. if `paii` supplied,
#' ignored and taken from length of `paii`.
#' @details Variable `obs_time` in `climdata` assumes times are in UTC. To derive
#' within canopy temperatures , the canopy is divided into `n` layers
#' and each element `i` of `paii` represents the plant area index within a layer at
#' height `(i / n) * hgt` where `n` is the total number of layers and `hgt` is canopy
#' height. If `paii` is not supplied, `vegp` must contain entries for `skew` and
#' `spread` as used by [PAIgeomtry()], which is used to generate a plausable
#' vector of `paii` values. Note that `sum(paii)` is assumed  to equal `vegp$pai`,
#' the later representing the total leaf area of the canopy.
#' @return if `reqhgt > 0` a Dataframe with the following columns:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{tair}{Air temperature (degrees C) at height `reqhgt`}
#'  \item{tleaf}{Leaf temperature (degrees C) at height `reqhgt` or `tcanopy` average canopy temperature if `reqhgt` above canopy}
#'  \item{relhum}{relative humidity (percentage) at height `reqhgt`}
#'  \item{windspeed}{wind speed (m/s) at height `reqhgt`}
#'  \item{Rdirdown}{Total downward flux of direct radiation (W / m^2) at height `reqhgt` perpendicular to solar beam}
#'  \item{Rdifdown}{Total downward flux of diffuse radiation (W / m^2) at height `reqhgt`}
#'  \item{Rlwdown}{Total downward flux of longwave radiation (W / m^2) at height `reqhgt`}
#'  \item{Rswup}{Total upward flux of shortwave radiation (W / m^2) at height `reqhgt` (assumed diffuse)}
#'  \item{Rlwup}{Total upward flux of longwave radiation (W / m^2) at height `reqhgt`}
#' }
#' if `reqhgt = 0` `tair` is replaced by `tground` representing ground surface
#' temperature (deg C), `relhum` by `soilm` (fractional soil moisture) and `tleaf`
#' and `windspeed` are not returned. If `reqhgt < 0` only `obs_time`, `tground` and
#'`soilm` are returned.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
#' @rdname runpointmodel
#' @examples
#' # run model using inbuilt datasets and defaults
#' mout<-runpointmodel(climdata, reqhgt = 1, forestparams, paii = NA, groundparams,
#'   lat =49.96807, long= -5.215668)
#' plot(mout$tair,type="l")
runpointmodel<-function(climdata,  reqhgt, vegp, paii = NA, groundp, lat, long, zref = 2, uref = zref, soilm = NA, surfwet = NA, dTmx = 25,
                        maxiter = 20, n = 20, runchecks = FALSE) {

  if (runchecks) {
    rcd<-runchecksfun(climdata, vegp, groundp, lat, long)
    climdata<-rcd$climdata
    vegp<-rcd$vegp
    groundp<-rcd$groundp
  }
  # generate paii if doesn't exist'
  if (class(paii) == "logical") {
    paii<-PAIgeometry(vegp$pai, vegp$skew, vegp$spread, n)
  }
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
  # Create Bigleaf Data.Frame
  bigleafvars<-data.frame(Tc=microp$Tc,Tg=microp$Tg,H=microp$H,psih=microp$psih,psim=microp$psim,
                          phih=microp$phih,uf=microp$uf,sm=soilm,albedo=microp$albedo)
  if (class(surfwet) == "logical") {
    trc<-exp(-vegp$pai)
    rmu<-ifelse(climdata$precip > 0,1,0.8)
    swet<-(soilm-groundp$Smin)/(groundp$Smax-groundp$Smin)
    bigleafvars$surfwet<-trc*swet+(1-trc)*rmu
  } else bigleafvars$surfwet<-surfwet
  # Run point model
  modp<-runmodel(reqhgt,zref,lat,long,obstime,climdata,bigleafvars,maxiter,vegpp,paii,groundpp)
  obs_time<-tme
  modp<-cbind(obs_time,modp)
  return(modp)
}
#' @title plot temperature or humidity profile
#' @description The function `plotprofile` plots a temperature or humidity height
#' profile and returns a list of heights or humidites.
#' above or below ground.
#' @param climdata a data.frame of hourly weather as in the example dataset
#' `climdata` (see details)
#' @param hr the hour in climdata for which the height profile is to be plotted
#' @param plotout of `tair`, `tleaf`, `ea` or `relhum`
#' @param vegp an object of class vegparams or forestparams formatted as for the inbuilt example
#' datasets `vegparams` or `forestparams` (see details).
#' @param paii optional vector of plant area index values for each canopy layer
#' as returned by [PAIgeometry()] (see details).
#' @param groundp an object of class groundparams formatted as for the inbuilt
#' example dataset `groundparams`.
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long longitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @param zref height above ground (m) of measurements in `weather` (see details)
#' @param uref height above gorund (m) of wind speed measurements in `climdata` (if different from zref)
#' @param soilm optional vector of hourly volumetric soil moisture fractions in
#' the top 10 cm of soil. Calculated if not supplied.
#' @param surfwet optional vector of canopy surface acting like a saturated
#' water surface. Calculated if not supplied.
#' @param dTmx = maximum by which vegetation or soil surface temperature can
#' exceed air temperature (deg C, set to ensure model convergence).
#' @param maxiter = maximum number of iterations over which to run the model to
#' achieve model convergence.
#' @param n number of canopy layers used to plot the model. if `paii` supplied,
#' ignored and taken from length of `paii`.
#' @details if `plotout` = `tleaf` the height propofile of leaf tmeperature is plotted
#' from  the ground surface to the top of the canopy. If `plotout` = `tair`, `ea` or `relhum`
#' the air temperature (tair), vapour pressure (ea) or relative humidity (relhum) profiles
#' are plotted from the ground to two metres above canopy.
#' @return a list of the following
#' \describe{
#'  \item{z}{height above gorund (m)}
#'  \item{var}{one of leaf temperature (degrees C, air temperature (degrees C),
#'  vapour pressure (kPa) or relative humidity (percentage) for eahc height `z`}
#' }
#'`soilm` are returned.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib micropoint, .registration = TRUE
#' @export
#' @rdname plotprofile
#' @examples
#' # plot air temperature for hottest hour
#' xx<-plotprofile(climdata, hr = 4094, "tair", forestparams, paii = NA, groundparams,
#'   lat =49.96807, long= -5.215668)
#' # plot leaf temperature for hottest hour
#' xx<-plotprofile(climdata, hr = 4094, "tleaf", forestparams, paii = NA, groundparams,
#'   lat =49.96807, long= -5.215668)
plotprofile<-function(climdata, hr, plotout = "tair", vegp, paii = NA, groundp, lat, long, zref = 2, uref = zref, soilm = NA, surfwet = NA, dTmx = 25,
                      maxiter = 20, n = 100) {
  # generate paii if doesn't exist'
  if (class(paii) == "logical") {
    paii<-PAIgeometry(vegp$pai, vegp$skew, vegp$spread, n)
  }
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
  # Run profile
  obsvars<-as.numeric(obstime[hr,])
  climvars<-c(climdata$temp[hr],climdata$relhum[hr],climdata$pres[hr],climdata$swdown[hr],climdata$difrad[hr],
              climdata$lwdown[hr])
  bigleafvarsone=c(microp$Tc[hr],microp$Tg[hr],microp$psih[hr],microp$psim[hr],microp$phih[hr],microp$uf[hr],soilm[hr])
  wc<-CanopyWindCpp(vegp$h,paii)
  tleaf<-rep(microp$Tc[hr],length(paii))
  tair<-rep(climdata$temp[hr],length(paii))
  ea<-.satvap(tleaf)*climdata$relhum[hr]/100
  if (class(surfwet) == "logical") {
    trc<-exp(-vegp$pai)
    rmu<-ifelse(climdata$precip > 0,1,0.8)
    swet<-(soilm-groundp$Smin)/(groundp$Smax-groundp$Smin)
    surfwet<-trc*swet+(1-trc)*rmu
  }
  reqhgt<-vegp$h/2
  n<-length(paii)
  zb<-(c(1:n)/n)*vegp$h
  below<-SmallLeafOne(reqhgt,zref,lat,long,obsvars,climvars,bigleafvarsone,maxiter,wc,vegpp,
                      paii,groundpp,tleaf,tair,ea,surfwet[hr],zb)
  if (plotout == "tleaf") {
    plot(zb~below$tleaf,type="l", , xlab = "Leaf temperature (deg C)",
         ylab = "Height (m)", col = "darkgreen", lwd = 2)
    zz<-zb
    out<-below$tleaf
  } else {
    nn<-(2/vegp$h)*n
    za<-(c(1:nn)/nn)*2+vegp$h
    d<-zeroplanedisCpp(vegp$h, vegp$pai)
    zm<-roughlengthCpp(vegp$h, vegp$pai, d, microp$psih[hr]);
    zh<-0.2*zm
    lnr<-log((za-d)/zh)/log((zref-d)/zh)
    zz<-c(zb,za)
    if (plotout == "tair") {
      ta<-climdata$temp[hr]+(microp$Tc[hr]-climdata$temp[hr])*(1-lnr)
      tair<-c(below$tair,ta)
      plot(zz~tair,type="l", xlab = "Air temperature (deg C)",
           ylab = "Height (m)", col = "red", lwd = 2)
      out<-tair
    } else {
      ea<-satvapCpp(climdata$temp[hr])*climdata$relhum[hr]/100;
      estl = satvapCpp(microp$Tc[hr]);
      eza<-ea+(estl-ea)*surfwet[hr]*(1-lnr)
      ezba<-c(below$ea,eza)
      if (plotout == "relhum") {
        ta<-climdata$temp[hr]+(microp$Tc[hr]-climdata$temp[hr])*(1-lnr)
        tair<-c(below$tair,ta)
        relhum<-(.satvap(tair)/ezba)*100
        relhum[relhum > 100]<-100
        plot(zz~relhum,type="l", xlab = "Relative humidity (%)",
             ylab = "Height (m)", col = "blue", lwd = 2)
        out<-relhum
      } else {
        plot(zz~ezba,type="l", xlab = "Vapour pressure (kPa)",
             ylab = "Height (m)", col = "blue", lwd = 2)
        out<-ezba
      }
    }
  }
  return(list(z=zz,var=out))
}
