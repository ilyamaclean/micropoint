#' @title Ground parameters from soil types
#' @description Generates an object of class `groundparams` from a specified
#' soil type.
#' @param soiltype a character vector of soil type - one of those listed in [soilparams()]
#' @return object of class `groundparams`
#' @rdname creategroundp
#' @export
#' @examples
#' creategroundp("Loam")
creategroundp<-function(soiltype) {
  ii<-which(soilparams$Soil.type==soiltype)
  gref<-0.15
  if (soiltype == "Sand") gref<-0.25
  if (soiltype == "Sandy loam") gref<-0.2
  if (soiltype == "Silt loam") gref<-0.12
  if (soiltype == "Sandy clay loam") gref<-0.2
  if (soiltype == "Silty clay loam") gref<-0.12
  if (soiltype == "Silty clay") gref<-0.12
  em<-0.97
  if (soiltype == "Sand") em<-0.92
  if (soiltype == "Sandy loam") em<-0.95
  if (soiltype == "Sandy clay loam") em<-0.95
  groundp<-list(gref=gref,slope=0,aspect=180,em=em,
                rho=soilparams$rho[ii],Vm=soilparams$Vm[ii],
                Vq=soilparams$Vq[ii],Mc=soilparams$Mc[ii],
                b=soilparams$b[ii],Psie=soilparams$psi_e[ii]*-1,
                Smax=soilparams$Smax[ii],Smin=soilparams$Smin[ii],
                alpha=soilparams$alpha[ii],n=soilparams$n[ii],
                Ksat=soilparams$Ksat[ii])
  class(groundp)<-"groundparams"
  groundp
}
#' Calculates moving average by default over 6 months
.ma <- function(x, n = 2190) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}
#' @title Computes ground heat flux
#' @description Computes ground heat flux from hourly surface temperature data
#' @param Tg a vector of ground surface temperatures (deg C)
#' @param soilm a vector of fractional soil moistures
#' @param rho Soil bulk density (Mg / m^3)
#' @param Vm Volumetric mineral fraction of soil
#' @param Vq Volumetric quartz fraction of soil
#' @param Mc Mass fraction of clay
#' @param Gmax Maximum allowable value of daily G (included to ensure model convergence - W/m^2)
#' @param Gmin Minumum allowable value of daily G (included to ensure model convergence - W/m^2)
#' @param i iteration of model
#' @param yearG optional logical indicating whether or not to calculate and incorporate yearly ground heat flux
#' @return a list of the following:
#' (1) G - ground heat fluxes (W/m^2)
#' (2) Gmax - newly compued values of Gmax if on first iteration of model
#' (2) Gmin - newly compued values of Gmax if on first iteration of model
#' @details The maximum allowable values of Gmax and Gmin are computed on the first
#' iteration of the model in which they are derived from the diurnal soil surface temperature
#' values with G set to 0. It is thus assumed that G can only dampen rather than enhance soil
#' surface temperature fluctuations,an assumption that increases the likelihood of iterative
#' convergence
#' @rdname GFlux
#' @export
GFlux<-function(Tg,soilm,rho,Vm,Vq,Mc,Gmax=200,Gmin=200,i=1,yearG=TRUE) {
  # Find soil diffusivity
  cs<-(2400*rho/2.64+4180*soilm) # specific heat of soil in J/kg/K
  ph<-(rho*(1-soilm)+soilm)*1000   # bulk density in kg/m3
  frs<-Vm+Vq
  c1<-(0.57+1.73*Vq+0.93*Vm)/(1-0.74*Vq-0.49*Vm)-2.8*frs*(1-frs)
  c2<-1.06*rho*soilm
  c3<-1+2.6*Mc^-0.5
  c4<-0.03+0.7*frs^2
  k<-c1+c2*soilm-(c1-c4)*exp(-(c3*soilm)^4) # Thermal conductivity   W/m/K
  ka<-k/(cs*ph)
  # Get T multiplier  for daily
  omdy<-(2*pi)/(24*3600)
  DD<-sqrt(2*ka/omdy)
  Gmu<-sqrt(2)*(k/DD)*0.5
  # Calculate T fluctuation from daily mean
  Td<-matrix(Tg,ncol=24,byrow=T)
  Td<-rep(apply(Td,1,mean),each=24)
  dT<-Tg-Td
  # Calculate 6 hour back rolling mean of Gmu and dT to get 3 hour lag
  # NB 1.1171 correction applied as average dampens sinusoidal cycle by this amount
  Gmud<-.ma(Gmu,6)
  Gday<-.ma(dT,6)*Gmud*1.1171
  n<-length(dT)
  if (i == 1) {
    Gd<-matrix(Gday,ncol=24,byrow=T)
    Gmax<-rep(apply(Gd,1,max),each=24)
    Gmin<-rep(apply(Gd,1,min),each=24)
  }
  Gday<-pmin(Gday,Gmax)
  Gday<-pmax(Gday,Gmin)
  # Test whether sequence is for a year or longer
  if (n >= 8760 & yearG == TRUE) {
    # Calculate 6 month back-rolling mean of Gmu and Td to get 3 month time lag
    # NB 1.1171 correction applied as average dampens sinusoidal cycle by this amount
    omyr<-(2*pi)/(n*3600)
    Gmuy<-sqrt(2)*.ma(k)/sqrt(2*.ma(ka)/omyr)
    dTy<-Td-mean(Td)
    Gyear<-.ma(dTy)*Gmuy*1.1171
  } else Gyear<-0
  G<-Gday+Gyear
  return(list(G=G,Gmax=Gmax,Gmin=Gmin))
}
#' Calculate soil effective relative humidity
.soilrh <- function(theta, b, Psie, Smax, tc = 11) {
  matric <- -Psie*(theta/Smax)^-b
  hr<-exp((0.018*matric)/(8.31*(tc+273.15)))
  hr[hr>1]<-1
  hr
}
#' Get soil type form gorund parameters
.soilmatch<-function(groundp) {
  # Calculate standard deviations away from each variable
  Smx<-abs(groundp$Smax-soilparams$Smax)/sd(soilparams$Smax)
  Smn<-abs(groundp$Smin-soilparams$Smin)/sd(soilparams$Smin)
  Rho<-abs(groundp$rho-soilparams$rho)/sd(soilparams$rho)
  Kst<-abs(groundp$Ksat-soilparams$Ksat)/sd(soilparams$Ksat)
  sm<-Smx+Smn+Rho+Kst
  i<-which.min(sm)[1]
  soiltype<-soilparams$Soil.type[i]
  return(soiltype)
}
#' @title Runs a simple two-layer soil model
#' @description Runs a simple two-layer soil model to derive soil moisture in upper 10 cm
#' @param climdata a data.frame of hourly weather formatted as [climdata()]
#' @param soiltype a character vector of soil type - one of those listed in [soilparams()]
#' @return a vector of hourly volumetric soil moisture fractions
#' @rdname soilmmodel
#' @export
#' @examples
#' soilm<-soilmmodel(climdata, "Loam")
#' plot(soilm,type="l")
soilmmodel<-function(climdata, soiltype) {
  # Create date data.frame of obstime
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour+tme$min/60+tme$sec/3600)
  climdata$obs_time<-NULL
  spa<-micropoint::soilparams
  ii<-which(spa$Soil.type==soiltype)
  soilm<-soilmCpp(climdata,spa$rmu[ii],spa$mult[ii],spa$pwr[ii],spa$Smax[ii],spa$Smin[ii],spa$Ksat[ii],spa$a[ii])
  soilm<-stats::spline(soilm,n=length(climdata$temp))$y
  return(soilm)
}

