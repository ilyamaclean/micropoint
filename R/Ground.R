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
                b=soilparams$b[ii],Psie=soilparams$Psie[ii],
                Smax=soilparams$Smax[ii],Smin=soilparams$Smin[ii],
                alpha=soilparams$alpha[ii],n=soilparams$n[ii],
                Ksat=soilparams$Ksat[ii])
  class(groundp)<-"groundparams"
  groundp
}
#' Calculates moving average by default over 6 months
.ma <- function(x, n = 4380) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}
#' @title Computes ground heat flux
#' @description Computes ground heat flux from hourly surface temperature data
#' @param Tg a vector of ground surface tmeperatures (deg C)
#' @param soilm a vector of fractional soil moistures
#' @param rho Soil bulk density (Mg / m^3)
#' @param Vm Volumetric mineral fraction of soil
#' @param Vq Volumetric quartz fraction of soil
#' @param Mc Mass fraction of clay
#' @param RnetG Net flux density of radiation at soil surface (W/m^2)
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
GFlux<-function(Tg,soilm,rho,Vm,Vq,Mc,RnetG,Gmax,Gmin,i,yearG=TRUE) {
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
  tt<-c(0:(length(Tg)-1))*3600
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
#' @title Runs a simple two-layer soil model
#' @description Runs a simple two-layer soil model to derive soil moisture in upper 10 cm
#' @param weather a data.frame of hourly weather formatted as [climdata()]
#' @param soiltype a character vector of soil type - one of those listed in [soilparams()]
#' @return a vector of hourly volumetrivc soil moisture fractions
#' @rdname soilmmodel
#' @export
#' @examples
#' soilm<-soilmmodel(climdata, "Loam")
#' plot(soilm,type="l")
soilmmodel<-function(weather, soiltype) {
  .onestep<-function(rnetd,rain,soilparams,si,si2,ti,ii) {
    s<-si+soilparams$rmu[ii]*rain[ti]-soilparams$mult[ii]*rnetd[ti]
    sav<-(si+si2)/2
    k<-soilparams$Ksat[ii]*(sav/soilparams$Smax[ii])^soilparams$pwr[ii]
    dif<-si2-si
    s<-s+soilparams$a[ii]*k*dif
    s2<-si2-((soilparams$a[ii]*k*dif)/10)
    s<-ifelse(s>soilparams$Smax[ii],soilparams$Smax[ii],s)
    s<-ifelse(s<soilparams$Smin[ii],soilparams$Smin[ii],s)
    s2<-ifelse(s2>soilparams$Smax[ii],soilparams$Smax[ii],s2)
    s2<-ifelse(s2<soilparams$Smin[ii],soilparams$Smin[ii],s2)
    return(list(s=s,s2=s2))
  }
  # Get soil moisture model coefficients
  ii<-which(soilparams$Soil.type==soiltype)
  # Calculate daily positive net radiation
  swrad<-(1-0.15)*weather$swdown
  lwout<-5.67e-8*0.95*(weather$temp+273.15)^4
  lwnet<-lwout-weather$lwdown
  rnet<-swrad-lwnet
  rnet[rnet<0]<-0
  rnetd<-matrix(rnet,ncol=24,byrow=TRUE)
  rnetd<-apply(rnetd,1,mean)
  # Get daily rainfall
  rain<-matrix(weather$precip,ncol=24,byrow=TRUE)
  rain<-apply(rain,1,sum)
  s<-soilparams$Smax[ii]
  s2<-s
  for (ti in 2:length(rain)) {
    ss<-.onestep(rnetd,rain,soilparams,s[ti-1],s2[ti-1],ti,ii)
    s[ti]<-ss$s
    s2[ti]<-ss$s2
  }
  if (length(rain)>364 & length(rain) < 367) {
    s<-s[length(s)]
    s2<-s2[length(s2)]
    for (ti in 2:length(rain)) {
      ss<-.onestep(rnetd,rain,soilparams,s[ti-1],s2[ti-1],ti,ii)
      s[ti]<-ss$s
      s2[ti]<-ss$s2
    }
  }
  soilm<-(s+s2)/2
  soilm<-stats::spline(soilm,n=length(rnet))$y
  soilm
}

