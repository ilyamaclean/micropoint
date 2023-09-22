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
  y<-(1-9*ze)^-0.5
  # stable
  s<-which(ze>=0)
  p<-(0.74+4.7*ze)/0.74
  y[s]<-p[s]
  y
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
  gs<-(gsmax * rpar) / (rpar + q50)
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
  # Fost point
  e0 <- 610.78/1000
  L <- 2.834*10^6
  T0 <- 273.16
  it <- 1/T0 - (Rv/L) * log(ea/e0)
  Tfrost <- 1/it - 273.15
  Tdew[Tdew < 0] <- Tfrost[Tdew < 0]
  Tdew
}
.PenmanMonteith<-function(Rabs,gHa,gV,tc,pk,ea,em=0.97,G=0,tms=3,erh=1) {
  .PMon<-function(Rabs,gHa,gV,tc,pk,ea,em,G,te,erh) {
    sb<-5.67*10^-8
    Rema<-em*sb*(tc+273.15)^4
    la<-(-42.575*te+44994)
    cp<-.cpair(te)
    Da<-.satvap(tc)-ea
    gR<-(4*em*sb*(te+273.15)^3)/cp
    De<-.satvap(te+0.5)-.satvap(te-0.5)
    Ts<-tc+((Rabs-Rema-la*(gV/pk)*Da*erh-G)/(cp*(gHa+gR)+la*(gV/pk)*De*erh))
    Ts
  }
  Ts<-.PMon(Rabs,gHa,gV,tc,pk,ea,em,G,tc,erh)
  if (tms>1) {
    for (i in 2:tms) {
      te<-(Ts+tc)/2
      Ts<-.PMon(Rabs,gHa,gV,tc,pk,ea,em,G,te,erh)
    }
  }
  Ts
}
#' Run big-leaf model: one iteration
.BigLeaf<-function(weather,zref,vegp,groundp,bigleafp,solar,twostreamp,soilm,lat,long,merid,dst,gmn,i,method,yearG) {
  # Calculate absorbed radiation
  Rads<-RadiationBigLeaf(weather,vegp,groundp,bigleafp,solar,twostreamp,lat,long,merid,dst)
  RabsG<-with(Rads,radGsw+radGlw)
  RabsC<-with(Rads,radCsw+radClw)
  # Calculate canopy temperature
  d<-with(vegp,.zeroplanedis(h,pai))
  zm<-with(vegp,.roughlength(h,pai,d,bigleafp$psih))
  uf<-(0.4*weather$windspeed)/(log((zref-d)/zm)+bigleafp$psim)
  uf[uf<0.0002]<-0.0002
  gmin<-with(vegp,.gfree(leafd,abs(bigleafp$H))*2*pai)
  ph<-.phair(bigleafp$tcc,weather$pres)
  gHa<-.gturb(uf,d,zm,zref,z0=NA,ph,bigleafp$psih,gmin=gmin)
  gV<-with(vegp,.stomcond(weather$swdown,gsmax*3,q50*3))
  ea<-.satvap(weather$temp)*weather$relhum/100
  Tc<-.PenmanMonteith(RabsC,gHa,gV,weather$temp,weather$pres,ea,vegp$em,bigleafp$G)
  tdew<-.dewpoint(weather$temp,ea)
  Tc[Tc<tdew]<-tdew[Tc<tdew]
  # Calculate effective soil relative humidity
  if (method == 0) {
    srh<-ifelse(weather$precip>0,0.9,0)
  } else if (method == 1) {
    srh<-with(groundp,.soilrh(soilm, b, Psie, Smax, bigleafp$Tg))
  } else {
    srh<-with(groundp,(soilm-Smin)/(Smax-Smin))
  }
  srh<-ifelse(weather$precip>0,0.9,0)
  #srh<-with(groundp,.soilrh(soilm, b, Psie, Smax, bigleafp$Tg))
  Tg<-.PenmanMonteith(RabsG,gHa,gHa,weather$temp,weather$pres,ea,groundp$em,bigleafp$G,3,srh)
  Tg[Tg<tdew]<-tdew[Tg<tdew]
  # Recalculate Big Leaf parameters
  # ** diabatic coefficients
  tcc<-(Tc+weather$temp)/2
  tcg<-(Tg+weather$temp)/2
  Tk<-273.15+tcc
  ph<-.phair(tcc,weather$pres)
  cp<-.cpair(tcc)
  H<-gHa*cp*(Tc-weather$temp)
  # ** set limits to H
  Rnet<-RabsC-5.67*10^-8*vegp$em*(Tc+273.15)^4
  s<-which(Rnet>0 & H>Rnet)
  H[s]<-Rnet[s]
  # Stability
  LL<-(ph*cp*uf^3*Tk)/(-0.4*9.81*H)
  psim<-.dpsim(zm/LL)-.dpsim((zref-d)/LL)
  psih<-.dpsih(zm/LL)-.dpsih((zref-d)/LL)
  phih<-.dphih((zref-d)/LL)
  # ** Set limits to diabatic coefficients
  ln1<-log((zref-d)/zm)
  ln2<-log((zref-d)/(0.2*zm))
  psim<-pmax(psim,-0.9*ln1)
  psih<-pmax(psih,-0.9*ln1)
  psim<-pmin(psim,0.9*ln1)
  psih<-pmin(psih,0.9*ln1)
  phih[phih>1.5]<-1.5
  phih[phih<0.75]<-0.75
  # ** G
  RnetG<-RabsG-5.67*10^-8*vegp$em*(weather$temp+273.15)^4 # deliberately capped at weather temp
  GG<-with(groundp,GFlux(Tg,soilm,rho,Vm,Vq,Mc,RnetG,bigleafp$Gmax,bigleafp$Gmin,i,yearG))
  # Assign values to bigelafp
  bigleafp<-list(psih=psih,psim=psim,H=H,G=GG$G,tcc=tcc,tcg=tcg,Tc=Tc,Tg=Tg,phih=phih,
                 RnetG=RnetG,LL=LL,uf=uf,RabsG=RabsG,Gmax=GG$Gmax,Gmin=GG$Gmin)
  return(bigleafp)
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
#' @param zref height above ground of measurements in `weather` (see details)
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
#' (1) tme - POSIXlt object of times and dates corresponding to model outputs
#' (2) Tc - a vector of canopy heat exchange surface temperatures (deg C).
#' (3) Tg - a vector of ground surface temperatures (deg C).
#' (4) H - a vector of sensible heat fluxes (W/m^2)
#' (5) G - a vector of ground heat fluxes (W/m^2)
#' (6) psih - diabatic correction factor for heat
#' (7) psim - diabatic correction factor for momentum
#' (8) phih - diabatic influencing factor for heat
#' (9) OL - Obukhov length
#' (10) error.mar - maximum temperature difference between ultimate and penultimate iteration of model

#' @details When running the model `zref` must be greater than he height of vegetation.
#' If running for a location with tall vegetation, function [WeatherHeight()] can
#' be used to adjust e.g. temperature and wind speed to a user-specified height.
#'
#' @rdname RunBigLeaf
#' @export
#'
#' @examples
#' # Run Big leaf model with inbuilt parameters assuming soil moisture = 0.2
#' bigleafp <- RunBigLeaf(climdata, vegparams, groundparams, soilm = 0.2,
#'                        lat = 49.96807, long = -5.215668)
#' # Plot ground and canopy temperature
#' tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
#' par(mar = c(6, 6, 3, 3))
#' # Ground temperature
#' plot(bigleafp$Tg ~ tme, type = "l", cex.axis = 2, cex.lab = 2,
#'      xlab = "Month", ylab = "Temperature", ylim = c(-5, 35),
#'      col = rgb(1, 0, 0, 0.5))
#' par(new = TRUE)
#' # Canopy temperature
#' plot(bigleafp$Tc ~ tme, type = "l", cex.axis = 2, cex.lab = 2,
#'      xlab = "", ylab = "", ylim = c(-5, 35), col = rgb(0, 0.5, 0, 0.5))
RunBigLeaf<-function(weather,vegp,groundp,soilm,lat,long,dTmx=25,zref=2,merid=0,dst=0,maxiter=100,bwgt=0.5,
                     tol=0.5,gmn=0.1,plotout=FALSE,swmethod=2,yearG=TRUE) {
  # Calculate two-stream parameters
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  if (length(soilm) == 1) soilm<-rep(soilm,length(tme))
  solar<-solarposition(lat,long,year=tme)
  twostreamp<-twostreamparams(vegp,groundp,solar)
  # Assign initial values to bigleafp
  H<-0.5*weather$swdown-vegp$em*5.67*10^-8*(weather$temp+273.15)^4
  bigleafp<-list(psih=0,psim=0,H=H,G=0,tcc=weather$temp,tcg=weather$temp,Tc=weather$temp,
                 Tg=weather$temp,phih=1,uf=0.1,LL=0,RabsG=0,Gmax=1000,Gmin=-1000)
  tst<-1
  i<-1
  while (tst > tol) {
    bigleafp1<-.BigLeaf(weather,zref,vegp,groundp,bigleafp,solar,twostreamp,soilm,lat,long,merid,dst,gmn,i,swmethod)
    # Cap values
    dTc<-bigleafp1$Tc-weather$temp
    dTg<-bigleafp1$Tg-weather$temp
    dTc[dTc>dTmx]<-dTmx
    dTg[dTg>dTmx]<-dTmx
    bigleafp1$Tc<-weather$temp+dTc
    bigleafp1$Tg<-weather$temp+dTg
    # Calculate differences
    m1<-max(abs(bigleafp$Tc-bigleafp1$Tc))
    m2<-max(abs(bigleafp$Tg-bigleafp1$Tg))
    tst<-max(m1,m2)
    # Add together
    bigleafp<-list(psih=bwgt*bigleafp$psih+(1-bwgt)*bigleafp1$psih,
                   psim=bwgt*bigleafp$psim+(1-bwgt)*bigleafp1$psim,
                   H=bwgt*bigleafp$H+(1-bwgt)*bigleafp1$H,
                   G=bwgt*bigleafp$G+(1-bwgt)*bigleafp1$G,
                   tcc=bwgt*bigleafp$tcc+(1-bwgt)*bigleafp1$tcc,
                   tcg=bwgt*bigleafp$tcg+(1-bwgt)*bigleafp1$tcg,
                   Tc=bwgt*bigleafp$Tc+(1-bwgt)*bigleafp1$Tc,
                   Tg=bwgt*bigleafp$Tg+(1-bwgt)*bigleafp1$Tg,
                   phih=bwgt*bigleafp$phih+(1-bwgt)*bigleafp1$phih,
                   uf=bwgt*bigleafp$uf+(1-bwgt)*bigleafp1$uf,
                   LL=bwgt*bigleafp$LL+(1-bwgt)*bigleafp1$LL,
                   RabsG=bigleafp1$RabsG,Gmax=bigleafp1$Gmax,Gmin=bigleafp1$Gmin)
    # Plot
    if (plotout) {
      plot(bigleafp$Tg~as.POSIXct(tme),type="l",main=round(tst,3))
    }
    if (i>maxiter) tst<-0
    i<-i+1
  }
  bigleafp<-list(tme=tme,Tc=bigleafp$Tc,Tg=bigleafp$Tg,H=bigleafp$H,G=bigleafp$G,
                 psih=bigleafp$psih,psim=bigleafp$psim,phih=bigleafp$phih,
                 OL=bigleafp$LL,uf=bigleafp$uf,RabsG=bigleafp$RabsG,error.mar=tst)
  class(bigleafp)<-"pointmicro"
  return(bigleafp)
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
weather_hgt<-function(weather, zin = 2, uzin = zin, zout = 10, lat, long) {
  bigleafp <- RunBigLeaf(weather,vegparams,groundparams,0.2,lat,long)
  d<-with(vegparams,.zeroplanedis(h,pai))
  zm<-with(vegparams,.roughlength(h,pai,d,0))
  zh<-0.2*zm
  lnr<-log((zout-d)/zh)/log((zin-d)/zh)
  # Temperature
  Tz<-(bigleafp$Tc-weather$temp)*(1-lnr)+weather$temp
  # Vapour pressure
  ea<-.satvap(weather$temp)*weather$relhum/100
  es<-.satvap(bigleafp$Tc)
  ez<-ea+(es-ea)*(1-lnr)
  es<-.satvap(Tz)
  rh<-(ez/es)*100
  # Wind speed
  lnr<-log((zout-d)/zm)/log((uzin-d)/zm)
  uz<-weather$windspeed*lnr
  weathero<-data.frame(obs_time=weather$obs_time,temp=Tz,relhum=rh,
                       pres=weather$pres,swdown=weather$swdown,
                       difrad=weather$difrad,lwdown=weather$lwdown,
                       windspeed=uz,winddir=weather$winddir,precip=weather$precip)
  return(weathero)
}
