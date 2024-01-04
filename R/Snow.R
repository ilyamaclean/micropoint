#' Calculate zero-plane displacement for snow
.szeroplanedis<-function(h,pai,snowdepth) {
  if (snowdepth > h) {
    d<-snowdepth
  } else {
    pai<-(h-snowdepth)/h
    pai[pai<0.001]<-0.001
    hs<-h-snowdepth
    d<-snowdepth+(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*hs
  }
  d
}
#' Calculate roughness length for snow
.sroughlength<-function(h,pai,snowdepth,d=NA,psi_h=0) {
  if (class(d)[1]=="logical") d<-.szeroplanedis(h,pai,snowdepth)
  if  (snowdepth > h) {
    zm<-0.0005
  } else {
    pai<-(h-snowdepth)/h
    Be<-sqrt(0.003+(0.2*pai)/2)
    zm<-(h-d)*exp(-0.4/Be)*exp(psi_h)
    zm[zm<0.0005]<-0.0005
  }
  zm
}
#' Calculates snow albedo
.snowalb <- function(age) {
  snowalb<-(-9.8740*log(age/24)+78.3434)/100
  snowalb[snowalb>0.95]<-0.95
  return(snowalb)
}
#' Calculates snow density parameters
.snowdens<-function(snowenv="Tundra") {
  densfun<-c(0.5975,0.2237,0.0012,0.0038)
  if (snowenv == "Maritime") densfun<-c(0.5979,0.2578,0.001,0.0038)
  if (snowenv == "Prairie") densfun<-c(0.594,0.2332,0.0016,0.0031)
  if (snowenv == "Tundra") densfun<-c(0.363,0.2425,0.0029,0.0049)
  if (snowenv == "Taiga") densfun<-c(0.217,0.217,0,0)
  densfun
}
#' Derives steady-state temperature and energy balance
.steadystate<-function(Rabs,gHa,tc,pk,rh,G=0) {
  .PMon2<-function(Rabs,gHa,tc,pk,ea,te) {
    sb<-5.67*10^-8
    Rema<-0.97*sb*(tc+273.15)^4
    la<-45068.7-42.8428*te
    sel<-which(te<0)
    la[sel]<-51078.69-4.338*te[sel]-0.06367*te[sel]^2
    cp<-.cpair(te)
    Da<-.satvap(tc)-ea
    gR<-(4*0.97*sb*(te+273.15)^3)/cp
    De<-.satvap(te+0.5)-.satvap(te-0.5)
    Ts<-tc+((Rabs-Rema-la*(gHa/pk)*Da-G)/(cp*(gHa+gR)+la*(gHa/pk)*De))
    Ts
  }
  ea<-.satvap(tc)*rh/100
  Ts<-.PMon2(Rabs,gHa,tc,pk,ea,tc)
  for (i in 2:3) {
    te<-(Ts+tc)/2
    Ts<-.PMon2(Rabs,gHa,tc,pk,ea,te)
  }
  tmn<-.dewpoint(tc, ea)
  Ts<-pmax(Ts,tmn)
  # Energy balance
  la<-45068.7-42.8428*te
  sel<-which(te<0)
  la[sel]<-51078.69-4.338*te[sel]-0.06367*te[sel]^2
  es<-.satvap(Ts)
  cp<-.cpair(te)
  return(list(Ts=Ts,Rem=0.97*5.67*10^-8*(Ts+273.15)^4,
              H=cp*gHa*(Ts-tc),L=la*gHa/pk*(es-ea),la=la))
}
#' @title RUns a simple snow model
#' @description Calculates snow depth and temperature
#' @param weather a data.frame of hourly weather as in the example dataset `climdata`
#' @param vegp an object of class vegparams formatted as for the inbuilt example dataset `vegparams`
#' @param initdepth snow depth (m) at start of model run
#' @param initdepth snow temperature (deg C) at start of model run
#' @param snowenv either one of 'Tundra', 'Taiga', 'Prairie', 'Alpine' or 'Maritime'
#' used for estimating snow density as function of snowpack age following the exponential
#' function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394. Alternatively a vector of
#' parameters for the Sturm et al function.
#' @param zref height above ground of measurements in `weather`. Adjust using [WeatherHeight()]
#' if `zref` < `vegp$h`.
#' @examples
#' snowm <- pointsnow(snowclimdata, vegparams)
#' plot(snowm$snowdepth ~ as.POSIXct(snowclimdata$obs_time), type = "l",
#'      xlab = "Date", ylab = "Depth (m)")
#' @rdname pointsnow
#' @return a list of the following:
#' \describe{
#' \item{tme}{POSIXlt object of times and dates corresponding to model outputs}
#'  \item{snowdepth}{a vector of snow depths (m)}
#'  \item{snowtemp}{a vector of snow temperatures (deg C)}
#'  \item{snowalb}{a vector of snow albedos (0-1)}
#'  \item{snowdens}{a vector of snow densities (kg/m^3)}
#'  \item{psim}{diabatic correction factors for heat}
#'  \item{psih}{diabatic correction factors for momentum}
#' }
#' @export
pointsnow<-function(weather, vegp, initdepth = 0.01, inittemp = 0, snowenv = "Taiga", zref = 2) {
  if (length(snowenv) == 1) {
    densfun<-.snowdens(snowenv) # Snow density function
  } else densfun<-snowenv

  # Initialise values
  snowdepth<-initdepth
  snowtemp<-inittemp
  snowage<-1
  psim<-0
  psih<-0
  snowalb<-0
  snowdens<-rep(1,length(weather$temp))
  for (i in 1:(length(weather$temp)-1)) {
    # Calculate Absorbed radiation
    snowalb[i]<-.snowalb(snowage)  # Snow albedo
    Rabs<-(1-snowalb[i])*weather$swdown[i]+0.97*weather$lwdown[i] # Absorbed radiation
    # Calculate conductivity
    d<-with(vegp,.szeroplanedis(h,pai,snowdepth[i]))
    zm<-with(vegp,.sroughlength(h,pai,snowdepth[i],d,psih[i]))
    uf<-(0.4*weather$windspeed[i])/(log((zref-d)/zm)+psim[i])
    ph<-with(weather,.phair(temp[i],pres[i]))
    gHa<-.gturb(uf,d,zm,zref,z0=NA,ph,psih[i],0.1)
    # Steady state model
    ss<-.steadystate(Rabs,gHa,weather$temp[i],weather$pres[i],weather$relhum[i])
    if (snowdepth[i] > 0) {
      # Calculate snow density
      pdensity<-((densfun[1]-densfun[2])*(1-exp(-densfun[3]*snowdepth[i]*100
                                                -densfun[4]*snowage))+densfun[2])*1000
      # Transient state
      Rem<-0.97*5.67*10^-8*(snowtemp[i]+273.15^4)
      # Calculate sensible heat
      cp<-.cpair(weather$temp[i])
      H<-cp*gHa*(snowtemp[i]-weather$temp[i])
      # Calculate latent heat
      if (snowtemp[i] > 0) {
        la<-45068.7-42.8428*snowtemp[i]
      } else la<-51078.69-4.338*snowtemp[i]-0.06367*snowtemp[i]^2
      es<-.satvap(snowtemp[i])
      ea<-with(weather,.satvap(temp[i])*relhum[i]/100)
      L<-la/weather$pres[i]*gHa*cp*(es-ea)
      # Transient energy budget
      EB<-Rabs-Rem-H-L
      # Determine whether to use steady-state or transient state model
      dTS<-ss$Ts-weather$temp[i]
      dTT<-(EB*3600)/(2108*pdensity*snowdepth[i])
      ssc<-0 # steady state check
      if (abs(dTT) > abs(dTS)) { # Steady state
        subl<-ss$L*3600*0.01802/ss$la # kg/m^2 (or mm SWE) # sublimation
        # Excess energy above zero
        EAZ<-snowdepth[i]*pdensity*2108*ss$Ts # J
        EAZ[EAZ<0]<-0
        melt<-EAZ/314000 # kg/m^2 (or mm SWE)
        # Temperature change due to melt
        dT<-EAZ/(2108*snowdepth[i]*pdensity)
        snowtemp[i+1]<-ss$Ts-dT
        tmn<-.dewpoint(snowtemp[i+1], ea)
        snowtemp[i+1]<-max(snowtemp[i+1],tmn)
        snowtemp[i+1]<-min(snowtemp[i+1],5)
        ssc[i]<-1
        LL<-(ph*cp*uf^3*(weather$temp[i]+273.15))/(-0.4*9.81*ss$H) # Stability
      } else { # Transient state
        subl<-L*3600*0.01802/la
        snowtemp[i+1]<-snowtemp[i]+dTT
        tmn<-.dewpoint(snowtemp[i+1], ea)
        snowtemp[i+1]<-max(snowtemp[i+1],tmn)
        EAZ<-Rem-H-L # Excess energy above zero
        EAZ[EAZ<0]<-0
        melt<-EAZ/314000 # kg/m^2 (or mm SWE)
        # Temperature change due to melt
        dT<-EAZ/(2108*snowdepth[i]*pdensity)
        snowtemp[i+1]<-snowtemp[i+1]-dT
        tmn<-.dewpoint(snowtemp[i+1], ea)
        snowtemp[i+1]<-max(snowtemp[i+1],tmn)
        snowtemp[i+1]<-min(snowtemp[i+1],5)
        ssc[i]<-0
        LL<-(ph*cp*uf^3*(weather$temp[i]+273.15))/(-0.4*9.81*H) # Stability
      }
      # Rain melt
      rainmelt<-with(weather,temp[i]*precip[i]*0.0125)
      rainmelt[rainmelt<0]<-0
      snowage<-snowage+1
    } else {
      subl<-0
      rainmelt<-0
      melt<-0
      snowtemp[i+1]<-0
      snowage<-1
    }
    # Snowfall
    snow<-weather$precip[i]
    snow[weather$temp[i]>1.5]<-0
    # Snow balance
    snowbalance<-snow-melt-rainmelt-subl
    snowdepth[i+1]<-snowdepth[i]+snowbalance/pdensity
    snowdepth[i+1]<-ifelse(snowdepth[i+1]<0,0,snowdepth[i+1])
    # Calculate diabatic coefficients
    psim[i+1]<-.dpsim(zm/LL)-.dpsim((zref-d)/LL)
    psih[i+1]<-.dpsih(zm/LL)-.dpsih((zref-d)/LL)
    snowdens[i]<-pdensity
  }
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  return(list(tme=tme,snowdepth=snowdepth,snowtemp=snowtemp,snowalb=snowalb,snowdens=snowdens,psim=psim,psih=psih))
}
