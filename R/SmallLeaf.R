#' Calculate resistance below canopy
.rhcanopy<-function(uf,h,pai,z,phih=1) {
  d<-.zeroplanedis(h,pai)
  dF<-d/h
  a2<-0.4*(1-dF)/(phih*1.25^2)
  int<-(2*h*((48*atan((sqrt(5)*sin((pi*z)/h))/(cos((pi*z)/h)+1)))/5^(3/2)+
               (32*sin((pi*z)/h))/((cos((pi*z)/h)+1)*((25*sin((pi*z)/h)^2)/
                                                        (cos((pi*z)/h)+1)^2+5))))/pi
  inth<-4.293251*h
  s<-which(z==h)
  int[s]<-inth
  mu<-((uf/(a2*h))*(1/uf^2))
  rHa<-int*mu
  rHa
}
#' @title Generates plant area index profile
#' @description Generates a vector of length `n` of plausible plant area index values
#' @param PAI Total plant area index of canopy
#' @param skew number between 0 and 10 indicating the degree of skew towards top of
#' canopy in foliage (see details)
#' @param spread positive non-zero number less than 100 indicating the degree of spread in
#' canopy foliage (see details)
#' @param n Number of plant area index values to generate. Default: 1000 (see details)
#' @return a vector of length `n` of plant area index values the sum of which equals `PAI`
#' @details when specifying `skew`, lower numbers indicate greater skew towards top of
#' canopy (5 = symmetrical). In specifying `spread` a value of one indicates almost
#' all the foliage in concentrated in one canopy layer, whereas a value of 100 indicates
#' completely evenly spread.
#' @examples
#' pai <- PAIgeometry(5, 7, 70)
#' z <- c(1:1000) / 100
#' # plant area within each layer
#' plot(z ~ pai, type = "l", main = paste("Total PAI:", sum(pai)))
#' # foliage density
#' fd <- pai * 1000 / 10 # 1000 layers, 10 m tall
#' plot(z ~ fd, type = "l", main = paste("Total PAI:", mean(fd) * 10))
#' @rdname PAIgeometry
#' @export
PAIgeometry <- function(PAI, skew, spread, n = 1000) {
  skew<-10-skew
  # Plant area index of canopy layer
  shape1<-100/spread
  x<-c(1:n)/(n+1)
  if (skew>5) {
    shape2<-(10-skew)/5+1
    shape2<-shape2/2*shape1
    y<-rev(stats::dbeta(x,shape1,shape2))
  } else {
    shape2<-(skew+5)/5
    shape2<-shape2/2*shape1
    y<-stats::dbeta(x,shape1,shape2)
  }
  y<-PAI/sum(y)*y
  y
}
#' @title Calculates canopy wind shelter coefficient
#' @description Generates a vector of the same length as `pai` of coefficients by
#' which to multiply wind speed at the top of the canopy
#' @param h canopy height (m)
#' @param pai a vector of plant area index values as returned by [PAIgeometry()]
#' @return a vector of of the same length as `pai` of canopy wind shelter coefficients.
#' @examples
#' # plant area within each layer
#' pai1 <- PAIgeometry(5, 7, 70)
#' pai2 <- PAIgeometry(5, 5, 50)
#' z <- c(1:1000) / 100
#' uh <- 3 # Wind speed at top of canopy
#' uz1 <- uh * CanopyWindS(10, pai1)
#' uz2 <- uh * CanopyWindS(10, pai2)
#' plot(z ~ uz1, type = "l", col = "blue",
#'      xlab = "Wind speed", xlim = c(0, 3), lwd = 2)
#' par(new = TRUE)
#' plot(z ~ uz2, type = "l", col = "red",
#'      xlab = "", xlim = c(0, 3), lwd = 2)
#' @rdname CanopyWindS
#' @export
CanopyWindS<-function(h,pai) {
  # Calculate whole canopy attenuation coefficient
  Be<-0.205*sum(pai)^0.445+0.1
  a<-sum(pai)/h
  Lc<-(0.25*a)^-1
  Lm<-2*Be^3*Lc
  at<-Be*h/Lm
  # Calculate variable attenuation coefficient
  Be<-0.205*pai^0.445+0.1
  a<-pai/h
  Lc<-(0.25*a)^-1
  Lm<-2*Be^3*Lc
  ati<-Be*h/Lm
  # Adjust attenuation coefficient
  ati<-(ati/sum(ati))*at
  # Calculate canopy wind shelter coefficient
  n<-length(pai)
  n2<-trunc(n/10)
  ui<-rep(1,n)
  for (i in n:n2) {
    ui[i-1]<-ui[i]*(1-ati[i])
  }
  # Calculate bottom 10%
  z2<-(c(1:n2)/n2)*h/10
  zm<-z2[1]/2
  uf<-(0.4*ui[n2])/log(h/(10*zm))
  ui[1:n2]<-(uf/0.4)*log(z2/zm)
  return(ui)
}
#' Canopy sensible heat flux
.CanopyHL<-function(weather,forestp,solar,twostreamp,smallleafp,groundem=0.97,i,surfwet=0.98) {
  # Calculate absorbed radiation
  radp<-RadiationSmallLeaf(weather,forestp,solar,twostreamp,twostreamp,smallleafp,groundem=0.97,i)
  Rabs<-with(radp,0.5*(swupper+lwupper)+0.5*(swunder+lwunder))
  # Calculate conductance
  # Heat
  gHa<-0.135*sqrt(smallleafp$uz/(0.71*forestp$leafd))
  gFo<-with(smallleafp,0.05*(abs(tleaf-tair)/(0.71*forestp$leafd))^0.25)
  gHa<-pmax(gHa,gFo)
  # Stomatal
  Rsw<-with(radp,0.5*PARupper+0.5*PARunder)
  gS<-with(forestp,.stomcond(Rsw,gsmax,q50))
  gV<-1/(1/gHa+1/gS)
  # Penmann-Monteith
  tleaf<-with(smallleafp,.PenmanMonteith(Rabs,gHa,gV,tair,weather$pres[i],ea,forestp$em,0,surfwet))
  # cap at dewpoint
  eaa<-with(weather,.satvap(temp[i])*relhum[i]/100)
  tdew<-with(weather,.dewpoint(temp[i], eaa))
  tleaf[tleaf<tdew]<-tdew
  # Calculate H
  te<-(smallleafp$tair+tleaf)/2
  cp<-.cpair(te)
  H<-cp*gHa*(tleaf-smallleafp$tair)
  # Calculate L
  la<-(-42.575*te+44994)
  es<-.satvap(tleaf)
  L<-la*gV/weather$pres[i]*(es-smallleafp$ea)*surfwet
  L[L<0]<-0
  return(list(H=H,L=L,tleaf=tleaf,radout=radp$radout))
}
#' Run one iteration of Langrangian below canopy model
.LangrangianOne<-function(i,weather,forestp,solar,twostreamp,smallleafp,groundem=0.97,
                          surfwet=1,theta=1,a0=0.25,a1=1.25) {
  # Calculate sensible heat flux
  Ht<-.CanopyHL(weather,forestp,solar,twostreamp,smallleafp,groundem,i,surfwet)
  tleaf<-Ht$tleaf
  # Calculate thermal diffusivity
  n<-length(smallleafp$pai_a)
  d<-with(forestp,.zeroplanedis(h,pai))
  zm<-with(forestp,.roughlength(h,pai,d,smallleafp$psih))
  uf<-(0.4*smallleafp$uz[n])/(log((forestp$h-d)/zm)+smallleafp$psim)
  a2<-0.4*(1-d/forestp$h)/(smallleafp$phih*a1^2)
  TL<-a2*forestp$h/uf
  z<-(c(1:n)/n)*forestp$h
  x<-pi*(1-z/forestp$h)
  ow<-uf*(0.5*(a1+a0)+0.5*(a1-a0)*cos(x))
  KH<-TL*ow^2
  RH<-1/KH
  rHa<-cumsum(RH)*(forestp$h/n)
  # Calculate ground heat flux
  ph<-.phair(smallleafp$tair,weather$pres[i])
  cp<-.cpair(smallleafp$tair)
  GT<-((ph*cp)/rHa)*(smallleafp$tground-smallleafp$tair)
  # Calculate ground evapotranspiration
  es<-.satvap(tleaf)
  te<-(tleaf+smallleafp$tair)/2
  la<-(-42.575*te+44994)
  GL<-(la/(rHa*weather$pres[i]))*(es-smallleafp$ea)*theta
  # Add to total flux
  H<-Ht$H+GT
  L<-Ht$L+GL
  # Calculate source concentration
  ST<-(smallleafp$pai/forestp$h)*Ht$H
  SL<-(smallleafp$pai/forestp$h)*Ht$L
  # Calculate mu for correcting near-field concentration with low sample size
  if (forestp$h < 10) {
    alpha<-0.756+0.3012*forestp$h
  } else alpha<-1.31
  btm<- -0.399*forestp$h*log(1-exp(-abs(forestp$h/2-z)))/n
  btm[is.infinite(btm)]<-NA
  mu<-alpha/sum(btm,na.rm=TRUE)
  # Compute reference height near-field and far-field
  Zeta<-abs((forestp$h-z)/(ow*TL))
  kn<- -0.39894*log(1-exp(-Zeta))-0.15623*exp(-Zeta)
  CnzrT<-with(forestp,sum((ST/ow)*(kn*((h-z)/(ow*TL))+kn*((h+z)/(ow*TL))),na.rm=T))*mu
  CnzrL<-with(forestp,sum((SL/ow)*(kn*((h-z)/(ow*TL))+kn*((h+z)/(ow*TL))),na.rm=T))*mu
  CT<-0
  CL<-0
  for (ii in 1:n) {
    # Near field
    Zeta<-abs((z[ii]-z)/(ow*TL))
    kn<- -0.39894*log(1-exp(-Zeta))-0.15623*exp(-Zeta)
    CnT<-sum((ST/ow)*(kn*((z[ii]-z)/(ow*TL))+kn*((z[ii]+z)/(ow*TL))),na.rm=T)*mu
    CnL<-sum((SL/ow)*(kn*((z[ii]-z)/(ow*TL))+kn*((z[ii]+z)/(ow*TL))),na.rm=T)*mu
    # Far field
    CfsT<-sum(H[ii:n]/KH[ii:n])*(forestp$h/n)
    CfsL<-sum(L[ii:n]/KH[ii:n])*(forestp$h/n)
    CfT<-(smallleafp$tair[n]*ph[n]*cp[n])-CnzrT+CfsT
    CfL<-(smallleafp$ea[n]*ph[n]*(la[n]/weather$pres[i]))-CnzrL+CfsL
    # Both
    CT[ii]<-CfT+CnT
    CL[ii]<-CfL+CnL
  }
  tair<-CT/(cp*ph)
  ea<-(CL*weather$pres[i])/(la*ph)
  # top of canopy
  tair[n]<-smallleafp$tair[n]
  ea[n]<-smallleafp$ea[n]
  # Limits
  mn<-min(smallleafp$tground,smallleafp$tair,tleaf)
  mx<-max(smallleafp$tground,smallleafp$tair,tleaf)
  # Apply Limits
  tair[tair<mn]<-mn
  tair[tair>mx]<-mx
  emx<-.satvap(tair)
  ea<-pmin(ea,emx)
  ea[ea<0.01]<-0.01
  # Return values
  smallleafp$tair<-tair
  smallleafp$ea<-ea
  smallleafp$tleaf<-tleaf
  return(list(smallleafp=smallleafp,radout=Ht$radout))
}
#' Initializes smallleafp
.SmallLeafInit <- function(th, eh, uh, forestp, bigleafp, n, i = 1) {
  pais<-with(forestparams,PAIgeometry(pai, skew, spread, n))
  d<-with(forestp,.zeroplanedis(h,pai))
  zm<-with(forestp,.roughlength(h,pai,d,bigleafp$psih[i]))
  zh<-0.2*zm
  smallleafp=list(pai=pais,
                  pai_a=rev(cumsum(rev(pais))),
                  lwgt=.paiwgt(pais),
                  tair=rep(th[i],n),
                  tleaf=rep(th[i],n),
                  uz=CanopyWindS(forestp$h,pais)*uh[i],
                  ea=rep(eh[i],n),
                  psih=.dpsih(zm/bigleafp$OL[i])-.dpsih((forestp$h-d)/bigleafp$OL[i]),
                  psim=.dpsim(zm/bigleafp$OL[i])-.dpsim((forestp$h-d)/bigleafp$OL[i]),
                  phih=.dphih((forestp$h-d)/bigleafp$OL[i]),
                  tground=bigleafp$Tg[i])
  class(smallleafp)<-"smallleafparams"
  return(smallleafp)
}
#' Assign values to smallleafp for timestep i
.SmallLeafAs <- function(smallleafp, th, eh, uh, h, d, zm, ws, bigleafp, n, i) {
  smallleafp$tair[n]<-th[i]
  smallleafp$ea[n]<-eh[i]
  smallleafp$uz<-ws*uh[i]
  smallleafp$psih<-.dpsih(zm/bigleafp$OL[i])-.dpsih((h-d)/bigleafp$OL[i])
  smallleafp$psim<-.dpsim(zm/bigleafp$OL[i])-.dpsim((h-d)/bigleafp$OL[i])
  smallleafp$phih<-.dphih((h-d)/bigleafp$OL[i])
  smallleafp$tground<-bigleafp$Tg[i]
return(smallleafp)
}

#' Run one step of small leaf model
.RunSmallLeafOne<-function(i,weather,forestp,solar,twostreamp,smallleafp,lat,long,groundem=0.97,
                          surfwet=1,theta=1,a0=0.25,a1=1.25,iter=10,bwgt=0.75) {
  for (ii in 1:iter) {
    sprad<-.LangrangianOne(i,weather,forestp,solar,twostreamp,smallleafp,groundem,surfwet,theta,a0,a1)
    smallleafp1<-sprad$smallleafp
    smallleafp$tair<-(1-bwgt)*smallleafp1$tair+(bwgt)*smallleafp$tair
    smallleafp$tleaf<-(1-bwgt)*smallleafp1$tleaf+(bwgt)*smallleafp$tleaf
    smallleafp$ea<-(1-sbwgt)*smallleafp1$ea+(sbwgt)*smallleafp$ea
  }
  return(list(smallleafp=smallleafp,radout=sprad$radout))
}
#' @title Runs small leaf model to derive below canopy temperatures
#' @description Runs small leaf model for each time-step in turn to derive below
#' canopy temperature, humidity and radiation fluxes.
#' @param weather a data.frame of hourly weather as in the example dataset `climdata`
#' @param forestp an object of class forestparams formatted as for the inbuilt example dataset `forestparams`
#' @param groundp an object of class vegparams formatted as for the inbuilt example dataset `groundparams`
#' @param reqhgt height (m) above ground for which microclimate variables are required
#' @param soilm a vector of volumetric soil moisture fractions in the top 10 cm of soil. Calculated if not supplied.
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long longitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @param n number of layers to divide canopy into (higher = more accurate but slower). Default 100.
#' @param dTmx maximum by which vegetation or soil surface temperature can exceed air temperature (deg C, set to ensure model convergence)
#' @param zref height above ground (m) of measurements in `weather`
#' @param merid Longitude of local time zone meridian (decimnal degrees), Default: 0 (Greenwich Mean Time)
#' @param dst Daylight saving hours. E.g. with `merid = 0` 0 for UTC or 1  for British Summer Time.
#' @param maxiter Maximum number of iterations over which to run the big leaf model
#' @param smalliter number of iterations over which to run the small leaf model (default 10)
#' @param bwgt backward weighting to apply when iteratively running the big leaf model model (default 0.5)
#' @param sbwgt backward weighting to apply when iteratively running the small leaf model model (default 0.5)
#' @param tol Error margin (deg C) deemed acceptable for big leaf model model convergence (default 0.5).
#' @param gmn optional minimum convective conductance value (mol/m^2/s) for big leaf model. Lower values are more
#' physically realistic, but can reduce likelihood of model convergence.
#' @param a0 thermal diffusivity scale parameter from Raupach 1989 Agr Forest Meteorol, 47: 85-108.
#' @param a1 thermal diffusivity shape parameter from Raupach 1989 Agr Forest Meteorol, 47: 85-108.
#' @param plotout optional logical indicating whether to plot vertical profile of air temperatures every 12 hours (default FALSE)
#' temperatures on each iteration
#' @param swmethod method used to calculate fraction of ground surface acting as free water surface in big leaf model
#' (0 = based on rainfall, 1 computed from soil effective relative humidity, 2 computed from soil moisture fraction)
#' @param yearG optional logical indicating whether or not to calculate and account for annual ground heat flux cycle in big leaf model
#' @param saveprofile optional logical indicating whether to return microclimate variables for all canopy layers (TRUE) or
#' or just for height `reqhgt` (FALSE). Default FALSE.
#' @return a a list of the following:
#' \describe{
#'  \item{tme}{POSIXlt object of times and dates corresponding to model outputs}
#'  \item{tair}{a vector or matrix of air temperatures (deg C)}
#'  \item{tleaf}{a vector or matrix of leaf temperatures (deg C)}
#'  \item{relhum}{a vector or matrix of relative humidities (Percentage)}
#'  \item{windspeed}{a vector or matrix of wind speeds (m/s)}
#'  \item{Rdni_down}{a vector or matrix of downward direct normal irradiance fluxes (W/m^2)}
#'  \item{Rdif_down}{a vector or matrix of downward diffuse radiation fluxes (W/m^2)}
#'  \item{Rsw_up}{a vector or matrix of upward shortwave radiation fluxes (W/m^2), assumed entirely diffuse}
#'  \item{Rlw_down}{a vector or matrix of downward longward radiation fluxes (W/m^2)}
#'  \item{Rlw_up}{a vector or matrix of upward longwave radiation fluxes (W/m^2)}
#' }
#' @details If `saveprofile = TRUE` the returned values are matrices, with rows corresponding
#' to each canopy layer and columns corresponding to each time step. If `saveprofile = FALSE`
#' returned values are vectors for height `reqhgt`, with elements corresponding to each time
#' step.
#' @examples
#' # Select first 5 days of May from inbuilt climate dataset
#' tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
#' sel <- which(tme$mon + 1 == 5 & tme$mday < 6)
#' weather <- climdata[sel,]
#' # Run model with inbuilt vegetation and ground parameters and default settings
#' smallleafp <- RunSmallLeaf(weather, forestparams, groundparams, reqhgt = 5,
#'                            lat = 49.96807, long = -5.215668, yearG = FALSE)
#' # Plot air temperature at reqhgt
#' plot(smallleafp$tair ~ as.POSIXct(tme[sel]), type = "l",
#'      xlab = "Day", ylab = "Temperature")
#' @rdname RunSmallLeaf
#' @export
RunSmallLeaf <- function(weather,forestp,groundp,reqhgt,lat,long,zref=2,n = 10,soilm=NA,dTmx=25,merid=0,dst=0,maxiter=100,
                         smalliter=10,bwgt=0.5,sbwgt=0.5,tol=0.5,gmn=0.1,a0=0.25,a1=1.25,plotout=FALSE,swmethod=2,yearG=TRUE,
                         saveprofile=FALSE) {
  # Run soil moisture model
  if (class(soilm) == "logical") {
    soiltype<-.soilmatch(groundp)
    soilm<-soilmmodel(weather, soiltype)
  }
  if (zref < forestp$h) {
    weather<-weather_hgt(weather,zref,zref,forestp$h,lat,long)
    zref<-forestp$h
  }
  # Calculate theta
  theta<-(soilm-groundp$Smin)/(groundp$Smax-groundp$Smin)
  bigleafp <- RunBigLeaf(weather,forestp,groundp,soilm,lat,long,dTmx,zref,merid,dst,maxiter,bwgt,
                         tol,gmn,plotout=FALSE,swmethod,yearG)
  # Get twostream parameters
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  solar <- solarposition(lat,long,year=tme,merid=merid,dst=dst)
  twostreamp<-twostreamparams(vegp,groundp,solar)
  # adjust to top of canopy
  d<-with(forestp,.zeroplanedis(h,pai))
  zm<-with(forestp,.roughlength(h,pai,d,bigleafp$psih))
  if (zref > forestp$h) {
    d<-with(forestp,.zeroplanedis(h,pai))
    zm<-with(forestp,.roughlength(h,pai,d,bigleafp$psih))
    zh<-0.2*zm
    lnr<-log((forestp$h-d)/zh)/log((zref-d)/zh)
    # Temperature
    th<-weather$temp+(bigleafp$Tc-weather$temp)*(1-lnr)
    # Vapour pressure
    ea<-with(weather,.satvap(temp)*(relhum/100))
    eh<-ea+(.satvap(bigleafp$Tc)-ea)*(1-lnr)
    # Wind speed
    psim<-.dpsim(zm/bigleafp$OL)-.dpsim((forestp$h-d)/bigleafp$OL)
    uh<-(uf/0.4)*(log((forestp$h-d)/zm)+psim)
  } else {
    th<-weather$temp
    eh<-with(weather,.satvap(temp)*(relhum/100))
    uh<-weather$windspeed
  }
  # Intialise
  smallleafp<-.SmallLeafInit(th, eh, uh, forestp, bigleafp, n, 1)
  # which layer = reqhgt
  z<-(c(1:n)/n)*vegp$h
  s<-which.min(abs(reqhgt-z))
  # values to save
  pais<-with(vegp,PAIgeometry(pai, skew, spread, n))
  ws<-CanopyWindS(vegp$h,pais)
  if (saveprofile) {
    tleaf<-matrix(0,ncol=length(tme),nrow=n)
    tair<-tleaf
    ea<-tleaf
    Rdnidn<-tleaf
    Rdifdn<-tleaf
    Rswup<-tleaf
    Rlwdn<-tleaf
    Rlwup<-tleaf
    uz<-tleaf
    for (i in 1:length(tme)) uz[,i]<-ws*uh[i]
  } else {
    tleaf<-0
    tair<-0
    ea<-0
    Rdnidn<-0
    Rdifdn<-0
    Rswup<-0
    Rlwdn<-0
    Rlwup<-0
    uz<-ws[s]*uh
  }
  for (i in 1:length(tme)) {
    sprad<-.RunSmallLeafOne(i,weather,forestp,solar,twostreamp,smallleafp,lat,long,groundp$em,
                            surfwet,theta[i],a0,a1,smalliter,sbwgt)
    # Extract values for height z
    smallleafp<-sprad$smallleafp
    radout<-sprad$radout
    if (saveprofile) {
      tleaf[,i]<-smallleafp$tleaf
      tair[,i]<-smallleafp$tair
      ea[,i]<-smallleafp$ea
      uz[,i]<-smallleafp$uz
      Rdnidn[,i]<-radout$Rdni_down
      Rdifdn[,i]<-radout$Rdif_down
      Rswup[,i]<-radout$Rdif_up
      Rlwdn[,i]<-radout$Rlw_down
      Rlwup[,i]<-radout$Rlw_up
    } else {
      tleaf[i]<-smallleafp$tleaf[s]
      tair[i]<-smallleafp$tair[s]
      ea[i]<-smallleafp$ea[s]
      uz[i]<-smallleafp$uz[s]
      Rdnidn[i]<-radout$Rdni_down[s]
      Rdifdn[i]<-radout$Rdni_down[s]
      Rswup[i]<-radout$Rdif_up[s]
      Rlwdn[i]<-radout$Rlw_down[s]
      Rlwup[i]<-radout$Rlw_up[s]
    }
    # Reassign values for next model run
    if (plotout & i%%12 == 0) {
      plot(z~smallleafp$tair, type = "l", xlab = "Air temperature", ylab = "Height (m)", main = tme[i])
    }
    smallleafp<-.SmallLeafAs(smallleafp, th, eh, uh, forestp$h, d, zm[i], ws, bigleafp, n, i+1)
  }
  # Calculate relative humidity
  relhum<-(ea/.satvap(tair))*100
  return(list(tair=tair,tleaf=tleaf,relhum=relhum,windspeed=uz,
              Rdni_down=Rdnidn,Rdif_down=Rdifdn,Rsw_up=Rswup,Rlw_down=Rlwdn,Rlw_up=Rlwup))
}
