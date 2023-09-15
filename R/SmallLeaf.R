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
.CanopyHL<-function(weather,forestp,solar,twostreamp,smallleafp,tground,groundem=0.97,i,surfwet=0.98) {
  # Calculate absorbed radiation
  radp<-RadiationSmallLeaf(weather,forestp,solar,twostreamp,twostreamp,smallleafp,tground,groundem=0.97,i)
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
  tleaf<-with(smallleafp,.PenmanMonteith(Rabs,gHa,gV,tair,weather$pres[i],ea,forestp$em,0,3,surfwet,TRUE))
  # cap at dewpoint
  ea<-with(weather,.satvap(temp[i])*relhum[i]/100)
  tdew<-with(weather,.dewpoint(temp[i], ea))
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
  return(list(H=H,L=L,tleaf=tleaf))
}
#' Run one iteration of Langrangian below canopy model
.LangrangianOne<-function(i,weather,forestp,solar,twostreamp,smallleafp,tground,groundem=0.97,
                          surfwet=1,theta=1,a0=0.25,a1=1.25) {
  # Calculate sensible heat flux
  Ht<-.CanopyHL(weather,forestp,solar,twostreamp,smallleafp,tground,groundem,i,surfwet)
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
  GT<-((ph*cp)/rHa)*(tground-smallleafp$tair)*0.5
  # Calculate ground evapotranspiration
  es<-.satvap(tleaf)
  te<-(tleaf+smallleafp$tair)/2
  la<-(-42.575*te+44994)
  GL<-(la/(rHa*weather$pres[i]))*(es-smallleafp$ea)*theta*0.5
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
  mn<-min(tground,smallleafp$tair[n])-5
  mx<-max(tground,smallleafp$tair[n])+5
  # Apply Limits
  tair[tair<mn]<-mn
  tair[tair>mx]<-mx
  ea[ea<0.01]<-0.01
  # Return values
  smallleafp$tair<-tair
  smallleafp$ea<-ea
  smallleafp$tleaf<-tleaf
  return(smallleafp)
}

#' @title Run one step of small leaf model
#' @description Runs a single step of the small leaf model
#' @param i time increment of model
#' @param weather a data.frame of hourly weather as in the example dataset `climdata`
#' @param forestp an object of class vegparams formatted as for the inbuilt example dataset `forestparams`
#' @param twostreamp parameters of two-stream radiation model as returned by [twostreamparams()]
#' @param smallleafp a list returned by [RunSmallLeafOne()] in previous time increment
#' @param tground ground surface temperature (deg C)
#' @param lat latitude (decimal degrees). Negative south of the equator.
#' @param long longitude (decimal degrees). Negative west of the Greenwich Meridian.
#' @param groundem ground emissivity temperature (deg C)
#' @param surfwet fraction of canopy surface acting like a saturated evaporating surface
#' @param theta fraction of ground surface acting like a saturated evaporating surface
#' @param a0 optional coefficient controlling vertical variation in eddy variance
#' @param a1 optional coefficient controlling vertical variation in eddy variance
#' @param iter number of iterations over which to run model to stablise outputs
#' @return A list of the following:
#' (1) tleaf - a vector representing vertical variation in leaf temperatures (deg C)
#' (2) tair - a vector representing vertical variation in air temperatures a vector representing vertical variation in leaf temperatures (deg C)
#' (3) pai - a vector representing vertical variation in plant area index values
#' (4) pai_a - a vector representing vertical variation in plant area index values above each point of interest
#' (5) uz - a vector representing vertical variation in wind speed (m/s)
#' (6) ea - a vector representing vertical variation in vapour pressure (kPa)
#' (7) psim - diabatic correction coefficient for momentum at canopy top
#' (8) psih - diabatic correction coefficient for heat at canopy top
#' (9) phih - diabatic influencing factor at canopy top
#' (10) lwgt - a list of weightings used for calculating longwave radiation
#' @rdname RunSmallLeafOne
#' @export
RunSmallLeafOne<-function(i,weather,forestp,twostreamp,smallleafp,tground,lat,long,groundem=0.97,
                          surfwet=1,theta=1,a0=0.25,a1=1.25,iter=10) {
  tme<-as.POSIXlt(weather$obs_time, tz = "UTC")
  solar<-solarposition(lat,long,year=tme)
  for (ii in 1:iter) {
    smallleafp<-.LangrangianOne(i,weather,forestp,solar,twostreamp,smallleafp,tground,groundem,surfwet,theta,a0,a1)
  }
  return(smallleafp)
}
