#' Calculate weights for below-canopy longwave radiation
.paiwgt<-function(pai) {
  n<-length(pai)
  # Tranmission from each node to every other node
  tr<-matrix(0,ncol=n,nrow=n)
  tr[1,]<-exp(-cumsum(pai))
  tr[n,]<-exp(-rev(cumsum(pai)))
  tr[n-1,]<-exp(-c(rev(cumsum(pai[1:(n-1)])),0))
  for (i in 2:(n-2)) {
    tr[i,]<-exp(-c(rev(cumsum(pai[1:(i-1)])),0,cumsum(pai[(i+1):n])))
  }
  # Calculate what total weighting from foliage should be
  paia<-rev(cumsum(pai))
  paib<-cumsum(pai)
  trg<-exp(-paib)
  trh<-exp(-paia)
  trf<-2-trg-trh
  # multiply transmission by source density  and adjust for canopy weighting
  wgt<-tr*0
  for (i in 1:n) {
    xx<-tr[i,]*pai
    wgt[i,]<-(xx/sum(xx))*trf[i]
  }
  return(list(wgt=wgt,trg=trg,trh=trh))
}
#' @title Calculates solar position
#' @description Calculates the solar zenith and azimuth angles
#' @param lat latitude (decimal degrees)
#' @param long longitude(decimal degrees, -ve west of Greenwich meridian)
#' @param year year or POSIXlt object of times
#' @param month month. Ignored if year is a POSIXlt object
#' @param day day of month. Ignored if year is a POSIXlt object
#' @param hour decimal hour. Ignored if year is a POSIXlt object
#' @param merid meridian of local time zone (0 for UTC)
#' @param dst daylight saving yours (e.g. one for British Summer Time)
#' @return a list of zenith and azimuth angles (decimal degrees)
#' @rdname solarposition
#' @export
solarposition<-function(lat,long,year,month,day,hour=12,merid=0,dst=0) {
  # Extract year, month, day and hour if year is a POSIXlt object
  if (class(year)[1] == "POSIXlt") {
    tme<-year
    year<-tme$year+1900
    month<-tme$mon+1
    day<-tme$mday
    hour<-with(tme,hour+(min+sec/60)/60)
  }
  # Calculate Astronomical Julian day
  dd<-day+0.5 # decimal day
  ma<-month+(month<3)*12 # adjusted month
  ya<-year+(month<3)*-1 # adjusted year
  jd<-trunc(365.25*(ya+4716))+trunc(30.6001*(ma+1))+dd-1524.5
  B<-(2-trunc(ya/100)+trunc(trunc(ya/100)/4))
  jd<-jd+(jd>2299160)*B
  # Calculate solar time
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-hour+(4*long+eot)/60
  # Calculate solar zenith
  latr<-lat*pi/180 # radians
  tt<-0.261799*(st-12)
  dec<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  coh<-sin(dec)*sin(latr)+cos(dec)*cos(latr)*cos(tt)
  z<-acos(coh)*(180/pi)
  # Calculate solar azimuth
  sh<-sin(dec)*sin(latr)+cos(dec)*cos(latr)*cos(tt)
  hh<-(atan(sh/sqrt(1-sh^2)))
  sazi<-cos(dec)*sin(tt)/cos(hh)
  cazi<-(sin(latr)*cos(dec)*cos(tt)-cos(latr)*sin(dec))/
    sqrt((cos(dec)*sin(tt))^2+(sin(latr)*cos(dec)*cos(tt)-cos(latr)*sin(dec))^2)
  sqt<-1-sazi^2
  sqt[sqt<0]<-0
  azi<-180+(180*atan(sazi/sqrt(sqt)))/pi
  azi[cazi<0 & sazi<0]<-180-azi[cazi<0 & sazi<0]
  azi[cazi<0 & sazi>=0]<-540-azi[cazi<0 & sazi>=0]
  solar<-list(zen=z,azi=azi)
  return(solar)
}
#' Internal function for calculating solar index
.solarindex<- function(slope,aspect,solar) {
  i<-with(solar,cos(zen*pi/180)*cos(slope*pi/180)+sin(zen*pi/180)*sin(slope*pi/180)*cos((azi-aspect)*pi/180))
  i[i<0]<-0
  i[solar$zen>90]<-0
  i
}
#' Internal function for calculating canopy extinction coefficient
.cank<-function(zen,x) {
  zen[zen>90]<-90
  Z<-zen*pi/180
  if (x==1) {
    k<-1/(2*cos(Z))
  } else if (x==0) {
    k<-2/(pi*tan(0.5*pi-Z))
  } else if (is.infinite(x)) {
    k<-1
  } else {
    k<-sqrt((x^2+(tan(Z)^2)))/(x+1.774*(x+1.182)^(-0.733))
  }
  k[k>6000]<-6000
  k
}
#' Internal function for calculating canopy extinction coefficients for sloped ground surfaces
.cank2 <- function(zen,x,si) {
  # Raw k
  k<-.cank(zen,x)
  k0<-sqrt(x^2)/(x+1.774*(x+1.182)^(-0.733))
  # k dash
  kd<-k*cos(zen*pi/180)/si
  sel<-which(si==0)
  kd[sel]<-1
  return(list(k=k,kd=kd,k0=k0))
}
#' @title Calculates two-stream radiation model parameters
#' @description Calculates parameters of Dickinson-Sellers two-stream
#' radiation model
#' @param vegp an object of class vegparams formatted as for the inbuilt example dataset `vegparams`
#' @param groundp an object of class vegparams formatted as for the inbuilt example dataset `groundparams`
#' @param solar a list of zenith and azimuth angles as returned by [solarposition()]
#' @return a list of parameters for two-stream approximation equations
#'
#' @details For details of the Dickinson-Sellers two-stream radiation
#' model see Sellers (1985) https://doi.org/10.1080/01431168508948283.
#' The function implements improvements to the model proposed by Yuan et
#' al (2017) https://doi.org/10.1002/2016MS000773
#'
#' @examples
#' # Run two-sream model with inbuilt datasets
#' tme <- as.POSIXlt(climdata$obs_time, tz="UTC")
#' # POSIXlt object used in place of year:
#' solar <- solarposition(49.96807, -5.215668, year=tme)
#' twostreamp <- twostreamparams(vegparams, groundparams, solar)
#' # Calculate white- and black-sky albedo
#' cl <- vegparams$clump
#' Kc <- with(twostreamp, kkd$kd / kkd$k0)
#' albd <- with(twostreamp, p1 + p2 + cl^2 * groundparams$gref)
#' albb <- with(twostreamp, p5 / sig + p6 + p7 + cl^Kc * gref2)
#' # Calculate and plot blue-sky albedo
#' alb <- with(climdata, (albd * difrad + albb * (swdown - difrad)) / swdown)
#' # Set blue-sky albedo to be white-sky albedo when direct radiation is zero
#' dni<-with(climdata, swdown - difrad)
#' alb[dni == 0] <- albd
#' par(mar = c(6, 6, 3, 3))
#' plot(alb ~ as.POSIXct(tme), type="l", cex.axis = 2, cex.lab = 2,
#'      xlab = "Month", ylab = "Albedo", ylim = c(0, 1))
#' @rdname twostreamparams
#' @export
twostreamparams<-function(vegp,groundp,solar) {
  # Calculate solar index
  si<-with(groundp,.solarindex(slope,aspect,solar))
  # === (1a) Calculate canopy k
  kkd<-.cank2(solar$zen,vegp$x,si)
  # === (1c) Adjust paramaters for gap fraction and inclined surface
  pai_t<-with(vegp,(pai/(1-clump)))
  gref2<-1-(((1-groundp$gref)*cos(solar$zen*pi/180))/si)
  gref2[is.infinite(gref2)]<-groundp$gref
  # === (1d) Calculate two stream base parameters
  om<-with(vegp,lref+ltra)
  a<-1-om
  del<-with(vegp,lref-ltra)
  mla<-(9.65*(3+vegp$x)^(-1.65))
  mla[mla>pi/2]<-pi/2
  J<-cos(mla)^2
  gma<-0.5*(om+J*del)
  s<-0.5*(om+J*del/kkd$k)*kkd$k
  sstr<-om*kkd$k-s
  # === (1e) Calculate two stream base parameters
  h<-sqrt(a^2+2*a*gma)
  sig<-kkd$kd^2+gma^2-(a+gma)^2
  S1<-exp(-h*pai_t)
  S2<-exp(-kkd$kd*pai_t)
  u1<-a+gma*(1-1/groundp$gref)
  u2<-a+gma*(1-groundp$gref)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  # === (1f) Calculate Diffuse radiation parameters
  p1<-(gma/(D1*S1))*(u1-h)
  p2<-(-gma*S1/D1)*(u1+h)
  p3<-(1/(D2*S1))*(u2+h)
  p4<-(-S1/D2)*(u2-h)
  # === (1f) Calculate Direct radiation parameters
  u1<-a+gma*(1-1/gref2)
  u2<-a+gma*(1-gref2)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  p5<- -s*(a+gma-kkd$kd)-gma*sstr
  v1<-s-(p5*(a+gma+kkd$kd))/sig
  v2<-s-gma-(p5/sig)*(u1+kkd$kd)
  p6<-(1/D1)*((v1/S1)*(u1-h)-(a+gma-h)*S2*v2)
  p7<-(-1/D1)*((v1*S1)*(u1+h)-(a+gma+h)*S2*v2)
  p8<-sstr*(a+gma+kkd$kd)-gma*s
  v3<-(sstr+gma*gref2-(p8/sig)*(u2-kkd$kd))*S2
  p9<-(-1/D2)*((p8/(sig*S1))*(u2+h)+v3)
  p10<-(1/D2)*(((p8*S1)/sig)*(u2-h)+v3)
  # Tests
  return(list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,p10=p10,
              h=h,sig=sig,gref2=gref2,kkd=kkd,si=si))
}

