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
#' @title Generates plant area index profile for grass habitat
#' @description Generates a vector of length `n` of plausible plant area index values
#' @param hgt grass height (m) - of tallest grass spike (~1.5x drop disk measured height)
#' @param n Number of plant area index values to generate. Default: 1000 (see details)
#' @param PAI Optionally total plant area index of canopy.
#' @param fraccover Optionally grass fractional cover
#' @param taper parameter controlling how tapered the plant area index
#' is towards the top. A value of one is typical of grass. Lower values give
#' a more even spread. Cannot be negative or zero.
#' @return a vector of length `n` of plant area index values the sum of which equals `PAI`
#' @details if `PAI` is not supplied it is estimated from `fraccover`, which must
#' be between 0.001 and 0.9999. If neither `PAI` or `fraccover` is supplied, `PAI`
#' is estimated using a simple allometric relationship with sward height,
#' derived from field measurements from a meadow in northern Spain.
#' @examples
#' paii <- PAIgrass(hgt = 0.25, n = 1000)
#' folden <- paii * 1000
#' z <- (c(1:1000) / 1000) * 0.25
#' plot(z ~ folden, type = "l", main = sum(paii)) # plots foliage density
PAIgrass<-function(hgt, n = 1000, PAI = NA, fraccover = NA, taper = 1) {
  if (taper <= 0) stop ("taper must be greater than 0")
  # Generate profile
  z<-c(0:50)/50
  const<-0.80304-0.09719*log(hgt*1000)
  const <- const * (1/taper)
  ifd<-  13.79233*z+const
  folden<-1/ifd
  folden<-spline(folden,n=n)$y
  # Calculate PAI
  if (is.na(PAI)) {
    if (is.na(fraccover) == FALSE) {
      if (fraccover<0.0001) stop("fraccover must be greater than 0.0001")
      if (fraccover>0.9999) stop("fraccover must be less than 0.9999")
      PAI<- -log(1-fraccover)/0.9867526
    }  else {
      lpai<- -8.9388+1.3898 *log(hgt*1000)
      PAI<-exp(lpai)
    }
  }
  # rescale foliage density profile
  paii<-(folden/sum(folden))*PAI
  return(paii)
}
#' @title Generates vegetation parameters for grass
#' @description Generates an object of class
#' `vegparams` and a vector of of length `n` of plausible plant
#' area index values for grassy vegetation.
#' @param hgt grass height (m) - of tallest grass spike (~1.5x drop disk measured height)
#' @param n Number of plant area index values to generate. Default: 20 (see details)
#' @param taper parameter controlling how tapered the plant area index
#' is towards the top. A value of one is typical of grass. Lower values give
#' a more even spread. Cannot be negative or zero.
#' @return a vector of length `n` of plant area index values the sum of which equals `PAI`
#' @details a list comrpising the following:
#' \describe{
#'  \item{vegp}{An object of vegetation parameters for running the microclimate model}
#'  \item{paii}{a vector of length `n` of plant area index values}
#' }
#' @examples
#' vp <- vegpforgrass(0.25, 20)
#' vegp <- vp$vegp
#' paii <- vp$paii
#' mout <- runpointmodel(climdata, reqhgt = 0.12, vegp, paii,  groundparams, lat = 49.96807, long= -5.215668)
#' tme <- as.POSIXct(climdata$obs_time)
#' plot(mout$tair ~ tme, type="l", xlab = "Month", ylab = "Air temperature",
#'      ylim = c(-5, 40), col = "red")
vegpforgrass <- function(hgt, n = 20, taper = 1) {
  lpai<- -8.9388+1.3898 *log(hgt*1000)
  PAI<-exp(lpai)
  vegp<-list(h=hgt,pai=PAI,x=0.161808,
             clump=0.3,lref=0.25,ltra=0.22,
             leafd=0.08*hgt,em=0.97,
             gsmax=0.38,q50=100)
  class(vegp)<-"vegparams"
  paii<-PAIgrass(hgt, n, taper = taper)
  return(list(vegp=vegp,paii=paii))
}



