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
