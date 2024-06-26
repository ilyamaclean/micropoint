% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Runfunctions.R
\name{runchecksfun}
\alias{runchecksfun}
\title{runs checks on microclimate input data}
\usage{
runchecksfun(climdata, vegp, groundp, lat, long)
}
\arguments{
\item{climdata}{a data.frame of hourly weather as in the example dataset
`climdata` (see details)}

\item{vegp}{an object of class vegparams or forestparams formatted as for the inbuilt example
datasets `vegparams` or `forestparams` (see details).}

\item{groundp}{an object of class groundparams formatted as for the inbuilt
example dataset `groundparams`.}

\item{latitude}{(decimal degrees) - used for radiation checks}

\item{longitude}{(decimal degrees) - used for radiation checks}
}
\description{
The function `runchecksfun` runs checks on the microclimate model
inputs
}
\details{
this function runs checks on provided data testing whether provided data lie within plausable ranges
and returning errors or warnings accordingly. The original data are returned as a list with
some values capped if checks fail.
}
