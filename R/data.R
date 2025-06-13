#' A data frame of hourly weather
#'
#' A data frame of hourly weather in 2017 at Caerthillean Cove, Lizard, Cornwall (49.96807N, 5.215668W)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swdown}{Total downward shortwave radiation (W / m^2)}
#'  \item{difrad}{Total downward diffuse radiation (W / m^2)}
#'  \item{lwdown}{Total downward longwave radiation (W / m^2)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{precip}{precipitation (mm)}
#' }
"climdata"
#' A data frame of hourly weather for a snowy environment
#'
#' A data frame of hourly weather (2nd Oct 2017 to 1st Apr 2019) for Sodankylä in Finland (67.367N, 26.629E)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swdown}{Total downward shortwave radiation (W / m^2)}
#'  \item{difrad}{Total downward diffuse radiation (W / m^2)}
#'  \item{lwdown}{Total downward shortwave radiation (W / m^2)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{precip}{precipitation (mm)}
#' }
"snowclimdata"
#'
#' A dataset of vegetation parameters for running the Small Leaf model.
#'
#' An object of class forestparams
#' @format A list of  the following objects:
#' \describe{
#'   \item{h}{vegetation height (m)}
#'   \item{pai}{plant area index}
#'   \item{x}{ratio of vertical to horizontal projections of leaf foliage}
#'   \item{clump}{canopy clumping factor quantified as fraction of direct radiation passing through larger gaps in the canopy when the sun is at its zenith}
#'   \item{lref}{leaf reflectance}
#'   \item{ltra}{leaf transmittance}
#'   \item{leafd}{leaf diameter (m)}
#'   \item{em}{emissvity of vegetation}
#'   \item{gsmax}{maximum stomatal conductances (mol / m^2 / s)}
#'   \item{q50}{value of PAR when stomatal conductances is at 50 percent of its maximum value (micro mol/m^2/s)}
#'   \item{skew}{degree of skew towards top of canopy in foliage}
#'   \item{spread}{degree of spread in foliage}
#' }
"forestparams"
#'
#' A dataset of ground parameters.
#'
#' An object of class groundparams
#'
#' @format A list of  the following objects:
#' \describe{
#'   \item{gref}{ground reflectance}
#'   \item{slope}{slope of ground surface (deg)}
#'   \item{aspect}{aspect of ground surface (deg from north)}
#'   \item{em}{emissivity of ground surface}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{Vm}{Volumetric mineral fraction of soil}
#'   \item{Vq}{Volumetric quartz fraction of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{b}{Shape parameter for Campbell soil moisture model}
#'   \item{Psie}{Matric potential (J / m^3)}
#'   \item{Smax}{Volumetric water content at saturation}
#'   \item{Smin}{Residual water content}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#' }
"groundparams"
#'
#' A dataset of soilparameters parameters.
#'
#' An object of class groundparams
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter (dimensionless, > 1)}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#'   \item{Vq}{Volumetric quartz fraction of soil}
#'   \item{Vm}{Volumetric mineral fraction of soil}
#'   \item{Vo}{Volumetric organic fraction of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{b}{Shape parameter for Campbell model (dimensionless, > 1)}
#'   \item{psi_e}{Matric potential (J / m^3)}
#'   \item{mult}{soil moisture model radiation coefficient}
#'   \item{rmu}{soil moisture model rainfall coefficient}
#'   \item{a}{soil moisture model deeper layer multiplier coefficient}
#'   \item{pwr}{soil moisture model deeper layer power coefficient}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparams"
#'
#' A dataset of vegetation parameters for running Big Leaf model.
#'
#' A data.frame of parameters for running the soil models
#'
#' @format A data.frame with the following elements:
#' \describe{
#'   \item{h}{vegetation height (m)}
#'   \item{pai}{plant area index}
#'   \item{x}{ratio of vertical to horizontal projections of leaf foliage}
#'   \item{clump}{canopy clumping factor quantified as fraction of direct radiation passing through larger gaps in the canopy when the sun is at its zenith}
#'   \item{lref}{leaf reflectance}
#'   \item{ltra}{leaf transmittance}
#'   \item{leafd}{leaf diameter (m)}
#'   \item{em}{emissvity of vegetation}
#'   \item{gsmax}{maximum stomatal conductances (mol / m^2 / s)}
#'   \item{q50}{value of PAR when stomatal conductances is at 50 percent of its maximum value (micro mol/m^2/s)}
#' }
"vegparams"
