knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(micropoint)

# require(devtools)
# install_github("ilyamaclean/micropoint")

head(climdata)

# climdata_adj <- weatherhgt_adjust(climdata, zin = 2, zout = 27,
#                                   lat = 49.96807, long = -5.215668)

# vegp <-  createvegp(vegtype = 'C3')

# soilp <- createsoilc(soiltype = "Clay loam", nlayers = 15, totalDepth = 3)

# Cafromyear(1980)
# Cafromyear(2030)

# # ------------------------- Forest example ----------------------------------- #
# n <- 100 # Number of layers
# # Create vegetation parameters from habitat type:
#   # BET-Te – Broadleaf Evergreen Tree (Temperate)
# vegp <-  createvegp(vegtype = "BET.Te")
# paii <- PAIgeometry(vegp$pai, 7, 70, n = n) # pai of each layer
# Lfrac <- seq(0.5, 0.85, length.out = length(paii)) # Living vegetation fraction
# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Select hottest hour
# hr <- which.max(climdata$temp) # hottest hour (4094)
# # return and plot profile
# mout <- return_profile(climdata, hr, vegp, soilc, paii, Lfrac, lat = 49.96807,
#                        long = -5.215668, varn = "temp")

# climdata_adj <- weatherhgt_adjust(climdata, zin = 2, zout = 35,
#                                   lat = 49.96807, long = -5.215668)
# mout <- return_profile(climdata_adj, hr = 4, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, zref = 35, varn = "temp")

# mout <- return_profile(climdata_adj, hr = 13, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, zref = 35,
#                        varn = "relhum")
# mout <- return_profile(climdata_adj, hr = 13, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, zref = 35,
#                        varn = "radiation")

# # ------------------------- Grassland example --------------------------------
# n <- 20
# vegp <-  createvegp(vegtype = "C3") # C3 grassland
# vegp$h <- 1 # adjust vegetation height to 1m
# paii <- PAIgrass(vegp$h, n, vegp$pai)
# Lfrac <- seq(0.9, 1, length.out = length(paii))
# par(mfrow =c (1, 2))
# mout <- return_profile(climdata, hr = 4094, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, varn = "temp")
# mout <- return_profile(climdata, hr = 4094, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, varn = "relhum")

# # ------------------------- Forest example ----------------------------------- #
# # NB assumes climdata_adj already created
# n <- 20 # Number of layers
# # Create vegetation parameters for Broadleaf Evergreen Tree (Temperate)
# vegp <-  createvegp(vegtype = "BET.Te")
# paii <- PAIgeometry(vegp$pai, 7, 70, n = n) # pai of each layer
# Lfrac <- seq(0.5, 0.85, length.out = length(paii)) # Living vegetation fraction
# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Run model for 25 cm above ground and compare to input forcing data
# mout <- RunMicro(climdata_adj, reqhgt = 0.25, vegp, soilc, paii, Lfrac,
#                  lat = 49.96807, long = -5.215668, zref = 35)
# # Plot results for temperature (grey = macroclimate, red = microclimate)
# tme <- as.POSIXct(climdata_adj$obs_time, tz = "UTC")
# par(mfrow =c (1, 1))
# plot(climdata_adj$temp ~ tme, type = "l", col = rgb(0, 0, 0, 0.5),
#      xlab = "Month", ylab = "Temperature", ylim = c(-2, 25))
# par(new = TRUE)
# plot(mout$tair ~ tme, type = "l", col = rgb(1, 0, 0, 0.5),
#      xlab = "", ylab = "", ylim = c(-2, 25))

# # ------------------------- Grassland example ---------------------------------- #
# n <- 20 # Number of layers
# # Create vegetation parameters for C3 grassland
# vegp <-  createvegp(vegtype = "C3") # C3 grassland
# vegp$h <- 1 # adjust vegetation height to 1 m
# vegp$pai <- 0.5 # reduce pai to 0.5
# paii <- PAIgrass(vegp$h, n, vegp$pai) # generate foliage profile
# Lfrac <- seq(0.9, 1, length.out = length(paii))
# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Run model for 5 cm above ground and compare to input forcing data
# mout <- RunMicro(climdata, reqhgt = 0.05, vegp, soilc, paii, Lfrac,
#                  lat = 49.96807, long = -5.215668)
# # Plot results for temperature
# tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
# plot(climdata$temp ~ tme, type = "l", col = rgb(0, 0, 0, 0.5),
#      xlab = "Month", ylab = "Temperature", ylim = c(-3, 32))
# par(new = TRUE)
# plot(mout$tair ~ tme, type = "l", col = rgb(1, 0, 0, 0.5),
#      xlab = "", ylab = "", ylim = c(-3, 32))

# # ------------------------- Forest example ----------------------------------- #
# # NB assumes climdata_adj already created
# n <- 20 # Number of layers
# # Create vegetation parameters for Broadleaf Evergreen Tree (Temperate)
# vegp <-  createvegp(vegtype = "BET.Te")
# paii <- PAIgeometry(vegp$pai, 7, 70, n = n) # pai of each layer
# Lfrac <- seq(0.5, 0.85, length.out = length(paii)) # Living vegetation fraction
# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Run model for all heights
# mout <- RunModelFull(climdata_adj, soilc, vegp, paii, Lfrac,
#                      lat = 49.96807, long = -5.215668, zref = 35)

# require(terra)
# # Heights of every column, from the base of the soil to the height of the weather data
# zall <- c(-rev(mout$zsoil), mout$z, mout$za)
# temp <- cbind(mout$tsoil[, ncol(mout$tsoil):1], mout$tair, mout$tair_above)
# # Interpolate each hour to every cm
# zout <- seq(min(zall), max(zall), by = 0.01)
# temp <- apply(temp, 1, function(x) approx(zall, x, zout)$y)
# # Convert to raster, with rows running from top to bottom
# temp <- rast(temp[nrow(temp):1, ])
# ext(temp) <- c(0, nrow(mout$tair), min(zall), max(zall))
# plot(temp, asp = NA)

# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Select hottest hour
# hr <- which.max(climdata$temp) # hottest hour (4094)
# # return and plot profile
# mout <- return_profile(climdata, hr, vegp = NA, soilc, paii = NA, Lfrac = NA,
#                        lat = 49.96807, long = -5.215668, varn = "temp")

# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Run model for 2 cm above ground
# mout <- RunMicro(climdata, reqhgt = 0.02, vegp = NA, soilc, paii = NA, Lfrac = NA,
#                  lat = 49.96807, long = -5.215668)
# # plot temperature time series
# tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
# plot(climdata$temp ~ tme, type = "l", col = rgb(0, 0, 0, 0.5),
#      xlab = "Month", ylab = "Temperature", ylim = c(-3, 36))
# par(new = TRUE)
# plot(mout$tair ~ tme, type = "l", col = rgb(1, 0, 0, 0.5),
#      xlab = "", ylab = "", ylim = c(-3, 36))

# # Run model for a year to initialize water profile
# ThetaIni <- InitailizeWater(climdata, vegp = NA, soilc,
#                             lat = 49.96807, long = -5.215668)
# # Run model to generate output
# mout <- RunMicro(climdata, reqhgt = -0.001, vegp = NA, soilc, paii = NA,
#                  Lfrac = NA, lat = 49.96807, long = -5.215668,
#                  ThetaIni = ThetaIni)
# plot(mout$soilwater ~ tme, type = "l", col = "blue",
#      xlab = "Month", ylab = "Soil water")

# # Run model for all heights
# mout <- RunModelFull(climdata, soilc, vegp = NA, paii = NA, Lfrac = NA,
#                      lat = 49.96807, long = -5.215668)

# ectop <- create_ectop(height = 4, width = 3, len = 8, refl = 0.15, adir = 0, position = "randomdir", type = "reptile")

# # ------------------------- Forest example ----------------------------------- #
# # Perform weather height adjustment
# climdata_adj <- weatherhgt_adjust(climdata, zin = 2, zout = 35,
#                                   lat = 49.96807, long = -5.215668)
# # Create vegetation parameters from Broadleaf Evergreen Tree (Temperate)
# n = 100
# vegp <-  createvegp(vegtype = "BET.Te")
# paii <- PAIgeometry(vegp$pai, 7, 70, n = n) # pai of each layer
# Lfrac <- seq(0.5, 0.85, length.out = length(paii)) # Living vegetation fraction
# # Create soil parameters
# soilc <- createsoilc(soiltype = "Clay loam")
# # Select hottest hour
# hr <- which.max(climdata$temp) # hottest hour (4094)
# # return and plot profile
# mout <- return_profile(climdata_adj, hr = hr, vegp, soilc, paii, Lfrac,
#                        lat = 49.96807, long = -5.215668, zref = 35, varn = "temp")
# # run ectotherm model using ectotherm variables previously created
# ecto_out <- profile_ecto(mout, climdata_adj, hr, vegp, soilc, ectop,
#                          lat = 49.96807, long = -5.215668, zref = 35,
#                          maxiter = 100, tolerance = 1e-2)
# # plot profile
# plot(z ~ Tbody, data = ecto_out, type = "l", lwd = 2,
#      col = "red", xlab = "Temperature", ylab = "Height (m)")
# abline(h = 0, lty = 2)
# abline(h = vegp$h, lty = 2)

# # Assume climdata_adj, ectop, vegp and soilc have already bveen created
# n = 20
# paii <- PAIgeometry(vegp$pai, 7, 70, n = n) # pai of each layer
# Lfrac <- seq(0.5, 0.85, length.out = length(paii)) # Living vegetation fraction
# # return time series of microclimate 1m above ground
# mout <- RunMicro(climdata_adj, reqhgt = 1, vegp, soilc, paii, Lfrac,
#                  lat = 49.96807, long = -5.215668, zref = 35)
# # run ectotherm modelwith default maxiter and tolerence
# ecto_out <- timeseries_ecto(mout, climdata_adj, reqhgt = 1, vegp,
#                             soilc, ectop, lat = 49.96807,
#                             long = -5.215668)
# # Plot results
# tme <- as.POSIXct(ecto_out$obs_time)
# plot(Tbody ~ tme, data = ecto_out, type = "l", col = "red",
#      xlab = "Month", ylab = "Temperature")

# # Run model for all heights and hours assuming inputs already created
# mout <- RunModelFull(climdata_adj, soilc, vegp, paii, Lfrac,
#                      lat = 49.96807, long = -5.215668, zref = 35)
# # Run ectotherm model with default maxiter and tolerence
# ecto_out <- full_ecto(mout, climdata_adj, vegp, soilc, ectop,
#                       lat = 49.96807, long = -5.215668)

# require(terra)
# # Heights of every column, from the base of the soil to the height of the weather data
# zall <- c(-rev(mout$zsoil), mout$z, mout$za)
# Tbody <- cbind(ecto_out$Tbodyb[, ncol(ecto_out$Tbodyb):1], ecto_out$Tbodya, ecto_out$Tbody_above)
# # Interpolate each hour to every cm
# zout <- seq(min(zall), max(zall), by = 0.01)
# Tbody <- apply(Tbody, 1, function(x) approx(zall, x, zout)$y)
# # Convert to raster, with rows running from top to bottom
# Tbody <- rast(Tbody[nrow(Tbody):1, ])
# ext(Tbody) <- c(0, nrow(mout$tair), min(zall), max(zall))
# plot(Tbody, asp = NA)

# # Create ectotherm inputs
# ectop <- create_ectop(height = 4, width = 3, len = 8, refl = 0.15, adir = 0,
#                       position = "randomdir", type = "reptile")
# # Calculate hieght above ground of bottm of organism half way a long length
# zBtm <-  ectop$len * 0.5 * sin(ectop$atilt * pi / 180) # in cm
# # Calculate centre and convert to m
# zCtr <- (zBtm + 0.5 * ectop$height) * 0.01
# # Assume sandy soil (more realistic for an ectotherm)
# soilc <- createsoilc(soiltype = "Loamy sand")
# # Run microclimate model with inbilt climate dataset
# mout <- RunMicro(climdata, reqhgt = zCtr, vegp = NA, soilc, paii = NA, Lfrac = NA,
#                  lat = 49.96807, long = -5.215668, zref = 2)
#  # Run ecotherm model time-series
# ecto_out <- timeseries_ecto(mout, climdata, zCtr, vegp = NA, soilc, ectop,
#                             lat = 49.96807, long = -5.215668)
# # Plot time-series
# tme <- as.POSIXct(ecto_out$obs_time)
# plot(Tbody ~ tme, data = ecto_out, type = "l", col = "red",
#      xlab = "Month", ylab = "Temperature")
