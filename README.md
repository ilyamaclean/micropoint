# micropoint

`micropoint` is an R package for mechanistic modelling of microclimate at a point location: above, within
and below a vegetation canopy and below ground. From standard hourly weather data and vegetation
and soil properties, it predicts temperature, humidity, wind speed, radiation and soil moisture
through the canopy and soil. It also includes a model of ectotherm body temperature that takes
the microclimate as its input.

Like most mechanistic microclimate models, it is based on surface energy balance. The system is
implemented as a multilayer canopy–soil model, allowing vertical variation in microclimate
conditions to be represented. For computational efficiency the core model is written in C++, with
an R interface that allows users to configure and run simulations within R.

## About this branch

- The `main` branch holds the current model, which improves the representation of soil
  temperature, soil water dynamics and plant transpiration.
- The `legacy` branch holds the earlier model. It is kept for reference and is no longer developed.
- **Snow is not represented.**
- This version is still under development. Functionality, parameter definitions and outputs may
  change without notice.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ilyamaclean/micropoint")
```

The package contains C++ code, so a compiler is needed: Rtools on Windows, Xcode command line
tools on macOS. It has few dependencies, but needs a recent version of `Rcpp`. The vignettes come
already built, so nothing extra is needed to read them.

To install the earlier model instead:

```r
remotes::install_github("ilyamaclean/micropoint", ref = "legacy")
```

## Quick start

```r
library(micropoint)
# Vegetation, foliage profile and soil for a temperate broadleaf evergreen forest
vegp <- createvegp(vegtype = "BET.Te")
paii <- PAIgeometry(vegp$pai, skew = 7, spread = 70, n = 20)
Lfrac <- rep(0.8, 20)
soilc <- createsoilc(soiltype = "Clay loam")
# Run the model for 1 m above ground, using the inbuilt weather dataset
mout <- RunMicro(climdata, reqhgt = 1, vegp, soilc, paii, Lfrac,
                 lat = 49.96807, long = -5.215668)
# Air temperature for every hour of 2017
tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
plot(mout$tair ~ tme, type = "l", xlab = "Month", ylab = "Temperature")
```

## Documentation

Two vignettes are included:

```r
vignette("new_model_tutorial")  # how to prepare inputs and run the model
vignette("modelequations")      # the equations the model solves
```

Each function also has its own help page, for example `?RunMicro`.

## Main functions

| Function | Purpose |
|---|---|
| `createvegp()`, `createsoilc()` | Create vegetation and soil parameters |
| `PAIgeometry()`, `PAIgrass()` | Generate the vertical profile of plant area |
| `return_profile()` | Vertical microclimate profile for a specified hour |
| `RunMicro()` | Time series of microclimate at a specified height |
| `RunModelFull()` | Microclimate for all heights and hours |
| `weatherhgt_adjust()` | Adjusts weather data to a height above the canopy |
| `create_ectop()` | Creates ectotherm parameters |
| `profile_ecto()`, `timeseries_ecto()`, `full_ecto()` | Ectotherm body temperature from the outputs of the three run functions |

## Problems and questions

Please report problems at <https://github.com/ilyamaclean/micropoint/issues>.
