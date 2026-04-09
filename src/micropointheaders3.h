// micropointheaders2.h
// Radiation model
#include <Rcpp.h>
using namespace Rcpp;

namespace newmodel {
    struct vegpstruct {
        // Constant within canopy
        double hgt; // vegetation height (m)
        double pai; // Total one-sided plant area index (m^2 / m^2)
        double x; // Campbell foliage angle coefficient (unitless)
        double lref; // Leaf reflectance (shortwave radiation, 0-1)
        double ltra; // Leaf transmittance (shortwave radiation, 0-1)
        double lrefp; // Leaf reflectance (PAR, 0-1)
        double ltrap; // Leaf transmittance (PAR, 0-1)
        double Vcmax25; // micromol / m ^ 2 / s ^ 1
        double Tup; // high temperature photosynthesis range (deg C)
        double Tlw; // low temperature photosynthesis range (deg C)
        double Dcrit; // Onset of strong stomatal limitation (kPa)
        double alpha; // photosynthesis quantum efficiency (mol CO2 / mol PAR)
        double Kxmx; // xylem maximum hydraulic conductance (mol / m^2 / s)
        double hv; // Huber value (m^2/m^2)
        double f0; // ratio of ci/(ci-Gma) under low D from Jacobs (1994) equation (unitless) 
        double fd; // reaction of Vcmax converted to leaf dark respiration (unitless)
        double psi50; // Water potential when plant loses 50% of krc_max (MPa) 
        double apsi; // xylem hydraulic conductance parameter. Calculated if set to less than zero (unitless)
        double len;  // Leaf length (m)
        double wid; // Leaf width (m)
        double vegem; // Vegetation emissivity (unitless)
        double mwft; // max water film thickness (mm)
        double pTAW; // Fraction of total available water plant can deplete before stress (unitless)
        double rootskew; // Skew towards top of root layer (unitless)
        std::vector<double> paii; // one-sided plant area index at each canopy node
        std::vector<double> pia; // one-sided plant area index above each canopy node
        std::vector<double> Lfrac; // Fraction of paii that is living plant material
    };
    struct soilpstruct {
        double slope; // slope (radians)
        double aspect; // aspect (radians)
        double gref; // Ground surface refletcance (shortwave radiation, 0-1)
        double grefPAR; // Ground PAR refletcance (shortwave radiation, 0-1)
        double groundem; // Ground surface emissivity (unitless)
        int nLayers; // number of layers
        bool deepSaturated; // whether deep layer is saturated or not
        std::vector<double> Vq; // Volumetric quartz fraction (m^3/m^3
        std::vector<double> Vm; // Volumetric mineral fraction (m^3/m^3)
        std::vector<double> Vo; // Volumetric organic fraction (m^3/m^3)
        std::vector<double> Mc; // Mass fraction of clay (kg/kg
        std::vector<double> psie; // Campbell air entry potential (J/kg)
        std::vector<double> b; // Water retention curve shapet parameter
        std::vector<double> thetaR; // Residual soil water fraction at wilting point 
        std::vector<double> thetaS; // Volumetric soil water fraction at saturation
        std::vector<double> Ksat; // Saturated hydrualic conductivity (kg s / m^3 - same as cm/s)
        std::vector<double> n; // Campbell n
        std::vector<double> psi_min; // Minimum psiw determined from thetaR
    };
    struct climstruct {
        double tref; // temperature (deg C)
        double relhum; // relhum (percentage)
        double pk; // atmospheric pressure (pk)
        double Rsw; // shortwave down (W/m^2)
        double Rdif;  // diffuse down (W/m^2)
        double Rlw; // longwave down (W/m^2)
        double uref; // wind speed (m/s)
        double winddir; // wind direction
        double precip; // precipitation (mm)
    };
    struct obsstruct {
        int year;
        int month;
        int day;
        double hour;
    };
    struct envstruct {
        double PARabs; // PAR absorbed by leaf W / m ^ 2
        double tair; // air temperature(deg C)
        double tleaf; // leaf temperature(deg C)
        double rh; // relative humidity(percentage)
        double pk; // atmospheric pressure(kPa)
        double psi_r; // mean water potential in root zone(MPa)
        double Ca; // CO2 ppm
        double precip; // precipitation (mm)
    };
    struct solmodel {
        double zenr;
        double azir;
    };
    struct kstruct {
        double k;
        double kd;
    };
    struct tsvegstruct {
        double om;
        double a;
        double del;
        double J;
        double gma;
        double h;
        double u1;
        double u2;
        double S1;
        double D1;
        double D2;
    };
    struct tsdifstruct {
        double p1;
        double p2;
        double p3;
        double p4;
    };
    struct tsdirstruct {
        double sig;
        double p5;
        double p6;
        double p7;
        double p8;
        double p9;
        double p10;
    };
    struct radmodel {
        std::vector<double> RswLsun; // Shortwave radiation absorbed by sunlit leaves
        std::vector<double> RswLshade; // Shortwave radiation absorbed by shaded leaves
        std::vector<double> RswLav; // Average shortwave radiation absorbed (used for woody vegetation)
        std::vector<double> RPARsun; // PAR absorbed by sunlit leaves
        std::vector<double> RPARshade; // PAR absorbed by shaded leaves
        std::vector<double> sunfrac; // Fraction of sunlit leaves
        std::vector<double> RPARLabs; // PAR absorbed by leaf
        std::vector<double> Rdirdown; // Downward direct stream
        std::vector<double> Rdifdown; // Downward diffuse stream
        std::vector<double> Rswup; // Upward shortwave stream (entirely diffuse)
        double RswGabs; // Shortwave radiation absorbed by ground surface
        double RswCabs; // Shortwave radiation absorbed by vegetated surface
    };
    struct LWweights {
        Rcpp::NumericMatrix wgts;
        std::vector<double> trh;
        std::vector<double> trg;
        std::vector<double> wgtg;
        double trsky;
    };
    struct radmodel2 {
        std::vector<double> RlwLabs; // Longwave radiation absorbed by leaf
        std::vector<double> Rlwdown; // Downward longwave stream
        std::vector<double> Rlwup; // Upward longwave stream
        double RlwGabs; // Longwave radiation absorbed by ground surface
        double RlwCabs; // Longwave radiation absorbed by vegetated surface
    };
    struct windmodel {
        std::vector<double> uz; // wind speed below canopy 
        double LL; // Monin–Obukhov length (m)
        double uf; // friction velocity
        double uh; // wind speed at top of canopy
        double a2; // Langrangian time-scale parameter
        double zm; // roughness length
        double psi_m;
        double psi_h;
        double phi_h;
    };
    struct rainmodel {
        std::vector<double> kd; // extinction coefficient
        std::vector<double> tr; // transmission of rain
    };
    struct Thomas {
        std::vector<double> bb;
        std::vector<double> cc;
        std::vector<double> dd;
        std::vector<double> x;
    };
    struct soilmod {
        int n; // number of layers
        std::vector<double> z; // node depths
        std::vector<double> dz; // layer thickness
        std::vector<double> zCenter; // centre between nodes
        std::vector<double> vol; // volume of layer
        std::vector<double> wc; // volumetric water fraction of layer
        std::vector<double> Te; // temperature of layer
        std::vector<double> oldTe; // old temperature of layer
        double Gflux; // Ground heat flux (W/m^2)
        int iters; // Number of iterations to convergence
    };
    struct climforwaterstruct {
        double Rabs;// ground absorbed radiation
        double Tair; // reference air temperature
        double relhum; // reference relative humidity
        double pk; // reference pressure
        double rHa; // reference resistance form gorund to zref
        double precip; // precipitation
        double Et; // Evapotranspiration
    };
    struct soilwatermod {
        int n; // number of layers
        std::vector<double> z; // node depths
        std::vector<double> dz; // layer thickness
        std::vector<double> zCenter; // centre between nodes
        std::vector<double> vol; // volume of layer (m3)
        std::vector<double> psiw; // water potential
        std::vector<double> k; // hydrualic conductivity
        std::vector<double> vapor; // vapour concentration (current time step)
        std::vector<double> oldvapor; // vapour concentration (previous time step)
        std::vector<double> theta; // volumetric water fraction (current time step)
        std::vector<double> oldtheta; // volumetric water fraction (previous time step)
        std::vector<double> Tc; // temperature (deg C) of soil
        std::vector<double> oldTc; // temperature (deg C) of soil
        std::vector<double> rootfrac; // root fraction in each soil layer
    };
    struct soilwaterout {
        soilwatermod swo; // main bit of soil water model
        bool success; // whether model converged
        int iterations; // how many iterations model run for
        double Evapmmhr; // surface evaporation
    };
    // Strucxtures used for convergence function
    struct WAitkenState {
        std::vector<double> r_prev;
        double omega = 0.3;
        bool have_prev = false;
    };
    struct WAitkenStateScalar {
        double r_prev = 0.0;
        double omega = 0.2;   // default: quite high backweight
        bool have_prev = false;
    };
    struct Hstruct {
        double Tsurf;
        double Htot;
    };
    struct cantop {
        double Th;
        double eh;
    };
    struct onestep {
        std::vector<double> Rdirdown;
        std::vector<double> Rdifdown;
        std::vector<double> Rswup;
        std::vector<double> Rlwdown;
        std::vector<double> Rlwup;
        std::vector<double> tair; // air temperature
        std::vector<double> tleaf; // avaerage temperature of foliage elements used for Langrangian model and returned by default
        std::vector<double> rh; // relative humidity
        std::vector<double> uz; // wind speed
        std::vector<double> rLB; // Leaf boundary layer conductance
        std::vector<double> swaterdepth; // Leaf surface water depth
        std::vector<double> Hz; // Sensible heat flux exchange with surrounding air
        std::vector<double> Lz;  // Sensible heat flux exchange with surrounding air
        std::vector<double> gs;  // Stomatal conductance
        double precipground; // Precipitation eaching the ground
        soilmod soilheatvars; // soil heat model input / output
        soilwatermod soilwatervars; // soil water model input / output;
        double H; // Total sensible heat (including from ground)
        double L; // Total latent heat (including from ground)
        double Et; // Total transpiration
        double Ev; // Bare soil evaporation
        double theta; // Volumetric soil moisture of surface layer
        double psim; // psi_m
        double psih; // psi_h
        double phih; // phi_h
        double LL; // // Monin–Obukhov length (m)
        int iters; // number of iterations for which model run
        int witers; // number of iterations for which water model run
        double error;
    };
    struct onestepbare {
        std::vector<double> tair; // air temperature
        std::vector<double> rh; // relative humidity
        std::vector<double> uz; // wind speed
        soilmod soilheatvars; // soil heat model input / output
        soilwatermod soilwatervars; // soil water model input / output;
        double Rb0;
        double Rswup;
        double Rlwup;
        double H; // Total sensible heat 
        double L; // Total latent heat
        double Ev; // Bare soil evaporation
        double theta; // Volumetric soil moisture of surface layer
        double psim; // psi_m
        double psih; // psi_h
        double LL; // // Monin–Obukhov length (m)
        int iters; // number of iterations for which model run
        double error;
    };
    struct bigsmallleafstruct {
        double Hsmall;
        double Lsmall;
        double Hbig;
        double Lbig;
        double zh;
        double gS;
        double T0;
    };
    struct silstruct {
        double A; // area (m^2)
        double V; // volume (m^3)
        double silA; // silhouette area (m^2)
    };
    struct albf {
        double albb; // direct
        double albd; // diffuse
        double wgtg; // weighting direct
        double wgti; // weighting diffuse
    };
    struct Aitken1DState {
        double r_prev = 0.0;
        double omega = 1.0;   // start with no damping
        bool have_prev = false;
    };
    struct bigleafone {
        soilmod soilheat;
        soilwaterout soilwater;
        double tcanopy;
        double uf;
        double LL;
        double Et;
        int iters;
    };
}