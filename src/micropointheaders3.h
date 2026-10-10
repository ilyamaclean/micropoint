// Data structures of the model.
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
        double alpha; // photosynthesis quantum efficiency (mol CO2 / mol PAR)
        double Kxmx; // xylem maximum hydraulic conductance (mol / m^2 / s)
        double hv; // Huber value (m^2/m^2)
        double fd; // reaction of Vcmax converted to leaf dark respiration (unitless)
        double psi50; // Water potential when plant loses 50% of krc_max (MPa) 
        double apsi; // xylem hydraulic conductance parameter. Calculated if set to less than zero (unitless)
        double rpmin = -1.0; // minimum hydraulic resistance of the plant; calculated on first use (negative until then)
        double len;  // Leaf length (m)
        double wid; // Leaf width (m)
        double vegem; // Vegetation emissivity (unitless)
        double mwft; // max water film thickness (mm)
        double pTAW; // Fraction of total available water plant can deplete before stress (unitless)
        double root50; // Depth above which half of all roots lie (m)
        double root95; // Depth above which 95% of all roots lie (m)
        std::vector<double> paii; // one-sided plant area index of each canopy layer
        std::vector<double> pia; // one-sided plant area index above the bottom of each canopy layer
        std::vector<double> Lfrac; // Fraction of paii that is living plant material
    };
    struct soilpstruct {
        double slope; // slope (radians)
        double aspect; // aspect (radians)
        double gref; // Ground surface reflectance (shortwave radiation, 0-1)
        double grefPAR; // Ground surface reflectance (PAR, 0-1)
        double groundem; // Ground surface emissivity (unitless)
        int nLayers; // number of layers
        bool FreeDrain; // whether deep layer is free draining or not
        double z0; // roughness length of the soil surface for momentum (m)
        std::vector<double> Vq; // Volumetric quartz fraction (m^3/m^3)
        std::vector<double> Vm; // Volumetric mineral fraction (m^3/m^3)
        std::vector<double> Vo; // Volumetric organic fraction (m^3/m^3)
        std::vector<double> Mc; // Mass fraction of clay (kg/kg)
        std::vector<double> psie; // Campbell air entry potential (J/kg)
        std::vector<double> b; // Shape parameter of the water retention curve
        std::vector<double> thetaR; // Residual volumetric soil water fraction
        std::vector<double> thetaS; // Volumetric soil water fraction at saturation
        std::vector<double> Ksat; // Saturated hydraulic conductivity (kg s / m^3)
        std::vector<double> psi_min; // lowest water potential allowed: oven dry (J/kg)
    };
    struct climstruct {
        double tref; // temperature (deg C)
        double relhum; // relhum (percentage)
        double pk; // atmospheric pressure (kPa)
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
        std::vector<double> uz; // wind speed in each canopy layer
        double LL; // Monin-Obukhov length (m)
        double uf; // friction velocity
        double uh; // wind speed at top of canopy
        double zm; // roughness length
        double psi_m;
        double psi_h;
        double phi_h;
        // Friction velocity per unit wind speed.
        double ufratio;
        double windscale; // actual wind over effective wind (which includes the convective addition)
        double fade; // weight given to canopy-scale turbulence (0 to 1)
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
        std::vector<double> vol; // thickness of soil each node stores heat for (m)
        std::vector<double> wc; // volumetric water fraction of layer
        std::vector<double> Te; // temperature of layer
        std::vector<double> oldTe; // temperature of layer at the start of the time step
        double Gflux; // Ground heat flux (W/m^2)
        int iters; // Number of iterations to convergence
    };
    struct climforwaterstruct {
        double Rabs;// ground absorbed radiation
        double Tair; // reference air temperature
        double relhum; // reference relative humidity
        double pk; // reference pressure
        double rHa; // resistance from the ground to the reference height
        double precip; // precipitation
        double Et; // transpiration (mm)
    };
    struct soilwatermod {
        int n; // number of layers
        std::vector<double> z; // node depths
        std::vector<double> dz; // layer thickness
        std::vector<double> vol; // thickness of soil each node stores water for (m)
        std::vector<double> psiw; // water potential
        std::vector<double> k; // hydraulic conductivity
        std::vector<double> vapor; // vapour concentration (current time step)
        std::vector<double> oldvapor; // vapour concentration (previous time step)
        std::vector<double> theta; // volumetric water fraction (current time step)
        std::vector<double> oldtheta; // volumetric water fraction (previous time step)
        std::vector<double> Tc; // temperature (deg C) of soil
        std::vector<double> oldTc; // temperature (deg C) of soil at the start of the time step
        std::vector<double> rootfrac; // root fraction in each soil layer
    };
    struct soilwaterout {
        soilwatermod swo; // state of the soil water model
        bool success; // whether model converged
        int iterations; // how many iterations model run for
        double Evapmmhr; // surface evaporation
        double surplus = 0.0; // water a saturated surface could not take, which runs off (mm)
        double drainage = 0.0; // water leaving through the base of the soil (mm), negative where it enters
    };
    // State of the damped (Aitken) updates
    struct WAitkenStateScalar {
        double r_prev = 0.0;
        double omega = 0.2; // starting weight given to the new value
        bool have_prev = false;
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
        std::vector<double> tleaf; // mean temperature of foliage, leaves and wood (deg C)
        std::vector<double> rh; // relative humidity
        std::vector<double> uz; // wind speed
        std::vector<double> rLB; // leaf boundary-layer resistance (s/m)
        std::vector<double> swaterdepth; // Leaf surface water depth
        std::vector<double> leafwet; // fraction of the leaf surface that is wet
        std::vector<double> Hz; // sensible heat from foliage to air, per unit plant area (W/m^2)
        std::vector<double> Lz; // latent heat from foliage to air, per unit ground area (W/m^2)
        std::vector<double> gs;  // Stomatal conductance
        std::vector<double> dSTdT; // change in each layer's sensible heat source per degree of its own air (W m-2 K-1)
        std::vector<double> dSLde; // change in each layer's latent heat source per kPa of its own air's vapour pressure (W m-2 kPa-1)
        std::vector<double> tfar, tnear; // far-field and near-field parts of each node's air temperature due to foliage (K)
        std::vector<double> efar, enear; // the same for vapour pressure (kPa)
        double precipground; // Precipitation reaching the ground (mm)
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
        double fade; // weight given to canopy-scale turbulence (0 to 1)
        double tcanopy; // air temperature at canopy top (deg C)
        double rhcanopy; // relative humidity at canopy top (%)
        double ucanopy; // effective wind speed at canopy top (m/s), which includes the convective addition
        double windscale; // actual wind over effective wind
        double LL; // Monin-Obukhov length (m)
        int iters; // passes the hour took
        int witers; // iterations of the soil water solve
        double runoff = 0.0; // water a saturated surface could not take (mm)
        double drainage = 0.0; // water leaving through the base of the soil (mm)
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
        double LL; // Monin-Obukhov length (m)
        int iters; // passes the hour took
        double error;
        double runoff = 0.0; // water a saturated surface could not take (mm)
        double drainage = 0.0; // water leaving through the base of the soil (mm)
        int witers = 0; // iterations of the soil water solve
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