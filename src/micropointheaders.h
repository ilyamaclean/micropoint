// micropointheaders.h
// Radiation model
#include <Rcpp.h>
struct radmodel {
    std::vector<double> ground; // ground absorbed radiation (W/m^2)
    std::vector<double> canopy; // canopy absorbed radiation (W/m^2)
    std::vector<double> albedo; // albedo
};
struct radmodel2 {
    std::vector<double> Rswabs;
    std::vector<double> PAR;
    std::vector<double> Rdirdown;
    std::vector<double> Rdifdown;
    std::vector<double> Rswup;
};
struct radmodel3 {
    std::vector<double> Rlwdown;
    std::vector<double> Rlwup;
};
struct Gmodel {
    std::vector<double> G;
    std::vector<double> Gmin;
    std::vector<double> Gmax;
};
struct LWweights {
    std::vector<double> trg;
    std::vector<double> trh;
    Rcpp::NumericMatrix wgt;
};
struct canHL {
    std::vector<double> H;
    std::vector<double> L;
    std::vector<double> tleaf;
    std::vector<double> uz;
};
struct Lang {
    std::vector<double> tleaf;
    std::vector<double> tair;
    std::vector<double> ea;
    std::vector<double> uz;
    std::vector<double> Rlwdown;
    std::vector<double> Rlwup;
    double mxdif;
};