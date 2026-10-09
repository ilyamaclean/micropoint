#include <Rcpp.h>
#include "micropointheaders3.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstddef>

// Optional timers and counters for profiling. Compiled in only with -DMP_PROFILE.
#ifdef MP_PROFILE
#include <chrono>
static double g_profSeconds[20] = { 0 };
static double g_profCount[24] = { 0 };
static std::chrono::steady_clock::time_point g_profLast = std::chrono::steady_clock::now();
static inline void profLap(int part) {
    const auto now = std::chrono::steady_clock::now();
    if (part >= 0) g_profSeconds[part] += std::chrono::duration<double>(now - g_profLast).count();
    g_profLast = now;
}
#define MP_LAP(part) profLap(part)
#define MP_COUNT(part) (g_profSeconds[part] += 1.0)
#define MP_TALLY(kind) (g_profCount[kind] += 1.0)
#define MP_TALLYN(kind, n) (g_profCount[kind] += static_cast<double>(n))
#else
#define MP_LAP(part)
#define MP_COUNT(part)
#define MP_TALLY(kind)
#define MP_TALLYN(kind, n)
#endif
using namespace Rcpp;
using namespace newmodel;

constexpr double pi = 3.14159265358979323846;
constexpr double sb = 5.67e-8;      // Stefan-Boltzmann constant (W/m2/K4)
constexpr double ka = 0.41;         // von Karman constant
constexpr double Mw = 0.018015;     // molar mass of water (kg/mol)
constexpr double RgasC = 8.314;     // gas constant (J/mol/K)
constexpr double g = 9.80665;       // gravity (m/s2)
constexpr double torad = 3.14159265358979323846 / 180.0;   // degrees to radians
// Step in soil water potential used to estimate how evaporation changes with
// surface wetness: a fraction of the potential, but never less than the minimum.
constexpr double PSI_SLOPE_FRAC = 1e-4;
constexpr double PSI_SLOPE_MIN = 1.0;
// Water potential of oven-dry soil (J/kg): the driest the soil is allowed to get.
constexpr double PSI_OVEN_DRY = -1.0e6;
// Water potentials (J/kg) between which roots can take up water: field capacity and wilting point.
constexpr double PSI_FIELD_CAPACITY = -33.0;
constexpr double PSI_WILTING_POINT = -1500.0;
inline double aitken1d(double oldv, double newv, Aitken1DState& st);
// Default stopping tolerance (K): an hour is solved when one more pass changes
// air temperature by less than this.
constexpr double OUTER_TOLERANCE = 0.1;
// Solver settings, fixed for a run by the accuracy level:
//   minPass            fewest passes an hour takes
//   leafSteps          refinements of each leaf temperature
//   soilHeatLoosening  factor by which the soil heat tolerance is relaxed
//   topTolerance       tolerance on canopy-top air temperature (K)
//   groundRow          0: ground held while canopy air is solved; 1 or 2: ground
//                      adjusted with the air (1 also re-solves soil heat)
//   nearLocal          near field from a layer's own foliage and its two neighbours only
//   leafClosed         leaf response to air temperature from a formula, not a second solve
struct solversettings { int minPass = 3; int leafSteps = 4; double soilHeatLoosening = 1.0; double topTolerance = 1e-8; int groundRow = 0;
    bool nearLocal = false; bool leafClosed = false; };
// Fewest iterations allowed to the wind, soil water and canopy-top solves.
constexpr int INNER_ITERATIONS_LEAST = 50;

static solversettings g_solver;
// [[Rcpp::export]]
void SolverSettingsCpp(int minpass = 3, int leafsteps = 4, double soilheatloosening = 1.0, double toptolerance = 1e-8, int groundrow = 0,
    bool nearlocal = false, bool leafclosed = false)
{
    g_solver.groundRow = groundrow;
    g_solver.nearLocal = nearlocal; g_solver.leafClosed = leafclosed;
    g_solver.minPass = std::max(1, minpass);
    g_solver.leafSteps = std::max(1, leafsteps);
    g_solver.soilHeatLoosening = std::max(1.0, soilheatloosening);
    g_solver.topTolerance = std::max(1e-8, toptolerance);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************** Solar model ***************************************** //
// Sun position and the angle of the beam on the ground.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Julian day from a calendar date.
static int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1.0)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// Solar time (hours) from local time and longitude.
static double soltimeCpp(int jd, double lt, double lond)
{
    double m = 6.24004077 + 0.01720197 * (static_cast<double>(jd) - 2451545.0);
    double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2.0 * m + 3.5932);
    double st = lt + (4.0 * lond + eot) / 60.0;
    return st;
}
// Solar zenith and azimuth (radians).
static solmodel solpositionCpp2(double latr, double lonr, int year, int month, int day, double lt)
{
    solmodel solpos{};
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lonr * 180.0 / pi);
    // Zenith
    double tt = 0.261799 * (st - 12.0);
    double dec = (pi * 23.5 / 180.0) * std::cos(2.0 * pi * ((jd - 159.5) / 365.25));
    double coh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
    solpos.zenr = std::acos(coh);
    // Azimuth
    double sh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
    double hh = std::atan(sh / std::sqrt(1.0 - sh * sh));
    double sazi = std::cos(dec) * std::sin(tt) / std::cos(hh);
    double cazi = (std::sin(latr) * std::cos(dec) * std::cos(tt) - std::cos(latr) * std::sin(dec)) /
        std::sqrt(std::pow(std::cos(dec) * std::sin(tt), 2) + std::pow(std::sin(latr) *
            std::cos(dec) * std::cos(tt) - std::cos(latr) * std::sin(dec), 2));
    double sqt = 1.0 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    solpos.azir = pi + std::atan(sazi / std::sqrt(sqt));
    if (cazi < 0.0) {
        solpos.azir = (sazi < 0.0) ? (pi - solpos.azir) : (3.0 * pi - solpos.azir);
    }
    return solpos;
}
// Solar index: cosine of the angle between the beam and the ground surface.
// Zero when the sun is below the horizon or behind the slope.
static double solarindexCpp2(double sloper, double aspectr, double zenr, double azir)
{
    double si;
    if (zenr > (pi/2.0)) {
        si = 0;
    }
    else {
        if (sloper == 0.0) {
            si = std::cos(zenr);
        }
        else {
            si = std::cos(zenr) * std::cos(sloper) + std::sin(zenr) * std::sin(sloper) * std::cos(azir - aspectr);
        }
    }
    if (si < 0.0) si = 0.0;
    return si;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Radiation model ****************************************************** //
// Shortwave radiation through the canopy, by a two-stream model.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Extinction coefficients for the direct beam: k over flat ground, kd for canopy
// layers lying parallel to sloping ground. x is the leaf-angle parameter.
static kstruct cankCpp(double zenr, double x, double si) {
    double k;
    if (zenr > (pi / 2.0)) zenr = pi / 2.0;
    if (si < 0.0) si = 0.0;
    // Flat ground
    if (x == 1.0) {
        k = 1.0 / (2.0 * std::cos(zenr));
    }
    else if (std::isinf(x)) {
        k = 1.0;
    }
    else if (x == 0.0) {
        k = (2.0 / pi) * std::tan(zenr);
    }
    else {
        k = std::sqrt(x * x + (std::tan(zenr) * std::tan(zenr))) / (x + 1.774 * std::pow((x + 1.182), -0.733));
    }
    if (k > 6000.0) k = 6000.0;
    // Sloping ground. Capped where the beam grazes the slope or is behind it.
    double kd = (si > 0.0) ? k * std::cos(zenr) / si : 6000.0;
    if (kd > 6000.0) kd = 6000.0;
    kstruct kparams{};
    kparams.k = k;
    kparams.kd = kd;
    return kparams;
}
// Optical properties of the canopy for the two-stream model (Dickinson 1983;
// Sellers 1985), from leaf reflectance and transmittance, leaf angle, plant area
// and ground reflectance.
static tsvegstruct twostreamvegCpp(double pai, double x, double lref, double ltra, double gref)
{
    tsvegstruct params{};
    params.om = lref + ltra;
    params.a = 1.0 - params.om;
    params.del = lref - ltra;
    params.J = 1.0 / 3.0;
    if (x != 1.0) {
        double mla = 9.65 * std::pow((3.0 + x), -1.65);
        if (mla > pi / 2.0) mla = pi / 2.0;
        params.J = std::cos(mla) * std::cos(mla);
    }
    params.gma = 0.5 * (params.om + params.J * params.del);
    params.h = std::sqrt(params.a * params.a + 2.0 * params.a * params.gma);
    params.S1 = std::exp(-params.h * pai);
    params.u1 = params.a + params.gma * (1.0 - 1.0 / gref);
    params.u2 = params.a + params.gma * (1.0 - gref);
    params.D1 = (params.a + params.gma + params.h) * (params.u1 - params.h) *
        1.0 / params.S1 - (params.a + params.gma - params.h) * (params.u1 + params.h) * params.S1;
    params.D2 = (params.u2 + params.h) * 1.0 / params.S1 - (params.u2 - params.h) * params.S1;
    return params;
}
// Two-stream coefficients for diffuse radiation.
static tsdifstruct twostreamdifCpp(tsvegstruct& tsvegp)
{
    tsdifstruct params{};
    params.p1 = (tsvegp.gma / (tsvegp.D1 * tsvegp.S1)) * (tsvegp.u1 - tsvegp.h);
    params.p2 = (-tsvegp.gma * tsvegp.S1 / tsvegp.D1) * (tsvegp.u1 + tsvegp.h);
    params.p3 = (1.0 / (tsvegp.D2 * tsvegp.S1)) * (tsvegp.u2 + tsvegp.h);
    params.p4 = (-tsvegp.S1 / tsvegp.D2) * (tsvegp.u2 - tsvegp.h);
    return params;
}
// Two-stream coefficients for the direct beam.
static tsdirstruct twostreamdirCpp(double pai, double kd, double gref, tsvegstruct tsvegp)
{
    tsdirstruct params{};
    double sig = kd * kd + tsvegp.gma * tsvegp.gma - std::pow((tsvegp.a + tsvegp.gma), 2.0);
    double ss = 0.5 * (tsvegp.om + tsvegp.J * tsvegp.del / kd) * kd;
    double sstr = tsvegp.om * kd - ss;
    double S2 = std::exp(-kd * pai);
    params.p5 = -ss * (tsvegp.a + tsvegp.gma - kd) - tsvegp.gma * sstr;
    double v1 = ss - (params.p5 * (tsvegp.a + tsvegp.gma + kd)) / sig;
    double v2 = ss - tsvegp.gma - (params.p5 / sig) * (tsvegp.u1 + kd);
    params.p6 = (1.0 / tsvegp.D1) * ((v1 / tsvegp.S1) * (tsvegp.u1 - tsvegp.h) - (tsvegp.a + tsvegp.gma - tsvegp.h) * S2 * v2);
    params.p7 = (-1.0 / tsvegp.D1) * ((v1 * tsvegp.S1) * (tsvegp.u1 + tsvegp.h) - (tsvegp.a + tsvegp.gma + tsvegp.h) * S2 * v2);
    params.sig = -sig;
    // Downward diffuse radiation produced by scattering of the direct beam.
    params.p8 = sstr * (tsvegp.a + tsvegp.gma + kd) + tsvegp.gma * ss;
    double v3 = (sstr + tsvegp.gma * gref - (params.p8 / params.sig) * (tsvegp.u2 - kd)) * S2;
    params.p9 = (-1.0 / tsvegp.D2) * ((params.p8 / (params.sig * tsvegp.S1)) * (tsvegp.u2 + tsvegp.h) + v3);
    params.p10 = (1.0 / tsvegp.D2) * (((params.p8 * tsvegp.S1) / params.sig) * (tsvegp.u2 - tsvegp.h) + v3);
    return params;
}
// Shortwave absorbed by the ground, the whole canopy and each layer, for total
// shortwave and for PAR. A layer absorbs the difference between the net downward
// flux at its top and at its bottom. Leaf values are per unit one-sided plant
// area; sunlit leaves receive the direct beam as well as diffuse radiation.
static radmodel shortwavemodelCpp(const std::vector<double>& pia, const std::vector<double>& paii, double pai, double gref, double grefPAR, double lref, double lrefp, double Rswdown, double Rdif, double si, const solmodel& solp, const kstruct& kp, const tsvegstruct& tspveg, const tsvegstruct& tspvegPAR, const tsdifstruct& tspdif, const tsdifstruct& tspdifPAR, const tsdirstruct& tspdir, const tsdirstruct& tspdirPAR) {
    const int n = static_cast<int>(pia.size());
    radmodel out;
    out.RswLsun.assign(n, 0.0); out.RswLshade.assign(n, 0.0); out.RswLav.assign(n, 0.0);
    out.RPARsun.assign(n, 0.0); out.RPARshade.assign(n, 0.0); out.sunfrac.assign(n, 0.0);
    out.Rdirdown.assign(n, 0.0); out.Rdifdown.assign(n, 0.0); out.Rswup.assign(n, 0.0);
    out.RswGabs = 0.0; out.RswCabs = 0.0;
    if (Rswdown <= 0.0) return out;

    // Direct beam normal to the sun (Rbeam) and onto the ground surface (Rb).
    const double cosz = std::cos(solp.zenr);
    const double Rbeam = (cosz > 0.0) ? (Rswdown - Rdif) / cosz : 0.0;
    const double Rb = Rbeam * si;

    // Upper limit on the fraction of the beam scattered: the larger of leaf and
    // ground reflectance.
    double amx = gref; if (amx < lref) amx = lref;
    double amxp = grefPAR; if (amxp < lrefp) amxp = lrefp;

    const double kd = kp.kd, h = tspveg.h, hp = tspvegPAR.h;
    const double sil = kp.k * cosz;
    const double a = tspveg.a, ap = tspvegPAR.a;

    // Downward and upward diffuse radiation at cumulative plant area p.
    auto streams = [&](double p, const tsvegstruct& tv, const tsdifstruct& td, const tsdirstruct& tb,
        double am, double& difd, double& up) {
        const double exp_mh = std::exp(-tv.h * p), exp_ph = std::exp(tv.h * p);
        double Rddu = td.p1 * exp_mh + td.p2 * exp_ph;
        if (Rddu > 1.0) Rddu = 1.0;
        if (Rddu < 0.0) Rddu = 0.0;
        double Rddd = td.p3 * exp_mh + td.p4 * exp_ph;
        if (Rddd > 1.0) Rddd = 1.0;
        if (Rddd < 0.0) Rddd = 0.0;
        double Rdbu = 0.0, Rdbd = 0.0;
        if (Rb > 0.0) {
            const double exp_kd = std::exp(-kd * p);
            Rdbu = (tb.p5 / -tb.sig) * exp_kd + tb.p6 * exp_mh + tb.p7 * exp_ph;
            if (Rdbu > am) Rdbu = am;
            if (Rdbu < 0.0) Rdbu = 0.0;
            Rdbd = (tb.p8 / tb.sig) * exp_kd + tb.p9 * exp_mh + tb.p10 * exp_ph;
            if (Rdbd > am) Rdbd = am;
            if (Rdbd < 0.0) Rdbd = 0.0;
        }
        difd = Rddd * Rdif + Rdbd * Rb;
        up = Rddu * Rdif + Rdbu * Rb;
    };
    // Net downward shortwave through the canopy plane at cumulative plant area p.
    auto netdown = [&](double p, bool par) {
        double difd, up;
        if (par) streams(p, tspvegPAR, tspdifPAR, tspdirPAR, amxp, difd, up);
        else streams(p, tspveg, tspdif, tspdir, amx, difd, up);
        return Rb * std::exp(-kd * p) + difd - up;
    };

    // Ground absorption: beam and diffuse radiation reaching the ground, less reflection.
    {
        double difd_g, up_g;
        streams(pai, tspveg, tspdif, tspdir, amx, difd_g, up_g);
        out.RswGabs = (1.0 - gref) * (difd_g + Rb * std::exp(-kd * pai));
    }
    // Whole-canopy absorption: what enters at the top and does not reach the ground.
    out.RswCabs = netdown(0.0, false) - netdown(pai, false);

    // Each layer: pia is cumulative plant area to its bottom. The radiation reported
    // for a layer is that entering its top.
    for (int i = 0; i < n; ++i) {
        const double pbot = pia[i];
        const double ptop = std::max(pbot - paii[i], 0.0);
        double difd, up;
        streams(ptop, tspveg, tspdif, tspdir, amx, difd, up);
        out.Rdirdown[i] = Rbeam * std::exp(-kd * ptop);
        out.Rdifdown[i] = difd;
        out.Rswup[i] = up;

        // Beam intercepted by the layer, and the layer's mean sunlit fraction.
        const double ebt = std::exp(-kd * ptop), ebb = std::exp(-kd * pbot);
        const double dE = std::max(ebt - ebb, 0.0);
        const double Bint = Rb * dE;
        double sf = ebt;
        if (paii[i] > 0.0 && kd * paii[i] > 1e-8) sf = dE / (kd * paii[i]);
        if (sf > 1.0) sf = 1.0;
        out.sunfrac[i] = sf;

        // Absorption per unit plant area: diffuse radiation shared by all leaves (Dsw,
        // Dpar), with the beam added for sunlit leaves.
        double Dsw, Dpar;
        if (paii[i] > 1e-9) {
            const double Asw = netdown(ptop, false) - netdown(pbot, false);
            const double Apar = netdown(ptop, true) - netdown(pbot, true);
            Dsw = std::max(Asw - a * Bint, 0.0) / paii[i];
            Dpar = std::max(Apar - ap * Bint, 0.0) / paii[i];
        }
        else {
            double difdp, upp;
            streams(ptop, tspvegPAR, tspdifPAR, tspdirPAR, amxp, difdp, upp);
            Dsw = a * (difd + up);
            Dpar = ap * (difdp + upp);
        }
        out.RswLsun[i] = a * sil * Rbeam + Dsw;
        out.RswLshade[i] = Dsw;
        out.RswLav[i] = sf * out.RswLsun[i] + (1.0 - sf) * out.RswLshade[i];
        out.RPARsun[i] = ap * sil * Rbeam + Dpar;
        out.RPARshade[i] = Dpar;
    }
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************** Longwave radiation model ***************************** //
// Longwave exchange between sky, ground and canopy layers.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Weights giving how much of each other layer, the ground and the sky each layer
// sees. They depend on plant area only, so are calculated once.
static LWweights lwradweights(const std::vector<double>& paii) {
    int n = static_cast<int>(paii.size());
    double pait = 0.0;
    for (int i = 0; i < n; ++i) pait += paii[i];
    // Plant area above and below the middle of each layer.
    std::vector<double> paia(n);
    std::vector<double> paib(n);
    double below = 0.0;
    for (int i = 0; i < n; ++i) {
        paib[i] = below + 0.5 * paii[i];
        paia[i] = pait - paib[i];
        below += paii[i];
    }
    // Weight between two layers: plant area of the other layer times transmission
    // through the foliage between them.
    NumericMatrix wgts(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double pij = std::abs(paia[i] - paia[j]);
                double tr = std::exp(-pij);
                wgts(i, j) = tr * paii[j];
            }
            else wgts(i, j) = 0.0;
        }
    }
    // Transmission from each layer to the sky (trh) and to the ground (trg).
    std::vector<double> trg(n);
    std::vector<double> trh(n);
    for (int i = 0; i < n; ++i) {
        trh[i] = std::exp(-paia[i]);
        trg[i] = std::exp(-paib[i]);
    }
    std::vector<double> wgt(n);
    for (int i = 0; i < n; ++i) {
        wgt[i] = 0.0;
        for (int j = 0; j < n; ++j) wgt[i] += wgts(i, j);
    }
    // Rescale so each layer's weights sum to 2: a leaf has two sides.
    std::vector<double> wsum(n);
    for (int i = 0; i < n; ++i) wsum[i] = wgt[i] + trh[i] + trg[i];
    std::vector<double> mu(n);
    for (int i = 0; i < n; ++i) {
        double mu1 = wsum[i] / 2.0;
        mu[i] = 1.0 + (2.0 - 2.0 * mu1) / wgt[i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) wgts(i, j) = wgts(i, j) * mu[i];
    }
    // The ground's weights, rescaled to sum to 1.
    std::vector<double> wgtg(n);
    for (int i = 0; i < n; ++i) {
        wgtg[i] = trg[i] * paii[i];
    }
    double gwsum = 0.0;
    for (int i = 0; i < n; ++i) gwsum += wgtg[i];
    double trsky = std::exp(-pait);
    double mug1 = (gwsum + trsky) / 1.0;
    double mug = 1.0 + (1.0 - mug1) / gwsum;
    for (int i = 0; i < n; ++i) wgtg[i] = wgtg[i] * mug;
    LWweights out;
    out.wgts = wgts; // weight matrix for each canopy element
    out.trh = trh; // transmission from sky
    out.trg = trg; // transmission from ground
    out.wgtg = wgtg; // weights for each canopy element for ground
    out.trsky = trsky;
    return out;
}
// Absolute temperature to the fourth power.
static double radem(double tc)
{
    return std::pow(tc + 273.15, 4.0);
}
// Longwave absorbed by each canopy layer and by the ground, given their temperatures.
static radmodel2 longwavemodelCpp(const LWweights& wgts, double lwdown, double tground, double groundem, double vegem,
    const std::vector<double>& tleaf)
{
    int n = static_cast<int>(wgts.trh.size());
    std::vector<double> RlwLabs(n);
    std::vector<double> Rlwdown(n);
    std::vector<double> Rlwup(n);
    std::vector<double> leafLW(n);
    // Longwave leaving each layer and the ground: black-body emission, which
    // includes the longwave a surface reflects as well as what it emits.
    for (int i = 0; i < n; ++i) leafLW[i] = sb * radem(tleaf[i]);
    double groundLW = sb * radem(tground);
    for (int i = 0; i < n; ++i) {
        double lwsky = lwdown * wgts.trh[i];
        double lwgro = groundLW * wgts.trg[i];
        double lwcand = 0.0;
        for (int j = i; j < n; ++j) {
            lwcand += wgts.wgts(i, j) * leafLW[j];
        }
        double lwcanu = 0.0;
        for (int j = 0; j <= i; ++j) {
            lwcanu += wgts.wgts(i, j) * leafLW[j];
        }
        Rlwdown[i] = lwsky + lwcand;
        Rlwup[i] = lwgro + lwcanu;
        // Absorbed by both sides of the leaves, per unit one-sided plant area.
        RlwLabs[i] = vegem * (Rlwdown[i] + Rlwup[i]);
    }
    double lwgrfromcan = 0.0;
    for (int i = 0; i < n; ++i) {
        lwgrfromcan += wgts.wgtg[i] * leafLW[i];
    }
    double lwgrfromsky = wgts.trsky * lwdown;
    double RlwGabs = groundem * (lwgrfromsky + lwgrfromcan);
    double em = wgts.trsky * groundem + (1.0 - wgts.trsky) * vegem;
    double RlwCabs = em * lwdown;
    radmodel2 out;
    out.Rlwdown = Rlwdown;
    out.Rlwup = Rlwup;
    out.RlwLabs = RlwLabs;
    out.RlwGabs = RlwGabs;
    out.RlwCabs = RlwCabs;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// *************************************************** Wind model ****************************************************** //
// Wind above and within the canopy. Friction velocity and atmospheric stability
// are solved together.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Zero-plane displacement height, from canopy height and plant area index.
// [[Rcpp::export]]
double zeroplanedisCpp2(double h, double pai)
{
    // Avoids division by zero when there is no foliage.
    if (pai < 1e-9) pai = 1e-9;
    double d = (1.0 - (1.0 - exp(-std::sqrt(7.5 * pai))) / std::sqrt(7.5 * pai)) * h;
    return d;
}
// Default roughness length of bare soil (m).
constexpr double Z0_GROUND = 0.004;
// Ratio of friction velocity to wind speed at canopy top (Raupach 1994, Eq. 7).
// It rises with plant area up to a maximum (UF_UH_MAX). With no foliage it is the
// ratio for bare soil at the same height.
constexpr double DRAG_ELEMENT = 0.3;
constexpr double UF_UH_MAX = 0.3;
// Bare soil
static double groundVelocityRatio(double h, double z0)
{
    const double hh = std::max(h, 1.001 * z0);
    return std::min(ka / std::log(hh / z0), UF_UH_MAX);
}
// Soil and foliage
static double frictionVelocityRatio(double h, double pai, double z0)
{
    const double b0 = groundVelocityRatio(h, z0);
    double be = std::sqrt(b0 * b0 + DRAG_ELEMENT * pai / 2.0);
    if (be > UF_UH_MAX) be = UF_UH_MAX;
    return be;
}
// Canopy closure: 0 for bare ground, 1 when more foliage adds no more drag.
static double canopyClosure(double h, double pai, double z0)
{
    const double b0 = groundVelocityRatio(h, z0);
    if (b0 >= UF_UH_MAX) return 1.0;
    const double f = (frictionVelocityRatio(h, pai, z0) - b0) / (UF_UH_MAX - b0);
    return std::min(std::max(f, 0.0), 1.0);
}
// Drag coefficient of a leaf.
constexpr double LEAF_DRAG = 0.25;
// Iteration limit, tolerance and relaxation for the mixing length of the canopy wind profile.
constexpr int WIND_MAX_ITER = 60;
constexpr double WIND_TOL = 1e-10;
constexpr double WIND_RELAX = 0.6;
// Strength of the roughness sublayer at canopy top (Raupach 1992):
// ln(cw) - 1 + 1/cw, where cw is 2 over a closed canopy and 1 over bare ground.
static double sublayerInfluence(double h, double pai, double z0, double fade)
{
    const double cw = 1.0 + canopyClosure(h, pai, z0) * fade;
    return std::log(cw) - 1.0 + 1.0 / cw;
}
// Smallest value allowed for the log term that gives wind at canopy top.
constexpr double LNTERM_MIN = 0.05;
// Roughness length for momentum (Raupach 1992, Eq. 27), from canopy geometry.
// [[Rcpp::export]]
double roughlengthCpp2(double h, double pai, double d, double z0 = 0.004, double fade = 1.0)
{
    double Be = frictionVelocityRatio(h, pai, z0);
    double zm = (h - d) * std::exp(-ka / Be + sublayerInfluence(h, pai, z0, fade));
    // No more than 0.9 of the height of canopy top above the displacement height.
    if (zm > (0.9 * (h - d))) zm = 0.9 * (h - d);
    // No less than the roughness length of the soil.
    if (zm < z0) zm = z0;
    return zm;
}
// Molar density of air (mol/m3).
static double phairCpp(double tc, double pk)
{
    double tk = tc + 273.15;
    double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
    return ph;
}
// Molar specific heat of air at constant pressure (J/mol/K).
static double cpairCpp(double tc)
{
    double cp = 2e-05 * tc * tc + 0.0002 * tc + 29.119;
    return cp;
}
// Integrated stability correction for momentum (ze is height over Obukhov length).
// Bounded in both stable and unstable air.
// [[Rcpp::export]]
double dpsimCpp2(double ze)
{
    double psim;
    // unstable
    if (ze < 0.0) {
        double x = std::pow((1.0 - 15.0 * ze), 0.25);
        psim = std::log(std::pow((1.0 + x) / 2.0, 2.0) * (1.0 + x * x) / 2.0) - 2.0 * std::atan(x) + pi / 2.0;
        if (psim > 3.0) psim = 3.0;
    }
    // stable
    else {
        const double zetaMaxM = 4.0 / 4.7;
        psim = -4.7 * (zetaMaxM * std::tanh(ze / zetaMaxM));
    }
    return psim;
}
// Integrated stability correction for heat. The unstable branch approaches its
// limit (PSIH_CAP) smoothly from PSIH_KNEE; the stable branch is bounded.
constexpr double PSIH_KNEE = 1.5;
constexpr double PSIH_CAP = 3.0;
constexpr double PSIH_STABLE_SLOPE = 4.7 / 0.74;
constexpr double PSIH_STABLE_ZMAX = 4.0 / PSIH_STABLE_SLOPE;
static double dpsihCpp2(double ze)
{
    double psih;
    // unstable
    if (ze < 0.0) {
        double y = std::sqrt(1.0 - 9.0 * ze);
        psih = 2.0 * std::log((1.0 + y) / 2.0);
        if (psih > PSIH_KNEE) {
            const double span = PSIH_CAP - PSIH_KNEE;
            psih = PSIH_KNEE + span * std::tanh((psih - PSIH_KNEE) / span);
        }
    }
    // stable
    else {
        psih = -PSIH_STABLE_SLOPE * (PSIH_STABLE_ZMAX * std::tanh(ze / PSIH_STABLE_ZMAX));
    }
    return psih;
}
// Dimensionless temperature gradient (phi_h) matching the correction above.
// Heat diffusivity in the surface layer is ka u* (z - d) / phi_h.
static double phihStar(double ze)
{
    if (ze < 0.0) {
        const double y = std::sqrt(1.0 - 9.0 * ze);
        const double phi = 1.0 / y;   // unmodified branch: Businger-Dyer form
        const double psiraw = 2.0 * std::log((1.0 + y) / 2.0);
        if (psiraw <= PSIH_KNEE) return phi;
        // Above the knee the gradient returns towards its neutral value.
        const double span = PSIH_CAP - PSIH_KNEE;
        const double c = std::cosh((psiraw - PSIH_KNEE) / span);
        return 1.0 + (phi - 1.0) / (c * c);
    }
    const double c = std::cosh(ze / PSIH_STABLE_ZMAX);
    return 1.0 + PSIH_STABLE_SLOPE * ze / (c * c);
}
// Obukhov length, in stable air, that gives a stated stability correction for
// momentum. The correction peaks at one length and falls either side: the length
// returned is at or above the peak, and a target above the peak returns the peak.
static double lStableFinalSearch(double psi_m, double zref, double d, double zm)
{
    if (psi_m <= 1e-12) return std::numeric_limits<double>::infinity();

    auto psiOfL = [&](double L) {
        return dpsimCpp2(zm / L) - dpsimCpp2((zref - d) / L);
    };

    // Find the length at which the correction peaks: coarse scan, then golden-section search.
    const int NSCAN = 60;
    double bestL = 1.0, bestVal = psiOfL(1.0);
    double logLo = std::log(1e-4), logHi = std::log(1e6);
    for (int i = 0; i <= NSCAN; ++i) {
        double L = std::exp(logLo + (logHi - logLo) * i / NSCAN);
        double v = psiOfL(L);
        if (v > bestVal) { bestVal = v; bestL = L; }
    }
    int bestIdx = static_cast<int>(std::round(NSCAN * (std::log(bestL) - logLo) / (logHi - logLo)));
    double loL = std::exp(logLo + (logHi - logLo) * std::max(0, bestIdx - 1) / NSCAN);
    double hiL = std::exp(logLo + (logHi - logLo) * std::min(NSCAN, bestIdx + 1) / NSCAN);

    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = loL, b = hiL;
    double c = b - gr * (b - a), e = a + gr * (b - a);
    for (int i = 0; i < 60; ++i) {
        if (psiOfL(c) > psiOfL(e)) { b = e; } else { a = c; }
        c = b - gr * (b - a); e = a + gr * (b - a);
        if (std::abs(b - a) < 1e-9 * b) break;
    }
    double Lpeak = 0.5 * (a + b);
    double psiPeak = psiOfL(Lpeak);

    if (psi_m >= psiPeak) return Lpeak;

    // Bisect above the peak for the length that gives the target.
    double lo = Lpeak, hi = std::max(Lpeak * 100.0, 1e6);
    while (psiOfL(hi) > psi_m && hi < 1e12) hi *= 10.0;
    for (int i = 0; i < 80; ++i) {
        double mid = 0.5 * (lo + hi);
        if (psiOfL(mid) > psi_m) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
// As above, reusing the last answer when the inputs are unchanged.
static double lStableFinalCpp(double psi_m, double zref, double d, double zm)
{
    static double lastPsi = -1.0, lastZref = -1.0, lastD = -1.0, lastZm = -1.0, lastL = 0.0;
    if (psi_m != lastPsi || zref != lastZref || d != lastD || zm != lastZm) {
        lastL = lStableFinalSearch(psi_m, zref, d, zm);
        lastPsi = psi_m; lastZref = zref; lastD = d; lastZm = zm;
    }
    return lastL;
}
// Limits the Obukhov length so that the stability correction for momentum is no
// more than a fraction (beta) of the log term ln((zref - d)/zm). This keeps
// resistances positive in very stable or very unstable air.
static double clipMOlength(double L, double zref, double d, double zm, double beta = 0.9)
{
    const double ln_z = std::log((zref - d) / zm);
    const double psim_min = -beta * ln_z;
    const double tol = 1e-4;
    double psim = dpsimCpp2(zm / L) - dpsimCpp2((zref - d) / L);
    if (L < 0.0 && psim < psim_min) {
        // Unstable: bisect for the length at which the correction equals the limit.
        double L_low = -500.0;
        double L_high = L;
        double L_mid;
        for (int i = 0; i < 30; ++i) {
            L_mid = 0.5 * (L_low + L_high);
            double psim_mid = dpsimCpp2(zm / L_mid) - dpsimCpp2((zref - d) / L_mid);
            if (psim_mid < psim_min)
                L_low = L_mid;
            else
                L_high = L_mid;
            if (std::abs(psim_mid - psim_min) < tol)
                break;
        }
        return L_high;  // corrected L that satisfies psim > min
    }
    if (L >= 0.0) {
        // Stable: the length itself is compared with its limit, because the correction
        // is not monotonic in the length.
        const double psim_max = beta * ln_z;
        const double L_bound = lStableFinalCpp(psim_max, zref, d, zm);
        if (L < L_bound) {
            return L_bound;
        }
    }
    return L;
}
// Wind in each canopy layer as a fraction of wind at canopy top. It depends on
// the canopy only, so is calculated once.
static std::vector<double> windprofileCpp(const vegpstruct& vegp, double z0) {
    int n = static_cast<int>(vegp.paii.size());
    const double dz = vegp.hgt / n;
    const double beta = frictionVelocityRatio(vegp.hgt, vegp.pai, z0);

    // Wind falls with the foliage passed on the way down: the cube root of foliage
    // density, summed from canopy top.
    std::vector<double> a(n), a13(n), cumA(n), ui(n);
    for (int i = 0; i < n; ++i) {
        a[i] = vegp.paii[i] / dz;
        a13[i] = std::cbrt(a[i]);
    }
    double running = 0.0;
    for (int i = n - 1; i >= 0; --i) {
        running += a13[i] * dz;
        cumA[i] = running;              // foliage encountered between here and canopy top
    }

    // The mixing length depends on the foliage density where momentum is absorbed,
    // which depends on the wind profile, so the two are iterated. For uniform foliage
    // the result is the exponential profile of Harman and Finnigan (2007).
    double aeff = vegp.pai / vegp.hgt;
    double lmix = 2.0 * beta * beta * beta / (LEAF_DRAG * aeff);
    for (int iter = 0; iter < WIND_MAX_ITER; ++iter) {
        const double kA = std::cbrt(LEAF_DRAG / (2.0 * lmix * lmix));
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; ++i) {
            ui[i] = std::exp(-kA * (cumA[i] - 0.5 * a13[i] * dz));
            double w = a[i] * ui[i] * ui[i];
            num += a[i] * w;
            den += w;
        }
        double anew = (den > 0.0) ? (num / den) : aeff;
        double lnew = 2.0 * beta * beta * beta / (LEAF_DRAG * anew);
        if (std::abs(lnew - lmix) < WIND_TOL * lmix) { lmix = lnew; break; }
        lmix += WIND_RELAX * (lnew - lmix);
    }

    // In the lowest tenth of the canopy the wind follows a logarithmic profile down
    // to the soil. Not applied if the lowest layer lies within the soil's roughness.
    int nground = n / 10;
    const double zLowest = 0.5 * dz;
    if (nground >= 1 && zLowest > z0) {
        const double zAnchor = (nground + 0.5) * dz;
        const double lnAnchor = std::log(zAnchor / z0);
        const double uAnchor = ui[nground];
        for (int i = 0; i < nground; ++i) {
            ui[i] = uAnchor * std::log((i + 0.5) * dz / z0) / lnAnchor;
        }
    }
    return ui;
}
// Vertical gust strength at the ground as a fraction of friction velocity. It
// decays with the drag area of the foliage above, reduced for sheltering
// (Massman and Weil 1999; Massman 1997).
constexpr double GUST_DECAY = 4.18;   // per unit of sheltered drag area
constexpr double SHELTER = 0.4;       // Massman's sheltering factor
static double floorGustRatio(const vegpstruct& vegp, double a1)
{
    const size_t n = vegp.paii.size();
    double zeta = 0.0;
    for (size_t i = 0; i < n; ++i) {
        // Drag area of the layer, reduced by sheltering where foliage is dense.
        zeta += LEAF_DRAG * vegp.paii[i] / (1.0 + SHELTER * static_cast<double>(n) * vegp.paii[i]);
    }
    return a1 * std::exp(-GUST_DECAY * zeta);
}
// Factor by which heat diffusivity at canopy top exceeds the surface-layer value:
// 3 over a closed canopy (Raupach et al. 1996), 1 over bare ground.
constexpr double CANOPY_TOP_ENHANCEMENT = 3.0;
static double canopyEnhancement(double h, double pai, double z0, double fade)
{
    return 1.0 + (CANOPY_TOP_ENHANCEMENT - 1.0) * canopyClosure(h, pai, z0) * fade;
}
// Lowest height at which the weather data may be measured: the height at which
// surface-layer diffusivity equals the diffusivity at canopy top.
// [[Rcpp::export]]
double referenceheightCpp(double h, double pai, double z0 = 0.004)
{
    if (h <= 0.0) return 0.0;
    const double d = zeroplanedisCpp2(h, pai);
    return d + canopyEnhancement(h, pai, z0, 1.0) * (h - d);
}
// Height above the displacement height used for stability at canopy top.
static double canopyTopStabilityHeight(double h, double d, double enhancement)
{
    return 0.5 * (1.0 + enhancement) * (h - d);
}
// Resistance to heat between the canopy treated as one big leaf (the spin-up
// model) and canopy top. Its roughness length for heat is a fifth of that for momentum.
static double bigLeafToCanopyTop(double h, double pai, double d, double uf, double LL, double z0)
{
    const double dt = h - d;
    double zm0 = dt * std::exp(-ka / frictionVelocityRatio(h, pai, z0));
    zm0 = std::max(std::min(zm0, 0.9 * dt), z0);
    const double zh0 = 0.2 * zm0;
    const double r = (std::log(dt / zh0) + dpsihCpp2(zh0 / LL) - dpsihCpp2(dt / LL)) / (ka * uf);
    return std::max(r, 0.0);
}
// Guard against still air; not physics. Friction velocity is never less than a
// wind of U_MIN_MIXING (m/s) would give over the same surface, the two combined
// in quadrature. Used for resistances, canopy mixing and stability.
constexpr double U_MIN_MIXING = 0.5;
static double guardedFrictionVelocity(double uf, double ufPerUnitWind)
{
    const double uf_calm = ufPerUnitWind * U_MIN_MIXING;
    return std::sqrt(uf * uf + uf_calm * uf_calm);
}
static double mixingFrictionVelocity(const windmodel& windvars)
{
    return guardedFrictionVelocity(windvars.uf, windvars.ufratio);
}
// Weight (0 to 1) given to canopy-scale turbulence, from the canopy's Reynolds
// number: wind at canopy top times canopy height over the viscosity of air. It is
// 0.5 at CANOPY_RE_TRANSITION, an order-of-magnitude value. It scales the
// canopy-top enhancement, the roughness sublayer and the near field.
constexpr double CANOPY_RE_TRANSITION = 1000.0;
constexpr double AIR_VISCOSITY = 1.5e-5;   // kinematic, m2/s
static double canopyTurbulence(double uDrive, double h, double pai, double d, double zref, double z0)
{
    const double zm1 = roughlengthCpp2(h, pai, d, z0, 1.0);
    double lnt = std::log((h - d) / zm1) + sublayerInfluence(h, pai, z0, 1.0);
    if (lnt < LNTERM_MIN) lnt = LNTERM_MIN;
    const double uh = uDrive * lnt / std::log((zref - d) / zm1);
    const double Re = uh * h / AIR_VISCOSITY;
    return Re * Re / (Re * Re + CANOPY_RE_TRANSITION * CANOPY_RE_TRANSITION);
}
// Friction velocity, Obukhov length, wind at canopy top and wind in each layer,
// from wind at the reference height and sensible heat flux (H). Friction velocity
// and Obukhov length depend on each other, so are iterated.
static windmodel windmodelCpp(const std::vector<double>& wc, double uref, double hgt, double pai,
    double zref, double H = 0.0, double tc = 15, double pk = 101.3, int maxiter = 100,
    double a1 = 1.25, double psi_h = 0.0, double psi_m = 0.0, double phi_h = 1.0,
    double zi = 1000.0, double beta = 1.0, double z0 = Z0_GROUND)
{
    if (zref < hgt) {
        Rcpp::stop("Height of wind speed measurement must be above canopy");
    }
    double Tk = tc + 273.15;
    double cp = cpairCpp(tc);
    double ph = phairCpp(tc, pk);
    double d = zeroplanedisCpp2(hgt, pai);
    double zm = roughlengthCpp2(hgt, pai, d, z0);
    double zh = 0.2 * zm;
    // With an upward heat flux, convection adds to the wind (Beljaars 1994).
    double Ueff = uref;
    if (H > 0.0) {
        double wstar = std::cbrt((g / Tk) * zi * (H / (ph * cp)));
        Ueff = std::sqrt(uref * uref + std::pow(beta * wstar, 2.0));
    }
    const double fade = canopyTurbulence(Ueff, hgt, pai, d, zref, z0);
    zm = roughlengthCpp2(hgt, pai, d, z0, fade);
    zh = 0.2 * zm;
    double uf = (ka * Ueff) / (std::log((zref - d) / zm) + psi_m);
    // Obukhov length: effectively infinite (neutral) when there is no heat flux.
    double LL = 1e99;
    // Iterate friction velocity and Obukhov length.
    double dif = 10.0;
    double olduf = -100.00;
    int iter = 1;
    // The friction velocity used for the Obukhov length is damped between iterations.
    Aitken1DState st_uf;
    double uf_iter = uf;
    if (H != 0.0) {
        MP_TALLY(11);
        while (dif > 0.00000001) {
            MP_TALLY(12);
            {
                const double ufL = guardedFrictionVelocity(uf_iter, ka / std::log((zref - d) / zm));
                LL = (ph * cp * std::pow(ufL, 3.0) * Tk) / (-ka * g * H);
            }
            LL = clipMOlength(LL, zref, d, zm);
            psi_m = dpsimCpp2(zm / LL) - dpsimCpp2((zref - d) / LL);
            psi_h = dpsihCpp2(zh / LL) - dpsihCpp2((zref - d) / LL);
            uf = (ka * Ueff) / (std::log((zref - d) / zm) + psi_m);
            dif = std::abs(olduf - uf);
            if (iter > maxiter) dif = 0.0;
            olduf = uf;
            uf_iter = aitken1d(uf_iter, uf, st_uf);
            iter += 1;
        }
    }
    // Friction velocity per unit wind speed.
    double ufratio = ka / (std::log((zref - d) / zm) + psi_m);
    // Stability correction for momentum between the roughness height and canopy top.
    if (H != 0.0) {
        double ln1 = std::log((zref - d) / zm);
        if (psi_m < -0.9 * ln1) psi_m = -0.9 * ln1;
        if (psi_m > 0.9 * ln1) psi_m = 0.9 * ln1;
        psi_m = dpsimCpp2(zm / LL) - dpsimCpp2((hgt - d) / LL);
    }
    // Wind at canopy top, including the effect of the roughness sublayer.
    double lnterm = std::log((hgt - d) / zm) + sublayerInfluence(hgt, pai, z0, fade) + psi_m;
    if (lnterm < LNTERM_MIN) lnterm = LNTERM_MIN;
    double uh = (uf / ka) * lnterm;
    int n = static_cast<int>(wc.size());
    std::vector<double> uz(n);
    for (int i = 0; i < n; ++i) uz[i] = wc[i] * uh;
    phi_h = phihStar(canopyTopStabilityHeight(hgt, d, canopyEnhancement(hgt, pai, z0, fade)) / LL);
    windmodel out;
    out.uz = uz;
    out.LL = LL;
    out.uf = uf;
    out.uh = uh;
    out.zm = zm;
    out.psi_m = psi_m;
    out.psi_h = psi_h;
    out.phi_h = phi_h;
    out.ufratio = ufratio;
    // Factor converting the effective wind (with convection) back to the actual wind, for output.
    out.windscale = (Ueff > 0.0) ? uref / Ueff : 1.0;
    out.fade = fade;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************** Plant model ****************************************************** //
// Leaf energy balance and stomatal conductance, and the heat and vapour each
// canopy layer gives to the air.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Saturated vapour pressure (kPa): over water above 0 degC, over ice below.
// [[Rcpp::export]]
double satvapCpp2(double tc)
{
    double es;
    if (tc > 0) {
        es = 0.61078 * exp(17.27 * tc / (tc + 237.3));
    }
    else {
        es = 0.61078 * exp(21.875 * tc / (tc + 265.5));
    }
    return es;
}
// Saturated vapour density (kg/m3) at temperature Tk (K).
static double satVapDensityCpp2(double Tk)
{
    double tc = Tk - 273.15;
    double es_Pa = 1000.0 * satvapCpp2(tc); // kPa -> Pa
    return Mw * es_Pa / (RgasC * Tk);
}
// Characteristic dimension of a leaf (m), from its length, width and leaf-angle parameter.
static double leafDimension(double len, double wid, double x)
{
    double y = 1.0 / x;
    double Gy = std::sqrt(y * y + (1.0 - y * y)) / (y + 1.774 * std::pow((y + 1.182), -0.733));
    return std::sqrt(len * wid) * Gy;
}
// Leaf boundary-layer resistance to heat (s/m): forced and free convection combined.
static double leafrHa(double tair, double dT, double uz, double d, double rHmax = 300.0)
{
    double Tk = tair + 273.15;
    double Kh = (1.6667e-10 * Tk * Tk + 2.9935e-8 * Tk - 1.7128e-6); // thermal diffusivity
    double v = 1.326 * std::pow(10.0, -5.0) * std::pow(Tk / 273.15, 1.5) * (393.55 / (Tk + 120.0)); // kinematic viscosity
    double Re = (uz * d) / v; // Reynolds number
    double Pr = v / Kh; // Prandtl number
    double Nuf; // Nusselt number, forced convection
    if (Re > 2e5) {
        Nuf = 0.37 * std::pow(Re, 0.6) * std::pow(Pr, 1.0 / 3.0) + 9.08;
    }
    else if (Re > 1000) {
        Nuf = 2.0 + (0.48 * std::pow(Re, 0.6) - 11.31) * std::pow(Pr, 1.0 / 3.0);
    }
    else {
        Nuf = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0 / 3.0);
    }
    double Gr = (g * std::pow(d, 3.0) * dT) / (Tk * v * v); // Grashof number
    double Nun = std::pow(0.825 + (0.387 * std::pow(Gr * Pr, 1.0 / 6.0)) /
        std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0), 2.0); // Nusselt number, free convection
    double Nu = std::pow(std::pow(Nuf, 3.0) + std::pow(Nun, 3.0), 1.0 / 3.0);
    double rHa = d / (Kh * Nu);
    if (rHa > rHmax) rHa = rHmax;
    return rHa;
}
// Minimum hydraulic resistance of the plant, from its height, Huber value (hv)
// and maximum xylem conductivity (Eller et al. 2020).
static double rpmin_calc(double h, double hv, double Kxmx)
{
    // Constants
    const double n_ext = 2.0;      // number of daughter branches per parent
    const double Lpet = 0.04;     // petiole length (m)
    const double r_intpet = 10.0;  // petiole conduit radius (micrometres)
    const double r_intref = 22.0;  // terminal branch conduit radius (micrometres)
    // Number of branching levels (at least 1)
    const double n_ext_13 = std::pow(n_ext, 1.0 / 3.0);
    const double arg = 1.0 - (h / Lpet) * (1.0 - n_ext_13);
    double Nh;
    if (arg <= 0.0 || !std::isfinite(arg)) {
        Nh = 1.0;
    }
    else {
        Nh = (3.0 * std::log(arg)) / std::log(n_ext);
        if (!std::isfinite(Nh) || Nh < 1.0) Nh = 1.0;
    }
    // Conduit taper factor
    const double Chi_1 = 2.888503;
    const double Chi_h = (6.6e-13 / 7.2e-13) * std::pow(Nh, 0.53);
    const double Chi_tap = Chi_h / Chi_1;
    // Maximum petiole hydraulic conductivity
    const double ratio = r_intpet / r_intref;
    const double Kpmx = Kxmx * (ratio * ratio);
    // Minimum resistance (infinite if any input is not positive)
    if (hv <= 0.0 || Kpmx <= 0.0 || Chi_tap <= 0.0 || !std::isfinite(Chi_tap)) {
        return std::numeric_limits<double>::infinity();
    }
    const double rpmin = h / (Kpmx * hv * Chi_tap);
    return rpmin;
}
// Everything stomatal conductance depends on except light (Eller et al. 2020):
// photosynthetic capacity at leaf temperature, vapour pressure deficit and water
// supply from the soil. Shared by the sunlit and shaded leaves of a layer.
struct stomatalState {
    double Vcmax, Rd, Oa, photocomp, DD, ca, ci;
    double Kc, Ko, Wc, We, Wc_ca;   // C3 only
    double zeta;
};
static stomatalState leafgsState(const envstruct& envdata, vegpstruct& vegp, double z, bool C3)
{
    stomatalState st{};
    // Vcmax temperature response
    double xx = 0.1 * (envdata.tleaf - 25.0);
    double Vcmax = vegp.Vcmax25 * std::exp2(xx) / ((1.0 + std::exp(0.3 * (envdata.tleaf - vegp.Tup))) *
        (1.0 + std::exp(0.3 * (vegp.Tlw - envdata.tleaf))));
    double Rd = vegp.fd * Vcmax;
    // Photocompensation
    double Q10rs = 0.57;
    double Oa = 0.2095 * envdata.pk * 1000;
    double photocomp = Oa / (2.0 * 2600.0 * std::pow(Q10rs, 0.1 * (envdata.tleaf - 25.0)));
    // Vapour pressure deficit between leaf and air; nil for a leaf at or below the dew point.
    double ea = satvapCpp2(envdata.tair) * (envdata.rh / 100.0);
    double es = satvapCpp2(envdata.tleaf);
    double DD = std::max(es - ea, 1e-9);
    double ca = envdata.Ca * envdata.pk / 1000.0;  // convert carbon concentration to Pa
    double ci = vegp.f0 * (1.0 - DD / vegp.Dcrit) * (ca - photocomp); // in Pa
    st.Vcmax = Vcmax; st.Rd = Rd; st.Oa = Oa; st.photocomp = photocomp; st.DD = DD; st.ca = ca; st.ci = ci;
    if (C3) {
        // Carbon-limited and transport-limited assimilation (Jacobs 1994; Eller et al. 2020, Eq. S2.1).
        double Q10Kc = 2.1;
        double Kc = 30.0 * std::pow(Q10Kc, 0.1 * (envdata.tair - 25.0));
        double Q10Ko = 1.2;
        double Ko = 30000.0 * std::pow(Q10Ko, 0.1 * (envdata.tair - 25.0));
        double Wc = Vcmax * ((ci - photocomp) /(ci + Kc * (1 + Oa / Ko))); // carbon-limited rate
        double We = 0.5 * Vcmax; // transport-limited rate
        if (Wc < 0.0) Wc = 0.0;
        if (We < 0.0) We = 0.0;
        double Wc_ca = Vcmax * ((ca - photocomp) / (ca + Kc * (1 + Oa / Ko)));
        if (Wc_ca < 0.0) Wc_ca = 0.0;
        st.Kc = Kc; st.Ko = Ko; st.Wc = Wc; st.We = We; st.Wc_ca = Wc_ca;
    }
    if (vegp.apsi < 0.0) {
        double stem_slope = 65.15 * pow(-vegp.psi50, -1.25);
        vegp.apsi = -4.0 * stem_slope / 100.0 * vegp.psi50;
    }
    double rhow = 1000.0 * (1 - (envdata.tleaf + 288.9414) * std::pow(envdata.tleaf - 3.9863, 2.0) /
        (508929.2 * (envdata.tleaf + 68.12963)));
    double psi_pd = envdata.psi_r - z * g * rhow * 1e-6; // pre-dawn water potential at this leaf's height
    double K_psi_pd = 1.0 / (1.0 + std::pow(psi_pd / vegp.psi50, vegp.apsi)); // fraction of max hydraulic conductance remaining (vulnerability curve)
    double K_50f = 0.5 / (1.0 + std::pow((psi_pd + vegp.psi50) / vegp.psi50, vegp.apsi));
    double psi_50f = (psi_pd + vegp.psi50) / 2.0;
    double dKdpKi = ((K_psi_pd - K_50f) / (psi_pd - psi_50f)) * (1.0 / K_psi_pd);
    // Minimum plant resistance is fixed for a run, so is calculated once.
    if (vegp.rpmin < 0.0) vegp.rpmin = rpmin_calc(vegp.hgt, vegp.hv, vegp.Kxmx);
    double rp = vegp.rpmin / K_psi_pd;
    double DDm = DD / envdata.pk; // vapour pressure deficit, mol/mol
    st.zeta = 2.0 / (dKdpKi * rp * 1.6 * DDm); // marginal carbon cost of water, Eller optimisation
    return st;
}
// Stomatal conductance of leaves absorbing PARabs, given the state above.
static double leafgsLight(const stomatalState& st, double PARabs, double pk, const vegpstruct& vegp, bool C3)
{
    if (!(PARabs > 0.0)) return 0.0;
    const double Vcmax = st.Vcmax, Rd = st.Rd, Oa = st.Oa, photocomp = st.photocomp, ca = st.ca, ci = st.ci;
    double IPAR = PARabs * (0.48 / 0.219) * 1e-6; // PAR in mol photons/m2/s
    double Aca;
    double Acol;
    double cicol;
    if (C3) { // C3 photosynthetic pathway
        const double Kc = st.Kc, Ko = st.Ko, We = st.We;
        // Light-limited assimilation (Jacobs 1994).
        double Wl = vegp.alpha * IPAR * ((ci - photocomp) / (ci + 2.0 * photocomp)); // light-limited rate
        if (Wl < 0.0) Wl = 0.0;
        double Wcol = ((We + Wl) - std::sqrt(std::pow(We + Wl, 2.0) - 4.0 * 0.93 * (We * Wl)))
            / (2.0 * 0.93); // co-limited rate
        Acol = Wcol - Rd;
        // Co-limiting CO2 concentration
        cicol = (-Vcmax * photocomp - Kc * (1.0 + Oa / Ko) * Wcol)
            / (Wcol - Vcmax);
        // Assimilation with stomata fully open: internal CO2 equal to ambient.
        double Wl_ca = vegp.alpha * IPAR * ((ca - photocomp) / (ca + 2.0 * photocomp));
        double We_ca = We;
        if (Wl_ca < 0.0) Wl_ca = 0.0;
        double Wca = st.Wc_ca;
        if (Wl_ca < Wca) Wca = Wl_ca;
        if (We_ca < Wca) Wca = We_ca;
        Aca = Wca - Rd;
    }
    else { // C4 pathway
        double k = 2.0e-4;
        double Wc = Vcmax;
        double Wl = vegp.alpha * IPAR;
        double We = k * Vcmax * (ci / (pk * 1000.0));
        double W = Wc;
        if (Wl < W) W = Wl;
        if (We < W) W = We;
        Aca = W - Rd; // net assimilation (can be negative)
        double Wcol = ((Wc + Wl) - std::sqrt(std::pow(Wc + Wl, 2.0) - 4.0 * 0.83 * (We * Wl)))
            / (2.0 * 0.83);
        Acol = Wcol - Rd;
        cicol = (Wcol * pk * 1000.0) / (k * Vcmax);
    }
    // Change in assimilation per unit CO2 between the co-limited point and fully
    // open stomata (Eq. S2.1).
    double dadc = (Aca - Acol) / (ca - cicol);
    if (dadc < 1e-99) dadc = 1e-99;
    double mu = 1.0 + (4.0 * st.zeta) / dadc;
    if (mu < 1.0) mu = 1.0;
    double gs = 0.5 * dadc * (std::sqrt(mu) - 1.0);
    double gsmax = vegp.Vcmax25 * 1e4;
    gs = std::fmin(gs, gsmax);
    return gs;
}
// Stomatal conductance in one call.
static double leafgs(const envstruct& envdata, vegpstruct& vegp, double z, bool C3 = true)
{
    if (!(envdata.PARabs > 0.0)) return 0.0;
    return leafgsLight(leafgsState(envdata, vegp, z, C3), envdata.PARabs, envdata.pk, vegp, C3);
}
// Fraction of a plant surface that is wet: all of it during rain, otherwise
// rising with the water held (mm). The wet part evaporates freely; the dry part
// of a leaf transpires.
static double wetFraction(double surfwater, double precip)
{
    if (precip > 0.0) return 1.0;
    return 1.0 - std::exp(-surfwater * 30.0);
}
// Latent heat of water vapour leaving a surface at Tc (J/mol): vaporisation above
// 0 degC, sublimation below -1 degC, and a linear blend between.
static inline double latentHeat(double Tc)
{
    const double vaporisation = 45068.7 - 42.8428 * Tc;
    if (Tc >= 0.0) return vaporisation;
    const double sublimation = 51078.69 - 4.338 * Tc - 0.06367 * Tc * Tc;
    return vaporisation + std::min(1.0, -Tc) * (sublimation - vaporisation);
}
// The air's part of a surface energy balance: longwave emitted at air temperature
// (Rema), heat capacity, molar density and vapour pressure deficit (Da).
struct surfaceAir { double Ta, pk, Rema, cp, ph, Da; };
static surfaceAir surfaceAirOf(double Ta, double pk, double rh, double em)
{
    surfaceAir air;
    air.Ta = Ta; air.pk = pk;
    air.Rema = em * sb * radem(Ta);
    air.cp = cpairCpp(Ta);
    air.ph = phairCpp(Ta, pk);
    air.Da = satvapCpp2(Ta) * (1.0 - rh / 100.0);
    return air;
}
// Surface temperature from the Penman-Monteith energy balance, given absorbed
// radiation (Rabs), resistances to heat (rHa) and vapour (rV) and ground heat
// flux (G). Iterated, as emission and latent heat depend on the answer.
static double PenmanMonteithCpp2(const surfaceAir& air, double Rabs, double em,
    double rHa, double rV, double G = 0.0, int iters = 4)
{
    const double Ta = air.Ta, pk = air.pk, Rema = air.Rema, cp = air.cp, ph = air.ph, Da = air.Da;
    double Ts = Ta;
    for (int i = 0; i < iters; ++i) {
        double Te = (Ts + Ta) / 2.0;
        double Rer = 4.0 * em * sb * std::pow(Te + 273.15, 3.0);
        const double la = latentHeat(Ts);
        double De = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5);
        Ts = Ta + ((Rabs - Rema - ((la * ph) / (pk * rV)) * Da - G) /
            (Rer + ph * (cp / rHa + ((la * De) / (pk * rV)))));
    }
    return (Ts);
}
// Fraction of rain reaching each canopy depth. Wind tilts the rain so that it
// passes more foliage; the beam extinction coefficient is applied along its path.
static rainmodel rainintercept(const std::vector<double>& wcm, const std::vector<double>& pia,
    double uh, double rain, double x)
{
    // Fall speed of rain (m/s)
    double vr = 3.78 * std::pow(rain, 0.067);
    int n = static_cast<int>(wcm.size());
    std::vector<double> kd(n);
    std::vector<double> raintr(n);
    for (int i = 0; i < n; ++i) {
        // calculate average wind speed above
        double uzm = wcm[i] * uh;
        // calculate rain angle from vertical (in radians)
        double rainZ = std::atan(uzm / vr);
        // calculate rain extinction coefficient
        const double si = std::cos(rainZ);
        kstruct kp = cankCpp(rainZ, x, si);
        kd[i] = kp.kd;
        // Fraction transmitted
        raintr[i] = std::exp(-kd[i] * pia[i]);
    }
    rainmodel out;
    out.tr = raintr;
    out.kd = kd;
    return out;
}
// Resistance of canopy mixing between the ground and height z: the integral of
// one over eddy diffusivity.
static double rhcanopy(double a2, double uf, double h, double z, double a0, double a1)
{
    // Gust strength rises from a0*uf at the ground to a1*uf at canopy top as a raised
    // cosine, which gives the integral in closed form.
    const double Aw = 0.5 * (a1 + a0);
    const double Bw = 0.5 * (a1 - a0);
    const double D = Aw * Aw - Bw * Bw; // positive unless gusts at the ground are nil
    const double mu = 1.0 / (a2 * h * uf);
    double inth;
    if (z == h) {
        // Canopy top
        inth = Aw * h / std::pow(D, 1.5);
    }
    else {
        const double x = pi * z / h;
        const double s = std::sin(x);
        const double c1 = std::cos(x) + 1.0;
        const double r = std::sqrt((Aw + Bw) / (Aw - Bw));
        const double sqrtD = std::sqrt(D);
        // Below canopy top. The half-angle form stays finite at mid-height.
        const double I1 = (2.0 / sqrtD) * std::atan(r * s / c1);
        const double I2 = Bw * s / (D * (Aw - Bw * std::cos(x))) + (Aw / D) * I1;
        inth = (h / pi) * I2;
    }
    return inth * mu;
}
// ~~~~~~~~~~~~~~ Diffusivity from the soil surface to the reference height ~~~~~~~~~~~~~~~ //
// One profile of eddy diffusivity runs from the soil to the reference height.
// Every resistance is the integral of its inverse between two heights, in closed
// form. It has three parts:
//   - inside the canopy: gust strength squared times a time scale;
//   - above the canopy: a surface layer with mixing enhanced near canopy top;
//   - near the soil: a surface layer driven by the gusts that reach the ground,
//     wherever that gives the lower diffusivity.
// With no foliage it is the surface layer over bare ground.
//
// Roughness height of the soil for heat: a fifth of that for momentum or, if
// larger, the height below which molecular diffusion exceeds turbulent diffusion.
constexpr double AIR_HEAT_DIFFUSIVITY = 2.1e-5;   // m2/s, near 20 deg C
static double heatRoughness(double z0, double ufloor)
{
    return std::max(0.2 * z0, AIR_HEAT_DIFFUSIVITY / (ka * std::max(ufloor, 1e-6)));
}
struct aerocolumn {
    double h, d, uf, LL, a0, a1;
    double zh;    // the soil's roughness height for heat: where the soil surface layer begins
    double a2;    // time a parcel keeps its vertical motion, as a fraction of h/uf
    double b;     // how far the air above the canopy is stretched (see stretchedHeight)
    double ug;    // velocity scale of the soil surface layer: floor gust strength as a friction velocity
    double Kc0;   // canopy diffusivity at the floor
    double zw;    // height at which the ground stops limiting eddy size
    double Rc0;   // canopy integral to the soil's roughness height
    double Rh;    // resistance from the soil to canopy top
};
// Resistance of the soil surface layer from the soil's roughness height to z.
static double soilLayerR(double z, double ug, double L, double zh)
{
    return (std::log(z / zh) + dpsihCpp2(zh / L) - dpsihCpp2(z / L)) / (ka * ug);
}
// Its rate of increase with height: one over diffusivity.
static double soilLayerSlope(double z, double ug, double L)
{
    return phihStar(z / L) / (ka * ug * z);
}
// Diffusivity of canopy mixing at height z.
static double canopyK(const aerocolumn& c, double z)
{
    const double sw = c.uf * (0.5 * (c.a1 + c.a0) - 0.5 * (c.a1 - c.a0) * std::cos(pi * z / c.h));
    return sw * sw * c.a2 * c.h / c.uf;
}
// Above the canopy, surface-layer theory is applied at a stretched height. Near
// canopy top this exceeds the true height above the displacement height, which
// gives the enhanced mixing there; far above, the two are equal.
static double stretchedHeight(const aerocolumn& c, double z)
{
    const double zp = z - c.d;
    return 0.5 * (zp + std::sqrt(zp * zp + c.b * c.b));
}
// Diffusivity above the canopy.
static double airK(const aerocolumn& c, double z)
{
    const double zp = z - c.d;
    return ka * c.uf * std::sqrt(zp * zp + c.b * c.b) / phihStar(stretchedHeight(c, z) / c.LL);
}
// Resistance between two heights above the canopy.
static double airR(const aerocolumn& c, double z1, double z2)
{
    const double s1 = stretchedHeight(c, z1), s2 = stretchedHeight(c, z2);
    return (std::log(s2 / s1) - dpsihCpp2(s2 / c.LL) + dpsihCpp2(s1 / c.LL)) / (ka * c.uf);
}
// Height up to which the ground limits eddy size: where the diffusivity of the
// soil surface layer equals canopy diffusivity at the ground. It lies above
// canopy top under sparse vegetation.
static double groundHandOver(double h, double ug, double L, double Kc0, double zh)
{
    const double target = 1.0 / Kc0;
    if (soilLayerSlope(zh * 1.0001, ug, L) <= target) return zh;
    double hi = h;
    while (soilLayerSlope(hi, ug, L) > target && hi < 1e4 * h) hi *= 2.0;
    double lo = zh;
    for (int k = 0; k < 80; ++k) {
        const double mid = std::sqrt(lo * hi);
        if (soilLayerSlope(mid, ug, L) > target) lo = mid; else hi = mid;
    }
    return std::sqrt(lo * hi);
}
// Resistance from the soil surface to height z.
static double soilToZ(const aerocolumn& c, double z)
{
    const double zc = std::max(z, c.zh);
    double r = (zc <= c.h) ? rhcanopy(c.a2, c.uf, c.h, zc, c.a0, c.a1) - c.Rc0
        : rhcanopy(c.a2, c.uf, c.h, c.h, c.a0, c.a1) - c.Rc0 + airR(c, c.h, zc);
    const double zz = std::min(zc, c.zw);
    // Add the soil surface layer, less the canopy mixing already counted over the
    // same depth.
    if (zz > c.zh) r += soilLayerR(zz, c.ug, c.LL, c.zh) - (zz - c.zh) / c.Kc0;
    return r;
}
// Diffusivity at height z.
static double columnK(const aerocolumn& c, double z)
{
    double inv = 1.0 / ((z <= c.h) ? canopyK(c, z) : airK(c, z));
    if (z < c.zw) inv += std::max(0.0, soilLayerSlope(z, c.ug, c.LL) - 1.0 / c.Kc0);
    return 1.0 / inv;
}
// Builds the diffusivity profile for the current friction velocity and stability.
static aerocolumn makeColumn(double h, double d, double uf, double LL, double a0, double a1,
    double enhancement, double z0)
{
    aerocolumn c{};
    c.h = h; c.d = d; c.uf = uf; c.LL = LL; c.a0 = a0; c.a1 = a1;
    // The soil surface layer is driven by the gusts that reach the ground.
    c.ug = (a0 / a1) * uf;
    c.zh = heatRoughness(z0, c.ug);
    const double dt = h - d;
    c.b = dt * std::sqrt(std::max(enhancement * enhancement - 1.0, 0.0));
    // Time scale of canopy mixing, set so that diffusivity at canopy top equals that
    // of the air just above. Stability enters through phi_h.
    c.a2 = ka * enhancement * (dt / h)
        / (a1 * a1 * phihStar(canopyTopStabilityHeight(h, d, enhancement) / LL));
    c.Kc0 = a0 * a0 * uf * c.a2 * h;
    c.zw = groundHandOver(h, c.ug, LL, c.Kc0, c.zh);
    c.Rc0 = rhcanopy(c.a2, uf, h, c.zh, a0, a1);
    c.Rh = soilToZ(c, h);
    return c;
}
// Energy and water balance of every canopy layer. Woody, sunlit and shaded
// surfaces are solved separately and combined to give leaf temperature, stomatal
// conductance, transpiration, evaporation of intercepted water, and the heat and
// vapour given to the air.
static void plantmodelCpp(onestep& onestepin, envstruct envdata, vegpstruct& vegp, const rainmodel& rainvars, const radmodel& swout,
    const radmodel2& lwout, const std::vector<double>& z, const std::vector<double>& dTs, double timestep = 3600.0, bool C3 = true)
{
    const int n = static_cast<int>(z.size());
    std::vector<double> tleafn(n);
    std::vector<double> swaterdepthn(n);
    std::vector<double> Ez(n);
    std::vector<double> Ezt(n);
    std::vector<double> Hz(n);
    std::vector<double> Lz(n);
    std::vector<double> rLB(n);
    std::vector<double> dSTdT(n), dSLde(n);
    onestepin.leafwet.assign(n, 0.0);
    // Rain is followed down from the top layer. Each layer intercepts part of what
    // arrives, and what it cannot hold at the end of the hour drips to the layer below.
    std::vector<double> Efilm(n);
    double arriving = envdata.precip; // water still falling, per unit ground area
    const double kdrip = cankCpp(0.0, vegp.x, 1.0).kd;
    const double leafSize = leafDimension(vegp.len, vegp.wid, vegp.x);
    // All values are per unit one-sided plant area. Both sides of a leaf exchange
    // heat and emit longwave, so heat conductance and emissivity are doubled.
    for (int i = n - 1; i >= 0; --i) {
        rLB[i] = leafrHa(onestepin.tair[i], dTs[i], onestepin.uz[i], leafSize);
        MP_TALLY(6);
        if (vegp.Lfrac[i] == 1.0) MP_TALLY(16);
        if (swout.sunfrac[i] == 0.0 || swout.sunfrac[i] == 1.0) MP_TALLY(17);
        envdata.tair = onestepin.tair[i];
        envdata.tleaf = onestepin.tleaf[i];
        envdata.rh = onestepin.rh[i];
        const double ph = phairCpp(onestepin.tair[i], envdata.pk);
        const double Tk = onestepin.tair[i] + 273.15;
        const double satvap_tair = satvapCpp2(onestepin.tair[i]);
        const double rH2 = 0.5 * rLB[i];      // both faces in parallel
        const double em2 = 2.0 * vegp.vegem;  // both faces emit
        // The wet part of a surface evaporates across the boundary layer; the dry part of
        // a leaf transpires through its stomata. A surface colder than the dew point
        // gains dew over both sides.
        const double fw = wetFraction(onestepin.swaterdepth[i], envdata.precip);
        onestepin.leafwet[i] = fw;
        double gwet = fw / rLB[i];
        double gtsun = 0.0, gtshade = 0.0;
        // Converts conductance times vapour pressure difference to mm of water over the hour.
        const double tomm = (Mw / (RgasC * Tk)) * timestep;
        // Each surface also returns how much its temperature changes per degree change in
        // air temperature (follow), which the canopy air solve uses.
        const double eair = satvap_tair * onestepin.rh[i] / 100.0;
        const double dTa = 0.1;
        const double rh_warm = 100.0 * eair / satvapCpp2(onestepin.tair[i] + dTa);
        const surfaceAir airNow = surfaceAirOf(onestepin.tair[i], envdata.pk, onestepin.rh[i], em2);
        const surfaceAir airWarm = surfaceAirOf(onestepin.tair[i] + dTa, envdata.pk, rh_warm, em2);
        // Change per degree of air temperature in emission and in vapour pressure deficit.
        const double emitSlope = (airWarm.Rema - airNow.Rema) / dTa;
        const double deficitSlope = (airWarm.Da - airNow.Da) / dTa;
        // Temperature (ts), evaporation (Ew) and transpiration (Et) of one surface.
        auto surface = [&](double Rabs, double gt, double& ts, double& Ew, double& Et,
            double& follow, double& gused) {
            const double g = gt + gwet;
            double rV = (g > 0.0) ? 1.0 / g : 1e9;
            ts = PenmanMonteithCpp2(airNow, Rabs, em2, rH2, rV, 0.0, g_solver.leafSteps);
            MP_TALLY(0);
            // vapour pressure deficit of the surface against the air (Pa)
            double DD = (satvapCpp2(ts) - eair) * 1000.0;
            if (DD >= 0.0) {
                Ew = tomm * DD * gwet;
                Et = tomm * DD * gt;
                gused = g;
            }
            else {
                // Below the dew point: condensation across the boundary layer.
                rV = rH2;
                ts = PenmanMonteithCpp2(airNow, Rabs, em2, rH2, rV, 0.0, g_solver.leafSteps);
                MP_TALLY(1);
                DD = (satvapCpp2(ts) - eair) * 1000.0;
                Ew = tomm * std::min(DD, 0.0) / rH2;
                Et = 0.0;
                gused = 1.0 / rH2;
            }
            if (!g_solver.leafClosed) {
                // Response to air temperature from a second solve in warmer air.
                const double tsw = PenmanMonteithCpp2(airWarm, Rabs, em2, rH2, rV, 0.0, g_solver.leafSteps);
                MP_TALLY(2);
                follow = (tsw - ts) / dTa;
            }
            else {
                // Response to air temperature from the linearised energy balance.
                const double latent = latentHeat(ts);
                const double vapourPath = (latent * airNow.ph) / (airNow.pk * rV);
                auto restoring = [&](double Tmid) {
                    return 4.0 * em2 * sb * std::pow(Tmid + 273.15, 3.0)
                        + airNow.ph * airNow.cp / rH2 + vapourPath * (satvapCpp2(Tmid + 0.5) - satvapCpp2(Tmid - 0.5));
                };
                const double Tmid = 0.5 * (ts + airNow.Ta);
                const double strength = restoring(Tmid);
                const double loadSlope = -emitSlope - vapourPath * deficitSlope;
                const double lean = (ts - airNow.Ta) * ((restoring(Tmid + 0.25) - strength) / 0.25) / (2.0 * strength);
                follow = (1.0 + loadSlope / strength - lean) / (1.0 + lean);
            }
        };
        // Stomatal conductance of sunlit and shaded leaves.
        double gssun = 0.0, gsshade = 0.0;
        if (swout.RPARsun[i] > 0.0 || swout.RPARshade[i] > 0.0) {
            const stomatalState stomata = leafgsState(envdata, vegp, vegp.hgt / 2.0, C3);
            MP_TALLY(5);
            gssun = leafgsLight(stomata, swout.RPARsun[i], envdata.pk, vegp, C3);
            gsshade = leafgsLight(stomata, swout.RPARshade[i], envdata.pk, vegp, C3);
        }
        const double sf = swout.sunfrac[i];
        onestepin.gs[i] = sf * gssun + (1.0 - sf) * gsshade; // average stomatal conductance
        // Water the layer can evaporate over the hour, per unit plant area: the film it
        // holds plus what it intercepts. Without rain only drip falls, vertically.
        const double kdfall = (envdata.precip > 0.0) ? rainvars.kd[i] : kdrip;
        const double intercepted = (vegp.paii[i] > 0.0 && arriving > 0.0) ? (1.0 - std::exp(-kdfall * vegp.paii[i])) * arriving : 0.0;
        const double avail = onestepin.swaterdepth[i] + ((vegp.paii[i] > 0.0) ? intercepted / vegp.paii[i] : 0.0);
        // If free evaporation would use more water than is available, the surface is wet
        // for only part of the hour (phi) and transpires for the rest. phi is refined a
        // few times, because a surface that evaporates less is warmer.
        double phi = 1.0;
        double twood, Ewwood, Etwood, fwood, gwood;
        double tsun, Ewsun, Etsun, fsun, gsun;
        double tshade, Ewshade, Etshade, fshade, gshade;
        for (int pass = 0; ; ++pass) {
            if (pass > 0) MP_TALLY(3);
            gwet = phi * fw / rLB[i];
            gtsun = (gssun > 0.0) ? (1.0 - phi * fw) / (rLB[i] + ph / gssun) : 0.0;
            gtshade = (gsshade > 0.0) ? (1.0 - phi * fw) / (rLB[i] + ph / gsshade) : 0.0;
            // A layer that is all leaf has no woody surface.
            if (vegp.Lfrac[i] == 1.0) { twood = 0.0; Ewwood = 0.0; Etwood = 0.0; fwood = 0.0; gwood = 0.0; }
            else surface(swout.RswLav[i] + lwout.RlwLabs[i], 0.0, twood, Ewwood, Etwood, fwood, gwood);
            surface(swout.RswLsun[i] + lwout.RlwLabs[i], gtsun, tsun, Ewsun, Etsun, fsun, gsun);
            // Sunlit and shaded leaves are the same in the dark and under wholly diffuse light.
            if (swout.RswLshade[i] == swout.RswLsun[i] && gtshade == gtsun) {
                MP_TALLY(4);
                tshade = tsun; Ewshade = Ewsun; Etshade = Etsun; fshade = fsun; gshade = gsun;
            }
            else surface(swout.RswLshade[i] + lwout.RlwLabs[i], gtshade, tshade, Ewshade, Etshade, fshade, gshade);
            // Layer average: leaves weighted by the leaf fraction, wood by the rest.
            Ez[i] = vegp.Lfrac[i] * (sf * Ewsun + (1.0 - sf) * Ewshade) + (1.0 - vegp.Lfrac[i]) * Ewwood; // per unit plant area
            if (Ez[i] <= 1.01 * avail || pass == 3) break;
            phi *= avail / Ez[i];
        }
        // Water left on the layer, up to the maximum film; the excess drips to the layer below.
        Efilm[i] = std::min(Ez[i], avail);
        const double remaining = avail - Efilm[i];
        swaterdepthn[i] = std::min(remaining, vegp.mwft); // leaf surface water depth
        arriving += (remaining - swaterdepthn[i]) * vegp.paii[i] - intercepted;
        const double tgreen = sf * tsun + (1.0 - sf) * tshade; // temperature of leaves
        tleafn[i] = vegp.Lfrac[i] * tgreen + (1.0 - vegp.Lfrac[i]) * twood; // temperature of foliage including woody
        Ezt[i] = vegp.Lfrac[i] * vegp.paii[i] * (sf * Etsun + (1.0 - sf) * Etshade); // per unit ground area
        const double cp = cpairCpp(onestepin.tair[i]);
        Hz[i] = ((ph * cp) / rH2) * (tleafn[i] - onestepin.tair[i]);
        // Change in the layer's heat and vapour sources per unit change in its air, per unit ground area.
        const double follow = vegp.Lfrac[i] * (sf * fsun + (1.0 - sf) * fshade) + (1.0 - vegp.Lfrac[i]) * fwood;
        const double gV = vegp.Lfrac[i] * (sf * gsun + (1.0 - sf) * gshade) + (1.0 - vegp.Lfrac[i]) * gwood;
        const double la_leaf = latentHeat(tleafn[i]);
        dSTdT[i] = vegp.paii[i] * ((ph * cp) / rH2) * (follow - 1.0);
        dSLde[i] = -vegp.paii[i] * la_leaf * gV * 1000.0 / (RgasC * Tk);
    }
    onestepin.precipground = arriving;
    // Latent heat given to the air by each layer, per unit ground area:
    // transpiration plus evaporation of intercepted water.
    for (int i = 0; i < n; ++i) {
        const double la = latentHeat(tleafn[i]);
        Lz[i] = (la / Mw) * (Ezt[i] + vegp.paii[i] * Efilm[i]) / timestep;
    }
    // Transpiration drawn from the soil through the roots
    double Et = 0.0;
    for (int i = 0; i < n; ++i) if (Ezt[i] > 0.0) Et += Ezt[i];
    onestepin.tleaf = std::move(tleafn);
    onestepin.swaterdepth = std::move(swaterdepthn);
    onestepin.Et = Et;
    onestepin.Hz = std::move(Hz);
    onestepin.Lz = std::move(Lz);
    onestepin.dSTdT = std::move(dSTdT);
    onestepin.dSLde = std::move(dSLde);
    onestepin.rLB = std::move(rLB);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Soil model *********************************************************** //
// Soil temperature and water content, carried from hour to hour.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Relative humidity of the air in the surface soil's pores, from its water
// potential (Kelvin equation).
static double soilrelhumCpp(const soilpstruct& soilp, double Tsoil, double theta)
{
    double psiw = soilp.psie[0] * std::pow(theta / soilp.thetaS[0], -soilp.b[0]);
    double Tk = Tsoil + 273.15;
    double hr = std::exp(Mw * psiw / (RgasC * Tk));
    return hr;
}
// Latent heat flux from the soil surface to the air, across resistance rHa.
static double groundLatentFlux(const soilpstruct& soilp, double Tsurface, double Tref,
    double relhum, double pk, double rHa, double theta)
{
    const double ph = phairCpp(Tref, pk);
    const double es = satvapCpp2(Tsurface) * soilrelhumCpp(soilp, Tsurface, theta);
    const double ea = satvapCpp2(Tref) * relhum / 100.0;
    const double la = latentHeat(Tsurface);
    return ((la * ph) / (rHa * pk)) * (es - ea);
}
// Energy available at the soil surface at temperature Tsurface: net radiation
// less sensible and latent heat. It is the heat flux into the soil.
static double soilsurfaceEB(const soilpstruct& soilp, double Rabs, double Tref,
    double Tsurface, double pk, double relhum, double rHa, double theta)
{
    double sb = 5.67e-8;
    double Rnet = Rabs - soilp.groundem * sb * radem(Tsurface);
    double cp = cpairCpp(Tref);
    double ph = phairCpp(Tref, pk);
    double H = ((ph * cp) / rHa) * (Tsurface - Tref);
    double L = groundLatentFlux(soilp, Tsurface, Tref, relhum, pk, rHa, theta);
    double Ba = Rnet - H - L;
    return Ba;
}
// Geometric or logarithmic mean of the conductivities of two adjacent soil nodes.
static double kMeanCpp(std::string meanType, double k1, double k2) {
    if (meanType == "GEOMETRIC") {
        return std::sqrt(k1 * k2);
    }
    else if (meanType == "LOGARITHMIC") {
        if (k1 == k2) {
            return k1;
        }
        else {
            return (k1 - k2) / std::log(k1 / k2);
        }
    }
    else {
        stop("Unknown mean type. Use 'GEOMETRIC' or 'LOGARITHMIC'.");
    }
}
// The parts of soil thermal conductivity that do not depend on temperature:
// conductivity of the solids, grain shape factors, pore space, and how far
// vapour contributes (wf).
struct soilConductivityFixed { double Vsolid, kSolid, Ga, Gc, porosity, gasPorosity, wf, wetShareSq; };
static soilConductivityFixed soilConductivityFixedOf(double Vq, double Vm, double Vo, double Vw, double Mc)
{
    soilConductivityFixed f;
    // Volume fraction of solids
    double Vsolid = Vq + Vm + Vo;
    // Conductivity of solids: quartz 8.8, other minerals 2.5, organic matter 0.25 W/m/K.
    f.kSolid = std::pow(8.8, Vq / Vsolid) * std::pow(2.5, Vm / Vsolid) * std::pow(0.25, Vo / Vsolid);
    // Shape factors
    f.Ga = 0.088 * ((Vq + Vm) / Vsolid) + 0.5 * (Vo / Vsolid);
    double Q = 7.25 * Mc + 2.52;
    double xwo = 0.33 * Mc + 0.078;
    double porosity = 1.0 - Vsolid;
    double gasPorosity = porosity - Vw;
    if (gasPorosity < 0.0)  gasPorosity = 0.0;
    if (Vw < 0.01 * xwo) {
        f.wf = 0.0;
    }
    else {
        f.wf = 1.0 / (1.0 + std::pow(Vw / xwo, -Q));
    }
    f.Gc = 1.0 - 2.0 * f.Ga;
    f.Vsolid = Vsolid; f.porosity = porosity; f.gasPorosity = gasPorosity;
    f.wetShareSq = pow(Vw / porosity, 2);
    return f;
}
// Soil thermal conductivity (W/m/K): water, ice, moist air and solids, each
// weighted by how it is arranged around the grains (de Vries).
static double thermalConductivityFrom(const soilConductivityFixed& f,
    double Vw, // volumetric water plus ice (m3/m3)
    double Vw_ice, // volumetric ice, part of Vw (m3/m3)
    double Tc, // temperature (deg C)
    double pk) // pressure (kPa)
{
    const double Vsolid = f.Vsolid, kSolid = f.kSolid, Ga = f.Ga, Gc = f.Gc, gasPorosity = f.gasPorosity, wf = f.wf;
    double Tk = Tc + 273.15;
    double Lv = 45144.0 - 48.0 * Tc;
    double svp = 0.611 * std::exp(17.502 * Tc / (Tc + 240.97));
    double slope = 17.502 * 240.97 * svp / std::pow(240.97 + Tc, 2);
    double Dv = 0.0000212 * (101.3 / pk) * std::pow(Tk / 273.15, 1.75);
    double rhoAir = 44.65 * (pk / 101.3) * (273.15 / Tk);
    double stcor = 1.0 - svp / pk;
    if (stcor < 0.3) stcor = 0.3;
    double kWater = 0.56 + 0.0018 * Tc;
    // Conductivity of water and of ice, combined by the fraction frozen.
    double kIce = 2.18 - 0.0068 * Tc;
    double kWaterPhase = (Vw > 1e-9) ? (kWater * (Vw - Vw_ice) + kIce * Vw_ice) / Vw : kWater;
    double kGas = 0.0242 + 0.00007 * Tc + wf * Lv * rhoAir * Dv * slope / (pk * stcor);
    double kFluid = kGas + (kWaterPhase - kGas) * f.wetShareSq;
    double ka = (2.0 / (1.0 + (kGas / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kGas / kFluid - 1.0) * Gc)) / 3.0;
    double kw = (2.0 / (1.0 + (kWaterPhase / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kWaterPhase / kFluid - 1.0) * Gc)) / 3.0;
    double ks = (2.0 / (1.0 + (kSolid / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kSolid / kFluid - 1.0) * Gc)) / 3.0;
    double out = (kw * kWaterPhase * Vw + ka * kGas * gasPorosity + ks * kSolid * Vsolid)
        / (kw * Vw + ka * gasPorosity + ks * Vsolid);
    return out;
}
// Soil thermal conductivity in one call.
static double thermalConductivityCpp(
    double Vq, // volumetric quartz (m3/m3)
    double Vm, // volumetric other minerals (m3/m3)
    double Vo, // volumetric organic matter (m3/m3)
    double Vw, // volumetric water plus ice (m3/m3)
    double Vw_ice, // volumetric ice, part of Vw (m3/m3)
    double Mc, // mass fraction of clay (kg/kg)
    double Tc, // temperature (deg C)
    double pk) // pressure (kPa)
{
    return thermalConductivityFrom(soilConductivityFixedOf(Vq, Vm, Vo, Vw, Mc), Vw, Vw_ice, Tc, pk);
}
// Solves a tridiagonal system between rows first and last (Thomas algorithm).
static Thomas ThomasBoundaryCondition(std::vector<double> aa, std::vector<double> bb,
    std::vector<double> cc, std::vector<double> dd, std::vector<double> x,
    int first, int last)
{
    for (int i = first; i < last; ++i) {
        cc[i] = cc[i] / bb[i];
        dd[i] = dd[i] / bb[i];
        bb[i + 1] -= aa[i + 1] * cc[i];
        dd[i + 1] -= aa[i + 1] * dd[i];
    }
    x[last] = dd[last] / bb[last];
    for (int i = (last - 1); i >= first; --i) {
        x[i] = dd[i] - cc[i] * x[i + 1];
    }
    Thomas out;
    out.bb = bb;
    out.cc = cc;
    out.dd = dd;
    out.x = x;
    return out;
}
// Degree of saturation from water potential (Campbell retention curve).
static double degreeOfSaturation(const soilpstruct& soilp, double psiw, int i) {
    if (psiw >= 0) return 1.0;
    double Se;
    if (psiw >= soilp.psie[i]) {
        Se = 1.0;
    }
    else {
        Se = std::pow(psiw / soilp.psie[i], -1.0 / soilp.b[i]);
    }
    return Se;
}
// Water potential (J/kg) from volumetric water content.
static double waterPotential(const soilpstruct& soilp, double theta, int i)
{
    double Se = 1.0;
    if (theta < soilp.thetaS[i]) Se = theta / soilp.thetaS[i];
    double psiw = soilp.psie[i] * std::pow(Se, -soilp.b[i]);
    return psiw;
}
// Volumetric water content from water potential.
static double thetaFromPsi(const soilpstruct& soilp, double psiw, int i)
{
    double Se = degreeOfSaturation(soilp, psiw, i);
    double theta = Se * soilp.thetaS[i];
    return theta;
}
// Water potential of the root zone (J/kg), which the stomata respond to: the
// potential at the mean liquid water content of the rooted layers, weighted by
// root fraction.
static double unfrozenTheta(const soilpstruct& soilp, double theta, double Tc, int i);
static double rootZonePotential(const soilpstruct& soilp, const std::vector<double>& theta,
    const std::vector<double>& rootfrac, const std::vector<double>& Tc)
{
    double sum = 0.0, wgt = 0.0;
    int first = -1;
    for (size_t i = 0; i < rootfrac.size() && i < theta.size(); ++i) {
        if (rootfrac[i] <= 0.0) continue;
        if (first < 0) first = static_cast<int>(i);
        // Only liquid water counts in a layer holding ice.
        const double liquid = (i < Tc.size() && Tc[i] < 0.0)
            ? std::min(theta[i], unfrozenTheta(soilp, theta[i], Tc[i], static_cast<int>(i))) : theta[i];
        sum += rootfrac[i] * liquid;
        wgt += rootfrac[i];
    }
    if (first < 0) return waterPotential(soilp, theta[0], 0);
    return waterPotential(soilp, sum / wgt, first);
}
// Liquid water content of soil below 0 degC. Freezing lowers the water potential
// in proportion to degrees below zero (Clapeyron), and the retention curve gives
// the water that stays liquid at that potential.
static double unfrozenTheta(const soilpstruct& soilp, double theta, double Tc, int i)
{
    if (Tc >= 0.0) return theta;
    constexpr double Lf = 334000.0; // J/kg, latent heat of fusion
    constexpr double T0 = 273.15;   // K
    double psi_ice = Lf * Tc / T0;  // J/kg, negative for Tc < 0
    double theta_u = thetaFromPsi(soilp, psi_ice, i);
    if (theta_u > theta) theta_u = theta;
    if (theta_u < 0.0) theta_u = 0.0;
    return theta_u;
}
// Hydraulic conductivity from water content (Campbell).
static double hydraulicConductivityFromTheta(const soilpstruct& soilp,
    double theta, int i) {
    double n = 2.0 * soilp.b[i] + 3.0;
    double k = soilp.Ksat[i] * std::pow(theta / soilp.thetaS[i], n);
    return k;
}
// Water vapour held in the pore air, from water potential and temperature.
static double vaporFromPsi(const soilpstruct& soilp, double psiw, double theta, double Tk, int i) {
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vapor = (soilp.thetaS[i] - theta) * satVapDensityCpp2(Tk) * humidity;
    return vapor;
}
// Slope of the retention curve: change in water content per unit water potential.
static double dTheta_dPsi(const soilpstruct& soilp, double psiw, int i)
{
    double psie = soilp.psie[i]; // assumed negative
    if (psiw >= psie) return 0.0;
    double Se = std::pow(psiw / psie, -1.0 / soilp.b[i]);
    double theta = soilp.thetaS[i] * Se;
    return -theta / (soilp.b[i] * psiw);
}
// Conductivity of the soil for water vapour.
static double vaporConductivityFromPsiTheta(const soilpstruct& soilp,
    double psiw, double theta, double Tk, int i)
{
    double dv = 0.000024;
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vp = satVapDensityCpp2(Tk);
    double k = 0.66 * (soilp.thetaS[i] - theta) * dv * vp * humidity * Mw / (RgasC * Tk);
    return k;
}
// Change in pore vapour per unit water potential.
static double dvapor_dPsi(const soilpstruct& soilp, double psiw,
    double theta, double Tk, int i)
{
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vp = satVapDensityCpp2(Tk);
    double capacity_vapor = (soilp.thetaS[i] - theta) * vp * humidity *
        (Mw / (RgasC * Tk)) - dTheta_dPsi(soilp, psiw, i) * vp * humidity;
    return capacity_vapor;
}
// Evaporation from the soil surface over time dT (mm).
static double evaporation_flux(const soilpstruct& soilp, double theta, double Tsurface,
    double Tair, double relhum, double rHa, double dT)
{
    double psiw = waterPotential(soilp, theta, 0); // J/kg
    double Tk = Tair + 273.15; // Temperature deg C -> K
    double hs = std::exp(Mw * psiw / (RgasC * Tk)); // relative humidity of soil air (0 to 1)
    double es = 1000.0 * satvapCpp2(Tsurface) * hs;            // kPa -> Pa
    double ea = 1000.0 * satvapCpp2(Tair) * (relhum / 100.0);  // kPa -> Pa
    double Ev = (Mw / (rHa * RgasC * Tk)) * (es - ea) * dT; // mm
    return Ev;
}
// Reduction of root uptake in wet soil: 0 at saturation, 1 beyond air entry.
static double alpha_wet(double psiw, double psie)
{
    if (psiw >= 0) return(0.0); // saturated
    if (psiw <= psie) return(1.0); // beyond air entry to well aerated
    return std::abs(psiw) / std::abs(psie);
}
// Reduction of root uptake in dry soil: 1 above psi_dry, 0 at the wilting point.
static double alpha_dry(double psiw, double psi_dry, double psi_wilt) {
    if (psiw <= psi_wilt) return(0.0);
    if (psiw >= psi_dry)  return(1.0);
    return (psiw - psi_wilt) / (psi_dry - psi_wilt); // linear ramp 0 to 1
}
// Shares transpiration among soil layers by root fraction, reduced where the
// soil is too wet or too dry.
static std::vector<double> transpiration_distribute(const soilpstruct& soilp, const std::vector<double>& rootfrac,
    double totalTransp_mm, double dT, const std::vector<double>& psiw, double p = 0.5)
{
    int n = static_cast<int>(psiw.size());
    std::vector<double> S(n);
    if (totalTransp_mm > 0.0) {
        std::vector<double> w(n);
        double Trate = totalTransp_mm / dT;
        for (int i = 0; i < n; ++i) {
            // Uptake falls once the fraction p of the water between field capacity and
            // the wilting point has been used, and stops at the wilting point.
            double theta_fc = thetaFromPsi(soilp, PSI_FIELD_CAPACITY, i);
            double theta_wilt = thetaFromPsi(soilp, PSI_WILTING_POINT, i);
            double psi_wilt = PSI_WILTING_POINT;
            double psi_dry = waterPotential(soilp, theta_fc - p * (theta_fc - theta_wilt), i);
            double aw = alpha_wet(psiw[i], soilp.psie[i]);
            double ad = alpha_dry(psiw[i], psi_dry, psi_wilt);
            double alpha = aw * ad;
            w[i] = rootfrac[i] * alpha;
        }
        double sumw = std::accumulate(w.begin(), w.end(), 0.0);
        // If every rooted layer is past wilting, uptake is shared by root fraction alone.
        if (!(sumw > 0.0)) {
            for (int i = 0; i < n; ++i) w[i] = rootfrac[i];
            sumw = std::accumulate(w.begin(), w.end(), 0.0);
        }
        for (int i = 0; i < n; ++i) {
            if (sumw > 0.0) {
                w[i] = w[i] / sumw;
            }
            else {
                w[i] = 0.0;
            }
            S[i] = w[i] * Trate;
        }
    }
    else {
        for (int i = 0; i < n; ++i) S[i] = 0.0;
    }
    return S;
}
// Depths of the soil nodes: thin layers near the surface, thickening with depth.
// Node n, the lower boundary, is at totalDepth. One further depth is returned
// below it.
// [[Rcpp::export]]
std::vector<double> geometricCpp(int n, double totalDepth) {
    std::vector<double> z(n + 2);
    double weightSum = 0.0;
    for (int i = 1; i <= n; ++i) {
        weightSum += static_cast<double>(i) * static_cast<double>(i);
    }
    double dz_unit = totalDepth / weightSum;
    z[0] = 0.0;
    z[1] = dz_unit; // thin top layer
    for (int i = 2; i <= n + 1; ++i) {
        z[i] = z[i - 1] + dz_unit * static_cast<double>(i) * static_cast<double>(i);
    }
    return z;
}
// Heat content of soil (J/m3), relative to the same soil at 0 degC with all its
// water liquid: sensible heat of solids, water, ice and air, less the latent heat
// of fusion of the ice.
static double soilHeatContent(const soilpstruct& soilp, double theta, double Tc, int i)
{
    constexpr double waterDensity = 1000.0;     // kg/m3
    constexpr double fusionHeat = 334000.0;     // J/kg
    constexpr double airHeatCapacity = 0.0013;  // MJ/m3/K
    double Vq = soilp.Vq[i], Vm = soilp.Vm[i], Vo = soilp.Vo[i], Vw = theta;
    double ice = theta - unfrozenTheta(soilp, theta, Tc, i);
    if (ice < 0.0) ice = 0.0;
    // Air fills the remaining volume. Fractions summing to more than 1 are rescaled.
    const double Vsum = Vq + Vm + Vo + Vw;
    double Va = 0.0;
    if (Vsum > 1.0) { Vq /= Vsum; Vm /= Vsum; Vo /= Vsum; Vw /= Vsum; ice /= Vsum; }
    else Va = 1.0 - Vsum;
    // Sensible heat (MJ/m3). The heat capacity of ice varies with temperature, so is
    // integrated from 0 degC.
    const double sensible = (Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + (Vw - ice) * 4.18 + Va * airHeatCapacity) * Tc
        + ice * (1.93 * Tc + 0.00335 * Tc * Tc);
    return sensible * 1e6 - waterDensity * fusionHeat * ice;
}
// Temperature of soil with a given heat content (J/m3): the inverse of
// soilHeatContent. Direct for unfrozen soil; by iteration where there is ice.
static double soilTemperatureFromHeat(const soilpstruct& soilp, double theta, double heat, int i)
{
    const double perDegree = soilHeatContent(soilp, theta, 1.0, i);
    const double unfrozen = heat / perDegree;
    if (unfrozen >= 0.0 || unfrozenTheta(soilp, theta, unfrozen, i) >= theta) return unfrozen;
    double cold = -100.0, warm = 0.0;
    double heatCold = soilHeatContent(soilp, theta, cold, i) - heat;
    double heatWarm = -heat; // heat content is nil at 0 degC
    if (heatCold >= 0.0) return cold;
    // False position, with the Illinois modification.
    int stuck = 0;
    double T = 0.5 * (cold + warm);
    for (int it = 0; it < 200; ++it) {
        T = (cold * heatWarm - warm * heatCold) / (heatWarm - heatCold);
        const double excess = soilHeatContent(soilp, theta, T, i) - heat;
        if (std::abs(excess) < 1e-3 || warm - cold < 1e-12) break;
        if (excess > 0.0) {
            warm = T; heatWarm = excess;
            if (stuck == 1) heatCold *= 0.5;
            stuck = 1;
        }
        else {
            cold = T; heatCold = excess;
            if (stuck == -1) heatWarm *= 0.5;
            stuck = -1;
        }
    }
    return T;
}
// Tolerance for the soil heat balance: the largest imbalance at any node (W/m2).
constexpr double SOIL_HEAT_TOL = 1e-3;
// Soil temperatures after one time step. The surface energy balance is the upper
// boundary; the deepest node is held at a fixed temperature. At each node the
// change in heat content over the step (which includes the latent heat of
// freezing) equals the heat conducted in at end-of-step temperatures. Solved by
// Newton iteration.
static soilmod SoilHeatCpp(soilmod state, const soilpstruct& soilp, double Rabs, double Tref, double relhum, double atmPressure,
    double rHa, double dT = 3600.0, double tolerance = SOIL_HEAT_TOL)
{
    // Iteration limit
    constexpr int maxNrIterations = 500;
    const int n = state.n;
    const std::vector<double>& wc = state.wc;
    const std::vector<double>& dz = state.dz;
    // Heat content at the start of the step, and a first guess for the end.
    std::vector<double> heatAtStart(n + 1), heat(n + 1), heatTrial(n + 1), T(n + 1), Ttrial(n + 1);
    for (int i = 0; i <= n; ++i) {
        heatAtStart[i] = soilHeatContent(soilp, wc[i], state.oldTe[i], i);
        T[i] = (i < n) ? state.Te[i] : state.oldTe[n];
        heat[i] = soilHeatContent(soilp, wc[i], T[i], i);
    }
    std::vector<double> ff(n + 1), ffTrial(n + 1), imbalance(n + 1), imbalanceTrial(n + 1);
    std::vector<double> aa(n + 1), bb(n + 1), cc(n + 1), dd(n + 1), step(n + 1, 0.0);
    std::vector<double> warmingPerHeat(n + 1, 0.0), coolingPerHeat(n + 1, 0.0);
    std::vector<double> conductivityWarming(n + 1, 0.0), conductivityCooling(n + 1, 0.0);
    std::vector<char> warming(n + 1, 1), byTemperature(n + 1, 0);
    std::vector<double> degreesStep(n + 1, 0.0);
    // Temperature at which each node's water begins to freeze: the freezing
    // point lowered by the suction the soil already holds its water at.
    std::vector<double> freezingOnset(n + 1, 0.0);
    for (int i = 0; i <= n; ++i) freezingOnset[i] = waterPotential(soilp, wc[i], i) * 273.15 / 334000.0;
    double restoringWarming = 0.0, restoringCooling = 0.0;

    // Conductivity of the soil at a node, given its temperature and so its ice.
    std::vector<soilConductivityFixed> conductivityFixed(n + 1);
    for (int i = 0; i <= n; ++i) conductivityFixed[i] = soilConductivityFixedOf(soilp.Vq[i], soilp.Vm[i], soilp.Vo[i], wc[i], soilp.Mc[i]);
    auto conductivityAt = [&](int i, double Ti) -> double {
        double ice = wc[i] - unfrozenTheta(soilp, wc[i], Ti, i);
        if (ice < 0.0) ice = 0.0;
        return thermalConductivityFrom(conductivityFixed[i], wc[i], ice, Ti, atmPressure);
    };
    // Imbalance at each node (W/m2): heat gained over the step less heat conducted
    // in, and at the top node less the surface energy balance. Returns the largest.
    auto heatImbalance = [&](const std::vector<double>& ht, const std::vector<double>& Tt, std::vector<double>& R, std::vector<double>& f,
        std::vector<double>& conductivity) -> double {
        double lambda0 = 0.0, lambda1 = 0.0;
        for (int i = 0; i < n; ++i) {
            const double lambda = conductivityAt(i, Tt[i]);
            conductivity[i] = lambda;
            if (i == 0) lambda0 = lambda;
            if (i == 1) lambda1 = lambda;
            f[i] = lambda / dz[i];
        }
        f[0] = kMeanCpp("LOGARITHMIC", lambda0, lambda1) / dz[0];
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double in = f[i] * (Tt[i + 1] - Tt[i]);
            if (i == 0) in += soilsurfaceEB(soilp, Rabs, Tref, Tt[0], atmPressure, relhum, rHa, wc[0]);
            else in += f[i - 1] * (Tt[i - 1] - Tt[i]);
            R[i] = state.vol[i] * (ht[i] - heatAtStart[i]) / dT - in;
            if (!std::isfinite(R[i])) worst = 1e300;
            else if (std::abs(R[i]) > worst) worst = std::abs(R[i]);
        }
        R[n] = 0.0;
        return worst;
    };

    std::vector<double> lambdaNow(n + 1, 0.0), lambdaTrial(n + 1, 0.0);
    double worst = heatImbalance(heat, T, imbalance, ff, lambdaNow);
    int nrIterations = 0;
    tolerance *= g_solver.soilHeatLoosening;
    while (worst > tolerance && nrIterations < maxNrIterations) {
        // Newton step in heat content. A node's temperature change per unit heat, and its
        // change in conductivity, differ greatly either side of the temperature at which
        // its water begins to freeze, so both are found separately for warming and
        // cooling. The step is re-solved if a node moves the opposite way to that assumed.
        constexpr double probe = 100.0; // J/m3
        for (int i = 0; i < n; ++i) {
            const double warmer = soilTemperatureFromHeat(soilp, wc[i], heat[i] + probe, i);
            const double cooler = soilTemperatureFromHeat(soilp, wc[i], heat[i] - probe, i);
            warmingPerHeat[i] = (warmer - T[i]) / probe;
            coolingPerHeat[i] = (T[i] - cooler) / probe;
            if (i == 0) {
                // Change in the surface energy balance per unit heat gained by the surface node.
                const double q = soilsurfaceEB(soilp, Rabs, Tref, T[0], atmPressure, relhum, rHa, wc[0]);
                restoringWarming = (q - soilsurfaceEB(soilp, Rabs, Tref, warmer, atmPressure, relhum, rHa, wc[0])) / probe;
                restoringCooling = (soilsurfaceEB(soilp, Rabs, Tref, cooler, atmPressure, relhum, rHa, wc[0]) - q) / probe;
            }
            conductivityWarming[i] = 0.0; conductivityCooling[i] = 0.0;
            if (cooler < 0.0) {
                const double lambda = lambdaNow[i];
                conductivityWarming[i] = (conductivityAt(i, warmer) - lambda) / probe;
                conductivityCooling[i] = (lambda - conductivityAt(i, cooler)) / probe;
            }
            warming[i] = (step[i] >= 0.0);
        }
        // Change in the mean conductivity of the top two nodes per unit change in each.
        const double lambda0 = lambdaNow[0], lambda1 = lambdaNow[1];
        const double meanPerTop = (kMeanCpp("LOGARITHMIC", lambda0 * (1.0 + 1e-6), lambda1) - kMeanCpp("LOGARITHMIC", lambda0, lambda1)) / (lambda0 * 1e-6);
        const double meanPerSecond = (kMeanCpp("LOGARITHMIC", lambda0, lambda1 * (1.0 + 1e-6)) - kMeanCpp("LOGARITHMIC", lambda0, lambda1)) / (lambda1 * 1e-6);
        for (int round = 0; round < 4; ++round) {
            auto perHeat = [&](int i) { return (i >= n) ? 0.0 : (warming[i] ? warmingPerHeat[i] : coolingPerHeat[i]); };
            // Change in the conductance below node i per unit heat it gains.
            auto conductanceSlope = [&](int i) { return (warming[i] ? conductivityWarming[i] : conductivityCooling[i]) / dz[i]; };
            // The conductance below the surface node depends on both of the top two nodes.
            const double topToOwn = meanPerTop * conductanceSlope(0);
            const double topToSecond = meanPerSecond * (warming[1] ? conductivityWarming[1] : conductivityCooling[1]) / dz[0];
            for (int i = 0; i < n; ++i) {
                double below = 0.0, own, above = -ff[i] * perHeat(i + 1);
                if (i == 0) {
                    own = ff[0] * perHeat(0) + (warming[0] ? restoringWarming : restoringCooling) - topToOwn * (T[1] - T[0]);
                    above -= topToSecond * (T[1] - T[0]);
                }
                else {
                    own = (ff[i - 1] + ff[i]) * perHeat(i) - conductanceSlope(i) * (T[i + 1] - T[i]);
                    if (i == 1) {
                        below = -ff[0] * perHeat(0) - topToOwn * (T[0] - T[1]);
                        own -= topToSecond * (T[0] - T[1]);
                    }
                    else below = -ff[i - 1] * perHeat(i - 1) - conductanceSlope(i - 1) * (T[i - 1] - T[i]);
                }
                aa[i] = below;
                bb[i] = state.vol[i] / dT + own;
                cc[i] = above;
                dd[i] = -imbalance[i];
            }
            aa[n] = 0.0; bb[n] = 1.0; cc[n] = 0.0; dd[n] = 0.0;
            step = ThomasBoundaryCondition(aa, bb, cc, dd, step, 0, n).x;
            bool sidesHold = true;
            for (int i = 0; i < n; ++i) {
                if ((step[i] >= 0.0) != warming[i]) { warming[i] = !warming[i]; sidesHold = false; }
            }
            if (sidesHold) break;
        }
        // The step is applied as a change in temperature where exchange with neighbours
        // dominates a node's balance (thin layers), and as a change in heat content
        // where storage dominates (thick layers).
        for (int i = 0; i < n; ++i) {
            const double perHeat = warming[i] ? warmingPerHeat[i] : coolingPerHeat[i];
            byTemperature[i] = (bb[i] - state.vol[i] / dT > state.vol[i] / dT);
            degreesStep[i] = perHeat * step[i];
        }
        // The step is shortened so that no node passes the temperature at which its
        // water begins to freeze: the first to reach it is placed exactly there.
        double share = 1.0;
        int reaching = -1;
        for (int i = 0; i < n; ++i) {
            const double fromOnset = T[i] - freezingOnset[i];
            if (std::abs(fromOnset) < 1e-10) continue;
            const double full = byTemperature[i] ? T[i] + degreesStep[i]
                : soilTemperatureFromHeat(soilp, wc[i], heat[i] + step[i], i);
            if ((full - freezingOnset[i]) * fromOnset >= 0.0) continue;
            const double part = byTemperature[i] ? -fromOnset / degreesStep[i]
                : (soilHeatContent(soilp, wc[i], freezingOnset[i], i) - heat[i]) / step[i];
            if (part < share) { share = part; reaching = i; }
        }
        // The step is halved until the imbalance falls.
        double worstTrial = worst;
        for (int halvings = 0; halvings < 30; ++halvings) {
            for (int i = 0; i < n; ++i) {
                if (i == reaching) {
                    Ttrial[i] = freezingOnset[i];
                    heatTrial[i] = soilHeatContent(soilp, wc[i], Ttrial[i], i);
                }
                else if (byTemperature[i]) {
                    Ttrial[i] = T[i] + share * degreesStep[i];
                    heatTrial[i] = soilHeatContent(soilp, wc[i], Ttrial[i], i);
                }
                else {
                    heatTrial[i] = heat[i] + share * step[i];
                    Ttrial[i] = soilTemperatureFromHeat(soilp, wc[i], heatTrial[i], i);
                }
            }
            reaching = -1; // on the first trial only
            heatTrial[n] = heat[n]; Ttrial[n] = T[n];
            worstTrial = heatImbalance(heatTrial, Ttrial, imbalanceTrial, ffTrial, lambdaTrial);
            if (worstTrial < worst) break;
            share *= 0.5;
        }
        // Stop if no step reduces the imbalance.
        if (!(worstTrial < worst)) { nrIterations = maxNrIterations; break; }
        heat = heatTrial; T = Ttrial; imbalance = imbalanceTrial; ff = ffTrial; lambdaNow = lambdaTrial; worst = worstTrial;
        ++nrIterations;
    }
    state.Te = T;
    state.iters = nrIterations;
    return state;
}
// Fraction of roots in each soil layer. The fraction above depth D is
// 1/(1 + (D/D50)^c) (Schenk and Jackson 2002, Eq. 1), where D50 and D95 are the
// depths above which 50% and 95% of roots lie. The deepest solved node takes all
// roots below it; the boundary node has none.
static std::vector<double> root_distribute(const std::vector<double>& z, double D50, double D95)
{
    int n = static_cast<int>(z.size()) - 2;   // nodes 0..n, plus one depth below the boundary
    std::vector<double> v(n + 1, 0.0);
    // The exponent c follows from D95: (D95/D50)^c = 19.
    double p = std::log(19.0) / std::log(D95 / D50);
    auto cumulative = [&](double D) {
        if (D <= 0.0) return 0.0;
        double x = std::pow(D / D50, p);
        return x / (1.0 + x);
    };
    double above = 0.0;
    for (int i = 0; i < n; ++i) {
        double below = (i < n - 1) ? cumulative(0.5 * (z[i] + z[i + 1])) : 1.0;
        v[i] = below - above;
        above = below;
    }
    return v;
}
// Tolerance for the soil water balance: water unaccounted for, summed over layers.
constexpr double SOIL_WATER_TOL = 1e-6;   // kg m-2 s-1
// Air-filled pore space (m3/m3) below which frozen soil takes in less and less
// water, and none once full. A chosen value.
constexpr double FROZEN_PORE_CLOSING = 0.01;
// Soil water after one time step: liquid flow and vapour diffusion between layers
// (Campbell 1985), with rain and evaporation at the surface, root uptake, and
// drainage at the base. Solved by Newton iteration.
// Where a node holds ice, its liquid water is at the potential set by its
// temperature, conducts as unfrozen soil would at that potential, and is
// obstructed by the ice.
static soilwaterout SoilWaterCpp(soilwatermod soilmod, const soilpstruct& soilp,
    const climforwaterstruct& climdata, double dT = 3600.0, double pTAW = 0.5,
    int maxNrIterations = 100, double tolerance = SOIL_WATER_TOL)
{
    const double rho = 1000.0;
    const int n = soilp.nLayers;
    // State at the start of the step
    const std::vector<double> oldtheta = soilmod.oldtheta;
    const std::vector<double> oldvapor = soilmod.oldvapor;
    // Working values
    std::vector<double> psiw = soilmod.psiw;
    std::vector<double> theta = soilmod.theta;
    std::vector<double> vapor = soilmod.vapor;
    std::vector<double> k(n + 1, 0.0);
    std::vector<double> aa(n, 0.0), bb(n, 0.0), cc(n, 0.0), dd(n, 0.0);
    std::vector<double> ff(n, 0.0), Ca(n, 0.0), dpsi(n, 0.0);
    std::vector<double> kh(n + 1, 0.0), kv(n + 1, 0.0), dkh(n + 1, 0.0), dkv(n + 1, 0.0);
    std::vector<double> dFaceA(n, 0.0), dFaceB(n, 0.0);
    std::vector<double> psiFlow(n + 1, 0.0), ownFlow(n + 1, 1.0), passes(n + 1, 1.0), dpasses(n + 1, 0.0), room(n + 1, 1.0), droom(n + 1, 0.0);
    std::vector<double> obstructA(n, 0.0), obstructB(n, 0.0);
    // The last full Newton step, kept so that it can be halved if it increases the imbalance.
    std::vector<double> stepOrigin(n, 0.0), stepTaken(n, 0.0);
    double originResidual = 0.0, stepFraction = 1.0;
    bool haveStep = false;
    double surplusRate = 0.0;   // water a saturated surface cannot take (kg m-2 s-1)
    int iter = 1;
    double massBalance = 1.0;
    double Evapmmhr = 0.0;
    while (iter < maxNrIterations) {
        // Evaporation less rain at the surface.
        Evapmmhr = evaporation_flux(soilp, theta[0], soilmod.Tc[0], climdata.Tair,
            climdata.relhum, climdata.rHa, dT);
        double surfaceFlux = (Evapmmhr - climdata.precip) / dT;
        // Change in evaporation per unit surface water potential, included in the
        // surface layer's equation.
        double dEvap_dpsi = 0.0;
        {
            double hp = PSI_SLOPE_FRAC * std::abs(psiw[0]);
            if (hp < PSI_SLOPE_MIN) hp = PSI_SLOPE_MIN;
            double theta_probe = thetaFromPsi(soilp, psiw[0] + hp, 0);
            double Evprobe = evaporation_flux(soilp, theta_probe, soilmod.Tc[0], climdata.Tair,
                climdata.relhum, climdata.rHa, dT);
            dEvap_dpsi = ((Evprobe - Evapmmhr) / hp) / dT;
            if (!(dEvap_dpsi > 0.0)) dEvap_dpsi = 0.0;
        }
        // Root uptake from each layer. In a layer holding ice, roots see the potential
        // of the liquid water, which temperature sets.
        std::vector<double> psiRoot = psiw;
        for (size_t i = 0; i < psiRoot.size() && i < soilmod.Tc.size(); ++i)
            if (soilmod.Tc[i] < 0.0) psiRoot[i] = std::min(psiw[i], 334000.0 * soilmod.Tc[i] / 273.15);
        std::vector<double> STr = transpiration_distribute(soilp, soilmod.rootfrac,
            climdata.Et, dT, psiRoot, pTAW);

        // Hydraulic properties of each node
        for (int i = 0; i < n; ++i) {
            double Tkelvin = soilmod.Tc[i] + 273.15;
            // Only liquid water conducts. Ice impedes the flow (passes; Lundin 1990), and a
            // frozen node whose pores are full takes in no more water (room).
            double theta_u = unfrozenTheta(soilp, theta[i], soilmod.Tc[i], i);
            double theta_ice = theta[i] - theta_u;
            if (theta_ice < 0.0) theta_ice = 0.0;
            const bool holdsIce = (theta_ice > 0.0 && theta_u < theta[i]);
            const double Cw = dTheta_dPsi(soilp, psiw[i], i);
            const double airFilled = soilp.thetaS[i] - theta[i];
            kh[i] = hydraulicConductivityFromTheta(soilp, theta_u, i);
            passes[i] = 1.0; dpasses[i] = 0.0; room[i] = 1.0; droom[i] = 0.0;
            if (holdsIce) {
                constexpr double Omega = 7.0;
                passes[i] = std::pow(10.0, -Omega * theta_ice);
                dpasses[i] = -passes[i] * Omega * std::log(10.0) * Cw;
                room[i] = std::min(1.0, std::max(0.0, airFilled / FROZEN_PORE_CLOSING));
                droom[i] = (airFilled > 0.0 && airFilled < FROZEN_PORE_CLOSING) ? -Cw / FROZEN_PORE_CLOSING : 0.0;
            }
            // Potential driving flow: that of the liquid water, set by temperature where
            // ice is present.
            psiFlow[i] = holdsIce ? 334000.0 * soilmod.Tc[i] / 273.15 : psiw[i];
            ownFlow[i] = holdsIce ? 0.0 : 1.0;
            kv[i] = vaporConductivityFromPsiTheta(soilp, psiFlow[i], theta[i], Tkelvin, i);
            k[i] = kh[i] * passes[i] + kv[i];
            double Cv = dvapor_dPsi(soilp, psiw[i], theta[i], Tkelvin, i);
            Ca[i] = soilmod.vol[i] * (rho * Cw + Cv) / dT;
            // Change in liquid and vapour conductivity per unit water potential.
            dkh[i] = holdsIce ? 0.0 : -(2.0 + 3.0 / soilp.b[i]) * kh[i] / psiw[i];
            dkv[i] = (holdsIce ? 0.0 : kv[i] * Mw / (RgasC * Tkelvin)) - ((airFilled > 0.0) ? kv[i] / airFilled * Cw : 0.0);
        }
        // Lower boundary: saturated at a fixed potential, or free drainage under gravity.
        if (!soilp.FreeDrain) {
            psiw[n] = soilp.psie[n - 1];
            theta[n] = soilp.thetaS[n - 1];
            kh[n] = soilp.Ksat[n - 1];
            kv[n] = 0.0;
            psiFlow[n] = psiw[n]; passes[n] = 1.0; room[n] = 1.0;
        }
        else {
            psiw[n] = psiw[n - 1];
            theta[n] = theta[n - 1];
            kh[n] = kh[n - 1];
            kv[n] = kv[n - 1];
            psiFlow[n] = psiFlow[n - 1]; passes[n] = passes[n - 1]; room[n] = 1.0;
        }
        ownFlow[n] = 0.0; dpasses[n] = 0.0; droom[n] = 0.0;
        k[n] = kh[n] * passes[n] + kv[n];
        // Upward flux across the face below each layer. Liquid: Campbell's integration
        // of conductivity between the two nodes, plus gravity, reduced by ice impedance
        // and by the room in the receiving node. Vapour: diffusion along the potential
        // gradient.
        for (int i = 0; i < n; ++i) {
            double nPsi = 2.0 + 3.0 / soilp.b[i];
            double sh = soilmod.dz[i] * (1.0 - nPsi);
            double dpsiFace = psiFlow[i + 1] - psiFlow[i];
            double kvFace = 0.5 * (kv[i] + kv[i + 1]);
            const bool upperLimits = (passes[i] <= passes[i + 1]);
            const double liquid = (psiFlow[i + 1] * kh[i + 1] - psiFlow[i] * kh[i]) / sh - g * kh[i];
            const bool upperReceives = (liquid > 0.0);   // the flux is counted upward
            const double throughIce = upperLimits ? passes[i] : passes[i + 1];
            const double hasRoom = upperReceives ? room[i] : room[i + 1];
            const double pass = throughIce * hasRoom;
            // Change in the fraction passed per unit potential of the node above and of the node below.
            const double dpassA = (upperLimits ? dpasses[i] * hasRoom : 0.0) + (upperReceives ? throughIce * droom[i] : 0.0);
            const double dpassB = (upperLimits ? 0.0 : dpasses[i + 1] * hasRoom) + (upperReceives ? 0.0 : throughIce * droom[i + 1]);
            ff[i] = pass * liquid + kvFace * dpsiFace / soilmod.dz[i];
            // Change in flux through the ice impedance alone, kept separate (see below).
            obstructA[i] = dpassA * liquid; obstructB[i] = (i + 1 < n) ? dpassB * liquid : 0.0;
            // Change in flux per unit potential of the node above the face (A) and below it (B).
            if (i == n - 1 && soilp.FreeDrain) {
                dFaceA[i] = pass * (-g * dkh[i]);
                obstructA[i] = dpasses[i] * room[i] * liquid; obstructB[i] = 0.0;
                dFaceB[i] = 0.0;
                continue;
            }
            dFaceA[i] = pass * (-(kh[i] * ownFlow[i] + psiFlow[i] * dkh[i]) / sh - g * dkh[i])
                - kvFace * ownFlow[i] / soilmod.dz[i] + 0.5 * dkv[i] * dpsiFace / soilmod.dz[i];
            dFaceB[i] = (i + 1 < n) ? pass * ((kh[i + 1] * ownFlow[i + 1] + psiFlow[i + 1] * dkh[i + 1]) / sh)
                + kvFace * ownFlow[i + 1] / soilmod.dz[i] + 0.5 * dkv[i + 1] * dpsiFace / soilmod.dz[i] : 0.0;
        }
        // Water balance of each layer: flow in from below, flow out above, root uptake
        // and change in storage. What remains is the imbalance.
        massBalance = 0.0;
        aa[0] = 0.0;
        cc[0] = -dFaceB[0] - obstructB[0];
        // The impedance terms enter the Newton step only where they are stabilising.
        // This changes the path to the solution, not the solution.
        auto steadying = [](double v) { return (v > 0.0) ? v : 0.0; };
        bb[0] = -dFaceA[0] + Ca[0] + dEvap_dpsi + steadying(-obstructA[0]);
        dd[0] = surfaceFlux + STr[0] - ff[0]
            + soilmod.vol[0] * (rho * (theta[0] - oldtheta[0]) + (vapor[0] - oldvapor[0])) / dT;
        massBalance += std::abs(dd[0]);
        for (int i = 1; i < n; ++i) {
            aa[i] = dFaceA[i - 1] + obstructA[i - 1];
            cc[i] = -dFaceB[i] - obstructB[i];
            bb[i] = dFaceB[i - 1] - dFaceA[i] + Ca[i] + steadying(obstructB[i - 1] - obstructA[i]);
            dd[i] = ff[i - 1] + STr[i] - ff[i]
                + soilmod.vol[i] * (rho * (theta[i] - oldtheta[i]) + (vapor[i] - oldvapor[i])) / dT;
            massBalance += std::abs(dd[i]);
        }
        // A saturated surface takes no more water: its potential is held and the excess runs off.
        surplusRate = 0.0;
        if (psiw[0] >= soilp.psie[0] - 2e-8 && dd[0] < 0.0) {
            surplusRate = -dd[0];
            massBalance -= std::abs(dd[0]);
            aa[0] = 0.0;
            bb[0] = 1.0;
            cc[0] = 0.0;
            dd[0] = 0.0;
        }
        // Stop when the balance is within tolerance.
        if (massBalance <= tolerance) break;
        // If the last full step increased the imbalance, retake it at half the length
        // (down to 1/128).
        if (haveStep && massBalance > originResidual && stepFraction > 1.0 / 128.0) {
            stepFraction *= 0.5;
            for (int i = 0; i < n; ++i) {
                psiw[i] = stepOrigin[i] - stepFraction * stepTaken[i];
                psiw[i] = std::min(psiw[i], soilp.psie[i] - 1e-8);
                psiw[i] = std::max(psiw[i], soilp.psi_min[i]);
                theta[i] = thetaFromPsi(soilp, psiw[i], i);
                vapor[i] = vaporFromPsi(soilp, psiw[i], theta[i], soilmod.Tc[i] + 273.15, i);
            }
            ++iter;
            continue;
        }
        // The full Newton step, kept between the dry bound and air entry.
        Thomas TBC = ThomasBoundaryCondition(aa, bb, cc, dd, dpsi, 0, n - 1);
        dpsi = TBC.x;
        for (int i = 0; i < n; ++i) {
            stepOrigin[i] = psiw[i];
            stepTaken[i] = dpsi[i];
            psiw[i] = psiw[i] - dpsi[i];
            psiw[i] = std::min(psiw[i], soilp.psie[i] - 1e-8);
            psiw[i] = std::max(psiw[i], soilp.psi_min[i]);
            theta[i] = thetaFromPsi(soilp, psiw[i], i);
            vapor[i] = vaporFromPsi(soilp, psiw[i], theta[i], soilmod.Tc[i] + 273.15, i);
        }
        originResidual = massBalance;
        stepFraction = 1.0;
        haveStep = true;
        ++iter;
    }
    soilmod.psiw = psiw;
    soilmod.theta = theta;
    soilmod.vapor = vapor;
    soilmod.k = k;
    soilwaterout out;
    out.swo = soilmod;
    out.success = (massBalance <= tolerance);
    out.iterations = iter;
    out.Evapmmhr = Evapmmhr;
    out.surplus = surplusRate * dT;
    // Water crossing the base of the deepest layer, counted downward.
    out.drainage = -ff[n - 1] * dT;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ******************************************* Below-canopy Langrangian model ********************************************** //
// Air temperature and humidity within the canopy from the heat and vapour sources
// of each layer, by the Lagrangian dispersion theory of Raupach (1989).
static constexpr double KN_C1 = -0.39894;   // fixes the kernel's logarithmic behaviour at zero separation
static constexpr double KN_C2 = -0.15623;   // fixes its integral, and hence total scalar conservation
static constexpr double PI_SQ_OVER_6 = 1.6449340668482264;
// Near-field kernel of Raupach (1989). zeta is the distance between source and
// observation height, in units of sigma_w times the Lagrangian time scale.
static inline double nearFieldKernel(double zeta)
{
    double e = std::exp(-std::abs(zeta));
    return KN_C1 * std::log(1.0 - e) + KN_C2 * e;
}
// Dilogarithm on (0,1), needed to integrate the kernel.
static inline double dilog01(double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return PI_SQ_OVER_6;
    const bool reflect = (x > 0.5);
    const double y = reflect ? (1.0 - x) : x;
    double sum = 0.0, term = y;
    for (int k = 1; k <= 60; ++k) { sum += term / (static_cast<double>(k) * k); term *= y; }
    return reflect ? (PI_SQ_OVER_6 - std::log(x) * std::log(1.0 - x) - sum) : sum;
}
// Mean of the near-field kernel across a layer centred on the observation
// height, where the kernel itself is infinite. halfthick is half the layer
// thickness in kernel units.
static inline double nearFieldLayerMean(double halfthick)
{
    const double e = std::exp(-halfthick);
    const double integral = KN_C1 * (dilog01(e) - PI_SQ_OVER_6) + KN_C2 * (1.0 - e);
    return integral / halfthick;
}
// Corrections for sampling the kernel at layer midpoints, which underestimates
// the layers next to the observation height. NF_NEIGHBOUR applies at a layer
// centre (halved for the top node, which has a neighbour on one side only).
// NF_EDGE applies at canopy top, the edge of the top layer.
constexpr double LN_PI = 1.1447298858494002;
constexpr double NF_NEIGHBOUR = -KN_C1 * (LN_PI - 1.0);
constexpr double NF_EDGE = -KN_C1 * (1.5 * 0.6931471805599453 - 1.0);
// Solves the dense linear system A x = b (A row-major, n by n) by Gaussian
// elimination with partial pivoting.
static std::vector<double> solveDense(std::vector<double> A, std::vector<double> b, size_t n)
{
    for (size_t k = 0; k < n; ++k) {
        size_t p = k;
        for (size_t r = k + 1; r < n; ++r) if (std::abs(A[r * n + k]) > std::abs(A[p * n + k])) p = r;
        if (p != k) {
            for (size_t c = 0; c < n; ++c) std::swap(A[k * n + c], A[p * n + c]);
            std::swap(b[k], b[p]);
        }
        for (size_t r = k + 1; r < n; ++r) {
            const double f = A[r * n + k] / A[k * n + k];
            if (f == 0.0) continue;
            for (size_t c = k; c < n; ++c) A[r * n + c] -= f * A[k * n + c];
            b[r] -= f * b[k];
        }
    }
    std::vector<double> x(n);
    for (size_t k = n; k-- > 0;) {
        double sum = b[k];
        for (size_t c = k + 1; c < n; ++c) sum -= A[k * n + c] * x[c];
        x[k] = sum / A[k * n + k];
    }
    return x;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// How the ground surface temperature responds to canopy air temperature, so that
// the ground can be adjusted within the canopy air solve.
struct groundfollow {
    double share = 0.0;      // ground's change per degree change in the air it exchanges with
    double rToRef = 0.0;     // resistance from the soil surface to the reference height (s/m)
    double topShare = 0.0;   // canopy-top air's change per degree change in the ground
    double rhocp = 1.0;      // volumetric heat capacity of the reference air (J m-3 K-1)
    double moved = 0.0;      // set by the solve: how far the ground moved (K)
};
// Air temperature and humidity at each canopy node (Raupach 1989). A node's air
// is the canopy-top value plus a far-field part (diffusion of heat and vapour
// from all sources) and a near-field part (the extra effect of nearby foliage).
// Air, and the sources that depend on it, are solved together as a linear system.
static void LangrangianOne(onestep& onestepin, double pk, double tground, double soilrelhum,
    double th, double eh, const vegpstruct& vegp,
    const std::vector<double>& z, const windmodel& windvars,
    double a0, double a1, double rTop, const aerocolumn& col, const std::vector<double>& Rz,
    groundfollow* ground = nullptr)
{
    auto& tair = onestepin.tair;
    const auto& tleaf = onestepin.tleaf;
    auto& rh = onestepin.rh;
    const size_t nn = rh.size();
    const double nnd = static_cast<double>(nn);
    // Friction velocity with the still-air guard applied.
    const double uf_mix = mixingFrictionVelocity(windvars);
    const double dz = vegp.hgt / nnd;
    // Limits on the vapour pressure of canopy air: the driest and moistest of the
    // ground, the leaves and the air above.
    double esg = satvapCpp2(tground) * soilrelhum;  // ground vapour pressure
    double emn = std::min(eh, esg);
    double emx = std::max(eh, esg);
    for (size_t i = 0; i < nn; ++i) {
        double el = satvapCpp2(tleaf[i]);
        emn = std::min(emn, el);
        emx = std::max(emx, el);
    }
    // For each layer: gust strength sigma_w (ow), rising from a0*uf at the ground to
    // a1*uf at canopy top; one over sigma_w times the Lagrangian time scale (inowTL),
    // which scales distance in the near-field kernel; and the sensible and latent
    // heat sources per unit ground area (ST, SL).
    std::vector<double> ow(nn), inowTL(nn), ST(nn), SL(nn);
    const double mu1 = (a1 + a0) * 0.5 * uf_mix;
    const double mu2 = (a1 - a0) * 0.5 * uf_mix;
    for (size_t i = 0; i < nn; ++i) {
        ow[i] = mu1 + mu2 * std::cos(pi * (1.0 - z[i] / vegp.hgt));
        inowTL[i] = ow[i] / columnK(col, z[i]);
        ST[i] = vegp.paii[i] * onestepin.Hz[i];
        SL[i] = onestepin.Lz[i];
    }
    // Near-field weight of each layer's source at canopy top.
    std::vector<double> cTop(nn, 0.0);
    // Number of layers either side of a node that contribute to its near field.
    const size_t band = g_solver.nearLocal ? 1 : nn;
    for (size_t j = 0; j < (nn - 1); ++j) {
        if (nn - 1 - j > band) continue;
        // Direct path, and the path reflected at the ground.
        double zeta1 = (vegp.hgt - z[j]) * inowTL[j];
        double zeta2 = (vegp.hgt + z[j]) * inowTL[j];
        cTop[j] = (nearFieldKernel(zeta1) + nearFieldKernel(zeta2)) / ow[j];
    }
    // The top layer ends at canopy top, so the kernel's exact mean across it is used.
    {
        const size_t t = nn - 1;
        cTop[t] = (nearFieldLayerMean(dz * inowTL[t])
            + nearFieldKernel((vegp.hgt + z[t]) * inowTL[t])
            + NF_EDGE) / ow[t];
    }
    // Canopy-top temperature times the volumetric heat capacity of air.
    const double phh = phairCpp(th, pk);
    const double CfTh = phh * cpairCpp(th) * th;
    // W[i][j]: effect of layer j's source on the air at node i, relative to canopy-top air.
    // Far field: a source below the node adds the resistance from the node to canopy
    // top, times the fraction of the ground-to-node resistance lying below the
    // source. A source above adds the resistance from the source to canopy top.
    // Near field: the kernel at the distance between the two layer centres, direct
    // and reflected at the ground. A node's own layer enters as the kernel's exact mean.
    MP_LAP(9);
    std::vector<double> W(nn * nn), Wnear(nn * nn);
    // Powers of the kernel's exponential decay across one layer, reused for every pair.
    const size_t npow = 2 * nn;
    std::vector<double> decay(nn * npow);
    for (size_t j = 0; j < nn; ++j) {
        const double q = std::exp(-dz * inowTL[j]);
        double pw = 1.0;
        for (size_t k = 0; k < npow; ++k) { decay[j * npow + k] = pw; pw *= q; }
    }
    for (size_t i = 0; i < nn; ++i) {
        const double Rup = col.Rh - Rz[i];
        for (size_t j = 0; j < nn; ++j) {
            double nf = 0.0;
            if (i != j && ((i > j) ? i - j : j - i) > band) {}
            else if (i != j) {
                const size_t apart = (i > j) ? i - j : j - i;
                const double direct = decay[j * npow + apart], reflected = decay[j * npow + i + j + 1];
                nf = (KN_C1 * std::log((1.0 - direct) * (1.0 - reflected)) + KN_C2 * (direct + reflected)) / ow[j];
            }
            else {
                double halfthick = 0.5 * dz * inowTL[i];
                const double neighbours = (i + 1 == nn) ? 0.5 * NF_NEIGHBOUR : NF_NEIGHBOUR;
                nf = (nearFieldLayerMean(halfthick) + nearFieldKernel(2.0 * z[i] * inowTL[i])
                    + neighbours) / ow[i];
            }
            // A node within the soil's roughness has no resistance below it.
            const double below = (Rz[i] > 0.0) ? Rz[j] / Rz[i] : 1.0;
            const double ff = (j <= i) ? Rup * below : (col.Rh - Rz[j]);
            // The near field is scaled by the weight given to canopy-scale turbulence.
            Wnear[i * nn + j] = windvars.fade * (nf - cTop[j]);
            W[i * nn + j] = Wnear[i * nn + j] + ff;
        }
    }
    // Linear systems for air temperature (AT, bT) and vapour pressure (AL, bL) at all
    // nodes. Each includes the ground's exchange with the node, the change in each
    // layer's source with its own air, and the response of canopy-top air (across
    // rTop). Vapour is counted in moles.
    MP_LAP(10);
    std::vector<double> AT(nn * nn), AL(nn * nn), bT(nn), bL(nn);
    // Vapour pressure, molar density and heat capacity of each node's air, and the
    // latent heat at its leaf temperature.
    std::vector<double> eAir(nn), phAir(nn), cpAir(nn), laLeaf(nn);
    for (size_t i = 0; i < nn; ++i) {
        eAir[i] = satvapCpp2(tair[i]) * rh[i] / 100.0;
        phAir[i] = phairCpp(tair[i], pk);
        cpAir[i] = cpairCpp(tair[i]);
        laLeaf[i] = latentHeat(tleaf[i]);
    }
    for (size_t i = 0; i < nn; ++i) {
        const double Rup = col.Rh - Rz[i];
        const double ph = phAir[i];
        const double cp = cpAir[i];
        double sT = CfTh, sL = 0.0;
        for (size_t j = 0; j < nn; ++j) {
            const double w = W[i * nn + j];
            const double t0 = tair[j];
            const double e0 = eAir[j];
            const double wr = w + rTop;
            sT += w * ST[j] - wr * onestepin.dSTdT[j] * t0;
            sL += (w * SL[j] - wr * onestepin.dSLde[j] * e0) / laLeaf[j];
            AT[i * nn + j] = -Rz[i] * wr * onestepin.dSTdT[j];
            AL[i * nn + j] = -pk * Rz[i] * wr * onestepin.dSLde[j] / laLeaf[j];
        }
        AT[i * nn + i] += ph * cp * col.Rh;
        AL[i * nn + i] += phh * col.Rh;
        bT[i] = Rz[i] * sT + ph * cp * Rup * tground;
        bL[i] = Rz[i] * (phh * eh + pk * sL) + phh * Rup * esg;
    }
    // Ground adjusted with the air: a change in a layer's air changes its source,
    // which changes the air the ground exchanges with, and the ground follows by its share.
    std::vector<double> groundPerAir;
    if (ground != nullptr) {
        groundPerAir.assign(nn, 0.0);
        double atStart = 0.0;
        for (size_t j = 0; j < nn; ++j) {
            groundPerAir[j] = (ground->share / ground->rhocp) * (ground->rToRef - Rz[j]) * onestepin.dSTdT[j];
            atStart += groundPerAir[j] * tair[j];
        }
        for (size_t i = 0; i < nn; ++i) {
            const double reach = phAir[i] * cpAir[i] * ((col.Rh - Rz[i]) + Rz[i] * ground->topShare);
            for (size_t j = 0; j < nn; ++j) AT[i * nn + j] -= reach * groundPerAir[j];
            bT[i] -= reach * atStart;
        }
    }
    MP_LAP(11);
    const std::vector<double> tnew = solveDense(AT, bT, nn);
    if (ground != nullptr) {
        ground->moved = 0.0;
        for (size_t j = 0; j < nn; ++j) ground->moved += groundPerAir[j] * (tnew[j] - tair[j]);
    }
    const std::vector<double> enew = solveDense(AL, bL, nn);
    MP_LAP(12);
    // Far-field and near-field parts of each node's departure from canopy-top air.
    // Diagnostic output only.
    onestepin.tfar.assign(nn, 0.0); onestepin.tnear.assign(nn, 0.0);
    onestepin.efar.assign(nn, 0.0); onestepin.enear.assign(nn, 0.0);
    for (size_t i = 0; i < nn; ++i) {
        const double ph = phAir[i];
        const double cp = cpAir[i];
        double tn = 0.0, tf = 0.0, en = 0.0, ef = 0.0;
        for (size_t j = 0; j < nn; ++j) {
            const double sTj = ST[j] + onestepin.dSTdT[j] * (tnew[j] - tair[j]);
            const double sLj = (SL[j] + onestepin.dSLde[j] * (enew[j] - eAir[j])) / laLeaf[j];
            const double wn = Wnear[i * nn + j], wf = W[i * nn + j] - wn;
            tn += wn * sTj; tf += wf * sTj; en += wn * sLj; ef += wf * sLj;
        }
        onestepin.tnear[i] = Rz[i] * tn / (ph * cp * col.Rh);
        onestepin.tfar[i] = Rz[i] * tf / (ph * cp * col.Rh);
        onestepin.enear[i] = pk * Rz[i] * en / (phh * col.Rh);
        onestepin.efar[i] = pk * Rz[i] * ef / (phh * col.Rh);
    }
    // Update air temperature and relative humidity, keeping vapour pressure within its limits.
    for (size_t i = 0; i < nn; ++i) {
        tair[i] = tnew[i];
        const double ean = std::min(std::max(enew[i], emn), emx);
        rh[i] = (ean / satvapCpp2(tair[i])) * 100.0;
        if (rh[i] > 100.0) rh[i] = 100.0;
        if (rh[i] < 0.0) rh[i] = 0.0;
    }
}
// The air the ground exchanges heat and vapour with: the reference air, raised
// by each foliage layer's source times the resistance from that layer to the
// reference height. The ground exchanges with it across the whole resistance from
// the soil to the reference height. This is the same as exchange with the air
// just above the ground, with that air eliminated from the equations.
struct groundair {
    double T;        // temperature of the air the ground exchanges with (deg C)
    double relhum;   // its relative humidity (%)
    double R;        // resistance from the soil surface to it (s/m)
    double Ttop;     // the foliage's raise of the air at the ground above canopy-top air (K)
    double etop;     // the same for vapour pressure (kPa)
    double Hfoliage; // sensible heat released by all foliage (W/m2 of ground)
    double Lfoliage; // latent heat released by all foliage (W/m2 of ground)
};
static groundair airSeenByGround(const climstruct& climdata, const vegpstruct& vegp,
    const std::vector<double>& Hz, const std::vector<double>& Lz,
    const std::vector<double>& rz_zref, double rHa, double rhz)
{
    groundair out{};
    double Trise = 0.0, erise = 0.0, Tsrc = 0.0, esrc = 0.0;
    const double phr = phairCpp(climdata.tref, climdata.pk), cpr = cpairCpp(climdata.tref);
    const double lar = latentHeat(climdata.tref);
    for (size_t i = 0; i < rz_zref.size(); ++i) {
        const double SH = vegp.paii[i] * Hz[i], SL = Lz[i];
        out.Hfoliage += SH; out.Lfoliage += SL;
        Trise += SH * rz_zref[i]; erise += SL * rz_zref[i];
        Tsrc += SH * (rz_zref[i] - rhz); esrc += SL * (rz_zref[i] - rhz);
    }
    Trise /= (phr * cpr); erise *= climdata.pk / (phr * lar);
    Tsrc /= (phr * cpr); esrc *= climdata.pk / (phr * lar);
    out.T = climdata.tref + Trise;
    out.relhum = 100.0 * (satvapCpp2(climdata.tref) * climdata.relhum / 100.0 + erise) / satvapCpp2(out.T);
    out.R = rHa;
    out.Ttop = Tsrc;
    out.etop = esrc;
    return out;
}
// Temperature and vapour pressure of the air at canopy top, where the heat and
// vapour arriving from the ground and the foliage equal what leaves to the air
// above. The ground's exchange is with canopy-top air raised by the foliage
// sources (Tsrc, esrc).
static cantop canopytop(vegpstruct& vegpc, windmodel& wind, climstruct climdata,
    std::vector<double>& Hz, std::vector<double>& Lz, double zref,
    double Th, double eh, double tground, double soilrh,
    double rH_g, double rH_h_zref, double Tsrc, double esrc, int maxIter)
{
    // Tolerance on canopy-top temperature and vapour pressure.
    const double tolerance = g_solver.topTolerance;
    size_t nb = vegpc.paii.size();
    double FcH = 0.0;
    double FcL = 0.0;
    // Sensible and latent heat from all foliage, per unit ground area.
    for (size_t i = 0; i < nb; ++i) {
        FcH += Hz[i] * vegpc.paii[i];
        FcL += Lz[i];
    }
    const double eground = satvapCpp2(tground) * soilrh - esrc;
    tground -= Tsrc;
    const double eref = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    double err = 1e99;
    int nrIterations = 0;
    MP_TALLY(13);
    while (err > tolerance && nrIterations < maxIter) {
        MP_TALLY(14);
        double ph = phairCpp(Th, climdata.pk);
        double cp = cpairCpp(Th);
        const double tg = tground + Tsrc;
        const double la = latentHeat(tg);
        // Balance of canopy-top air between the reference air (across rH_h_zref), the
        // ground (across rH_g) and the foliage sources. Iterated only because air
        // density and heat capacity depend on temperature.
        const double rr = rH_h_zref / rH_g;
        double Th_new = (climdata.tref + FcH * rH_h_zref / (ph * cp) + rr * tground) / (1.0 + rr);
        double eh_new = (eref + FcL * rH_h_zref * climdata.pk / (ph * la) + rr * eground) / (1.0 + rr);
        err = std::abs(Th_new - Th);
        double err2 = std::abs(eh_new - eh);
        if (err2 > err) err = err2;
        Th = Th_new;
        eh = eh_new;
        ++nrIterations;
    }
    // Vapour pressure no higher than saturation.
    double rh = (eh / satvapCpp2(Th));
    if (rh > 1.0) eh = satvapCpp2(Th);
    cantop out;
    out.Th = Th;
    out.eh = eh;
    return out;
 }
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************** Run model for one step ********************************************** //
// One hour of the model. Wind, leaves, canopy air, ground and soil depend on one
// another and are iterated until they agree.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Values carried from one pass of the hour's iteration to the next.
struct outercarry {
    double tground = 0.0, Th = 0.0, eh = 0.0, psir = 0.0, H_iter = 0.0, G_iter = 0.0;
    double tdif = 1e99;        // how far the last pass moved the canopy air (K)
    double uf = 0.0;           // friction velocity the last pass found
    std::vector<double> dTs;   // leaf-to-air temperature difference seeding leaf convection
    int pass = 0;
};
// One pass of the hour: wind and stability, leaves, ground and soil, then canopy
// air, each solved with the others as they stand. Nothing is damped.
static void outerPass(onestep& onestepin, outercarry& c, const climstruct& climdata,
    vegpstruct& vegpc, const soilpstruct& soilpc, const std::vector<double>& z, const LWweights& wgts,
    const std::vector<double>& wc, const radmodel& swrad, const std::vector<double>& swater0,
    double Ca, double zref, int maxIter, double tolerance, double a0v, double a1, double dcan, bool C3)
{
    MP_LAP(14);
    const int innerIter = std::max(maxIter, INNER_ITERATIONS_LEAST);
    const size_t na = onestepin.tair.size();
    std::vector<double> rz_zref(na); // resistance from each layer to zref
    std::vector<double> rz_soil(na); // resistance from the soil to each layer
    std::vector<double> oldTe_fixed = onestepin.soilheatvars.oldTe; // soil temperatures at the start of the hour
    // Water on the leaves is reset to its value at the start of the hour.
    onestepin.swaterdepth = swater0;
    // Longwave radiation.
    radmodel2 lwrad = longwavemodelCpp(wgts, climdata.Rlw, c.tground, soilpc.groundem, vegpc.vegem, onestepin.tleaf);
    onestepin.Rlwdown = lwrad.Rlwdown;
    onestepin.Rlwup = lwrad.Rlwup;
    MP_LAP(2);
    // Wind and stability.
    windmodel wind = windmodelCpp(wc, climdata.uref, vegpc.hgt, vegpc.pai, zref, c.H_iter, climdata.tref, climdata.pk, innerIter, a1,
        onestepin.psih, onestepin.psim, onestepin.phih, 1000.0, 1.0, soilpc.z0);
    onestepin.uz = wind.uz;
    onestepin.ucanopy = wind.uh;
    onestepin.windscale = wind.windscale;
    onestepin.psih = wind.psi_h;
    onestepin.psim = wind.psi_m;
    onestepin.phih = wind.phi_h;
    onestepin.fade = wind.fade;
    onestepin.LL = wind.LL;
    c.uf = wind.uf;
    MP_LAP(3);
    // Rain interception.
    rainmodel rainvars = rainintercept(wc, vegpc.pia, wind.uh, climdata.precip, vegpc.x);
    envstruct envdata{};
    envdata.pk = climdata.pk;
    envdata.psi_r = c.psir;
    envdata.Ca = Ca;
    envdata.precip = climdata.precip;
    for (size_t i = 0; i < na; ++i) c.dTs[i] = std::abs(onestepin.tleaf[i] - onestepin.tair[i]);
    // Diffusivity profile, and resistances from the soil and from each layer to the
    // reference height.
    const aerocolumn col = makeColumn(vegpc.hgt, dcan, mixingFrictionVelocity(wind), wind.LL, a0v, a1,
        canopyEnhancement(vegpc.hgt, vegpc.pai, soilpc.z0, wind.fade), soilpc.z0);
    double rHa = soilToZ(col, zref); // soil surface to zref
    double rhg = col.Rh;             // soil surface to canopy top
    double rhz = rHa - rhg;          // canopy top to zref
    for (size_t i = 0; i < na; ++i) {
        rz_soil[i] = soilToZ(col, z[i]);
        rz_zref[i] = rHa - rz_soil[i];
    }
    MP_LAP(4);
    // Leaf energy balance.
    plantmodelCpp(onestepin, envdata, vegpc, rainvars, swrad, lwrad, z, c.dTs, 3600.0, C3); // updates onestepin in place
    MP_LAP(5);
    // Ground: the air it exchanges with, then soil heat and soil water.
    const groundair ga = airSeenByGround(climdata, vegpc, onestepin.Hz, onestepin.Lz, rz_zref, rHa, rhz);
    double Rabs = swrad.RswGabs + lwrad.RlwGabs;
    soilmod soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, ga.T, ga.relhum, climdata.pk, ga.R, 3600);
    MP_TALLY(7); MP_TALLYN(8, soilheat.iters);
    MP_LAP(6);
    climforwaterstruct cfw = {};
    cfw.Rabs = Rabs; cfw.Tair = ga.T; cfw.relhum = ga.relhum; cfw.pk = climdata.pk; cfw.rHa = ga.R;
    cfw.precip = onestepin.precipground; cfw.Et = onestepin.Et;
    onestepin.soilwatervars.oldTc = soilheat.oldTe;
    onestepin.soilwatervars.Tc = soilheat.Te;
    soilwaterout soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, vegpc.pTAW, innerIter, SOIL_WATER_TOL);
    MP_LAP(7);
    MP_TALLY(9); MP_TALLYN(10, soilwater.iterations);
    onestepin.witers = soilwater.iterations;
    onestepin.runoff = soilwater.surplus;
    onestepin.drainage = soilwater.drainage;
    c.tground = soilheat.Te[0];
    double soilrh = soilrelhumCpp(soilpc, c.tground, soilwater.swo.theta[0]);
    soilheat.wc = soilwater.swo.theta;
    double G_raw = soilsurfaceEB(soilpc, Rabs, ga.T, soilheat.Te[0], climdata.pk,
        ga.relhum, ga.R, soilwater.swo.theta[0]);
    c.G_iter = G_raw;
    soilheat.Gflux = c.G_iter;
    onestepin.soilheatvars = soilheat;
    onestepin.soilwatervars = soilwater.swo;
    onestepin.Ev = soilwater.Evapmmhr;
    onestepin.theta = soilwater.swo.theta[0];
    // Sensible and latent heat from the whole surface: foliage plus ground.
    {
        const double phr = phairCpp(climdata.tref, climdata.pk), cpr = cpairCpp(climdata.tref);
        onestepin.H = ga.Hfoliage + (phr * cpr / ga.R) * (c.tground - ga.T);
    }
    c.H_iter = onestepin.H;
    onestepin.L = ga.Lfoliage + groundLatentFlux(soilpc, c.tground, ga.T, ga.relhum,
        climdata.pk, ga.R, soilwater.swo.theta[0]);
    // Canopy-top air.
    if (zref > vegpc.hgt) {
        cantop Theh = canopytop(vegpc, wind, climdata, onestepin.Hz, onestepin.Lz, zref, c.Th, c.eh, c.tground, soilrh,
            rhg, rhz, ga.Ttop, ga.etop, innerIter);
        c.Th = Theh.Th;
        c.eh = Theh.eh;
    }
    else {
        c.Th = climdata.tref;
        c.eh = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    }
    MP_LAP(8);
    std::vector<double> tair = onestepin.tair;
    groundfollow follower;
    const bool groundFollows = (g_solver.groundRow > 0);
    std::vector<double> eBefore;
    if (groundFollows) {
        // How far the ground follows the air: its exchange with the air as a share of
        // everything restoring surface temperature (the slope of the surface energy
        // balance plus heat uptake by the soil over the hour).
        const double theta0 = soilwater.swo.theta[0];
        const double balanceSlope = (soilsurfaceEB(soilpc, Rabs, ga.T, c.tground - 0.05, climdata.pk, ga.relhum, ga.R, theta0)
            - soilsurfaceEB(soilpc, Rabs, ga.T, c.tground + 0.05, climdata.pk, ga.relhum, ga.R, theta0)) / 0.1;
        const double heatNow = soilHeatContent(soilpc, theta0, c.tground, 0);
        const double span = soilTemperatureFromHeat(soilpc, theta0, heatNow + 100.0, 0) - soilTemperatureFromHeat(soilpc, theta0, heatNow - 100.0, 0);
        const double capacity = (span > 0.0) ? 200.0 / span : 1e12;
        double ice = theta0 - unfrozenTheta(soilpc, theta0, c.tground, 0);
        if (ice < 0.0) ice = 0.0;
        const double conductivity = thermalConductivityFrom(soilConductivityFixedOf(soilpc.Vq[0], soilpc.Vm[0], soilpc.Vo[0], theta0, soilpc.Mc[0]),
            theta0, ice, c.tground, climdata.pk);
        const double soilUptake = std::sqrt(conductivity * capacity / 3600.0);
        follower.rhocp = phairCpp(climdata.tref, climdata.pk) * cpairCpp(climdata.tref);
        const double toAir = follower.rhocp / ga.R;
        follower.share = toAir / std::max(balanceSlope + soilUptake, toAir);
        follower.rToRef = rHa;
        follower.topShare = (zref > vegpc.hgt) ? (rhz / rhg) / (1.0 + rhz / rhg) : 0.0;
        eBefore.resize(na);
        for (size_t i = 0; i < na; ++i) eBefore[i] = satvapCpp2(onestepin.tair[i]) * onestepin.rh[i] / 100.0;
    }
    // Canopy air.
    LangrangianOne(onestepin, climdata.pk, c.tground, soilrh, c.Th, c.eh, vegpc, z, wind, a0v, a1,
        (zref > vegpc.hgt) ? rhz * rhg / (rhz + rhg) : 0.0, col, rz_soil, groundFollows ? &follower : nullptr); // updates onestepin.tair/rh in place
    if (groundFollows) {
        // Leaf temperatures and sources updated to the solved air.
        for (size_t i = 0; i < na; ++i) {
            if (!(vegpc.paii[i] > 0.0)) continue;
            const double dT = onestepin.tair[i] - tair[i];
            const double de = satvapCpp2(onestepin.tair[i]) * onestepin.rh[i] / 100.0 - eBefore[i];
            const double perLeaf = onestepin.dSTdT[i] / vegpc.paii[i];
            onestepin.Hz[i] += perLeaf * dT;
            onestepin.Lz[i] += onestepin.dSLde[i] * de;
            const double toLeaf = phairCpp(tair[i], climdata.pk) * cpairCpp(tair[i]) / (0.5 * onestepin.rLB[i]);
            onestepin.tleaf[i] += (1.0 + perLeaf / toLeaf) * dT;
        }
        // Ground updated to match: by re-solving soil heat (groundRow 1) or by the
        // change found in the canopy air solve (groundRow 2).
        const groundair gaNow = airSeenByGround(climdata, vegpc, onestepin.Hz, onestepin.Lz, rz_zref, rHa, rhz);
        const double theta0 = soilwater.swo.theta[0];
        if (g_solver.groundRow == 1) {
            onestepin.soilheatvars.oldTe = oldTe_fixed;
            soilmod again = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, gaNow.T, gaNow.relhum, climdata.pk, gaNow.R, 3600);
            MP_TALLY(7); MP_TALLYN(8, again.iters);
            again.wc = soilwater.swo.theta;
            c.tground = again.Te[0];
            onestepin.soilheatvars = again;
        }
        else {
            c.tground += follower.moved;
            onestepin.soilheatvars.Te[0] = c.tground;
        }
        c.G_iter = soilsurfaceEB(soilpc, Rabs, gaNow.T, c.tground, climdata.pk, gaNow.relhum, gaNow.R, theta0);
        onestepin.soilheatvars.Gflux = c.G_iter;
        onestepin.H = gaNow.Hfoliage + (follower.rhocp / gaNow.R) * (c.tground - gaNow.T);
        c.H_iter = onestepin.H;
        onestepin.L = gaNow.Lfoliage + groundLatentFlux(soilpc, c.tground, gaNow.T, gaNow.relhum, climdata.pk, gaNow.R, theta0);
    }
    MP_LAP(13);
    // Largest change in canopy air temperature over the pass.
    c.tdif = 0.0;
    for (size_t i = 0; i < na; ++i) {
        double dif = std::abs(tair[i] - onestepin.tair[i]);
        if (dif > c.tdif) c.tdif = dif;
    }
    c.psir = rootZonePotential(soilpc, soilwater.swo.theta, soilwater.swo.rootfrac, soilwater.swo.Tc) / 1000.0; // J/kg -> MPa
    onestepin.soilheatvars.oldTe = oldTe_fixed; // restore soil temperatures at the start of the hour
    ++c.pass;
    MP_LAP(8); MP_COUNT(15);
}
// The hour's state as one vector: air temperature, vapour pressure and leaf
// temperature at each node, then ground temperature, canopy-top temperature and
// vapour pressure, and the sensible and ground heat fluxes.
static std::vector<double> packState(const onestep& s, const outercarry& c)
{
    const size_t n = s.tair.size();
    std::vector<double> x(3 * n + 5);
    for (size_t i = 0; i < n; ++i) {
        x[i] = s.tair[i];
        x[n + i] = satvapCpp2(s.tair[i]) * s.rh[i] / 100.0;
        x[2 * n + i] = s.tleaf[i];
    }
    x[3 * n] = c.tground; x[3 * n + 1] = c.Th; x[3 * n + 2] = c.eh; x[3 * n + 3] = c.H_iter; x[3 * n + 4] = c.G_iter;
    return x;
}
// The reverse: writes the vector back to the state.
static void injectState(const std::vector<double>& x, onestep& s, outercarry& c)
{
    const size_t n = s.tair.size();
    for (size_t i = 0; i < n; ++i) {
        s.tair[i] = x[i];
        s.rh[i] = std::min(100.0, std::max(0.0, 100.0 * x[n + i] / satvapCpp2(x[i])));
        s.tleaf[i] = x[2 * n + i];
    }
    c.tground = x[3 * n]; c.Th = x[3 * n + 1]; c.eh = x[3 * n + 2]; c.H_iter = x[3 * n + 3]; c.G_iter = x[3 * n + 4];
}
// Number of passes remembered by the acceleration of the hour's iteration.
constexpr int OUTER_MEMORY = 10;
// Greatest distance (K or equivalent) allowed between the accelerated state and
// the latest pass's own result.
constexpr double OUTER_REACH = 10.0;
// Recent states (X) and how far a pass moved each (F), in comparable units.
struct outermemory {
    std::vector<std::vector<double>> X, F;
    std::vector<double> passStart, passResult;
    double lastMove = 1e99;
    int depth = OUTER_MEMORY; // zero gives plain passes with no acceleration
};
// One step of the hour's iteration: a pass, then Anderson (1965) acceleration to
// choose where the next pass starts: the combination of recent states that the
// passes together moved least. The memory is cleared if a pass moves the state
// three times further than the one before, or if the combination is not finite
// or lies too far from the latest result.
static void outerStep(onestep& s, outercarry& c, outermemory& mem, const climstruct& climdata,
    vegpstruct& vegpc, const soilpstruct& soilpc, const std::vector<double>& z, const LWweights& wgts,
    const std::vector<double>& wc, const radmodel& swrad, const std::vector<double>& swater0,
    double Ca, double zref, int maxIter, double tolerance, double a0v, double a1, double dcan, bool C3)
{
    const size_t n = s.tair.size(), m = 3 * n + 5;
    const std::vector<double> x = packState(s, c);
    outerPass(s, c, climdata, vegpc, soilpc, z, wgts, wc, swrad, swater0, Ca, zref, maxIter, tolerance, a0v, a1, dcan, C3);
    mem.passStart = x;
    mem.passResult = packState(s, c);
    if ((c.tdif <= tolerance && c.pass >= g_solver.minPass) || mem.depth <= 0) return;   // the pass's own result stands
    // Scale each variable so that one unit is comparable to 1 K.
    std::vector<double> sc(m, 1.0);
    for (size_t i = 0; i < n; ++i) sc[n + i] = 0.1;         // kPa
    sc[3 * n + 2] = 0.1; sc[3 * n + 3] = 50.0; sc[3 * n + 4] = 50.0;   // kPa, W m-2, W m-2
    std::vector<double> xs(m), f(m);
    double move = 0.0;
    for (size_t i = 0; i < m; ++i) {
        xs[i] = x[i] / sc[i];
        f[i] = (mem.passResult[i] - x[i]) / sc[i];
        move = std::max(move, std::abs(f[i]));
    }
    if (move > 3.0 * mem.lastMove) { mem.X.clear(); mem.F.clear(); }
    mem.lastMove = move;
    mem.X.push_back(xs); mem.F.push_back(f);
    if (mem.X.size() > static_cast<size_t>(mem.depth) + 1) { mem.X.erase(mem.X.begin()); mem.F.erase(mem.F.begin()); }
    const size_t k = mem.X.size() - 1;                      // differences between remembered passes
    if (k == 0) return;                                     // nothing to combine yet: the plain pass stands
    // Weights that make the combined move as small as possible, by least squares.
    std::vector<double> A(k * k, 0.0), b(k, 0.0);
    for (size_t p = 0; p < k; ++p) {
        for (size_t q = 0; q < k; ++q) {
            double acc = 0.0;
            for (size_t i = 0; i < m; ++i) acc += (mem.F[p + 1][i] - mem.F[p][i]) * (mem.F[q + 1][i] - mem.F[q][i]);
            A[p * k + q] = acc;
        }
        double acc = 0.0;
        for (size_t i = 0; i < m; ++i) acc += (mem.F[p + 1][i] - mem.F[p][i]) * f[i];
        b[p] = acc;
    }
    double tr = 0.0;
    for (size_t p = 0; p < k; ++p) tr += A[p * k + p];
    for (size_t p = 0; p < k; ++p) A[p * k + p] += 1e-10 * tr + 1e-300;   // keeps the fit defined when two moves coincide
    const std::vector<double> gam = solveDense(A, b, k);
    // The accelerated state.
    std::vector<double> xn(m);
    bool sound = true;
    for (size_t i = 0; i < m; ++i) {
        double v = xs[i] + f[i];
        for (size_t p = 0; p < k; ++p)
            v -= gam[p] * ((mem.X[p + 1][i] - mem.X[p][i]) + (mem.F[p + 1][i] - mem.F[p][i]));
        xn[i] = v * sc[i];
        if (!std::isfinite(v) || std::abs(v - (xs[i] + f[i])) > OUTER_REACH) sound = false;
    }
    if (sound) injectState(xn, s, c);
    else { mem.X.clear(); mem.F.clear(); }
}
// One hour for a vegetated surface. Shortwave radiation is calculated once;
// everything else is iterated until a pass changes canopy air temperature by
// less than the tolerance.
static onestep OneStepBelow(onestep onestepin, const obsstruct& obsdata, const climstruct& climdata, vegpstruct& vegpc, const soilpstruct& soilpc,
    const std::vector<double>& z, const tsvegstruct& tspveg, const tsvegstruct& tspvegPAR, const tsdifstruct& tspdif,
    const tsdifstruct& tspdifPAR, const LWweights& wgts, const std::vector<double>& wc,
    double Ca, double latr, double lonr, double zref, int maxIter = 100, double  tolerance = OUTER_TOLERANCE,
    double a0 = -1.0, double a1 = 1.25, bool C3 = true)
{
    MP_LAP(-1); MP_COUNT(16);
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double si = solarindexCpp2(soilpc.slope, soilpc.aspect, solp.zenr, solp.azir);
    kstruct kp = cankCpp(solp.zenr, vegpc.x, si);
    tsdirstruct tspdir = twostreamdirCpp(vegpc.pai, kp.kd, soilpc.gref, tspveg);
    tsdirstruct tspdirPAR = twostreamdirCpp(vegpc.pai, kp.kd, soilpc.grefPAR, tspvegPAR);
    radmodel swrad = shortwavemodelCpp(vegpc.pia, vegpc.paii, vegpc.pai, soilpc.gref, soilpc.grefPAR, vegpc.lref, vegpc.lrefp,
        climdata.Rsw, climdata.Rdif, si, solp, kp, tspveg, tspvegPAR, tspdif, tspdifPAR, tspdir, tspdirPAR);
    onestepin.Rdirdown = swrad.Rdirdown;
    onestepin.Rdifdown = swrad.Rdifdown;
    onestepin.Rswup = swrad.Rswup;
    size_t na = onestepin.tair.size();
    outercarry c;
    MP_LAP(1);
    c.tground = onestepin.soilheatvars.Te[0];
    // Canopy-top air from the previous hour.
    c.Th = onestepin.tcanopy;
    c.eh = satvapCpp2(c.Th) * (onestepin.rhcanopy / 100.0);
    // root-zone water potential, J/kg -> MPa
    c.psir = rootZonePotential(soilpc, onestepin.soilwatervars.theta, onestepin.soilwatervars.rootfrac, onestepin.soilwatervars.Tc) / 1000.0;
    c.dTs.assign(na, 0.0);
    // Sensible and ground heat fluxes from the previous hour.
    c.H_iter = onestepin.H;
    c.G_iter = onestepin.soilheatvars.Gflux;
    // Gust strength at the ground: as given, or from the foliage.
    const double a0v = (a0 > 0.0) ? a0 : floorGustRatio(vegpc, a1);
    const double dcan = zeroplanedisCpp2(vegpc.hgt, vegpc.pai);
    // Water on the foliage at the start of the hour, restored at the start of every pass.
    const std::vector<double> swater0 = onestepin.swaterdepth;
    outermemory mem;
    std::vector<double> steadiestStart;
    double steadiestMove = 1e99;
    while ((c.pass < g_solver.minPass || c.tdif > tolerance) && c.pass < maxIter) {
        outerStep(onestepin, c, mem, climdata, vegpc, soilpc, z, wgts, wc, swrad, swater0, Ca, zref, maxIter, tolerance,
            a0v, a1, dcan, C3);
        if (c.pass >= g_solver.minPass && c.tdif < steadiestMove) { steadiestMove = c.tdif; steadiestStart = mem.passStart; }
    }
    // An hour that has not converged returns the pass that changed canopy air least.
    if (c.tdif > tolerance && steadiestMove < c.tdif) {
        injectState(steadiestStart, onestepin, c);
        outerPass(onestepin, c, climdata, vegpc, soilpc, z, wgts, wc, swrad, swater0, Ca, zref, maxIter, tolerance, a0v, a1, dcan, C3);
        c.pass = maxIter;
    }
    MP_LAP(14);
    onestepin.iters = c.pass;
    onestepin.error = c.tdif;
    onestepin.tcanopy = c.Th;
    onestepin.rhcanopy = std::min(100.0, 100.0 * c.eh / satvapCpp2(c.Th));
    return onestepin;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************************** Above canopy ************************************************************** //
// Temperature, humidity and wind between canopy top and the reference height.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Fraction of the resistance between canopy top and the reference height that
// lies below height za. Temperature and vapour pressure at za lie this fraction
// of the way from their canopy-top values to their reference values.
static double aboveCanopyFraction(double za, double zref, double hgt, double pai, double LL, double z0, double fade)
{
    const double a1 = 1.25;
    const double d = zeroplanedisCpp2(hgt, pai);
    const double a0 = a1 * std::exp(-GUST_DECAY * LEAF_DRAG * pai / (1.0 + SHELTER * pai));
    const aerocolumn col = makeColumn(hgt, d, 1.0, LL, a0, a1, canopyEnhancement(hgt, pai, z0, fade), z0);
    return (soilToZ(col, za) - col.Rh) / (soilToZ(col, zref) - col.Rh);
}
// Air temperature at height za above the canopy, between canopy-top air (th)
// and the reference temperature (tref). LL defaults to neutral.
// [[Rcpp::export]]
double Tabove(double za, double zref, double th, double tref, double hgt, double pai, double LL = 1e99, double z0 = 0.004, double fade = 1.0)
{
    return th + (tref - th) * aboveCanopyFraction(za, zref, hgt, pai, LL, z0, fade);
}
// Relative humidity at height za above the canopy. Vapour pressure follows the
// same resistance as temperature, and is converted back at the temperature there.
// [[Rcpp::export]]
double RHabove(double za, double zref, double rh, double th, double tref, double tz, double relhum, double hgt, double pai, double LL = 1e99, double z0 = 0.004, double fade = 1.0)
{
    const double eh = satvapCpp2(th) * (rh / 100.0);
    const double eref = satvapCpp2(tref) * (relhum / 100.0);
    const double ez = eh + (eref - eh) * aboveCanopyFraction(za, zref, hgt, pai, LL, z0, fade);
    double rz = (ez / satvapCpp2(tz)) * 100.0;
    if (rz > 100.0) rz = 100.0;
    return rz;
}
// Wind speed at height za above the canopy, from the log profile with its
// stability correction.
// [[Rcpp::export]]
double Uabove(double za, double zref, double uh, double uref, double hgt, double pai, double LL) {
    double d = 0.0;
    if (hgt > 0.0) d = zeroplanedisCpp2(hgt, pai);
    double psitop = dpsimCpp2((hgt - d) / LL) - dpsimCpp2((za - d) / LL);
    double psibtm = dpsimCpp2((hgt - d) / LL) - dpsimCpp2((zref - d) / LL);
    double uz = uh + (uref - uh) * ((std::log((za - d) / (hgt - d)) + psitop)
        / (std::log((zref - d) / (hgt - d)) + psibtm));
    return uz;
}
// One hour for bare ground. Stability, soil heat and soil water are iterated
// until ground temperature stops changing. Air temperature, humidity and wind
// are then found at each height.
static onestepbare OneStepBare(onestepbare onestepin, const obsstruct& obsdata, const climstruct& climdata, const soilpstruct& soilpc,
    const std::vector<double>& z, double latr, double lonr, double zref, double zm = 0.004, int maxIter = 100,
    double  tolerance = OUTER_TOLERANCE)
{
    double Rswin = 0.0;   // shortwave arriving on the (possibly sloping) ground
    double Rb0 = 0.0;
    if (climdata.Rsw > 0.0) {
        solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
        if (solp.zenr < pi / 2.0) {
            Rb0 = (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr);
        }
        double si = solarindexCpp2(soilpc.slope, soilpc.aspect, solp.zenr, solp.azir);
        if (si < 0.0) si = 0.0;
        Rswin = climdata.Rdif + si * Rb0;
    }
    const double Rswabs = (1.0 - soilpc.gref) * Rswin;
    double Rlwabs = soilpc.groundem * climdata.Rlw;
    double Rabs = Rswabs + Rlwabs;
    double Tk = climdata.tref + 273.15;
    double cp = cpairCpp(climdata.tref);
    double ph = phairCpp(climdata.tref, climdata.pk);
    double cpph = cp * ph;
    double Ts = onestepin.soilheatvars.Te[0];
    double psi_m = onestepin.psim;
    double psi_h = onestepin.psih;
    double H = onestepin.H;
    // Wind driving exchange: with an upward heat flux, convection adds to the wind
    // (Beljaars 1994), as over vegetation.
    double zh = 0.2 * zm;
    auto driveWind = [&](double Hnow) {
        double Ueff = climdata.uref;
        if (Hnow > 0.0) {
            const double wstar = std::cbrt((g / Tk) * 1000.0 * (Hnow / cpph));
            Ueff = std::sqrt(climdata.uref * climdata.uref + wstar * wstar);
        }
        return Ueff;
    };
    double uf = (ka * driveWind(H)) / (std::log(zref / zm) + psi_m);
    // Obukhov length.
    double LL = (H != 0.0) ? (cpph * std::pow(guardedFrictionVelocity(uf, ka / std::log(zref / zm)), 3.0) * Tk) / (-ka * g * H) : 1e99;
    double dif = 1e99;
    int nrIterations = 0;
    soilmod soilheat;
    soilwaterout soilwater;
    // Resistance from the ground to the reference height, with the still-air guard
    // on friction velocity.
    double ufr = guardedFrictionVelocity(uf, ka / (std::log(zref / zm) + psi_m));
    zh = heatRoughness(zm, ufr);
    double rHa = (std::log(zref / zh) + psi_h) / (ka * ufr);
    // Ground heat flux, damped between passes.
    Aitken1DState st_G;
    double G_iter = onestepin.soilheatvars.Gflux;
    // At least two passes.
    while ((nrIterations < 2 || dif > tolerance) && nrIterations < maxIter) {
        if (H != 0.0) {
            LL = (cpph * std::pow(guardedFrictionVelocity(uf, ka / std::log(zref / zm)), 3.0) * Tk) / (-ka * g * H);
            LL = clipMOlength(LL, zref, 0.0, zm);
            psi_m = dpsimCpp2(zm / LL) - dpsimCpp2(zref / LL);
            psi_h = dpsihCpp2(zh / LL) - dpsihCpp2(zref / LL);
        }
        else {
            LL = 1e99;
            psi_m = 0.0;
            psi_h = 0.0;
        }
        uf = (ka * driveWind(H)) / (std::log(zref / zm) + psi_m);
        ufr = guardedFrictionVelocity(uf, ka / (std::log(zref / zm) + psi_m));
        zh = heatRoughness(zm, ufr);
        rHa = (std::log(zref / zh) + psi_h) / (ka * ufr);
        soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600);
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = climdata.precip; cfw.Et = 0.0;
        onestepin.soilwatervars.oldTc = soilheat.oldTe;
        onestepin.soilwatervars.Tc = soilheat.Te;
        // No transpiration from bare ground.
        soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, 0.5, std::max(maxIter, INNER_ITERATIONS_LEAST), SOIL_WATER_TOL);
        onestepin.runoff = soilwater.surplus;
        onestepin.drainage = soilwater.drainage;
        onestepin.witers = soilwater.iterations;
        // The next pass starts the water solve from this pass's result.
        onestepin.soilwatervars = soilwater.swo;
        onestepin.soilheatvars.wc = soilwater.swo.theta;
        double G_raw = soilsurfaceEB(soilpc, Rabs, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rHa, soilwater.swo.theta[0]);
        G_iter = aitken1d(G_iter, G_raw, st_G);
        soilheat.Gflux = G_iter;
        dif = std::abs(Ts - soilheat.Te[0]);
        Ts = soilheat.Te[0];
        H = (cpph / rHa) * (Ts - climdata.tref);
        ++nrIterations;
    }
    double soilrh = soilrelhumCpp(soilpc, Ts, soilwater.swo.theta[0]);
    double eg = satvapCpp2(Ts) * soilrh;
    double la = latentHeat(Ts);
    onestepin.L = groundLatentFlux(soilpc, Ts, climdata.tref, climdata.relhum, climdata.pk, rHa,
        soilwater.swo.theta[0]);
    // Air temperature, humidity and wind at each height.
    size_t n = z.size();
    std::vector<double> uz(n); // wind speed
    for (size_t i = 0; i < n; ++i) {
        if (z[i] > zm) {
            double psimz = dpsimCpp2(zm / LL) - dpsimCpp2(z[i] / LL);
            // Actual wind: the convective addition is removed.
            const double Udrive = driveWind(H);
            uz[i] = (uf / ka) * (std::log(z[i] / zm) + psimz) * ((Udrive > 0.0) ? climdata.uref / Udrive : 1.0);
        }
        else {
            uz[i] = 0.0;
        }
    }
    std::vector<double> tair(n); // air temperature
    std::vector<double> rh(n); // relative humidity
    for (size_t i = 0; i < n; ++i) {
        if (z[i] > zh) {
            double psihz = dpsihCpp2(zh / LL) - dpsihCpp2(z[i] / LL);
            // Resistance from the roughness height to z.
            double rHz = (std::log(z[i] / zh) + psihz) / (ka * ufr);
            tair[i] = Ts - (rHz / cpph) * H;
            double ez = eg - onestepin.L * (climdata.pk * rHz) / (la * ph);
            rh[i] = (ez / satvapCpp2(tair[i])) * 100.0;
        }
        else {
            tair[i] = Ts;
            rh[i] = (eg / satvapCpp2(Ts)) * 100.0;
        }
    }
    onestepin.Rb0 = Rb0;
    onestepin.Rswup = soilpc.gref * Rswin;   // reflected by the ground
    onestepin.Rlwup = soilpc.groundem * sb * radem(Ts);
    onestepin.uz = uz;
    onestepin.tair = tair;
    onestepin.rh = rh;
    onestepin.soilheatvars = soilheat;
    onestepin.soilwatervars = soilwater.swo;
    onestepin.H = H;
    onestepin.Ev = soilwater.Evapmmhr;
    onestepin.theta = soilwater.swo.theta[0];
    onestepin.psih = psi_h;
    onestepin.psim = psi_m;
    onestepin.LL = LL;
    onestepin.iters = nrIterations;
    onestepin.error = dif;
    return onestepin;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ******************* Bigleaf model currently only used for soil initialization ****************************************** //
// A simpler model that treats the canopy as a single layer. Used only to give
// the soil realistic starting temperatures and water contents.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Canopy albedo for diffuse (albd) and direct (albb) radiation, and the fractions
// of direct and diffuse radiation that reach the ground.
albf albedo(double kd, double pai, const tsdifstruct& tspdif, const tsdirstruct& tspdir)
{
    albf out;
    out.albd = tspdif.p1 + tspdif.p2; // white sky albedo
    out.albb = tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7; // black sky albedo
    out.wgtg = std::exp(-kd * pai); // ground weighting(direct)
    out.wgti = std::exp(-pai); // ground weighting(diffuse)
    return out;
}
// Stomatal resistance of the whole canopy: conductance of sunlit and shaded
// leaves, weighted by their leaf areas.
double bulkstomatalresist(double tair, double tcan, double relhum, double pk, double psiw,
    double Ca, double precip, double Rsw, double Rdif, double k, double LAIfrac,
    vegpstruct vegp, const solmodel& solp, bool C3 = true)
{
    double omp = vegp.lrefp + vegp.ltrap;
    double Lsun = (1.0 - std::exp(-k * vegp.pai * LAIfrac)) / k;
    double Lshade = (vegp.pai * LAIfrac) - Lsun;
    double Rshade_abs = Rdif * ((1.0 - std::exp(-vegp.pai)) / vegp.pai) * (1.0 - omp);
    double Rsun_abs = (Rsw - Rdif) * k * (1.0 - omp) + Rshade_abs;
    envstruct envdata;
    envdata.PARabs = Rsun_abs; // PAR absorbed
    envdata.tair = tair; // air temperature(deg C)
    envdata.tleaf = tcan; // canopy temperature(deg C)
    envdata.rh = relhum; // relative humidity(percentage)
    envdata.pk = pk; // atmospheric pressure(kPa)
    envdata.psi_r = psiw * 0.001; // mean water potential in root zone(MPa)
    envdata.Ca = Ca; // CO2 ppm
    envdata.precip = precip; // precipitation (mm)
    double z = vegp.hgt / 2.0;
    const double ph = phairCpp(tair, pk);
    double gs_sun = leafgs(envdata, vegp, z, C3);
    envdata.PARabs = Rshade_abs;
    double gs_shade = leafgs(envdata, vegp, z, C3);
    double Gs = gs_sun * Lsun + gs_shade * Lshade;
    // Convert to bulk surface resistance
    return ph / Gs;
}
// Dew-point temperature from air temperature and relative humidity.
double dewpointCpp2(double tc, double rh) {
    // actual vapour pressure (kPa)
    double ea = satvapCpp2(tc) * rh / 100.0;
    // Magnus constants
    const double a_w = 17.27;
    const double b_w = 237.7;   // water
    const double a_i = 21.875;
    const double b_i = 265.5;   // ice
    // first guess assuming water
    double gamma = std::log(ea / 0.6108);
    double td = (b_w * gamma) / (a_w - gamma);
    // correct for ice if dewpoint < 0
    if (td < 0.0) {
        double gamma_i = std::log(ea / 0.6108);
        td = (b_i * gamma_i) / (a_i - gamma_i);
    }
    return td;
}
// Damped update of an iterated value (Aitken): the weight given to the new value
// adapts to how the change between iterations is shrinking.
inline double aitken1d(double oldv, double newv, Aitken1DState& st)
{
    double r = newv - oldv;
    if (!st.have_prev) {
        st.r_prev = r;
        st.have_prev = true;
        return newv;  // no damping on first step
    }
    double dr = r - st.r_prev;
    if (dr != 0.0) {
        st.omega = -st.omega * st.r_prev / dr;
    }
    if (st.omega < 0.05) st.omega = 0.05;
    if (st.omega > 0.9)  st.omega = 0.9;
    st.r_prev = r;
    return oldv + st.omega * r;
}
// One hour of the big-leaf model for a vegetated surface. Canopy temperature,
// transpiration, soil heat and water and stability are iterated together.
bigleafone solveonestep(const obsstruct& obsdata, const climstruct& climdata, const vegpstruct& vegp,
    const tsvegstruct& tsvegp, const tsdifstruct& tspdif, const soilpstruct soilp,
    soilmod soilheat, soilwaterout soilwater,
    double latr, double lonr, double Ca = 430, double zref = 2.0,
    int maxiter = 100, bool C3 = true)
{
    if (zref < vegp.hgt) {
        Rcpp::stop("Height of climdate measurement must be above canopy");
    }
    double LAIfrac = 0.0;
    for (size_t i = 0; i < vegp.Lfrac.size(); ++i) {
        LAIfrac += vegp.Lfrac[i];
    }
    LAIfrac = LAIfrac / static_cast<double>(vegp.Lfrac.size());
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double si = solarindexCpp2(soilp.slope, soilp.aspect, solp.zenr, solp.azir);
    kstruct kk = cankCpp(solp.zenr, vegp.x, si);
    tsdirstruct tspdir = twostreamdirCpp(vegp.pai, kk.kd, soilp.gref, tsvegp);
    albf albs = albedo(kk.kd, vegp.pai, tspdif, tspdir);
    double Rabs_sw = 0.0;
    double RabsG_sw = 0.0;
    double amx = soilp.gref; if (amx < vegp.lref) amx = vegp.lref;
    if (climdata.Rsw > 0.0) {
        // beam normal to the sun; none with the sun at or below the horizon
        double Rb0 = (std::cos(solp.zenr) > 0.0) ? (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr) : 0.0;
        double Rdirdowng = 0.0;
        double Rdbdg = 0.0;
        double kksi = 0.0;
        if (Rb0 > 0.0) {
            double albc = (albs.albb - albs.wgtg * soilp.gref) / (1.0 - albs.wgtg); // canopy albedo
            kksi = kk.k * si;
            if (kksi > 1.0) kksi = 1.0;
            Rabs_sw = (1.0 - albs.albd) * climdata.Rdif +
                (1.0 - soilp.gref) * albs.wgtg * Rb0 * si +
                (1.0 - albc) * (1.0 - albs.wgtg) * Rb0 * kksi;
            Rdirdowng = Rb0 * std::exp(-kk.kd * vegp.pai);
            Rdbdg = (tspdir.p8 / tspdir.sig) * std::exp(-kk.kd * vegp.pai) +
                tspdir.p9 * std::exp(-tsvegp.h * vegp.pai) +
                tspdir.p10 * std::exp(tsvegp.h * vegp.pai);
            if (Rdbdg > amx) Rdbdg = amx;
            if (Rdbdg < 0.0) Rdbdg = 0.0;
        }
        else {
            Rabs_sw = (1.0 - albs.albd) * climdata.Rdif;
        }
        double Rdddg = tspdif.p3 * std::exp(-tsvegp.h * vegp.pai) + tspdif.p4 * std::exp(tsvegp.h * vegp.pai);
        if (Rdddg > 1.0) Rdddg = 1.0;
        if (Rdddg < 0.0) Rdddg = 0.0;
        // Diffuse radiation reaching the ground: from the sky, and scattered from the beam.
        double Rdifdowng = Rdddg * climdata.Rdif + Rdbdg * Rb0 * si;
        RabsG_sw = (1.0 - soilp.gref) * (Rdifdowng + Rdirdowng * si);
    }
    double Rabs_lw = climdata.Rlw * (albs.wgti * soilp.groundem + (1.0 - albs.wgti) * vegp.vegem);
    double Rabs = Rabs_sw + Rabs_lw; // total radiation absorbed by canopy
    double d = zeroplanedisCpp2(vegp.hgt, vegp.pai);
    double cp = cpairCpp(climdata.tref);
    double ph = phairCpp(climdata.tref, climdata.pk);
    double tr = std::exp(-vegp.pai);
    double em = tr * soilp.groundem + (1.0 - tr) * vegp.vegem;
    double Rema = vegp.vegem * sb * radem(climdata.tref);
    double Da = satvapCpp2(climdata.tref) * (1.0 - climdata.relhum / 100.0);
    double Tk = climdata.tref + 273.15;
    double tdew = dewpointCpp2(climdata.tref, climdata.relhum);
    double ea = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    // Starting values
    double psim = 0.0;
    double H = 0.01;
    double G = 0.0;
    double tcanopy = climdata.tref + 0.01 * Rabs;
    double error = 1e99;
    double Tr = 0.0;
    int itr = 0.0;
    double uf = 0.0;
    double LL = 1e10;
    Aitken1DState st;
    // At least two passes.
    while ((itr < 2 || error > OUTER_TOLERANCE) && itr < maxiter) {
        const double fade = canopyTurbulence(climdata.uref, vegp.hgt, vegp.pai, d, zref, soilp.z0);
        double zm = roughlengthCpp2(vegp.hgt, vegp.pai, d, soilp.z0, fade);
        uf = (ka * climdata.uref) / (std::log((zref - d) / zm) + psim);
        // Resistances to the reference height: for the ground, across the whole
        // diffusivity profile (rGz); for the big leaf, from its exchange height to canopy
        // top and then through the profile above (rHa).
        const double a1s = 1.25;
        const aerocolumn col = makeColumn(vegp.hgt, d, uf, LL, floorGustRatio(vegp, a1s), a1s,
            canopyEnhancement(vegp.hgt, vegp.pai, soilp.z0, fade), soilp.z0);
        const double rGz = soilToZ(col, zref);
        const double rHa = bigLeafToCanopyTop(vegp.hgt, vegp.pai, d, uf, LL, soilp.z0) + rGz - col.Rh;
        // Stomatal resistance of the canopy, from the water potential of the root zone.
        const double psiw = rootZonePotential(soilp, soilwater.swo.theta, soilwater.swo.rootfrac, soilwater.swo.Tc);
        double rS = bulkstomatalresist(climdata.tref, tcanopy, climdata.relhum, climdata.pk, psiw,
            Ca, climdata.precip, climdata.Rsw, climdata.Rdif, kk.k, LAIfrac, vegp, solp, C3);
        rS = std::min(rS, 1e10);
        double rV = rS + rHa;
        double RabsG_lw = (tr * climdata.Rlw + (1.0 - tr) * sb * radem(tcanopy)) * soilp.groundem;
        double RabsG = RabsG_sw + RabsG_lw;
        double Tkc = tcanopy + 273.15;
        double es = satvapCpp2(tcanopy);
        Tr = Mw / (RgasC * Tkc * rV) * (es - ea) * 1000.0 * 3600.0 * vegp.pai * LAIfrac; // transpiration
        if (Tr < 0.0) Tr = 0.0;
        soilheat.wc = soilwater.swo.theta;
        soilheat = SoilHeatCpp(soilheat, soilp, RabsG, climdata.tref, climdata.relhum, climdata.pk, rGz, 3600);
        climforwaterstruct cfw = {};
        cfw.Rabs = RabsG; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rGz;
        cfw.precip = climdata.precip; cfw.Et = Tr;
        soilwater.swo.oldTc = soilheat.oldTe;
        soilwater.swo.Tc = soilheat.Te;
        soilwater = SoilWaterCpp(soilwater.swo, soilp, cfw, 3600, vegp.pTAW, maxiter, SOIL_WATER_TOL);
        double G_raw = soilsurfaceEB(soilp, RabsG, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rGz, soilwater.swo.theta[0]); // ground heat flux
        G = aitken1d(G, G_raw, st);
        soilheat.Gflux = G;
        // Canopy temperature from the Penman-Monteith energy balance.
        double Te = (tcanopy + climdata.tref) / 2.0;
        double Rer = 4.0 * em * sb * std::pow(Te + 273.15, 3.0);
        double la = latentHeat(tcanopy);
        double De = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5);
        double oldtcanopy = tcanopy;
        tcanopy = climdata.tref + ((Rabs - Rema - ((la * ph) / (climdata.pk * rV)) * Da - G) /
            (Rer + ph * (cp / rHa + ((la * De) / (climdata.pk * rV)))));
        if (tcanopy < tdew) tcanopy = tdew;
        // Sensible heat flux and stability.
        H = (ph * cp / rHa) * (tcanopy - climdata.tref);
        LL = (ph * cp * std::pow(guardedFrictionVelocity(uf, ka / std::log((zref - d) / zm)), 3.0) * Tk) / (-ka * g * H);
        LL = clipMOlength(LL, zref, d, zm);
        psim = dpsimCpp2(zm / LL) - dpsimCpp2((zref - d) / LL);
        error = std::abs(tcanopy - oldtcanopy);
        ++itr;
    }
    bigleafone out;
    out.soilheat = soilheat;
    out.soilwater = soilwater;
    out.tcanopy = tcanopy;
    out.uf = uf;
    out.LL = LL;
    out.Et = Tr;
    out.iters = itr;
    return out;
}
// One hour of the big-leaf model for bare ground.
bigleafone solveonestepbare(const obsstruct& obsdata, const climstruct& climdata, const soilpstruct soilp,
    soilmod soilheat, soilwaterout soilwater, double zmr, double latr, double lonr, double zref = 2.0,
    int maxiter = 20)
{
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double si = solarindexCpp2(soilp.slope, soilp.aspect, solp.zenr, solp.azir);
    double Rabs_sw = 0.0;
    if (climdata.Rsw > 0.0) {
        // beam normal to the sun; none with the sun at or below the horizon
        double Rb0 = (std::cos(solp.zenr) > 0.0) ? (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr) : 0.0;
        Rabs_sw = (1.0 - soilp.gref) * (climdata.Rdif + si * Rb0);
    }
    double Rabs = Rabs_sw + soilp.groundem * climdata.Rlw;
    double cp = cpairCpp(climdata.tref);
    double ph = phairCpp(climdata.tref, climdata.pk);
    double Tk = climdata.tref + 273.15;
    double psih = 0.0;
    double psim = 0.0;
    double H = 0.01;
    double error = 1e99;
    int itr = 0.0;
    double uf = 0.0;
    double LL = 1e10;
    double oldtsoil = soilheat.Te[0];
    double newtsoil = soilheat.Te[0];
    Aitken1DState st;
    // Ground heat flux, damped between passes.
    Aitken1DState st_G;
    double G = soilheat.Gflux;
    // At least two passes.
    while ((itr < 2 || error > OUTER_TOLERANCE) && itr < maxiter) {
        // Friction velocity and resistance, with convection added to the wind when the
        // heat flux is upward.
        const double zm = zmr;
        double Ueff = climdata.uref;
        if (H > 0.0) {
            const double wstar = std::cbrt((g / Tk) * 1000.0 * (H / (ph * cp)));
            Ueff = std::sqrt(climdata.uref * climdata.uref + wstar * wstar);
        }
        uf = (ka * Ueff) / (std::log(zref / zm) + psim);
        const double ufr = guardedFrictionVelocity(uf, ka / (std::log(zref / zm) + psim));
        double zh = heatRoughness(zm, ufr);
        double rHa = (std::log(zref / zh) + psih) / (ka * ufr);
        soilheat.wc = soilwater.swo.theta;
        soilheat = SoilHeatCpp(soilheat, soilp, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600);
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = climdata.precip; cfw.Et = 0.0;
        soilwater.swo.oldTc = soilheat.oldTe;
        soilwater.swo.Tc = soilheat.Te;
        // No transpiration from bare ground.
        soilwater = SoilWaterCpp(soilwater.swo, soilp, cfw, 3600, 0.5, maxiter, SOIL_WATER_TOL);
        double G_raw = soilsurfaceEB(soilp, Rabs, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rHa, soilwater.swo.theta[0]);
        G = aitken1d(G, G_raw, st_G);
        soilheat.Gflux = G;
        // Change in surface temperature over the pass, before damping.
        error = std::abs(soilheat.Te[0] - oldtsoil);
        newtsoil = aitken1d(newtsoil, soilheat.Te[0], st);
        soilheat.Te[0] = newtsoil;
        H = (ph * cp / rHa) * (newtsoil - climdata.tref);
        LL = (ph * cp * std::pow(guardedFrictionVelocity(uf, ka / std::log(zref / zm)), 3.0) * Tk) / (-ka * g * H);
        LL = clipMOlength(LL, zref, 0.0, zm);
        psim = dpsimCpp2(zm / LL) - dpsimCpp2(zref / LL);
        psih = dpsihCpp2(zh / LL) - dpsihCpp2(zref / LL);
        oldtsoil = newtsoil;
        ++itr;
    }
    bigleafone out;
    out.soilheat = soilheat;
    out.soilwater = soilwater;
    out.tcanopy = newtsoil;
    out.uf = uf;
    out.LL = LL;
    out.Et = 0.0;
    out.iters = itr;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ****************************************** Ecotherm model *********************************************************** //
// Body temperature of an ectotherm from its energy balance, given the
// microclimate around it. The animal is treated as an ellipsoid.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Surface area and volume of the animal, and the area it presents to the direct
// beam (silA), for its shape and orientation.
static silstruct silhouette(double zenr, double azir, double height, double width, double length,
    double adir = 0.0, double atilt = 0.0, std::string position = "fixed")
{
    // Compute semi-axis (converted from cm to m)
    double A = length / 200.0;
    double B = width / 200.0;
    double C = height / 200.0;
    // Convert directions to radians
    double adirr = adir * torad;
    double atiltr = atilt * torad;
    silstruct out;
    // Surface area by Knud Thomsen's formula for an ellipsoid
    double P = 1.6075;  // Approximation constant
    out.A = 4 * pi * std::pow((std::pow(A, P) * std::pow(B, P) +
        std::pow(A, P) * std::pow(C, P) +
        std::pow(B, P) * std::pow(C, P)) / 3, 1.0 / P);
    out.V = (4.0 / 3.0) * pi * A * B * C;
    if (position == "max") { // assume solar radiation maximised
        out.silA = pi * std::max({ B * C, A * C, A * B });
    }
    else if (position == "min") { // assume solar radiation minimised
        out.silA = pi * std::min({ B * C, A * C, A * B });
    }
    else if (position == "randomdir") { // direction (horizontal rotation) assumed random
        double phi = std::acos(std::cos(zenr) * std::cos(-atiltr) +
            std::sin(zenr) * std::sin(-atiltr) * std::cos(azir));
        double M = std::sin(phi);
        double N = std::cos(phi);
        out.silA = pi * std::sqrt(
            0.5 * M * M * (B * B * C * C + C * C * A * A) +
            N * N * A * A * B * B);
    }
    else if (position == "random") { // direction and tilt assumed random
        out.silA = pi * std::sqrt((A * A * B * B + A * A * C * C + B * B * C * C) / 3.0);
    }
    else if (position == "fixed") {
        // Calculate silhouette area
        double theta = (adirr - azir);
        double phi = std::acos(std::cos(zenr) * std::cos(-atiltr) +
            std::sin(zenr) * std::sin(-atiltr) * cos(theta));
        double L = cos(theta) * sin(phi);
        double M = sin(theta) * sin(phi);
        double N = cos(phi);
        out.silA = pi * std::sqrt(L * L * B * B * C * C +
            M * M * C * C * A * A +
            N * N * A * A * B * B);
    }
    else Rcpp::stop("Position not recongised");
    return out;
}
// Characteristic dimension of the animal across the wind, for convective heat exchange.
static double chardim(double wdir, double zenr, double azir, double height, double width, double length,
    double adir = 0.0, double atilt = 0.0, std::string position = "fixed")
{
    // Compute projected area in direction of wind flow (horizontal)
    if (position == "max") {
        if (length <= width && length <= height) {
            // Beam along A-axis
            adir = azir * 180.0 / pi;
            atilt = (pi / 2.0 - zenr) * 180.0 / pi;
        }
        else if (width <= length && width <= height) {
            // Beam along B-axis
            adir = (azir + pi / 2.0) * 180.0 / pi;
            atilt = 90.0;
        }
        else {
            // Beam along C-axis
            adir = azir * 180.0 / pi;
            atilt = -zenr * 180 / pi;
        }
        position = "fixed";
    }
    if (position == "min") {
        if (length >= width && length >= height) {
            // Beam along A-axis
            adir = azir * 180.0 / pi;
            atilt = (pi / 2.0 - zenr) * 180.0 / pi;
        }
        else if (width >= length && width >= height) {
            // Beam along B-axis
            adir = (azir + pi / 2.0) * 180.0 / pi;
            atilt = 90.0;
        }
        else {
            // Beam along C-axis
            adir = azir * 180.0 / pi;
            atilt = -zenr * 180.0 / pi;
        }
        position = "fixed";
    }
    silstruct sa = silhouette(pi / 2.0, wdir, height, width, length, adir, atilt, position);
    // Volume over the area presented to the wind
    double d = sa.V / sa.silA;
    return d;
}
// Boundary-layer resistance of the animal to heat (s/m): forced and free
// convection combined.
static double animalrHa(double tair, double dT, double uz, double d, double rHmax = 300.0)
{
    double Tk = tair + 273.15;
    double Kh = (1.6667e-10 * Tk * Tk + 2.9935e-8 * Tk - 1.7128e-6); // thermal diffusivity
    double v = 1.326 * std::pow(10.0, -5.0) * std::pow(Tk / 273.15, 1.5) * (393.55 / (Tk + 120.0)); // kinematic viscosity
    double Re = (uz * d) / v; // Reynolds number
    double Pr = v / Kh; // Prandtl number
    double Nuf; // Nusselt number, forced convection
    if (Re > 2e5) {
        Nuf = 0.37 * std::pow(Re, 0.6) * std::pow(Pr, 1.0 / 3.0) + 9.08;
    }
    else if (Re > 1000) {
        Nuf = 2.0 + (0.48 * std::pow(Re, 0.6) - 11.31) * std::pow(Pr, 1.0 / 3.0);
    }
    else {
        Nuf = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0 / 3.0);
    }
    double Gr = (g * std::pow(d, 3.0) * dT) / (Tk * v * v); // Grashof number
    double Nun = std::pow(0.825 + (0.387 * std::pow(Gr * Pr, 1.0 / 6.0)) /
        std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0), 2.0); // Nusselt number, free convection
    double Nu = std::pow(std::pow(Nuf, 3.0) + std::pow(Nun, 3.0), 1.0 / 3.0);
    double rHa = d / (Kh * Nu);
    if (rHa > rHmax) rHa = rHmax;
    return rHa;
}
// Diffusivity of water vapour in air (m2/s).
double Dw_waterVapour(double Tf, double pk)
{
    const double D0 = 2.26e-5;   // m^2 s^-1 at 273.15 K and 101325 Pa
    const double T0 = 273.15;    // reference temperature (K)
    const double P0 = 101325.0;  // reference pressure (Pa)
    double Tfk = Tf + 273.15;
    double Pa = pk * 1000.0;
    return D0 * std::pow(Tfk / T0, 1.75) * (P0 / Pa);
}
// Body temperature from the Penman-Monteith energy balance: absorbed radiation
// and metabolic heat (M) against emission, convection and evaporation to the air,
// and conduction and evaporation through the fraction of the body in contact
// with the surface (confrac).
// [[Rcpp::export]]
double PenmanMonteith_animal(double Rabs, double Ta, double Ts, double Te, double Tf, double pk, double rh,
    double rHa, double height, double rc, double confrac, double M, double em = 0.97, double k = 0.5,
    double surfrh = 1.0, double surfwet = 1.0)
{
    double Tbody = Ta;
    double cp = cpairCpp(Ta);
    double ph = phairCpp(Ta, pk);
    double dT = Ts - Ta; // contact surface minus air temperature
    double cd = 1.0 - confrac;
    double ea = satvapCpp2(Ta);
    double Da = ea * (1.0 - rh / 100.0); // vapour pressure deficit, air
    double Dac = ea - surfrh * satvapCpp2(Ts); // vapour pressure deficit against the contact surface
    // Conduction from the body core to the contact surface, across half the body
    // height (given in cm).
    const double halfheight = 0.5 * height / 100.0; // m
    double keff = k / halfheight;
    double mr = cd * em * sb; // radiative exchange coefficient, air-facing
    double mrc = confrac * em * sb; // radiative exchange coefficient, contact-facing
    double mc = keff * confrac; // conductive exchange coefficient, contact-facing
    double mv = ph * cp * cd / rHa; // sensible heat exchange coefficient, air-facing
    double la = (Tbody >= 0.0) ? 45068.7 - 42.8428 * Tbody : 51078.69 - 4.338 * Tbody - 0.06367 * Tbody * Tbody;
    double rv = rHa + rc; // total vapour resistance (s/m)
    double ml = (ph * la * cd) / (rv * pk); // latent heat exchange coefficient, air-facing
    double Dw = Dw_waterVapour(Tf, pk);
    // Vapour passes between the skin and the wet share of the contact surface, across the same half height.
    double rv_contact = halfheight / Dw + rc;
    double mlc = (ph * la * confrac * surfwet) / (rv_contact * pk); // latent heat exchange coefficient, contact-facing
    double DeV = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5); // slope of the saturation vapour pressure curve, air-facing
    double top = (Rabs + M - mr * radem(Ta) - mrc * radem(Ta) + mc * dT - ml * Da - mlc * Dac);
    double btm = 4.0 * mr * std::pow(Te + 273.15, 3.0) + 4.0 * mrc * std::pow(Tf + 273.15, 3.0) +
        mc + mv + ml * DeV + mlc * DeV;
    Tbody = Ta + top / btm;
    return (Tbody);
}
// Damped update of an iterated value (Aitken), with limits on the weight given
// to the old value.
static inline void aitkin_weightdif_scalar(
    double oldv,
    double& newv,   // raw new value on input, updated in place
    WAitkenStateScalar& st,
    double backweight_min = 0.10,
    double backweight_max = 0.98
)
{
    double r = newv - oldv;
    double omega_max = 1.0 - backweight_min;
    double omega_min = 1.0 - backweight_max;
    // ---- First iteration ----
    if (!st.have_prev) {
        double omega = std::max(omega_min, std::min(st.omega, omega_max));
        st.r_prev = r;
        newv = oldv + omega * r;
        st.omega = omega;
        st.have_prev = true;
        return;
    }
    // ---- Learn omega (scalar Aitken) ----
    double dr = r - st.r_prev;
    double omega = st.omega;

    if (dr != 0.0) {
        omega = -st.omega * (st.r_prev / dr);
    }
    // Clamp through backweight limits
    omega = std::max(omega_min, std::min(omega, omega_max));
    // ---- Apply update ----
    st.r_prev = r;
    newv = oldv + omega * r;
    st.omega = omega;
}
// Metabolic rate (W): allometric scaling with body mass and a Q10 response to
// body temperature.
double metabolic_rate(double volume, double rho, double Q10, double a0, double b, double Tref, double Tbody)
{
    double m = volume * rho;
    double M = a0 * std::pow(m, b) * std::pow(Q10, (Tbody - Tref) / 10.0);
    return M;
}
// Body temperature at each time step for one microclimate series.
// position is one of: fixed (adir and atilt as given), max (radiation absorbed is
// maximised), min (minimised), randomdir (tilt as given, direction random),
// random (tilt and direction random).
// [[Rcpp::export]]
Rcpp::NumericVector Ectotherm(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List animal,
    double lat, // Latitude (decimal degrees)
    double lon, // longitude decimal degrees
    double leaft, // leaf transmittance
    int maxIter = 100, // max number of iterations
    double tolerance = 1e-2) // error tolerence for convergence
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Extract from climdata data.frame
    std::vector<double> Rdirdown = climdata["Rdirdown"]; // Direct downward radiation flux perpendicular to solar beam (W/m^2)
    std::vector<double> Rdifdown = climdata["Rdifdown"]; // Diffuse downward radiation flux (W/m^2)
    std::vector<double> Rswup = climdata["Rswup"]; // Upward shortwave radiation flux  (W/m^2)
    std::vector<double> Rlwdown = climdata["Rlwdown"]; // Downward lognwave flux (W/m^2)
    std::vector<double> Rlwup = climdata["Rlwup"]; // Upward longwave flux (W/m^2)
    std::vector<double> uz = climdata["uz"]; // Wind speed at height of animal (m/s
    std::vector<double> wdir = climdata["wdir"]; // wind direction (only needed if position = fixed
    std::vector<double> Ta = climdata["Ta"]; // Air temperature (deg C)
    std::vector<double> Ts = climdata["Ts"]; // Surface temperature (only needed if confrac > 0)
    std::vector<double> rh = climdata["rh"]; // Air relative humidity (percentage)
    std::vector<double> pk = climdata["pk"]; // Atmospheric pressure (kPa)
    std::vector<double> surfrh = climdata["surfrh"]; // Surface relative humidity (0-1)
    std::vector<double> surfwet(surfrh.size(), 1.0); // Fraction of the contact surface that is wet or porous (0-1)
    if (climdata.containsElementNamed("surfwet")) surfwet = Rcpp::as<std::vector<double>>(climdata["surfwet"]);
    // Extract animal parameters
    double height = animal["height"]; // height of animal (cm)
    double width = animal["width"]; // width of animal (cm)
    double length = animal["len"]; // length of animal (cm)
    double refl = animal["refl"]; // reflectance of animal (0-1)
    double confrac = animal["confrac"]; // fraction of animal in direct contact with surface
    double rc = animal["rc"]; // cutaneous resistance (s / m)
    double em = animal["em"]; // emissvity of animal 
    double rho = animal["rho"]; // animal density (kg / m^3)
    double volume = animal["volume"]; // animal volume (m^3)
    double area = animal["area"]; // animal surface area (m^2)
    double Q10 = animal["Q10"]; // factor by which metabolic rate changes for a 10 degrees C temperature increase
    double a0 = animal["a0"]; // normalization constant at reference temperature for calculating metabolic rate
    double b = animal["b"]; // mass scaling exponent for calculating metabolic rate
    double Tref = animal["Tref"]; // reference metabolic calibration temperature (deg C)
    double adir = animal["adir"]; // direction animal is facing relative to north (ignored if position != fixed)
    double atilt = animal["atilt"]; // direction of tilt of longest axis of animal relative to horizontal (ignored if position = random)
    double k = animal["k"]; // animal heat conductance W/m
    std::string position = animal["position"]; // see details above
    size_t n = Ta.size();
    double latr = lat * torad;
    double lonr = lon * torad;
    Rcpp::NumericVector Tbody(n);
    for (size_t i = 0; i < n; ++i) {
        Tbody[i] = Ta[i] + 3.0;
        double dT = Tbody[i] - Ta[i]; // initial animal air temperature difference
        solmodel solp = solpositionCpp2(latr, lonr, year[i], month[i], day[i], hour[i]);
        double Rsw_flux = 0.0;
        if (solp.zenr < pi / 2.0) {
            silstruct sa = silhouette(solp.zenr, solp.azir, height, width, length, adir, atilt, position);
            double sc = sa.silA / sa.A; // solar coefficient
            if (atilt < 0.0) { // animal suspended below leaf
                Rsw_flux = leaft * (sc * Rdirdown[i] + 0.5 * Rdifdown[i]) + 0.5 * Rswup[i];
            }
            else {
                Rsw_flux = sc * Rdirdown[i] + 0.5 * Rdifdown[i] + 0.5 * Rswup[i];
            }
        }
        double Rabs = (1.0 - refl) * Rsw_flux + 0.5 * em * (Rlwdown[i] + Rlwup[i]);
        double danim = chardim(wdir[i], solp.zenr, solp.azir, height, width, length, adir, atilt, position);
        double tdif = 1e99;
        WAitkenStateScalar st;
        int nrIterations = 0;
        while (tdif > tolerance && nrIterations < maxIter) {
            double Te = (Tbody[i] + Ta[i]) / 2.0;
            double Tf = (Tbody[i] + Ts[i]) / 2.0;
            double rHa = animalrHa(Ta[i], std::abs(dT), uz[i], danim, 300.0);
            double M = metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody[i]) / area;
            double Told = Tbody[i];
            Tbody[i] = PenmanMonteith_animal(Rabs, Ta[i], Ts[i], Te, Tf, pk[i], rh[i], rHa, height,
                rc, confrac, M, em, k, surfrh[i], surfwet[i]);
            aitkin_weightdif_scalar(Told, Tbody[i], st);
            double dTold = dT;
            dT = Tbody[i] - Ta[i];
            tdif = std::abs(dTold - dT);
            ++nrIterations;
        }
    }
    return Tbody;
}
// As Ectotherm, for microclimate given at several heights: returns body
// temperature at each time step and height.
// [[Rcpp::export]]
Rcpp::NumericVector EctothermM(Rcpp::DataFrame obstime, Rcpp::List climdata, Rcpp::List animal,
    double lat, // Latitude (decimal degrees)
    double lon, // longitude decimal degrees
    double leaft, // leaf transmittance
    int maxIter = 100, // max number of iterations
    double tolerance = 1e-2) // error tolerence for convergence
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    Rcpp::NumericMatrix Rdirdown = climdata["Rdirdown"]; // Direct downward radiation flux perpendicular to solar beam (W/m^2)
    Rcpp::NumericMatrix Rdifdown = climdata["Rdifdown"]; // Diffuse downward radiation flux (W/m^2)
    Rcpp::NumericMatrix Rswup = climdata["Rswup"]; // Upward shortwave radiation flux  (W/m^2)
    Rcpp::NumericMatrix Rlwdown = climdata["Rlwdown"]; // Downward lognwave flux (W/m^2)
    Rcpp::NumericMatrix Rlwup = climdata["Rlwup"]; // Upward longwave flux (W/m^2)
    Rcpp::NumericMatrix uz = climdata["uz"]; // Wind speed at height of animal (m/s
    Rcpp::NumericMatrix Ta = climdata["Ta"]; // Air temperature (deg C)
    Rcpp::NumericMatrix Ts = climdata["Ts"]; // Surface temperature (only needed if confrac > 0)
    Rcpp::NumericMatrix rh = climdata["rh"]; // Air relative humidity (percentage)
    Rcpp::NumericVector pk = climdata["pk"]; // Atmospheric pressure (kPa)
    const bool hasWet = climdata.containsElementNamed("surfwet"); // fraction of the contact surface that is wet (0-1)
    Rcpp::NumericMatrix surfwet = hasWet ? Rcpp::as<Rcpp::NumericMatrix>(climdata["surfwet"]) : Rcpp::NumericMatrix(1, 1);
    Rcpp::NumericVector wdir = climdata["wdir"]; // wind direction (only needed if position = fixed
    double height = animal["height"]; // height of animal (cm)
    double width = animal["width"]; // width of animal (cm)
    double length = animal["len"]; // length of animal (cm)
    double refl = animal["refl"]; // reflectance of animal (0-1)
    double confrac = animal["confrac"]; // fraction of animal in direct contact with surface
    double rc = animal["rc"]; // cutaneous resistance (s / m)
    double em = animal["em"]; // emissvity of animal 
    double rho = animal["rho"]; // animal density (kg / m^3)
    double volume = animal["volume"]; // animal volume (m^3)
    double area = animal["area"]; // animal surface area (m^2)
    double Q10 = animal["Q10"]; // factor by which metabolic rate changes for a 10 degrees C temperature increase
    double a0 = animal["a0"]; // normalization constant at reference temperature for calculating metabolic rate
    double b = animal["b"]; // mass scaling exponent for calculating metabolic rate
    double Tref = animal["Tref"]; // reference metabolic calibration temperature (deg C)
    double adir = animal["adir"]; // direction animal is facing relative to north (ignored if position != fixed)
    double atilt = animal["atilt"]; // direction of tilt of longest axis of animal relative to horizontal (ignored if position = random)
    double k = animal["k"]; // animal heat conductance W/m
    std::string position = animal["position"]; // see details above
    int tsteps = Ta.nrow();
    int nlayers = Ta.ncol();
    Rcpp::NumericMatrix Tbody(tsteps, nlayers);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (int i = 0; i < tsteps; ++i) {
        solmodel solp = solpositionCpp2(latr, lonr, year[i], month[i], day[i], hour[i]);
        double sc = 1.0; // solar coefficient
        if (solp.zenr < pi / 2.0) {
            silstruct sa = silhouette(solp.zenr, solp.azir, height, width, length, adir, atilt, position);
            sc = sa.silA / sa.A;
        }
        double danim = chardim(wdir[i], solp.zenr, solp.azir, height, width, length, adir, atilt, position);
        for (int j = 0; j < nlayers; ++j) {
            double Rsw_flux = 0.0;
            if (solp.zenr < pi / 2.0) {
                if (atilt < 0.0) { // animal suspended below leaf
                    Rsw_flux = leaft * (sc * Rdirdown(i,j) + 0.5 * Rdifdown(i, j)) + 0.5 * Rswup(i, j);
                }
                else {
                    Rsw_flux = sc * Rdirdown(i, j) + 0.5 * Rdifdown(i, j) + 0.5 * Rswup(i, j);
                }
            }
            double Rabs = (1.0 - refl) * Rsw_flux + 0.5 * em * (Rlwdown(i, j) + Rlwup(i, j));
            Tbody(i, j) = Ta(i, j) + 3.0;
            double dT = Tbody(i, j) - Ta(i, j); // initial animal air temperature difference
            double tdif = 1e99;
            WAitkenStateScalar st;
            int nrIterations = 0;
            while (tdif > tolerance && nrIterations < maxIter) {
                double Te = (Tbody(i, j) + Ta(i, j)) / 2.0;
                double Tf = (Tbody(i, j) + Ts(i, j)) / 2.0;
                double rHa = animalrHa(Ta(i, j), std::abs(dT), uz(i, j), danim, 300.0);
                double M = metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody(i, j)) / area;
                double Told = Tbody(i, j);
                Tbody(i, j) = PenmanMonteith_animal(Rabs, Ta(i, j), Ts(i, j), Te, Tf, pk[i], rh(i, j), rHa, height,
                    rc, confrac, M, em, k, 1.0, hasWet ? surfwet(i, j) : 1.0);
                aitkin_weightdif_scalar(Told, Tbody(i, j), st);
                double dTold = dT;
                dT = Tbody(i, j) - Ta(i, j);
                tdif = std::abs(dTold - dT);
                ++nrIterations;
            } // end iterative while loop
        } // end j
    } // end i
    return Tbody;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** R wrappers *********************************************************** //
// Functions called from R. They convert the R inputs, run the model hour by hour
// and return the results.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Plant area index above the bottom of each canopy layer (cumulative from the top down).
// [[Rcpp::export]]
std::vector<double> reverseCumsum(const std::vector<double>& paii) {
    std::vector<double> rev_paii = paii;
    std::reverse(rev_paii.begin(), rev_paii.end());
    std::vector<double> cum(rev_paii.size());
    cum[0] = rev_paii[0];
    for (size_t i = 1; i < rev_paii.size(); ++i) {
        cum[i] = cum[i - 1] + rev_paii[i];
    }
    std::reverse(cum.begin(), cum.end());
    return cum;
}
// Converts the R vegetation parameters to the C++ structure.
static vegpstruct tovegpstruct(List vegp, std::vector<double> paii, std::vector<double> Lfrac) {
    // Constant within canopy
    vegpstruct out;
    out.hgt = vegp["h"]; // vegetation height (m)
    out.x = vegp["x"]; // Campbell foliage angle coefficient (unitless)
    out.lref = vegp["lref"]; // Leaf reflectance (shortwave radiation, 0-1)
    out.ltra = vegp["ltra"]; // Leaf transmittance (shortwave radiation, 0-1)
    out.lrefp = vegp["lrefp"]; // Leaf reflectance (PAR, 0-1)
    out.ltrap = vegp["ltrap"]; // Leaf transmittance (PAR, 0-1)
    out.Vcmax25 = vegp["Vcmx25"]; // micromol to mol / m ^ 2 / s ^ 1
    out.Tup = vegp["Tup"]; // high temperature photosynthesis range (deg C)
    out.Tlw = vegp["Tlow"]; // low temperature photosynthesis range (deg C)
    out.Dcrit = vegp["Dcrit"]; // Onset of strong stomatal limitation (kPa)
    out.alpha = vegp["alpha"]; // photosynthesis quantum efficiency (mol CO2 / mol PAR)
    out.Kxmx = vegp["Kxmx"]; // xylem maximum hydraulic conductance (mol to micromol / m^2 / s)
    out.hv = vegp["hv"]; // Huber value (cm^2 to m^2/m^2)
    out.f0 = vegp["f0"]; // ratio of ci/(ci-Gma) under low D from Jacobs (1994) equation (unitless) 
    out.fd = vegp["fd"]; // reaction of Vcmax converted to leaf dark respiration (unitless)
    out.psi50 = vegp["psi50"]; // Water potential when plant loses 50% of krc_max (MPa) 
    out.apsi = vegp["apsi"]; // xylem hydraulic conductance parameter. Calculated if set to less than zero (unitless)
    out.len = vegp["len"];  // Leaf length (m)
    out.wid = vegp["wid"]; // Leaf width (m)
    out.vegem = vegp["vegem"]; // Vegetation emissivity (unitless)
    out.mwft = vegp["mwft"]; // max water film thickness (mm)
    out.pTAW = vegp["pTAW"]; // Fraction of total available water plant can deplete before stress (unitless)
    out.root50 = vegp["root50"]; // Depth above which half of all roots lie (m)
    out.root95 = vegp["root95"]; // Depth above which 95% of all roots lie (m)
    int n = static_cast<int>(paii.size());
    out.pai = 0.0; // Total one-sided plant area index (m^2 / m^2)
    for (int i = 0; i < n; ++i) out.pai += paii[i];
    out.paii = paii; // one-sided plant area index at each canopy node
    out.pia = reverseCumsum(paii); // one-sided plant area index above each canopy node
    out.Lfrac = Lfrac; // Fraction of paii that is living plant material
    return out;
}
// Converts the R soil parameters to the C++ structure.
static soilpstruct tosoilpstruct(List soilp) {
    soilpstruct out;
    // Scalars
    out.slope = Rcpp::as<double>(soilp["slope"]) * torad; // radians
    out.aspect = Rcpp::as<double>(soilp["aspect"]) * torad; // radians
    out.gref = Rcpp::as<double>(soilp["gref"]);
    out.grefPAR = Rcpp::as<double>(soilp["grefPAR"]);
    out.groundem = Rcpp::as<double>(soilp["groundem"]);
    out.nLayers = Rcpp::as<int>(soilp["nLayers"]);
    // Vectors
    out.Vq = Rcpp::as<std::vector<double>>(soilp["Vq"]);
    out.Vm = Rcpp::as<std::vector<double>>(soilp["Vm"]);
    out.Vo = Rcpp::as<std::vector<double>>(soilp["Vo"]);
    out.Mc = Rcpp::as<std::vector<double>>(soilp["Mc"]);
    out.psie = Rcpp::as<std::vector<double>>(soilp["psi_e"]); // key is "psi_e"
    out.b = Rcpp::as<std::vector<double>>(soilp["b"]);
    out.thetaR = Rcpp::as<std::vector<double>>(soilp["Smin"]);
    out.thetaS = Rcpp::as<std::vector<double>>(soilp["Smax"]);
    out.Ksat = Rcpp::as<std::vector<double>>(soilp["Ksat"]);
    out.FreeDrain = Rcpp::as<bool>(soilp["FreeDrain"]);
    // Roughness length of the soil for momentum, if given.
    out.z0 = soilp.containsElementNamed("zm") ? Rcpp::as<double>(soilp["zm"]) : Z0_GROUND;
    // Air entry potential is negative.
    for (double& v : out.psie) v = -std::abs(v);
    // Lower limit of soil water potential: oven dry.
    out.psi_min = std::vector<double>(out.psie.size(), PSI_OVEN_DRY);
    return out;
}
// Initial state of the soil heat model, from profiles of temperature and water content.
static soilmod toSoilheatmod(List soilp, std::vector<double> Te, std::vector<double> wc)
{
   
    int n = Rcpp::as<int>(soilp["nLayers"]);
    double totalDepth = Rcpp::as<double>(soilp["totalDepth"]);
    std::vector<double> z = geometricCpp(n, totalDepth);
    std::vector<double> dz(n + 1);
    for (int i = 0; i <= n; ++i) dz[i] = z[i + 1] - z[i];
    soilmod state;
    state.n = n;
    state.z = z; // node depths
    state.dz = dz; // layer thickness
    // Thickness of soil each node stores heat for: from halfway to the node above to
    // halfway to the node below (Campbell 1985, Eq. 4.7).
    std::vector<double> vol(n + 1);
    vol[0] = 0.5 * (z[1] - z[0]);
    for (int i = 1; i <= n; ++i) vol[i] = 0.5 * (z[i + 1] - z[i - 1]);
    state.vol = vol;
    state.wc = wc; // volumetric water fraction of layer
    state.Te = Te; // temperature of layer
    state.oldTe = Te; // old temperature of layer
    state.Gflux = 0.0;
    state.iters = 0;
    return state;
}
// Initial state of the soil water model: water potential, vapour and root fractions.
static soilwatermod toSoilwatermod(soilpstruct soilpc, soilmod state, double root50, double root95)
{
    soilwatermod out;
    out.n = state.n; // number of layers
    out.z = state.z; // node depths
    out.dz = state.dz; // layer thickness
    out.Tc = state.Te; // temperature of soil profile
    out.oldTc = state.oldTe; // old temperature of soil profile
    out.theta = state.wc; // volumetric water content of soil profile
    // Initialize
    std::vector<double> vol(out.n + 1);
    std::vector<double> psiw(out.n + 1);
    std::vector<double> k(out.n + 1);
    std::vector<double> vapor(out.n + 1);
    // Thickness of soil each node stores water for: half the distance to each
    // neighbour, and for the surface node half the top layer.
    vol[0] = (out.z[1] - out.z[0]) / 2.0;
    for (int i = 0; i <= out.n; ++i) {
        if (i > 0) vol[i] = (out.z[i + 1] - out.z[i - 1]) / 2.0;
        psiw[i] = waterPotential(soilpc, out.theta[i], i);
        k[i] = hydraulicConductivityFromTheta(soilpc, out.theta[i], i);
        vapor[i] = vaporFromPsi(soilpc, psiw[i], out.theta[i], out.Tc[i] + 273.15, i);
    }
    out.vol = vol;
    out.psiw = psiw;
    out.vapor = vapor;
    out.k = k;
    out.oldvapor = vapor;
    out.oldtheta = out.theta;
    out.rootfrac = root_distribute(out.z, root50, root95);
    return out;
}
// Heights of the canopy nodes: the centres of n layers of equal depth. Values
// for a layer are means over the layer. Canopy top is not a node.
static std::vector<double> canopyNodeHeights(size_t n, double hgt)
{
    std::vector<double> z(n);
    for (size_t i = 0; i < n; ++i) z[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(n) * hgt;
    return z;
}
// Initial state of the multilayer model: every layer at the same temperature and humidity.
static onestep CreateOneStep(soilmod soilheatvars, soilwatermod soilwatervars, double temp0, double rh0, double Rlw0, int na) {
    std::vector<double> v0(na, 0.0);
    std::vector<double> v1(na, 1.0);
    std::vector<double> Rlwv(na, Rlw0);
    std::vector<double> tv(na, temp0);
    std::vector<double> rv(na, rh0);
    onestep onestepin;
    onestepin.soilheatvars = std::move(soilheatvars);
    onestepin.soilwatervars = std::move(soilwatervars);
    onestepin.Rdirdown = v0;
    onestepin.Rdifdown = v0;
    onestepin.Rswup = v0;
    onestepin.Rlwdown = Rlwv;
    onestepin.Rlwup = Rlwv;
    onestepin.tair = tv;
    onestepin.tleaf = tv;
    onestepin.rh = rv;
    onestepin.uz = v1;
    onestepin.rLB = v1;
    onestepin.swaterdepth = v0;
    onestepin.leafwet = v0;
    onestepin.Hz = v0;
    onestepin.Lz = v0;
    onestepin.gs = v0;
    onestepin.dSTdT = v0;
    onestepin.dSLde = v0;
    onestepin.tfar = v0; onestepin.tnear = v0; onestepin.efar = v0; onestepin.enear = v0;
    onestepin.precipground = 0.0;
    onestepin.H = 0.0;
    onestepin.L = 0.0;
    onestepin.Et = 0.0;
    onestepin.Ev = 0.0;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    onestepin.psim = 0.0;
    onestepin.psih = 0.0;
    onestepin.phih = 1.0;
    onestepin.fade = 1.0;
    onestepin.tcanopy = temp0;
    onestepin.ucanopy = 1.0;
    onestepin.windscale = 1.0;
    onestepin.rhcanopy = rh0;
    onestepin.LL = 999.99;
    onestepin.iters = 0;
    onestepin.witers = 0;
    onestepin.error = 0.0;
    return onestepin;
}
// Converts the model state to an R list.
static Rcpp::List OneStepCpptoList(onestep onestepin, std::vector<double> z)
{
    Rcpp::List out;
    // Soil heat variables
    out["nb"] = Rcpp::wrap(onestepin.soilheatvars.n);
    out["zb"] = Rcpp::wrap(onestepin.soilheatvars.z);
    out["dz"] = Rcpp::wrap(onestepin.soilheatvars.dz);
    out["thetas"] = Rcpp::wrap(onestepin.soilwatervars.theta);
    out["Soiltemp"] = Rcpp::wrap(onestepin.soilheatvars.Te);
    out["G"] = Rcpp::wrap(onestepin.soilheatvars.Gflux);
    out["soilhiters"] = Rcpp::wrap(onestepin.soilheatvars.iters);
    // Soil water variables
    out["vol"] = Rcpp::wrap(onestepin.soilwatervars.vol);
    out["psiw"] = Rcpp::wrap(onestepin.soilwatervars.psiw);
    out["k"] = Rcpp::wrap(onestepin.soilwatervars.k);
    out["vapor"] = Rcpp::wrap(onestepin.soilwatervars.vapor);
    out["rootfrac"] = Rcpp::wrap(onestepin.soilwatervars.rootfrac);
    // Radiation streams
    out["Rdirdown"] = Rcpp::wrap(onestepin.Rdirdown);
    out["Rdifdown"] = Rcpp::wrap(onestepin.Rdifdown);
    out["Rswup"] = Rcpp::wrap(onestepin.Rswup);
    out["Rlwdown"] = Rcpp::wrap(onestepin.Rlwdown);
    out["Rlwup"] = Rcpp::wrap(onestepin.Rlwup);
    // Below canopy profiles
    out["tair"] = Rcpp::wrap(onestepin.tair);
    out["tleaf"] = Rcpp::wrap(onestepin.tleaf);
    out["rh"] = Rcpp::wrap(onestepin.rh);
    out["tcanopy"] = onestepin.tcanopy;
    out["rhcanopy"] = onestepin.rhcanopy;
    // Actual wind speeds: the convective addition is removed.
    std::vector<double> uzActual = onestepin.uz;
    for (double& u : uzActual) u *= onestepin.windscale;
    out["ucanopy"] = onestepin.ucanopy * onestepin.windscale;
    out["uz"] = Rcpp::wrap(uzActual);
    out["runoff"] = onestepin.runoff;
    out["witers"] = onestepin.witers;
    out["rLB"] = Rcpp::wrap(onestepin.rLB);
    out["swaterdepth"] = Rcpp::wrap(onestepin.swaterdepth);
    out["leafwet"] = Rcpp::wrap(onestepin.leafwet);
    out["Hz"] = Rcpp::wrap(onestepin.Hz);
    out["tfar"] = Rcpp::wrap(onestepin.tfar);
    out["tnear"] = Rcpp::wrap(onestepin.tnear);
    out["efar"] = Rcpp::wrap(onestepin.efar);
    out["enear"] = Rcpp::wrap(onestepin.enear);
    out["Lz"] = Rcpp::wrap(onestepin.Lz);
    out["gs"] = Rcpp::wrap(onestepin.gs);
    // Additional variables
    out["precipground"] = onestepin.precipground;
    out["H"] = onestepin.H;
    out["L"] = onestepin.L;
    out["Et"] = onestepin.Et;
    out["Ev"] = onestepin.Ev;
    out["theta"] = onestepin.theta;
    out["psim"] = onestepin.psim;
    out["psih"] = onestepin.psih;
    out["phih"] = onestepin.phih;
    out["fade"] = onestepin.fade;
    out["LL"] = onestepin.LL;
    out["iters"] = onestepin.iters;
    out["witers"] = Rcpp::wrap(onestepin.witers);
    out["error"] = onestepin.error;
    // Above ground z
    out["z"] = Rcpp::wrap(z);
    return out;
}
// Profile of air temperature, humidity and wind over bare ground at one chosen
// hour. The model is run up to that hour, which is then solved at the heights
// requested.
// [[Rcpp::export]]
List profilebareR(size_t hourtoplot, DataFrame obstime, DataFrame climdata, List soilc, std::vector<double> z,
    double zref, double lat, double lon, std::vector<double> SoilTempIni,
    std::vector<double> SoilThetaIni, double zm = 0.004, int maxNrIterations = 100, double tolerance = 0.1)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    soilpstruct soilpc = tosoilpstruct(soilc);
    std::vector<double> z2 = { zm, zref / 2.0 };
    onestepbare onestepin;
    onestepin.tair = { temp[0], temp[0] };
    onestepin.rh = { relhum[0], relhum[0] };
    onestepin.uz = { wspeed[0], wspeed[0] };
    onestepin.soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 0.05, 0.4);
    onestepin.H = 0.0;
    onestepin.L = 0.0;
    onestepin.psim = 0.0;
    onestepin.psih = 0.0;
    onestepin.Rswup = 0.0;
    onestepin.Rb0 = 0.0;
    onestepin.Rlwup = 0.0;
    onestepin.LL = 999.99;
    onestepin.iters = 1;
    onestepin.error = 1e99;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    double latr = torad * lat;
    double lonr = torad * lon;
    if (hourtoplot > 0) {
        for (size_t hr = 0; hr < hourtoplot; ++hr) {
            obsstruct obsdata;
            climstruct climin;
            obsdata.year = year[hr];
            obsdata.month = month[hr];
            obsdata.day = day[hr];
            obsdata.hour = hour[hr];
            climin.tref = temp[hr];
            climin.relhum = relhum[hr];
            climin.pk = pres[hr];
            climin.Rsw = Rsw[hr];
            climin.Rdif = Rdif[hr];
            climin.Rlw = Rlw[hr];
            climin.uref = wspeed[hr];
            climin.winddir = wdir[hr];
            climin.precip = precip[hr];
            onestepin = OneStepBare(onestepin, obsdata, climin, soilpc, z2, latr, lonr, zref, zm, maxNrIterations, tolerance);
            onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
            onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
            onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
        }
    }
    size_t hr = hourtoplot;
    obsstruct obsdata;
    climstruct climin;
    obsdata.year = year[hr];
    obsdata.month = month[hr];
    obsdata.day = day[hr];
    obsdata.hour = hour[hr];
    climin.tref = temp[hr];
    climin.relhum = relhum[hr];
    climin.pk = pres[hr];
    climin.Rsw = Rsw[hr];
    climin.Rdif = Rdif[hr];
    climin.Rlw = Rlw[hr];
    climin.uref = wspeed[hr];
    climin.winddir = wdir[hr];
    climin.precip = precip[hr];
    onestepin = OneStepBare(onestepin, obsdata, climin, soilpc, z, latr, lonr, zref, zm, maxNrIterations, tolerance);
    size_t nn = z.size();
    Rcpp::List out;
    out["nb"] = onestepin.soilheatvars.n;
    out["zb"] = onestepin.soilheatvars.z;
    out["dz"] = onestepin.soilheatvars.dz;
    out["thetas"] = onestepin.soilwatervars.theta;
    out["Soiltemp"] = onestepin.soilheatvars.Te;
    out["G"] = onestepin.soilheatvars.Gflux;
    out["soilhiters"] = onestepin.soilheatvars.iters;
    out["vol"] = onestepin.soilwatervars.vol;
    out["psiw"] = onestepin.soilwatervars.psiw;
    out["k"] = onestepin.soilwatervars.k;
    out["vapor"] = onestepin.soilwatervars.vapor;
    out["rootfrac"] = onestepin.soilwatervars.rootfrac;
    out["Rdirdown"] = NumericVector(nn, onestepin.Rb0);
    out["Rdifdown"] = NumericVector(nn, Rdif[hr]);
    out["Rswup"] = NumericVector(nn, onestepin.Rswup);
    out["Rlwdown"] = NumericVector(nn, Rlw[hr]);
    out["Rlwup"] = NumericVector(nn, onestepin.Rlwup);
    out["tair"] = onestepin.tair;
    out["rh"] = onestepin.rh;
    out["uz"] = onestepin.uz;
    out["H"] = onestepin.H;
    out["L"] = onestepin.L;
    out["Ev"] = onestepin.Ev;
    out["theta"] = onestepin.theta;
    out["psim"] = onestepin.psim;
    out["psih"] = onestepin.psih;
    out["LL"] = onestepin.LL;
    out["iters"] = onestepin.iters;
    out["error"] = onestepin.error;
    out["z"] = z;
    return out;
}
// Full state of the vegetated model at one chosen hour.
// [[Rcpp::export]]
List profileR(size_t hourtoplot, DataFrame obstime, DataFrame climdata, List soilc, List vegp,
    std::vector<double> paii20, std::vector<double> paii, std::vector<double> Lfrac20, std::vector<double> Lfrac, 
    double zref, double Ca, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double tolerance = 0.1, double a0 = -1.0, double a1 = 1.25, bool C3 = true)
{
    // Run the model up to the chosen hour
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    vegpstruct vegpc = tovegpstruct(vegp, paii20, Lfrac20);
    soilpstruct soilpc = tosoilpstruct(soilc);
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc, soilpc.z0);
    size_t na = paii20.size();
    const std::vector<double> z = canopyNodeHeights(na, vegpc.hgt);
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.root50, vegpc.root95);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    double latr = lat * torad;
    double lonr = lon * torad;
    if (hourtoplot > 0) {
        for (size_t hr = 0; hr < hourtoplot; ++hr) {
            obsstruct obsdata;
            climstruct climin;
            obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
            climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
            climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
            climin.winddir = wdir[hr]; climin.precip = precip[hr];
            onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR, tspdif,
                tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
            onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
            onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
            onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
        }
        soilheatvars = onestepin.soilheatvars;
        soilwatervars = onestepin.soilwatervars;
    }
    // Solve the chosen hour and return its full state
    vegpc = tovegpstruct(vegp, paii, Lfrac);
    wgts = lwradweights(vegpc.paii);
    wc = windprofileCpp(vegpc, soilpc.z0);
    size_t hr = hourtoplot;
    obsstruct obsdata;
    climstruct climin;
    obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
    climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
    climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
    climin.winddir = wdir[hr]; climin.precip = precip[hr];
    na = paii.size();
    const std::vector<double> z2 = canopyNodeHeights(na, vegpc.hgt);
    onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[hr], relhum[hr], Rlw[hr], na);
    onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z2, tspveg, tspvegPAR, tspdif,
        tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
    return OneStepCpptoList(onestepin, z2);
}
// Bare ground: the whole time series at one height (below ground if reqhgt is negative).
// [[Rcpp::export]]
DataFrame RunBareR(double reqhgt, DataFrame obstime, DataFrame climdata, List soilc, double zref,
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    double zm = 0.004, int maxNrIterations = 100, double  tolerance = 0.1)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    soilpstruct soilpc = tosoilpstruct(soilc);
    // ** Derive additional variables
    size_t n = temp.size();
    std::vector<double> z(2);
    if (reqhgt > 0.0) {
        z = { zm, reqhgt };
    }
    else {
        z = { zm, zref / 2.0 };
    }
    // Initialize onstepin
    onestepbare onestepin;
    onestepin.tair = { temp[0], temp[0] };
    onestepin.rh = { relhum[0], relhum[0] };
    onestepin.uz = { wspeed[0], wspeed[0] };
    onestepin.soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 0.05, 0.4);
    onestepin.H = 0.0; onestepin.L = 0.0; onestepin.psim = 0.0; onestepin.psih = 0.0; 
    onestepin.Rswup = 0.0; onestepin.Rb0 = 0.0; onestepin.Rlwup = 0.0;
    onestepin.LL = 999.99; onestepin.iters = 1; onestepin.error = 1e99;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    // Soil node nearest the depth requested, if below ground
    size_t sel = 0;
    if (reqhgt < 0.0) {
        // The last depth is below the boundary and has no node.
        std::vector<double> zz = onestepin.soilheatvars.z;
        double mn = std::abs(reqhgt + zz[0]);
        for (size_t i = 1, n = zz.size() - 1; i < n; ++i) {
            double dif = std::abs(reqhgt + zz[i]);
            if (dif < mn) {
                mn = dif;
                sel = i;
            }
        }
    }
    // Create output variables for storing
    std::vector<double> Rdirdown(n);
    std::vector<double> Rdifdown(n);
    std::vector<double> Rswup(n);
    std::vector<double> Rlwdown(n);
    std::vector<double> Rlwup(n);
    std::vector<double> tair(n);
    std::vector<double> tground(n);
    std::vector<double> rz(n);
    std::vector<double> uz(n);
    std::vector<double> H(n);
    std::vector<double> L(n);
    std::vector<double> G(n);
    std::vector<double> Ev(n);
    std::vector<double> surfwater(n); // water content of the soil surface
    std::vector<double> runoff(n), drainage(n);
    std::vector<int> iters(n), witers(n);
    std::vector<double> error(n);
    // Run model for every hr
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
        climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
        climin.winddir = wdir[hr]; climin.precip = precip[hr];
        onestepin = OneStepBare(onestepin, obsdata, climin, soilpc, z, latr, lonr, zref, zm,
            maxNrIterations, tolerance);
        // Extract values
        tground[hr] = onestepin.soilheatvars.Te[0];
        H[hr] = onestepin.H;
        L[hr] = onestepin.L;
        G[hr] = onestepin.soilheatvars.Gflux;
        Ev[hr] = onestepin.Ev;
        surfwater[hr] = onestepin.theta;
        runoff[hr] = onestepin.runoff;
        drainage[hr] = onestepin.drainage;
        witers[hr] = onestepin.witers;
        iters[hr] = onestepin.iters;
        error[hr] = onestepin.error;
        if (reqhgt >= 0.0) {
            Rdirdown[hr] = onestepin.Rb0;
            Rdifdown[hr] = Rdif[hr];
            Rswup[hr] = onestepin.Rswup;
            Rlwdown[hr] = Rlw[hr];
            Rlwup[hr] = onestepin.Rlwup;
            if (reqhgt > 0.0) {
                tair[hr] = onestepin.tair[1];
                rz[hr] = onestepin.rh[1];
                uz[hr] = onestepin.uz[1];
            }
            else {
                tair[hr] = onestepin.soilheatvars.Te[0];
                rz[hr] = onestepin.soilwatervars.theta[0];
                uz[hr] = 0.0;
            }
        }
        else {
            Rdirdown[hr] = 0.0;
            Rdifdown[hr] = 0.0;
            Rswup[hr] = 0.0;
            Rlwdown[hr] = 0.0;
            Rlwup[hr] = 0.0;
            tair[hr] = onestepin.soilheatvars.Te[sel];
            rz[hr] = onestepin.soilwatervars.theta[sel];
            uz[hr] = 0.0;
        }
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    return DataFrame::create(
        // Obstime
        Named("year") = Rcpp::wrap(year),
        Named("month") = Rcpp::wrap(month),
        Named("day") = Rcpp::wrap(day),
        Named("hour") = Rcpp::wrap(hour),
        // Radiation streams
        Named("Rdirdown") = Rcpp::wrap(Rdirdown),
        Named("Rdifdown") = Rcpp::wrap(Rdifdown),
        Named("Rswup") = Rcpp::wrap(Rswup),
        Named("Rlwdown") = Rcpp::wrap(Rlwdown),
        Named("Rlwup") = Rcpp::wrap(Rlwup),
        // Other climate variables
        Named("tair") = Rcpp::wrap(tair),
        Named("tground") = Rcpp::wrap(tground),
        Named("relhum") = Rcpp::wrap(rz),
        Named("windspeed") = Rcpp::wrap(uz),
        // Fluxes
        Named("H") = Rcpp::wrap(H),
        Named("L") = Rcpp::wrap(L),
        Named("G") = Rcpp::wrap(G),
        Named("Evaporation") = Rcpp::wrap(Ev),
        Named("surfwater") = Rcpp::wrap(surfwater),
        Named("runoff") = Rcpp::wrap(runoff),
        Named("drainage") = Rcpp::wrap(drainage),
        // Error variables
        Named("iters") = Rcpp::wrap(iters),
        Named("witers") = Rcpp::wrap(witers),
        Named("error") = Rcpp::wrap(error)
    );
}
// Vegetated surface: the whole time series at one height.
// [[Rcpp::export]]
Rcpp::List RunModelR(double reqhgt, Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List soilc,
    Rcpp::List vegp, std::vector<double> paii, std::vector<double> Lfrac, double zref, double Ca,
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni,
    int maxNrIterations = 100, double tolerance = 0.1, double a0 = -1.0,
    double a1 = 1.25, bool C3 = true)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    tsvegstruct tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc, soilpc.z0);
    size_t n = temp.size();
    size_t na = paii.size();
    const std::vector<double> z = canopyNodeHeights(na, vegpc.hgt);
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.root50, vegpc.root95);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    double latr = lat * torad;
    double lonr = lon * torad;
    std::vector<double> zz = z;
    if (reqhgt < 0.0) {
        // The last depth is below the boundary and has no node.
        zz = soilheatvars.z;
        zz.pop_back();
        for (size_t i = 0; i < zz.size(); ++i)
            zz[i] = -zz[i];
    }
    size_t sel = 0;
    double mn = std::abs(reqhgt - zz[0]);
    for (size_t i = 1; i < zz.size(); ++i) {
        double dif = std::abs(reqhgt - zz[i]);
        if (dif < mn) {
            mn = dif;
            sel = i;
        }
    }
    std::vector<double> Rdirdown(n), Rdifdown(n), Rswup(n), Rlwdown(n), Rlwup(n);
    std::vector<double> tair(n), tground(n), tleaf(n), leafwet(n, -999.99), rz(n), uz(n);
    std::vector<double> H(n), L(n), G(n), Ev(n), Et(n), error(n), runoff(n), drainage(n);
    std::vector<double> surfwater(n); // water content of the soil surface
    std::vector<int> iters(n), witers(n);
    for (size_t hr = 0; hr < n; ++hr) {
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr];
        obsdata.month = month[hr];
        obsdata.day = day[hr];
        obsdata.hour = hour[hr];
        climin.tref = temp[hr];
        climin.relhum = relhum[hr];
        climin.pk = pres[hr];
        climin.Rsw = Rsw[hr];
        climin.Rdif = Rdif[hr];
        climin.Rlw = Rlw[hr];
        climin.uref = wspeed[hr];
        climin.winddir = wdir[hr];
        climin.precip = precip[hr];
        onestepin = OneStepBelow(
            onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR,
            tspdif, tspdifPAR, wgts, wc, Ca, latr, lonr, zref,
            maxNrIterations, tolerance, a0, a1, C3
        );
        tground[hr] = onestepin.soilheatvars.Te[0];
        H[hr] = onestepin.H;
        L[hr] = onestepin.L;
        G[hr] = onestepin.soilheatvars.Gflux;
        Ev[hr] = onestepin.Ev;
        surfwater[hr] = onestepin.theta;
        Et[hr] = onestepin.Et;
        iters[hr] = onestepin.iters;
        error[hr] = onestepin.error;
        runoff[hr] = onestepin.runoff;
        drainage[hr] = onestepin.drainage;
        witers[hr] = onestepin.witers;
        if (reqhgt < 0.0) {
            Rdirdown[hr] = Rdifdown[hr] = Rswup[hr] = Rlwdown[hr] = Rlwup[hr] = 0.0;
            tair[hr] = onestepin.soilheatvars.Te[sel];
            tleaf[hr] = -999.99;
            rz[hr] = onestepin.soilwatervars.theta[sel];
            uz[hr] = 0.0;
        }
        else {
            Rdirdown[hr] = onestepin.Rdirdown[sel];
            Rswup[hr] = onestepin.Rswup[sel];
            Rlwup[hr] = onestepin.Rlwup[sel];
            if (reqhgt >= vegpc.hgt) {
                Rdifdown[hr] = Rdif[hr];
                Rlwdown[hr] = Rlw[hr];
                const double Th = onestepin.tcanopy;
                tair[hr] = Tabove(reqhgt, zref, Th, temp[hr], vegpc.hgt, vegpc.pai, onestepin.LL, soilpc.z0, onestepin.fade);
                tleaf[hr] = -999.99;
                rz[hr] = RHabove(reqhgt, zref, onestepin.rhcanopy, Th, temp[hr], tair[hr], relhum[hr], vegpc.hgt, vegpc.pai,
                    onestepin.LL, soilpc.z0, onestepin.fade);
                uz[hr] = Uabove(reqhgt, zref, onestepin.ucanopy * onestepin.windscale, wspeed[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
            }
            else {
                Rdifdown[hr] = onestepin.Rdifdown[sel];
                Rlwdown[hr] = onestepin.Rlwdown[sel];
                if (reqhgt == 0.0) {
                    tair[hr] = onestepin.soilheatvars.Te[0];
                    tleaf[hr] = -999.99;
                    rz[hr] = onestepin.soilwatervars.theta[0];
                    uz[hr] = 0.0;
                }
                else {
                    tair[hr] = onestepin.tair[sel];
                    tleaf[hr] = onestepin.tleaf[sel];
                    leafwet[hr] = onestepin.leafwet[sel];
                    rz[hr] = onestepin.rh[sel];
                    uz[hr] = onestepin.uz[sel] * onestepin.windscale;
                }
            }
        }
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    Rcpp::List out;
    out["year"] = year;
    out["month"] = month;
    out["day"] = day;
    out["hour"] = hour;
    out["Rdirdown"] = Rdirdown;
    out["Rdifdown"] = Rdifdown;
    out["Rswup"] = Rswup;
    out["Rlwdown"] = Rlwdown;
    out["Rlwup"] = Rlwup;
    out["tair"] = tair;
    out["tground"] = tground;
    out["tleaf"] = tleaf;
    out["leafwet"] = leafwet;
    out["relhum"] = rz;
    out["windspeed"] = uz;
    out["H"] = H;
    out["L"] = L;
    out["G"] = G;
    out["Evaporation"] = Ev;
    out["surfwater"] = surfwater;
    out["Transpiration"] = Et;
    out["runoff"] = runoff;
    out["drainage"] = drainage;
    out["iters"] = iters;
    out["witers"] = witers;
    out["error"] = error;
    out.attr("class") = "data.frame";
    out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(n));
    return out;
}
// Converts weather measured at height zin to the equivalent at height zout above
// the canopy. The model is run each hour with the weather at zin; temperature,
// humidity, wind and pressure are then taken along the above-canopy profile to zout.
// [[Rcpp::export]]
List WeatherhgtCpp2(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zin, double zout, double lat, double lon, 
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, double CO2ppm = 430.0,
    int maxNrIterations = 100, double tolerance = 0.1)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    size_t n = temp.size();
    size_t na = vegpc.paii.size();
    const std::vector<double> z = canopyNodeHeights(na, vegpc.hgt);
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.root50, vegpc.root95);
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc, soilpc.z0);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    std::vector<double> temp_new(n);
    std::vector<double> relhum_new(n);
    std::vector<double> windspeed_new(n);
    std::vector<double> pres_new(n);
    double d = zeroplanedisCpp2(vegpc.hgt, vegpc.pai);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
        climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
        climin.winddir = wdir[hr]; climin.precip = precip[hr];
        onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR, tspdif,
            tspdifPAR, wgts, wc, CO2ppm, latr, lonr, zin, maxNrIterations, tolerance);
        const double th = onestepin.tcanopy;
        temp_new[hr] = Tabove(zout, zin, th, temp[hr], vegpc.hgt, vegpc.pai, onestepin.LL, soilpc.z0, onestepin.fade);
        relhum_new[hr] = RHabove(zout, zin, onestepin.rhcanopy, th, temp[hr], temp_new[hr], relhum[hr], vegpc.hgt, vegpc.pai,
            onestepin.LL, soilpc.z0, onestepin.fade);
        // Friction velocity from the wind at zin and the stability correction
        // over the same span, then the wind at zout from the same profile.
        double zm = roughlengthCpp2(vegpc.hgt, vegpc.pai, d, soilpc.z0);
        double psi_in = dpsimCpp2(zm / onestepin.LL) - dpsimCpp2((zin - d) / onestepin.LL);
        double uf = (ka * wspeed[hr]) / (std::log((zin - d) / zm) + psi_in);
        double psi_m = dpsimCpp2(zm / onestepin.LL) - dpsimCpp2((zout - d) / onestepin.LL);
        windspeed_new[hr] = (uf / ka) * (std::log((zout - d) / zm) + psi_m);
        double Tv = ((temp_new[hr] + climin.tref) / 2.0) + 273.15;
        pres_new[hr] = pres[hr] * std::exp(-(g * (zout - zin)) / (287.05 * Tv)); // hypsometric pressure correction
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    return List::create(
        Named("new_temp") = Rcpp::wrap(temp_new),
        Named("new_relhum") = Rcpp::wrap(relhum_new),
        Named("new_windspeed") = Rcpp::wrap(windspeed_new),
        Named("new_pressure") = Rcpp::wrap(pres_new)
    );
}

// Copies a vector into row hr of a matrix.
void FillMyMatrix(NumericMatrix& mat, std::vector<double>& v, size_t hr)
{
    for (size_t i = 0; i < v.size(); ++i) {
        mat(hr, i) = v[i];
    }
}
// Bare ground: the whole time series at every height in z and every soil node.
// [[Rcpp::export]]
List RunBelowFullBare(DataFrame obstime, DataFrame climdata, List soilc, std::vector<double> z,
    double zref, double zm, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double  tolerance = 0.1)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    soilpstruct soilpc = tosoilpstruct(soilc);
    size_t n = temp.size();
    size_t na = z.size();
    onestepbare onestepin;
    onestepin.tair = std::vector<double>(na, temp[0]);
    onestepin.rh = std::vector<double>(na, relhum[0]);
    onestepin.uz = std::vector<double>(na, wspeed[0]);
    onestepin.soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 0.05, 0.4);
    onestepin.H = 0.0; onestepin.L = 0.0; onestepin.psim = 0.0; onestepin.psih = 0.0;
    onestepin.Rswup = 0.0; onestepin.Rb0 = 0.0; onestepin.Rlwup = 0.0;
    onestepin.LL = 999.99; onestepin.iters = 1; onestepin.error = 1e99;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    std::vector<double> runoff(n, 0.0); // rain the saturated surface could not take (mm)
    std::vector<double> drainage(n, 0.0); // water leaving through the base of the soil (mm)
    std::vector<int> witers(n, 0); // iterations of the soil water solve
    std::vector<int> iters(n, 0); // passes the hour took
    std::vector<double> error(n, 0.0); // change in air temperature over the last pass (K)
    // Create Numeric Matrices for above ground
    NumericMatrix tair(n, na);
    NumericMatrix rh(n, na);
    NumericMatrix uz(n, na);
    // Create Numeric Matrices for below ground
    size_t nb = onestepin.soilheatvars.Te.size();
    NumericMatrix tsoil(n, nb);
    NumericMatrix theta(n, nb);
    // Run model for every hr
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
        climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
        climin.winddir = wdir[hr]; climin.precip = precip[hr];
        onestepin = OneStepBare(onestepin, obsdata, climin, soilpc, z, latr, lonr, zref, zm,
            maxNrIterations, tolerance);
        // Extract values
        runoff[hr] = onestepin.runoff;
        drainage[hr] = onestepin.drainage;
        witers[hr] = onestepin.witers;
        iters[hr] = onestepin.iters;
        error[hr] = onestepin.error;
        FillMyMatrix(tair, onestepin.tair, hr);
        FillMyMatrix(rh, onestepin.rh, hr);
        FillMyMatrix(uz, onestepin.uz, hr);
        FillMyMatrix(tsoil, onestepin.soilheatvars.Te, hr);
        FillMyMatrix(theta, onestepin.soilwatervars.theta, hr);
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    return List::create(
        // Microclimate variables
        Named("tair") = Rcpp::wrap(tair),
        Named("relhum") = Rcpp::wrap(rh),
        Named("windspeed") = Rcpp::wrap(uz),
        // Below ground temp and theta
        Named("tsoil") = Rcpp::wrap(tsoil),
        Named("theta") = Rcpp::wrap(theta),
        // Surface runoff, drainage from the base, and the water solve's iterations
        Named("runoff") = Rcpp::wrap(runoff),
        Named("drainage") = Rcpp::wrap(drainage),
        Named("witers") = Rcpp::wrap(witers),
        // Passes each hour took, and the change over the last pass
        Named("iters") = Rcpp::wrap(iters),
        Named("error") = Rcpp::wrap(error)
    );
}
// Vegetated surface: the whole time series at every canopy node and soil node.
// [[Rcpp::export]]
List RunBelowFull(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zref, double Ca, double lat, double lon,
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, int maxNrIterations = 100,
    double  tolerance = 0.1, double a0 = -1.0, double a1 = 1.25, bool C3 = true)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc, soilpc.z0);
    size_t n = temp.size();
    size_t na = paii.size();
    const std::vector<double> z = canopyNodeHeights(na, vegpc.hgt);
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.root50, vegpc.root95);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    // Create Numeric Matrices for above ground
    NumericMatrix Rdirdown(n, na);
    NumericMatrix Rdifdown(n, na);
    NumericMatrix Rswup(n, na);
    NumericMatrix Rlwdown(n, na);
    NumericMatrix Rlwup(n, na);
    NumericMatrix tair(n, na);
    NumericMatrix tleaf(n, na);
    NumericMatrix leafwet(n, na);
    NumericMatrix rh(n, na);
    NumericMatrix uz(n, na);
    std::vector<double> runoff(n, 0.0); // rain the saturated surface could not take (mm)
    std::vector<double> drainage(n, 0.0); // water leaving through the base of the soil (mm)
    std::vector<int> witers(n, 0); // iterations of the soil water solve
    std::vector<int> iters(n, 0); // passes the hour took
    std::vector<double> error(n, 0.0); // change in air temperature over the last pass (K)
    std::vector<double> tcanopy(n), rhcanopy(n), ucanopy(n); // air temperature, humidity and wind at canopy top
    std::vector<double> LL(n), fade(n); // Obukhov length and canopy turbulence weight
    // Create Numeric Matrices for below ground
    size_t nb = onestepin.soilheatvars.Te.size();
    NumericMatrix tsoil(n, nb);
    NumericMatrix theta(n, nb);
    // Run model for every hr
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
        climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
        climin.winddir = wdir[hr]; climin.precip = precip[hr];
        onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR, tspdif,
            tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
        // Extract values
        runoff[hr] = onestepin.runoff;
        drainage[hr] = onestepin.drainage;
        witers[hr] = onestepin.witers;
        iters[hr] = onestepin.iters;
        error[hr] = onestepin.error;
        tcanopy[hr] = onestepin.tcanopy;
        rhcanopy[hr] = onestepin.rhcanopy;
        ucanopy[hr] = onestepin.ucanopy * onestepin.windscale;
        LL[hr] = onestepin.LL;
        fade[hr] = onestepin.fade;
        FillMyMatrix(Rdirdown, onestepin.Rdirdown, hr);
        FillMyMatrix(Rdifdown, onestepin.Rdifdown, hr);
        FillMyMatrix(Rswup, onestepin.Rswup, hr);
        FillMyMatrix(Rlwdown, onestepin.Rlwdown, hr);
        FillMyMatrix(Rlwup, onestepin.Rlwup, hr);
        FillMyMatrix(tair, onestepin.tair, hr);
        FillMyMatrix(tleaf, onestepin.tleaf, hr);
        FillMyMatrix(leafwet, onestepin.leafwet, hr);
        FillMyMatrix(rh, onestepin.rh, hr);
        {
            std::vector<double> uzActual = onestepin.uz; // actual wind, without the convective addition
            for (double& u : uzActual) u *= onestepin.windscale;
            FillMyMatrix(uz, uzActual, hr);
        }
        FillMyMatrix(tsoil, onestepin.soilheatvars.Te, hr);
        FillMyMatrix(theta, onestepin.soilwatervars.theta, hr);
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    List out = List::create(
        // Radiation streams
        Named("Rdirdown") = Rcpp::wrap(Rdirdown),
        Named("Rdifdown") = Rcpp::wrap(Rdifdown),
        Named("Rswup") = Rcpp::wrap(Rswup),
        Named("Rlwdown") = Rcpp::wrap(Rlwdown),
        Named("Rlwup") = Rcpp::wrap(Rlwup),
        // Other above ground climate variables
        Named("tair") = Rcpp::wrap(tair),
        Named("tleaf") = Rcpp::wrap(tleaf),
        Named("relhum") = Rcpp::wrap(rh),
        Named("windspeed") = Rcpp::wrap(uz),
        // Air at canopy top
        Named("tcanopy") = Rcpp::wrap(tcanopy),
        Named("rhcanopy") = Rcpp::wrap(rhcanopy),
        Named("ucanopy") = Rcpp::wrap(ucanopy),
        // Below ground temp and theta
        Named("tsoil") = Rcpp::wrap(tsoil),
        Named("theta") = Rcpp::wrap(theta),
        // Surface runoff, drainage from the base, and the water solve's iterations
        Named("runoff") = Rcpp::wrap(runoff),
        Named("drainage") = Rcpp::wrap(drainage),
        Named("witers") = Rcpp::wrap(witers),
        // Passes each hour took, and the change over the last pass
        Named("iters") = Rcpp::wrap(iters),
        Named("error") = Rcpp::wrap(error)
    );
    out["LL"] = Rcpp::wrap(LL);
    out["fade"] = Rcpp::wrap(fade);
    out["leafwet"] = Rcpp::wrap(leafwet);
    return out;
}
// Air temperature, relative humidity and wind speed at heights za between canopy top and the
// height of the weather data, for every hour of a run.
// [[Rcpp::export]]
List AboveCanopyCpp(std::vector<double> za, double zref, std::vector<double> tcanopy, std::vector<double> rhcanopy,
    std::vector<double> ucanopy, std::vector<double> tref, std::vector<double> rhref, std::vector<double> uref,
    double hgt, double pai, std::vector<double> LL, double z0, std::vector<double> fade)
{
    const size_t n = tcanopy.size(), na = za.size();
    NumericMatrix tair(n, na), rh(n, na), uz(n, na);
    for (size_t hr = 0; hr < n; ++hr) {
        const double eh = satvapCpp2(tcanopy[hr]) * rhcanopy[hr] / 100.0;
        const double eref = satvapCpp2(tref[hr]) * rhref[hr] / 100.0;
        for (size_t i = 0; i < na; ++i) {
            const double f = aboveCanopyFraction(za[i], zref, hgt, pai, LL[hr], z0, fade[hr]);
            tair(hr, i) = tcanopy[hr] + (tref[hr] - tcanopy[hr]) * f;
            rh(hr, i) = std::min(100.0, 100.0 * (eh + (eref - eh) * f) / satvapCpp2(tair(hr, i)));
            uz(hr, i) = Uabove(za[i], zref, ucanopy[hr], uref[hr], hgt, pai, LL[hr]);
        }
    }
    return List::create(Named("tair") = tair, Named("relhum") = rh, Named("windspeed") = uz);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Other useful functions *********************************************** //
// Solar utilities, and the big-leaf runs that give the soil its starting state.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Clear-sky shortwave radiation, allowing for scattering and absorption by
// gases, water vapour and aerosols.
// [[Rcpp::export]]
std::vector<double> clearskyradCpp2(DataFrame obstime, DataFrame climdata, double lat, double lon)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    double latr = lat * torad;
    double lonr = lon * torad;
    size_t n = hour.size();
    std::vector<double> csr(n, 0.0);
    for (size_t hr = 0; hr < n; ++hr) {
        solmodel solp = solpositionCpp2(latr, lonr, year[hr], month[hr], day[hr], hour[hr]);
        if (solp.zenr <= pi / 2.0) {
            double cosz = std::cos(solp.zenr);
            double m = 35.0 / std::sqrt(1224.0 * cosz * cosz + 1.0);
            double TrTpg = 1.021 - 0.084 * std::sqrt(m * 0.00949 * pres[hr] + 0.051);
            double xx = std::log(relhum[hr] / 100.0) + ((17.27 * temp[hr]) / (237.3 + temp[hr]));
            double Td = (237.3 * xx) / (17.27 - xx);
            double u = std::exp(0.1133 - std::log(3.78) + 0.0393 * Td);
            double Tw = 1.0 - 0.077 * std::pow(u * m, 0.3);
            double Ta = pow(0.935, m);
            double od = TrTpg * Tw * Ta;
            csr[hr] = 1352.778 * std::cos(solp.zenr) * od;
        }
    }
    return csr;
}
// Fraction of shortwave radiation that is diffuse, from the clearness index and
// solar altitude.
// [[Rcpp::export]]
std::vector<double> difpropCpp(DataFrame obstime, std::vector<double> swrad, double lat, double lon)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    size_t n = hour.size();
    std::vector<double> dp(n, 1.0);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        solmodel solp = solpositionCpp2(latr, lonr, year[hr], month[hr], day[hr], hour[hr]);
        if (solp.zenr < pi / 2.0) {
            double zd = solp.zenr * 180.0 / pi;
            double k1 = 0.83 - 0.56 * std::exp(-0.06 * (90 - zd));
            double si = std::cos(solp.zenr);
            double k = swrad[hr] / (1352.0 * si);
            if (k > k1) k = k1;
            if (k < 0.0) k = 0.0;
            double rho = k / k1;
            double sigma3 = 0;
            if (rho > 1.04) {
                sigma3 = 0.12 + 0.65 * (rho - 1.04);
            }
            else {
                sigma3 = 0.021 + 0.397 * rho - 0.231 * rho * rho - 0.13 * std::exp(-1.0 * std::pow((rho - 0.931) / 0.134, 2) * 0.834);
            }
            double k2 = 0.95 * k1;
            double d1 = 1.0;
            if (zd < 88.6) d1 = 0.07 + 0.046 * zd / (93 - zd);
            double K = 0.5 * (1.0 + std::sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2.0));
            double d2 = 1 - ((1.0 - d1) * (0.11 * std::sqrt(K) + 0.15 * K + 0.74 * K * K));
            double d3 = (d2 * k2) * (1.0 - k) / (k * (1.0 - k2));
            double alpha = std::pow(1.0 / std::cos(solp.zenr), 0.6);
            double kbmax = std::pow(0.81, alpha);
            double kmax = (kbmax + d2 * k2 / (1.0 - k2)) / (1.0 + d2 * k2 / (1.0 - k2));
            double dmax = (d2 * k2) * (1.0 - kmax) / (kmax * (1.0 - k2));
            dp[hr] = 1.0 - kmax * (1.0 - dmax) / k;
            if (k <= kmax) dp[hr] = d3;
            if (k <= k2) dp[hr] = d2;
            if (k <= 0.22) dp[hr] = 1.0;
            double kX = 0.56 - 0.32 * std::exp(-0.06 * (90 - zd));
            double kL = (k - 0.14) / (kX - 0.14);
            double kR = (k - kX) / 0.71;
            double delta = (k >= 0.14 && k < kX) ? (-3.0 * kL * kL * (1.0 - kL) * std::pow(sigma3, 1.3)) : 0;
            if (k >= kX && k < (kX + 0.71)) delta = 3.0 * kR * std::pow((1 - kR), 2) * std::pow(sigma3, 0.6);
            if (sigma3 > 0.01) dp[hr] = dp[hr] + delta;
        }
    }
    return dp;
}
// Solar altitude (degrees above the horizon).
// [[Rcpp::export]]
std::vector<double> solaltCpp(DataFrame obstime, double lat, double lon)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    size_t n = hour.size();
    std::vector<double> sa(n, 1.0);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        solmodel solp = solpositionCpp2(latr, lonr, year[hr], month[hr], day[hr], hour[hr]);
        double zd = solp.zenr * 180.0 / pi;
        sa[hr] = 90.0 - zd;
    }
    return sa;
}
// Linearly interpolates each column of a matrix to nout rows.
// [[Rcpp::export]]
NumericMatrix expand_outputCpp(const NumericMatrix& mat, int nout) {
    int nin = mat.nrow();
    int ncol = mat.ncol();
    if (nout < 2) stop("nout must be at least 2");
    if (nin < 2) stop("mat must have at least 2 rows");
    NumericMatrix out(nout, ncol);
    double scale = static_cast<double>(nin - 1) / static_cast<double>(nout - 1);
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nout; i++) {
            double x = i * scale;
            int i0 = static_cast<int>(x);
            int i1 = i0 + 1;
            if (i1 >= nin) {
                out(i, j) = mat(nin - 1, j);
            }
            else {
                double w = x - i0;
                out(i, j) = (1.0 - w) * mat(i0, j) + w * mat(i1, j);
            }
        }
    }
    return out;
}
// Runs the big-leaf model through the time series to give soil temperatures and
// water contents for the start of a multilayer run.
// [[Rcpp::export]]
List BigLeafCpp2(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List soilc,
    Rcpp::List vegp, std::vector<double> Lfrac, double zref, double Ca, double lat, double lon,
    double boundaryT, int maxiter = 100, bool C3 = true)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> precip = climdata["precip"];
    double pai = Rcpp::as<double>(vegp["pai"]);
    std::vector<double> paii(Lfrac.size());
    double val = pai / static_cast<double>(Lfrac.size());
    for (size_t i = 0; i < Lfrac.size(); ++i) paii[i] = val;
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    tsvegstruct tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    size_t n = soilpc.thetaS.size();
    std::vector<double> SoilTempIni(n);
    std::vector<double> SoilThetaIni(n);
    // Starting soil profile: temperature from the first air temperature at the
    // surface to the boundary temperature at depth; water content from 60% to 100%
    // of saturation.
    double end1 = boundaryT;
    double end2 = soilpc.thetaS[n - 1];
    double start1 = temp[0];
    double start2 = 0.6 * end2;
    for (size_t i = 0; i < n; ++i) {
        SoilTempIni[i] = start1 + (end1 - start1) * static_cast<double>(i) / static_cast<double>(n - 1);
        SoilThetaIni[i] = start2 + (end2 - start2) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    soilmod soilheat = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwater = toSoilwatermod(soilpc, soilheat, vegpc.root50, vegpc.root95);
    soilwaterout water;
    water.swo = soilwater;
    water.success = false;
    water.iterations = 0; // how many iterations model run for
    water.Evapmmhr = 0.0; // surface evaporation
    size_t tsteps = hour.size();
    NumericVector tcanopy(tsteps);
    NumericVector uf(tsteps);
    NumericVector LL(tsteps);
    NumericVector Et(tsteps);
    NumericVector soilt(tsteps);
    NumericVector theta(tsteps);
    IntegerVector iters(tsteps);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < tsteps; ++hr) {
        obsstruct obsdata;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climstruct climin;
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr]; climin.Rsw = Rsw[hr];
        climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr]; climin.winddir = 0.0;
        climin.precip = precip[hr];
        if (climin.uref < 0.5) climin.uref = 0.5;
        bigleafone blo = solveonestep(obsdata, climin, vegpc, tspveg, tspdif, soilpc,
            soilheat, water, latr, lonr, Ca, zref, maxiter, C3);
        soilheat = blo.soilheat;
        water = blo.soilwater;
        soilheat.oldTe = soilheat.Te;
        water.swo.oldvapor = water.swo.vapor;
        water.swo.oldtheta = water.swo.theta;
        tcanopy[hr] = blo.tcanopy;
        uf[hr] = blo.uf;
        LL[hr] = blo.LL;
        Et[hr] = blo.Et;
        iters[hr] = blo.iters;
        soilt[hr] = soilheat.Te[0];
        theta[hr] = water.swo.theta[0];
    }
    return List::create(
        Named("tcanopy") = tcanopy,
        Named("uf") = uf,
        Named("LL") = LL,
        Named("Et") = Et,
        Named("soilt") = soilt,
        Named("theta") = theta,
        Named("SoilTempIni") = Rcpp::wrap(soilheat.Te),
        Named("SoilThetaIni") = Rcpp::wrap(water.swo.theta)
    );
}
// As BigLeafCpp2, for bare ground.
// [[Rcpp::export]]
List BigLeafBareCpp(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List soilc,
    double zref, double zmr, double lat, double lon,
    double boundaryT, int maxiter = 100, bool C3 = true)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> precip = climdata["precip"];
    soilpstruct soilpc = tosoilpstruct(soilc);
    size_t n = soilpc.thetaS.size();
    std::vector<double> SoilTempIni(n);
    std::vector<double> SoilThetaIni(n);
    // Starting soil profile, as in BigLeafCpp2.
    double end1 = boundaryT;
    double end2 = soilpc.thetaS[n - 1];
    double start1 = temp[0];
    double start2 = 0.6 * end2;
    for (size_t i = 0; i < n; ++i) {
        SoilTempIni[i] = start1 + (end1 - start1) * static_cast<double>(i) / static_cast<double>(n - 1);
        SoilThetaIni[i] = start2 + (end2 - start2) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    soilmod soilheat = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwater = toSoilwatermod(soilpc, soilheat, 0.05, 0.4);
    soilwaterout water;
    water.swo = soilwater;
    water.success = false;
    water.iterations = 0; // how many iterations model run for
    water.Evapmmhr = 0.0; // surface evaporation
    size_t tsteps = hour.size();
    NumericVector uf(tsteps);
    NumericVector LL(tsteps);
    NumericVector Et(tsteps);
    NumericVector soilt(tsteps);
    NumericVector theta(tsteps);
    IntegerVector iters(tsteps);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < tsteps; ++hr) {
        obsstruct obsdata;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climstruct climin;
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr]; climin.Rsw = Rsw[hr];
        climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr]; climin.winddir = 0.0;
        climin.precip = precip[hr];
        if (climin.uref < 0.5) climin.uref = 0.5;
        bigleafone blo = solveonestepbare(obsdata, climin, soilpc, soilheat, water, zmr, latr, lonr, zref, maxiter);
        soilheat = blo.soilheat;
        water = blo.soilwater;
        soilheat.oldTe = soilheat.Te;
        water.swo.oldvapor = water.swo.vapor;
        water.swo.oldtheta = water.swo.theta;
        uf[hr] = blo.uf;
        LL[hr] = blo.LL;
        Et[hr] = blo.Et;
        iters[hr] = blo.iters;
        soilt[hr] = soilheat.Te[0];
        theta[hr] = water.swo.theta[0];
    }
    return List::create(
        Named("iters") = iters,
        Named("uf") = uf,
        Named("LL") = LL,
        Named("Et") = Et,
        Named("soilt") = soilt,
        Named("theta") = theta,
        // Other above ground climate variables
        Named("SoilTempIni") = Rcpp::wrap(soilheat.Te),
        Named("SoilThetaIni") = Rcpp::wrap(water.swo.theta)
    );
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** Examining the hour's iteration ********************************************** //
// Internal functions for examining how an hour is solved. They run one chosen
// hour pass by pass and record each pass.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Everything the chosen hour's iteration starts from.
struct diagcontext {
    onestep start, base;
    outercarry basecarry;
    bool hasBase = false;
    climstruct clim;
    vegpstruct vegp;
    soilpstruct soilp;
    std::vector<double> z, wc, swater0;
    LWweights wgts;
    radmodel swrad;
    double Ca = 0.0, zref = 0.0, a0v = 0.0, a1 = 1.25, dcan = 0.0;
    int maxIter = 100;
    bool C3 = true;
};
static outercarry diagStartCarry(const diagcontext& ctx, const onestep& s)
{
    outercarry c;
    c.tground = s.soilheatvars.Te[0];
    c.Th = s.tcanopy;
    c.eh = satvapCpp2(c.Th) * (s.rhcanopy / 100.0);
    c.psir = rootZonePotential(ctx.soilp, s.soilwatervars.theta, s.soilwatervars.rootfrac, s.soilwatervars.Tc) / 1000.0;
    c.dTs.assign(s.tair.size(), 0.0);
    c.H_iter = s.H;
    c.G_iter = s.soilheatvars.Gflux;
    return c;
}
// Runs the model as usual up to the chosen hour (counted from 0) and keeps
// everything that hour's iteration starts from.
// [[Rcpp::export]]
SEXP DiagSetupCpp(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List soilc,
    Rcpp::List vegp, std::vector<double> paii, std::vector<double> Lfrac, double zref, double Ca,
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, int hourIndex,
    int maxNrIterations = 100, double tolerance = 0.1, double a0 = -1.0, double a1 = 1.25, bool C3 = true)
{
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    diagcontext* ctx = new diagcontext();
    ctx->vegp = tovegpstruct(vegp, paii, Lfrac);
    ctx->soilp = tosoilpstruct(soilc);
    tsvegstruct tspveg = twostreamvegCpp(ctx->vegp.pai, ctx->vegp.x, ctx->vegp.lref, ctx->vegp.ltra, ctx->soilp.gref);
    tsvegstruct tspvegPAR = twostreamvegCpp(ctx->vegp.pai, ctx->vegp.x, ctx->vegp.lrefp, ctx->vegp.ltrap, ctx->soilp.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    ctx->wgts = lwradweights(ctx->vegp.paii);
    ctx->wc = windprofileCpp(ctx->vegp, ctx->soilp.z0);
    const size_t na = paii.size();
    ctx->z = canopyNodeHeights(na, ctx->vegp.hgt);
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(ctx->soilp, soilheatvars, ctx->vegp.root50, ctx->vegp.root95);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    const double latr = lat * torad, lonr = lon * torad;
    auto weather = [&](size_t hr, obsstruct& ob, climstruct& cl) {
        ob.year = year[hr]; ob.month = month[hr]; ob.day = day[hr]; ob.hour = hour[hr];
        cl.tref = temp[hr]; cl.relhum = relhum[hr]; cl.pk = pres[hr];
        cl.Rsw = Rsw[hr]; cl.Rdif = Rdif[hr]; cl.Rlw = Rlw[hr]; cl.uref = wspeed[hr];
        cl.winddir = wdir[hr]; cl.precip = precip[hr];
    };
    for (int hr = 0; hr < hourIndex; ++hr) {
        obsstruct obsdata; climstruct climin;
        weather(static_cast<size_t>(hr), obsdata, climin);
        onestepin = OneStepBelow(onestepin, obsdata, climin, ctx->vegp, ctx->soilp, ctx->z, tspveg, tspvegPAR, tspdif,
            tspdifPAR, ctx->wgts, ctx->wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    obsstruct obsdata;
    weather(static_cast<size_t>(hourIndex), obsdata, ctx->clim);
    // the hour's shortwave, as the hourly step forms it
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double si = solarindexCpp2(ctx->soilp.slope, ctx->soilp.aspect, solp.zenr, solp.azir);
    kstruct kp = cankCpp(solp.zenr, ctx->vegp.x, si);
    tsdirstruct tspdir = twostreamdirCpp(ctx->vegp.pai, kp.kd, ctx->soilp.gref, tspveg);
    tsdirstruct tspdirPAR = twostreamdirCpp(ctx->vegp.pai, kp.kd, ctx->soilp.grefPAR, tspvegPAR);
    ctx->swrad = shortwavemodelCpp(ctx->vegp.pia, ctx->vegp.paii, ctx->vegp.pai, ctx->soilp.gref, ctx->soilp.grefPAR,
        ctx->vegp.lref, ctx->vegp.lrefp, ctx->clim.Rsw, ctx->clim.Rdif, si, solp, kp, tspveg, tspvegPAR, tspdif, tspdifPAR,
        tspdir, tspdirPAR);
    onestepin.Rdirdown = ctx->swrad.Rdirdown;
    onestepin.Rdifdown = ctx->swrad.Rdifdown;
    onestepin.Rswup = ctx->swrad.Rswup;
    ctx->start = onestepin;
    ctx->swater0 = onestepin.swaterdepth;
    ctx->Ca = Ca; ctx->zref = zref; ctx->a1 = a1; ctx->C3 = C3; ctx->maxIter = maxNrIterations;
    ctx->a0v = (a0 > 0.0) ? a0 : floorGustRatio(ctx->vegp, a1);
    ctx->dcan = zeroplanedisCpp2(ctx->vegp.hgt, ctx->vegp.pai);
    return Rcpp::XPtr<diagcontext>(ctx, true);
}
// Runs the chosen hour step by step and records the state after each step, the
// result of each pass, and the change in canopy air. memory is the number of
// past passes combined (the model's own if negative, none if zero). The run can
// start from a given state (xstart).
// [[Rcpp::export]]
Rcpp::List DiagRunCpp(SEXP ctxptr, int memory, int maxPass, double tolerance, Rcpp::NumericVector xstart, bool setBase = true)
{
    Rcpp::XPtr<diagcontext> ctx(ctxptr);
    onestep s = ctx->start;
    outercarry c = diagStartCarry(*ctx, s);
    const size_t n = s.tair.size(), m = 3 * n + 5;
    if (xstart.size() == static_cast<R_xlen_t>(m)) injectState(Rcpp::as<std::vector<double>>(xstart), s, c);
    outermemory mem;
    if (memory >= 0) mem.depth = memory;
    std::vector<std::vector<double>> X, R, D;
    X.push_back(packState(s, c));
    while ((c.pass < 3 || c.tdif > tolerance) && c.pass < maxPass) {
        outerStep(s, c, mem, ctx->clim, ctx->vegp, ctx->soilp, ctx->z, ctx->wgts, ctx->wc, ctx->swrad, ctx->swater0,
            ctx->Ca, ctx->zref, ctx->maxIter, tolerance, ctx->a0v, ctx->a1, ctx->dcan, ctx->C3);
        X.push_back(packState(s, c));
        R.push_back(mem.passResult);
        D.push_back({ c.tdif, c.uf, s.LL, s.L, s.fade, static_cast<double>(s.witers), static_cast<double>(mem.X.size()),
            static_cast<double>(s.soilheatvars.iters), s.soilheatvars.Te[1], s.soilheatvars.Te[2], s.soilheatvars.Te[3],
            s.soilheatvars.oldTe[0], s.soilheatvars.oldTe[1], s.soilheatvars.wc[0] });
    }
    if (setBase) { ctx->base = s; ctx->basecarry = c; ctx->hasBase = true; }
    Rcpp::NumericMatrix Xm(X.size(), m), Rm(R.size(), m), Dm(D.size(), 14);
    for (size_t i = 0; i < X.size(); ++i) for (size_t j = 0; j < m; ++j) Xm(i, j) = X[i][j];
    for (size_t i = 0; i < R.size(); ++i) { for (size_t j = 0; j < m; ++j) Rm(i, j) = R[i][j]; for (size_t j = 0; j < 14; ++j) Dm(i, j) = D[i][j]; }
    Rcpp::colnames(Dm) = Rcpp::CharacterVector::create("tdif", "uf", "LL", "L", "fade", "witers", "remembered",
        "soiliters", "soil1", "soil2", "soil3", "old0", "old1", "wc0");
    return Rcpp::List::create(Rcpp::Named("state") = Xm, Rcpp::Named("target") = Rm, Rcpp::Named("diag") = Dm,
        Rcpp::Named("passes") = c.pass, Rcpp::Named("error") = c.tdif, Rcpp::Named("n") = static_cast<int>(n),
        Rcpp::Named("z") = ctx->z, Rcpp::Named("tref") = ctx->clim.tref, Rcpp::Named("uref") = ctx->clim.uref,
        Rcpp::Named("Rsw") = ctx->clim.Rsw, Rcpp::Named("precip") = ctx->clim.precip, Rcpp::Named("relhum") = ctx->clim.relhum);
}
// One pass applied to a given state.
// [[Rcpp::export]]
Rcpp::NumericVector DiagMapCpp(SEXP ctxptr, Rcpp::NumericVector x)
{
    Rcpp::XPtr<diagcontext> ctx(ctxptr);
    if (!ctx->hasBase) Rcpp::stop("Run DiagRunCpp first");
    onestep s = ctx->base;
    outercarry c = ctx->basecarry;
    injectState(Rcpp::as<std::vector<double>>(x), s, c);
    if (c.pass < 3) c.pass = 3;
    outerPass(s, c, ctx->clim, ctx->vegp, ctx->soilp, ctx->z, ctx->wgts, ctx->wc, ctx->swrad, ctx->swater0,
        ctx->Ca, ctx->zref, ctx->maxIter, OUTER_TOLERANCE, ctx->a0v, ctx->a1, ctx->dcan, ctx->C3);
    return Rcpp::wrap(packState(s, c));
}
// Seconds spent in each part of the vegetated hour since the last reset, with
// the numbers of passes and hours (the first 20 values), then the counts of
// each kind of repeated work (24 more); all zero unless built with MP_PROFILE.
// [[Rcpp::export]]
Rcpp::NumericVector ProfileReadCpp(bool reset = true)
{
    Rcpp::NumericVector out(44);
#ifdef MP_PROFILE
    for (int i = 0; i < 20; ++i) { out[i] = g_profSeconds[i]; if (reset) g_profSeconds[i] = 0.0; }
    for (int i = 0; i < 24; ++i) { out[20 + i] = g_profCount[i]; if (reset) g_profCount[i] = 0.0; }
#endif
    return out;
}
