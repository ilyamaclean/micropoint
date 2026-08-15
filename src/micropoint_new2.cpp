#include <Rcpp.h>
#include "micropointheaders3.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstddef>

using namespace Rcpp;
using namespace newmodel;

constexpr double pi = 3.14159265358979323846;
constexpr double sb = 5.67e-8;
constexpr double ka = 0.41;
constexpr double Mw = 0.018015; // kg/mol
constexpr double RgasC = 8.314; // J/mol/K
constexpr double g = 9.80665;
constexpr double torad = 3.14159265358979323846 / 180.0;
// Forward declaration: defined below (Bigleaf section), needed by several
// earlier functions' Aitken-damped iterates (windmodelCpp's uf_iter,
// SoilHeatCpp's Tsurf_iter, OneStepBelow's H_iter).
inline double aitken1d(double oldv, double newv, Aitken1DState& st);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************** Solar model ***************************************** //
// This section establishes solar geometry only.  It converts date/time and location into
// sun position and local beam incidence; canopy interception itself is handled by the radiation model below.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculates Astronomical Julian day ** //
// Converts a calendar date to Julian day for the solar-geometry calculations.
// This is purely astronomical bookkeeping; the meteorological timestep remains in local clock time.
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
// ** Calculates solar time ** //
// Converts local clock time to apparent solar time using longitude and the equation of time.
// Solar time is then used to place the sun correctly within the day.
static double soltimeCpp(int jd, double lt, double lond)
{
    double m = 6.24004077 + 0.01720197 * (static_cast<double>(jd) - 2451545.0);
    double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2.0 * m + 3.5932);
    double st = lt + (4.0 * lond + eot) / 60.0;
    return st;
}
// Returns solar zenith and azimuth for the requested place, date and local time.
// These angles drive slope illumination, direct-beam extinction and the animal radiation geometry.
static solmodel solpositionCpp2(double latr, double lonr, int year, int month, int day, double lt)
{
    solmodel solpos{};
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lonr * 180.0 / pi);
    // Calculate solar zenith (degrees)
    double tt = 0.261799 * (st - 12.0);
    double dec = (pi * 23.5 / 180.0) * std::cos(2.0 * pi * ((jd - 159.5) / 365.25));
    double coh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
    solpos.zenr = std::acos(coh);
    // Calculate solar azimuth (degrees)
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
// ** Calculates solar index ** //
// Projects the direct solar beam onto the local ground plane.
// The result is zero when the sun is below the horizon or lies behind the slope.
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
// Shortwave radiation is represented with a two-stream canopy model.  Static optical
// coefficients are separated from the per-timestep evaluation so canopy geometry can be reused efficiently.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculate canopy extinction coefficient for sloped ground surfaces** //
// Converts leaf-angle distribution and solar geometry into canopy beam-extinction coefficients.
// k describes extinction normal to the beam; kd adjusts it for the illuminated ground-plane geometry used by the two-stream model.
static kstruct cankCpp(double zenr, double x, double si) {
    double k;
    if (zenr > (pi / 2.0)) zenr = pi / 2.0;
    if (si < 0.0) si = 0.0;
    // Calculate normal canopy extinction coefficient
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
    // Calculate adjusted k
    double kd = k * std::cos(zenr) / si;
    kstruct kparams{};
    kparams.k = k;
    kparams.kd = kd;
    return kparams;
}
// Two-stream radiative transfer (Dickinson 1983 / Sellers 1985): canopy
// optical parameters -- single-scattering albedo (om), backscatter (gma),
// extinction (h) -- from leaf reflectance/transmittance and the leaf-angle
// distribution (x). Feeds twostreamdifCpp/twostreamdirCpp below, which use
// these to build the diffuse and direct-beam flux profiles through the canopy.
// Builds the canopy optical state used by the two-stream shortwave solution.
// Leaf reflectance/transmittance, PAI and leaf-angle distribution are reduced to the scattering, absorption and attenuation coefficients reused below.
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
// Closed-form two-stream coefficients (p1-p4) for the diffuse shortwave
// flux profile through the canopy, from twostreamvegCpp's optical
// parameters and the ground reflectance (gref).
// Solves the boundary coefficients for diffuse shortwave transport through the canopy.
// These coefficients are geometry-independent for a given canopy state and are subsequently evaluated at any cumulative PAI.
static tsdifstruct twostreamdifCpp(tsvegstruct& tsvegp)
{
    tsdifstruct params{};
    params.p1 = (tsvegp.gma / (tsvegp.D1 * tsvegp.S1)) * (tsvegp.u1 - tsvegp.h);
    params.p2 = (-tsvegp.gma * tsvegp.S1 / tsvegp.D1) * (tsvegp.u1 + tsvegp.h);
    params.p3 = (1.0 / (tsvegp.D2 * tsvegp.S1)) * (tsvegp.u2 + tsvegp.h);
    params.p4 = (-tsvegp.S1 / tsvegp.D2) * (tsvegp.u2 - tsvegp.h);
    return params;
}
// Closed-form two-stream coefficients (p5-p10) for the direct-beam
// shortwave flux profile through the canopy, analogous to
// twostreamdifCpp but driven by the beam extinction coefficient (kd)
// instead of diffuse penetration.
// Solves the direct-beam contribution to the two-stream canopy radiation field.
// It combines beam extinction with scattering by foliage and reflection from the ground.
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
    params.p8 = sstr * (tsvegp.a + tsvegp.gma + kd) - tsvegp.gma * ss;
    double v3 = (sstr + tsvegp.gma * gref - (params.p8 / params.sig) * (tsvegp.u2 - kd)) * S2;
    params.p9 = (-1.0 / tsvegp.D2) * ((params.p8 / (params.sig * tsvegp.S1)) * (tsvegp.u2 + tsvegp.h) + v3);
    params.p10 = (1.0 / tsvegp.D2) * (((params.p8 * tsvegp.S1) / params.sig) * (tsvegp.u2 - tsvegp.h) + v3);
    return params;
}
// Shortwave radiation absorbed by the ground, whole canopy and individual
// canopy layers. The two-stream solution is evaluated separately for total
// shortwave and PAR, and separates direct-beam from diffuse radiation.
// Layer outputs distinguish sunlit and shaded leaves for the leaf energy
// balance and photosynthesis calculations.
static radmodel shortwavemodelCpp(const std::vector<double>& pia, double pai, double gref, double grefPAR, double lref, double lrefp, double Rswdown, double Rdif, double si, const solmodel& solp, const kstruct& kp, const tsvegstruct& tspveg, const tsvegstruct& tspvegPAR, const tsdifstruct& tspdif, const tsdifstruct& tspdifPAR, const tsdirstruct& tspdir, const tsdirstruct& tspdirPAR) {
    const int n = static_cast<int>(pia.size());
    radmodel out;
    out.RswLsun.assign(n, 0.0); out.RswLshade.assign(n, 0.0); out.RswLav.assign(n, 0.0);
    out.RPARsun.assign(n, 0.0); out.RPARshade.assign(n, 0.0); out.sunfrac.assign(n, 0.0);
    out.Rdirdown.assign(n, 0.0); out.Rdifdown.assign(n, 0.0); out.Rswup.assign(n, 0.0);
    out.RswGabs = 0.0; out.RswCabs = 0.0;
    if (Rswdown <= 0.0) return out;

    // Separate incoming shortwave into direct-beam and diffuse components.
    // Rbeam is irradiance normal to the solar beam; Rb is its horizontal
    // component. The extraterrestrial-scale cap prevents unrealistically
    // large beam estimates when the sun is close to the horizon.
    const double cosz = std::cos(solp.zenr);
    double Rbeam = (Rswdown - Rdif) / cosz;
    if (Rbeam > 1352.0) Rbeam = 1352.0;
    const double Rb = Rbeam * cosz;

    // Bounds used below for beam-generated diffuse fluxes. These cannot
    // physically exceed the largest reflectance in the vegetation-ground
    // system, calculated separately for total shortwave and PAR.
    double amx = gref; if (amx < lref)amx = lref;
    double amxp = grefPAR; if (amxp < lrefp)amxp = lrefp;

    // kd controls attenuation of the direct beam; h and hp control the
    // diffuse two-stream solution for total shortwave and PAR respectively.
    // sil converts beam irradiance to interception per unit leaf area.
    const double kd = kp.kd, h = tspveg.h, hp = tspvegPAR.h;
    const double sil = kp.k * cosz;

    // Evaluate transmission through the entire canopy to obtain radiation
    // reaching the ground.
    const double exp_kd_pai = std::exp(-kd * pai);
    const double exp_mh_pai = std::exp(-h * pai);
    const double exp_ph_pai = std::exp(h * pai);

    // Fraction of incident diffuse radiation transmitted downward through
    // the canopy according to the diffuse two-stream solution.
    double Rdddg = tspdif.p3 * exp_mh_pai + tspdif.p4 * exp_ph_pai;
    if (Rdddg > 1.0) Rdddg = 1.0;
    if (Rdddg < 0.0) Rdddg = 1.0;

    // Additional downward diffuse radiation generated by scattering of the
    // direct beam within the canopy.
    double Rdbdg = 0.0;
    if (Rb > 0.0) {
        Rdbdg = (tspdir.p8 / tspdir.sig) * exp_kd_pai + tspdir.p9 * exp_mh_pai + tspdir.p10 * exp_ph_pai;
        if (Rdbdg > amx) Rdbdg = amx;
        if (Rdbdg < 0.0) Rdbdg = 0.0;
    }

    // Ground absorption combines the attenuated direct beam (adjusted for
    // slope/aspect through si) with diffuse radiation reaching the ground.
    const double Rdirdowng = Rbeam * exp_kd_pai;
    const double Rdifdowng = Rdddg * Rdif + Rdbdg * Rbeam;
    out.RswGabs = (1.0 - gref) * (Rdifdowng + si * Rdirdowng);

    // Whole-canopy absorption is obtained from the canopy-top reflectance
    // (albedo) of the diffuse and direct-beam two-stream solutions. For the
    // direct component, transmitted beam is treated according to whether it
    // reaches the sloping ground or is intercepted within the canopy.
    const double albd = tspdif.p1 + tspdif.p2;
    const double albb = tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7;
    const double trg = exp_kd_pai;
    if (Rb > 0.0) {
        const double Rbc = (trg * si + (1.0 - trg) * cosz) * Rbeam;
        out.RswCabs = (1.0 - albd) * Rdif + (1.0 - albb) * Rbc;
    }
    else out.RswCabs = (1.0 - albd) * Rdif;

    // Evaluate the radiation field at each canopy layer. pia is cumulative
    // PAI above the layer, so increasing p represents progressively deeper
    // positions within the canopy.
    for (int i = 0; i < n; ++i) {
        const double p = pia[i];
        const double exp_kd = std::exp(-kd * p);
        const double exp_mh = std::exp(-h * p);
        const double exp_ph = std::exp(h * p);

        // Upward and downward diffuse flux produced by diffuse illumination.
        double Rddu = tspdif.p1 * exp_mh + tspdif.p2 * exp_ph;
        if (Rddu > 1.0) Rddu = 1.0;
        if (Rddu < 0.0) Rddu = 0.0;
        double Rddd = tspdif.p3 * exp_mh + tspdif.p4 * exp_ph;
        if (Rddd > 1.0) Rddd = 1.0;
        if (Rddd < 0.0) Rddd = 1.0;

        // Upward and downward diffuse flux generated by scattering of the
        // direct solar beam.
        double Rdbu = 0.0, Rdbd = 0.0;
        if (Rb > 0.0) {
            Rdbu = (tspdir.p5 / -tspdir.sig) * exp_kd + tspdir.p6 * exp_mh + tspdir.p7 * exp_ph;
            if (Rdbu > amx) Rdbu = amx;
            if (Rdbu < 0.0) Rdbu = 0.0;
            Rdbd = (tspdir.p8 / tspdir.sig) * exp_kd + tspdir.p9 * exp_mh + tspdir.p10 * exp_ph;
            if (Rdbd > amx) Rdbd = amx;
            if (Rdbd < 0.0) Rdbd = 0.0;
        }

        // Actual direct, downward-diffuse and upward-diffuse shortwave fluxes
        // at this canopy depth.
        out.Rdirdown[i] = Rbeam * exp_kd;
        out.Rdifdown[i] = Rddd * Rdif + Rdbd * Rb;
        out.Rswup[i] = Rddu * Rdif + Rdbu * Rb;

        // Absorbed radiation per unit leaf area. Shaded leaves receive only
        // the diffuse field; sunlit leaves additionally intercept the direct
        // beam. RswLav is the average over all leaves at this canopy depth.
        out.RswLsun[i] = tspveg.a * 0.5 * (sil * Rbeam + out.Rdifdown[i] + out.Rswup[i]);
        out.RswLshade[i] = tspveg.a * 0.5 * (out.Rdifdown[i] + out.Rswup[i]);
        out.RswLav[i] = tspveg.a * 0.5 * (sil * out.Rdirdown[i] + out.Rdifdown[i] + out.Rswup[i]);

        // Repeat the two-stream calculation using PAR-specific optical
        // properties. These coefficients describe the PAR radiation field
        // used to determine absorbed PAR by sunlit and shaded leaves.
        const double exp_mhp = std::exp(-hp * p);
        const double exp_php = std::exp(hp * p);
        double Rddup = tspdifPAR.p1 * exp_mhp + tspdifPAR.p2 * exp_php;
        if (Rddup > 1.0) Rddup = 1.0;
        if (Rddup < 0.0) Rddup = 0.0;
        double Rdddp = tspdifPAR.p3 * exp_mhp + tspdifPAR.p4 * exp_php;
        if (Rdddp > 1.0) Rdddp = 1.0;
        if (Rdddp < 0.0) Rdddp = 1.0;
        double Rdbup = (tspdirPAR.p5 / -tspdirPAR.sig) * exp_kd + tspdirPAR.p6 * exp_mhp + tspdirPAR.p7 * exp_php;
        if (Rdbup > amxp) Rdbup = amxp;
        if (Rdbup < 0.0) Rdbup = 0.0;
        double Rdbdp = (tspdirPAR.p8 / tspdirPAR.sig) * exp_kd + tspdirPAR.p9 * exp_mhp + tspdirPAR.p10 * exp_php;
        if (Rdbdp > amxp) Rdbdp = amxp;
        if (Rdbdp < 0.0) Rdbdp = 0.0;

        // Absorbed PAR supplied to the photosynthesis/stomatal model.
        // The direct-beam contribution occurs only for the sunlit fraction.
        out.RPARshade[i] = tspvegPAR.a * (out.Rdifdown[i] + out.Rswup[i]);
        out.RPARsun[i] = tspvegPAR.a * (sil * Rbeam + out.Rdifdown[i] + out.Rswup[i]);

        // Beer-Lambert probability that a leaf at this depth is directly
        // illuminated; also provides the weighting between sunlit and shaded
        // leaf calculations elsewhere in the canopy model.
        out.sunfrac[i] = exp_kd;
    }
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************** Longwave radiation model ***************************** //
// Longwave exchange is treated as emission and absorption among sky, ground and canopy layers.
// Geometry-dependent view weights are precomputed; temperatures then determine the changing emitted fluxes.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Longwave view-factor weights between every canopy layer, the ground and
// the sky, from Beer-Lambert attenuation through the cumulative PAI
// between two points, then rescaled per layer so the total view factor
// (weights to other layers + transmission to ground/sky) matches the
// physical two-sided exchange each layer should see.
// Precomputes the geometric coupling among canopy layers, ground and sky for longwave exchange.
// The weights depend only on the vertical PAI distribution, so they can be reused while temperatures and emitted longwave change each timestep.
static LWweights lwradweights(const std::vector<double>& paii) {
    int n = static_cast<int>(paii.size());
    double pait = 0.0;
    for (int i = 0; i < n; ++i) pait += paii[i];
    // Cumulative PAI above (paia) and below (paib) each layer.
    std::vector<double> paia(n);
    std::vector<double> paib(n);
    paia[0] = pait;
    paib[0] = 0.0;
    for (int i = 1; i < n; ++i) {
        paib[i] = paib[i - 1] + paii[i];
        paia[i] = pait - paib[i];
    }
    // Longwave exchange weight between layers i and j: leaf area of j,
    // attenuated by Beer-Lambert transmission through the canopy gap
    // between them.
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
    // Calculate transmission from ground and canopy
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
    // A leaf has two sides, so layer i's total view factor (weights to
    // other layers, ground and sky) should sum to 2; rescale to correct
    // for the exponential-attenuation approximation's discretisation error.
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
    // Same rescaling for the ground's view of the canopy, which is
    // one-sided: weights plus transmission to sky should sum to 1.
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
// Stefan-Boltzmann helper: returns absolute temperature to the fourth power from Celsius input.
static double radem(double tc)
{
    return std::pow(tc + 273.15, 4.0);
}
// Longwave balance for each canopy layer and the ground, from
// lwradweights' view-factor weights: downward flux at a layer is sky
// emission attenuated to it plus emission from layers above; upward is
// ground emission attenuated to it plus emission from layers below;
// absorption is the emissivity-weighted mean of the two (both leaf faces).
// Computes longwave exchange for the current ground and leaf temperatures.
// Sky, ground and canopy emission are propagated with the precomputed view weights to give layer and ground absorption used by their energy balances.
static radmodel2 longwavemodelCpp(const LWweights& wgts, double lwdown, double tground, double groundem, double vegem,
    const std::vector<double>& tleaf)
{
    int n = static_cast<int>(wgts.trh.size());
    std::vector<double> RlwLabs(n);
    std::vector<double> Rlwdown(n);
    std::vector<double> Rlwup(n);
    std::vector<double> leafLW(n);
    for (int i = 0; i < n; ++i) leafLW[i] = vegem * sb * radem(tleaf[i]);
    double groundLW = groundem * sb * radem(tground);
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
        RlwLabs[i] = 0.5 * vegem * (Rlwdown[i] + Rlwup[i]);
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
// Wind is split into an above-canopy MOST solution and a dimensionless within-canopy attenuation profile.
// Sensible heat feeds back on atmospheric stability, so friction velocity and Obukhov length are solved iteratively.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Worker functions **
// Estimates canopy zero-plane displacement from canopy height and total PAI.
// This locates the effective momentum sink used by the above-canopy similarity profiles.
// [[Rcpp::export]]
double zeroplanedisCpp2(double h, double pai)
{
    if (pai < 0.001) pai = 0.001;
    double d = (1.0 - (1.0 - exp(-std::sqrt(7.5 * pai))) / std::sqrt(7.5 * pai)) * h;
    return d;
}
// ** Calculate roughness length ** //
// Follows Raupach (1994) Eq. 4: z0/h = (1-d/h)*exp(-ka*Uh/u* - PsiH), using
// Be (sqrt(0.003+0.1*pai)) as a stand-in for u*/Uh and a fixed
// roughness-sublayer influence constant PsiH = ln(cw)-1+1/cw ~ 0.193
// (Raupach's Eq. 5, cw = 2) -- a canopy-geometry constant, not a function
// of atmospheric stability, so no diabatic correction is passed in here.
// Estimates aerodynamic roughness length from canopy geometry using the Raupach formulation.
// The result sets the momentum-profile origin above the canopy and is bounded only to prevent degenerate profiles.
// [[Rcpp::export]]
double roughlengthCpp2(double h, double pai, double d)
{
    double Be = std::sqrt(0.003 + (0.2 * pai) / 2.0);
    const double PsiH = 0.193;
    double zm = (h - d) * std::exp(-ka / Be - PsiH);
    if (zm < 0.0005) zm = 0.0005;
    // safety check to stop the roughness-sublayer correction reversing profile
    if (zm > (0.9 * (h - d))) zm = 0.9 * (h - d);
    return zm;
}
// ** Calculate molar density of air ** //
// Returns molar density of air from temperature and pressure.
// The model works with molar heat capacity/conductance in several energy-balance calculations.
static double phairCpp(double tc, double pk)
{
    double tk = tc + 273.15;
    double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
    return ph;
}
// ** Calculate specific heat of air at constant pressure ** //
// Returns molar heat capacity of air at constant pressure as a weak function of temperature.
static double cpairCpp(double tc)
{
    double cp = 2e-05 * tc * tc + 0.0002 * tc + 29.119;
    return cp;
}
// Computes the integrated diabatic stability correction for momentum
// (Monin-Obukhov similarity theory) and adjusts the wind profile
// accordingly.
// Integrated Monin-Obukhov stability correction for momentum.
// It modifies the logarithmic wind profile under unstable and stable stratification, with strong-stability behaviour deliberately bounded.
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
// **  Calculate integrated diabatic correction coefficient for heat ** //
// Stable branch tapers the same way as dpsimCpp2 above, but with slope
// 4.7/0.74 (heat) instead of 4.7 (momentum).
// Integrated Monin-Obukhov stability correction for heat.
// This is the heat-transfer counterpart of dpsimCpp2 and enters aerodynamic resistance and temperature/humidity profiles.
static double dpsihCpp2(double ze)
{
    double psih;
    // unstable
    if (ze < 0.0) {
        double y = std::sqrt(1.0 - 9.0 * ze);
        psih = std::log(std::pow((1.0 + y) / 2.0, 2.0));
        if (psih > 3.0) psih = 3.0;
    }
    // stable
    else {
        const double zetaMaxH = 4.0 / (4.7 / 0.74);
        psih = -(4.7 / 0.74) * (zetaMaxH * std::tanh(ze / zetaMaxH));
    }
    return psih;
}
// **  Calculate diabatic influencing factor for heat ** //  
// Local stability factor for heat transfer at a specified z/L.
// It is used in the within-canopy exchange scaling, linking above-canopy stability to canopy turbulent transport.
static double dphihCpp2(double ze)
{
    double phih;
    // unstable
    if (ze < 0.0) {
        double phim = 1.0 / std::pow((1.0 - 16.0 * ze), 0.25);
        phih = std::pow(phim, 2.0);
    }
    // stable
    else {
        phih = 1.0 + ((6.0 * ze) / (1.0 + ze));
    }
    if (phih > 1.5) phih = 1.5;
    if (phih < 0.5) phih = 0.5;
    return phih;
}
// Inverts psi_m(L) on the stable branch: finds the Obukhov length L whose
// stability correction equals a given target. psi_m(L) is 0 at both
// L->0+ and L->infinity with an interior peak, so it is not one-to-one --
// a target beyond the peak has no exact solution and is taken as the
// peak's L (the most stable case the profile can represent).
// Maps a requested stable momentum correction back to an Obukhov length.
// This is needed when the model caps the stable MOST correction: because the bounded correction is not one-to-one in L, the physically milder branch is selected and unattainable targets saturate at the peak.
static double lStableFinalCpp(double psi_m, double zref, double d, double zm)
{
    if (psi_m <= 1e-12) return std::numeric_limits<double>::infinity();

    auto psiOfL = [&](double L) {
        return dpsimCpp2(zm / L) - dpsimCpp2((zref - d) / L);
    };

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

    double lo = Lpeak, hi = std::max(Lpeak * 100.0, 1e6);
    while (psiOfL(hi) > psi_m && hi < 1e12) hi *= 10.0;
    for (int i = 0; i < 80; ++i) {
        double mid = 0.5 * (lo + hi);
        if (psiOfL(mid) > psi_m) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
// Caps the stability correction psi_m(L) at a fraction (beta) of the
// neutral log-law term ln((zref-d)/zm), on both the unstable (L<0) and
// stable (L>=0) branches. Aerodynamic resistance uses ln(...) - psi_m, so
// an unbounded psi_m under strong instability/stability would drive
// resistance toward zero or negative -- an unphysical result.
// Keeps the diagnosed Obukhov length within the range for which the aerodynamic profile remains usable.
// Rather than clipping fluxes afterwards, it limits the corresponding MOST correction and returns a consistent L for subsequent wind and heat calculations.
static double clipMOlength(double L, double zref, double d, double zm, double beta = 0.9)
{
    const double ln_z = std::log((zref - d) / zm);
    const double psim_min = -beta * ln_z;
    const double tol = 1e-4;
    double psim = dpsimCpp2(zm / L) - dpsimCpp2((zref - d) / L);
    if (L < 0.0 && psim < psim_min) {
        // Bisect for the L giving psim = psim_min, between -500 and the original L.
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
        // Compares L itself against a bound, not psi_m(L) against a bound:
        // psi_m(L) is not monotonic on this branch (see lStableFinalCpp),
        // so an over-stable L can still show a small psi_m.
        const double psim_max = beta * ln_z;
        const double L_bound = lStableFinalCpp(psim_max, zref, d, zm);
        if (L < L_bound) {
            return L_bound;
        }
    }
    return L;
}
// ** Calculate scaled wind profile ** //
// Calculate canopy wind profile
// Within-canopy wind attenuation profile (relative to canopy-top speed):
// an exponential decay shaped by each layer's own leaf-area density,
// normalised so the whole-canopy decay matches an attenuation coefficient
// fitted from total PAI (Be/Lc/Lm follow the mixing-length form used
// throughout this section). The bottom 10% of layers switch to a neutral
// log-law profile instead, since the exponential decay becomes
// unrealistic (and can invert) very close to the ground.
// Constructs the dimensionless vertical wind profile inside the canopy.
// The profile is scaled later by canopy-top wind speed, allowing canopy geometry to determine attenuation while atmospheric forcing determines its magnitude.
static std::vector<double> windprofileCpp(const vegpstruct& vegp) {
    int n = static_cast<int>(vegp.paii.size());
    if (n < 10) Rcpp::stop("Wind profile requires at least 10 layers");
    double Be = std::sqrt(0.003 + 0.1 * vegp.pai);
    double a = vegp.pai / vegp.hgt;
    double Lc = std::pow(0.25 * a, -1.0);
    double Lm = 2.0 * std::pow(Be, 3.0) * Lc;
    double at = Be * vegp.hgt / Lm;
    std::vector<double> ati(n);
    double sati = 0;
    for (int i = 0; i < n; ++i) {
        double Bei = std::sqrt(0.003 + 0.1 * vegp.paii[i]);
        double ai = vegp.paii[i] / vegp.hgt;
        double Lci = std::pow(0.25 * ai, -1.0);
        double Lmi = 2.0 * std::pow(Bei, 3.0) * Lci;
        ati[i] = Bei * vegp.hgt / Lmi;
        sati += ati[i];
    }
    for (int i = 0; i < n; ++i) ati[i] = (ati[i] / sati) * at;
    int n2 = std::trunc(n / 10);
    std::vector<double> ui(n, 1.0);
    for (int i = n - 1; i >= n2; --i) {
        ui[i - 1] = ui[i] * (1.0 - ati[i - 1]);
    }
    double zm = vegp.hgt / (20.0 * n2);
    for (int i = 0; i < n2; ++i) {
        double z2 = (i + 1) * vegp.hgt / (10 * n2);
        double uf = (ka * ui[n2 - 1]) / std::log(vegp.hgt / (10.0 * zm));
        ui[i] = (uf / ka) * std::log(z2 / zm);
    }
    return ui;
}
// Friction velocity (uf) and Monin-Obukhov length (LL) above the canopy,
// solved jointly since each depends on the other (uf sets LL via the
// sensible heat flux H, LL's stability correction feeds back into uf).
// Links reference-height wind and sensible heat flux to the canopy wind field.
// It jointly solves friction velocity and Obukhov length, derives canopy-top wind, and scales the within-canopy profile; these quantities control both turbulent transport and leaf boundary layers.
static windmodel windmodelCpp(const std::vector<double>& wc, double uref, double hgt, double pai,
    double zref, double H = 0.0, double tc = 15, double pk = 101.3, int maxiter = 100,
    double a1 = 1.25, double psi_h = 0.0, double psi_m = 0.0, double phi_h = 1.0,
    double zi = 1000.0, double beta = 1.0)
{
    if (zref < hgt) {
        Rcpp::stop("Height of wind speed measurement must be above canopy");
    }
    double Tk = tc + 273.15;
    double cp = cpairCpp(tc);
    double ph = phairCpp(tc, pk);
    double d = zeroplanedisCpp2(hgt, pai);
    double zm = roughlengthCpp2(hgt, pai, d);
    double zh = 0.2 * zm;
    // Under free convection (H > 0), shear alone can't sustain the implied
    // friction velocity: combine uref with a Beljaars (1994) free-convective
    // velocity scale in quadrature so uf stays bounded as uref -> 0. Covers
    // the unstable/daytime case; LangrangianOne's TL floor covers
    // stable/calm-night instead.
    double Ueff = uref;
    if (H > 0.0) {
        double wstar = std::cbrt((g / Tk) * zi * (H / (ph * cp)));
        Ueff = std::sqrt(uref * uref + std::pow(beta * wstar, 2.0));
    }
    double uf = (ka * Ueff) / (std::log((zref - d) / zm) + psi_m);
    // Monin-Obukhov length; stays at the neutral sentinel until the loop
    // below (H != 0) computes a real value.
    double LL = 1e99;
    // Derive diabatic coefficients iteratively
    double dif = 10.0;
    double olduf = -100.00;
    int iter = 1;
    // Aitken-damps uf feeding into LL to stabilise the self-referential
    // uf<->LL coupling. The raw, undamped uf (not uf_iter) is still what's
    // returned in windmodel.uf and used for the convergence check below.
    Aitken1DState st_uf;
    double uf_iter = uf;
    if (H != 0.0) {
        while (dif > 0.00000001) {
            zm = roughlengthCpp2(hgt, pai, d);
            zh = 0.2 * zm;
            LL = (ph * cp * std::pow(uf_iter, 3.0) * Tk) / (-ka * g * H);
            double Lsafe = clipMOlength(LL, zref, d, zm);
            if (H > 0) {
                if (LL < Lsafe) LL = Lsafe;
            }
            else {
                if (LL > Lsafe) LL = Lsafe;
            }
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
    if (H != 0.0) {
        double ln1 = std::log((zref - d) / zm);
        if (psi_m < -0.9 * ln1) psi_m = -0.9 * ln1;
        if (psi_m > 0.9 * ln1) psi_m = 0.9 * ln1;
        psi_m = dpsimCpp2(zm / LL) - dpsimCpp2((hgt - d) / LL);
    }
    double uh = (uf / ka) * (std::log((hgt - d) / zm) + psi_m);
    int n = static_cast<int>(wc.size());
    std::vector<double> uz(n);
    for (int i = 0; i < n; ++i) uz[i] = wc[i] * uh;
    // a2 describes within-canopy exchange (feeds rhcanopy for the
    // ground-to-canopy-top leg), so phi_h reflects stability at canopy
    // top, not at the reference height zref.
    phi_h = dphihCpp2((hgt - d) / LL);
    double a2 = (phi_h * ka * (1.0 - d / hgt)) / (a1 * a1);
    windmodel out;
    out.uz = uz;
    out.LL = LL;
    out.uf = uf;
    out.uh = uh;
    out.a2 = a2;
    out.zm = zm;
    out.psi_m = psi_m;
    out.psi_h = psi_h;
    out.phi_h = phi_h;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************** Plant model ****************************************************** //
// The plant model closes leaf-scale energy, water and carbon exchange within each canopy layer.
// Local radiation and canopy air conditions determine leaf temperature, stomatal conductance and heat/moisture source terms.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculate saturated vapour pressure * * //  
// Saturation vapour pressure over water or ice, used throughout humidity and evaporation calculations.
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
// ** Saturated water vapour density (kg/m^3) at temperature Tk (Kelvin),
// via the ideal gas law applied to satvapCpp2's saturation vapour
// pressure: rho_vs = Mw * es / (RgasC * Tk). * //
// Converts saturation vapour pressure to saturation water-vapour density.
// This is used where the soil/water model represents vapour in mass rather than pressure units.
static double satVapDensityCpp2(double Tk)
{
    double tc = Tk - 273.15;
    double es_Pa = 1000.0 * satvapCpp2(tc); // kPa -> Pa
    return Mw * es_Pa / (RgasC * Tk);
}
// Leaf boundary-layer resistance to heat transfer, from forced and free
// convection Nusselt-number correlations (mixed by cube-summing, Nu =
// (Nuf^3 + Nun^3)^(1/3)) over the leaf's characteristic dimension.
// Boundary-layer resistance between a leaf surface and its surrounding canopy air.
// Forced convection from local wind and free convection from leaf-air temperature difference are combined, making this the leaf-scale aerodynamic link in the energy balance.
static double leafrHa(double tair, double dT, double uz, double len, double wid,
    double x, double rHmax = 300.0)
{
    double Tk = tair + 273.15;
    double Kh = (1.6667e-10 * Tk * Tk + 2.9935e-8 * Tk - 1.7128e-6); // thermal diffusivity
    double v = 1.326 * std::pow(10.0, -5.0) * std::pow(Tk / 273.15, 1.5) * (393.55 / (Tk + 120.0)); // kinematic viscosity
    // Leaf area projected in the wind direction, from the leaf-angle
    // distribution parameter x (same G-function form as cankCpp).
    double y = 1.0 / x;
    double Gy = std::sqrt(y * y + (1.0 - y * y)) / (y + 1.774 * std::pow((y + 1.182), -0.733));
    double d = std::sqrt(len * wid) * Gy; // characteristic dimension
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
// Compute minimum resistances for stomatal model
// Derives the minimum whole-plant hydraulic resistance used by the stomatal optimisation model.
// Plant height and xylem conductivity are translated through the branching/taper assumptions into the hydraulic cost seen by a leaf.
static double rpmin_calc(double h, double hv, double Kxmx)
{
    // Parameters (constants for this model)
    const double n_ext = 2.0;      // number of daughter branches per parent
    const double Lpet = 0.04;     // petiole length (m)
    const double r_intpet = 10.0;  // petiole conduit radius (micrometres)
    const double r_intref = 22.0;  // terminal branch conduit radius (micrometres)
    // ---- Nh calculation (matches: Nh <- (3*log(arg))/log(n_ext); Nh[Nh<1] <- 1 ) ----
    const double n_ext_13 = std::pow(n_ext, 1.0 / 3.0);
    const double arg = 1.0 - (h / Lpet) * (1.0 - n_ext_13);
    double Nh;
    if (arg <= 0.0 || !std::isfinite(arg)) {
        // R would give NaN here; choose a stable fallback consistent with Nh[Nh<1] <- 1
        Nh = 1.0;
        // Alternative (stricter): return NAN to expose the issue
        // return std::numeric_limits<double>::quiet_NaN();
    }
    else {
        Nh = (3.0 * std::log(arg)) / std::log(n_ext);
        if (!std::isfinite(Nh) || Nh < 1.0) Nh = 1.0;  // R clamp
    }
    // ---- Tapering factor ----
    // R: Chi_h <- (6.6e-13*Nh^1.85)/(7.2e-13*Nh^1.32) = (6.6/7.2)*Nh^0.53
    const double Chi_1 = 2.888503;
    const double Chi_h = (6.6e-13 / 7.2e-13) * std::pow(Nh, 0.53);
    const double Chi_tap = Chi_h / Chi_1;
    // ---- Maximum petiole hydraulic conductivity ----
    const double ratio = r_intpet / r_intref;
    const double Kpmx = Kxmx * (ratio * ratio);
    // ---- rpmin ----
    // Guard against zeros to avoid blow-ups
    if (hv <= 0.0 || Kpmx <= 0.0 || Chi_tap <= 0.0 || !std::isfinite(Chi_tap)) {
        return std::numeric_limits<double>::infinity();
    }
    const double rpmin = h / (Kpmx * hv * Chi_tap);
    return rpmin;
}
// Leaf stomatal conductance from the Eller et al. (2020) hydraulic-
// optimisation scheme: photosynthesis (Jacobs 1994 co-limitation of
// light/carbon/transport-limited assimilation, C3 or C4) traded off
// against the whole-plant hydraulic cost of transpiration at the leaf's
// height (z), via the plant's vulnerability curve (psi50-based).
// Predicts stomatal conductance by balancing carbon gain against hydraulic cost.
// Photosynthetic demand, VPD, plant water status and height-dependent hydraulic vulnerability combine to determine how open stomata can profitably remain.
static double leafgs(const envstruct& envdata, vegpstruct& vegp, double z, bool C3 = true)
{
    double gs = 0.0;
    if (envdata.PARabs > 0.0) {
        double IPAR = envdata.PARabs * (0.48 / 0.219) * 1e-6; // PAR converstion to(mol photons / m ^ 2 / s)
        // Vcmax temperature response
        double xx = 0.1 * (envdata.tleaf - 25.0);
        double Vcmax = vegp.Vcmax25 * std::exp2(xx) / ((1.0 + std::exp(0.3 * (envdata.tleaf - vegp.Tup))) *
            (1.0 + std::exp(0.3 * (vegp.Tlw - envdata.tleaf))));
        double Rd = vegp.fd * Vcmax;
        // Photocompensation
        double Q10rs = 0.57;
        double Oa = 0.2095 * envdata.pk * 1000;
        double photocomp = Oa / (2.0 * 2600.0 * std::pow(Q10rs, 0.1 * (envdata.tleaf - 25.0)));
        // Carbon conversion
        double ea = satvapCpp2(envdata.tair) * (envdata.rh / 100.0);
        double es = satvapCpp2(envdata.tleaf);
        double DD = es - ea; // vapour pressure deficit
        double ca = envdata.Ca * envdata.pk / 1000.0;  // convert carbon concentration to Pa
        double ci = vegp.f0 * (1.0 - DD / vegp.Dcrit) * (ca - photocomp); // in Pa
        // Initialize variables
        double Aca;
        double Acol;
        double cicol;
        if (C3) { // C3 photosynthetic pathway
            // Calculate Kc and Ko
            double Q10Kc = 2.1;
            double Kc = 30.0 * std::pow(Q10Kc, 0.1 * (envdata.tair - 25.0));
            double Q10Ko = 1.2;
            double Ko = 30000.0 * std::pow(Q10Ko, 0.1 * (envdata.tair - 25.0));
            // Gross assimilation evaluated at ci (Jacobs 1994): retains
            // the ci-response that shapes the light-response curve,
            // feeding the co-limitation smoothing (Wcol/cicol) below.
            double Wc = Vcmax * ((ci - photocomp) /(ci + Kc * (1 + Oa / Ko))); // Limitating rate due to carbon
            double Wl = vegp.alpha * IPAR * ((ci - photocomp) / (ci + 2.0 * photocomp)); // Limitating rate due to light
            double We = 0.5 * Vcmax; // Limiting rate due to transport
            if (Wc < 0.0) Wc = 0.0;
            if (Wl < 0.0) Wl = 0.0;
            if (We < 0.0) We = 0.0;
            // Co-limiting assimilation, derived from the ci-evaluated We/Wl above.
            double Wcol = ((We + Wl) - std::sqrt(std::pow(We + Wl, 2.0) - 4.0 * 0.93 * (We * Wl)))
                / (2.0 * 0.93);  // Point at which no longer limited by ci
            Acol = Wcol - Rd;
            // Co-limiting CO2 concentration
            cicol = (-Vcmax * photocomp - Kc * (1.0 + Oa / Ko) * Wcol)
                / (Wcol - Vcmax);
            // A(ca): Wc/Wl evaluated at ca instead of ci (Eqn S2.1, Eller
            // et al. 2020, New Phytologist 226:1622-1637), representing
            // assimilation with stomata fully open. Feeds only dadc's
            // numerator below -- substituting into Wcol/cicol would remove
            // Wl's light-response shape and risk a ca==cicol singularity
            // near ambient CO2.
            double Wc_ca = Vcmax * ((ca - photocomp) / (ca + Kc * (1 + Oa / Ko)));
            double Wl_ca = vegp.alpha * IPAR * ((ca - photocomp) / (ca + 2.0 * photocomp));
            double We_ca = We;
            if (Wc_ca < 0.0) Wc_ca = 0.0;
            if (Wl_ca < 0.0) Wl_ca = 0.0;
            double Wca = Wc_ca;
            if (Wl_ca < Wca) Wca = Wl_ca;
            if (We_ca < Wca) Wca = We_ca;
            Aca = Wca - Rd;
        }
        else { // C4 pathway
            double k = 2.0e-4;
            double Wc = Vcmax;
            double Wl = vegp.alpha * IPAR;
            double We = k * Vcmax * (ci / (envdata.pk * 1000.0));
            double W = Wc;
            if (Wl < W) W = Wl;
            if (We < W) W = We;
            Aca = W - Rd; // net assimilation (can be negative)
            // Co - limiting assimilation
            double Wcol = ((Wc + Wl) - std::sqrt(std::pow(Wc + Wl, 2.0) - 4.0 * 0.83 * (We * Wl)))
                / (2.0 * 0.83);
            Acol = Wcol - Rd;
            // Co-limiting CO2 concentration
            cicol = (Wcol * envdata.pk * 1000.0) / (k * Vcmax);
        }
        // dadc: finite difference of assimilation w.r.t. ci between A(ca)
        // and A(ci,col) (Eqn S2.1). Acol is already A at ci,col by
        // construction, so only the A(ca) term is evaluated separately.
        double dadc = (Aca - Acol) / (ca - cicol);
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
        // rpmin depends only on hgt/hv/Kxmx (fixed for a run); cached on
        // vegp, like apsi above, rather than recomputed every call.
        if (vegp.rpmin < 0.0) vegp.rpmin = rpmin_calc(vegp.hgt, vegp.hv, vegp.Kxmx);
        double rp = vegp.rpmin / K_psi_pd;
        double DDm = (es - ea) / envdata.pk; // vapour pressure deficit, mol/mol
        double zeta = 2.0 / (dKdpKi * rp * 1.6 * DDm); // marginal carbon cost of water, Eller optimisation
        if (dadc < 1e-99) dadc = 1e-99;
        double mu = 1.0 + (4.0 * zeta) / dadc;
        if (mu < 1.0) mu = 1.0;
        gs = 0.5 * dadc * (std::sqrt(mu) - 1.0);
        // Apply gsmax cap
        double gsmax = vegp.Vcmax25 * 1e4;
        gs = std::fmin(gs, gsmax);
    }
    return gs;
}
// Compute leaf vapour resistance
// Combines leaf boundary-layer and stomatal resistances into an effective vapour-transfer resistance.
// Wet intercepted water bypasses stomatal control, whereas a dry leaf exchanges vapour through stomata in series with its boundary layer.
static double leafrV(double rHa, double gs, double Lfrac, double ph,
    double surfwater = 0.0, double precip = 0.0)
{
    // compute leaf stomatal resistance
    double rs = 1e9;
    if (gs > 0.0) {
        rs = ph / gs;
    }
    double rVwet = rHa;
    double rVdry = 1e9;
    if (Lfrac > 0.0) rVdry = (rHa + rs) / Lfrac;
    double drywgt = std::exp(-surfwater * 30.0);
    double rV = rVwet;
    if (precip == 0.0) {
        rV = drywgt * rVdry + (1.0 - drywgt) * rVwet;
    }
    return rV;
}
// Solves surface temperature from a Penman-Monteith-style energy balance.
// Given absorbed radiation and heat/vapour resistances, it iterates radiative and latent-heat terms to obtain the temperature at which the surface energy budget closes.
static double PenmanMonteithCpp2(double Rabs, double Ta, double pk, double rh, double em,
    double rHa, double rV, double G = 0.0, int iters = 4)
{
    double Ts = Ta;
    double Rema = em * sb * radem(Ta);
    double cp = cpairCpp(Ta);
    double ph = phairCpp(Ta, pk);
    for (int i = 0; i < iters; ++i) {
        double Te = (Ts + Ta) / 2.0;
        double Rer = 4.0 * em * sb * std::pow(Te + 273.15, 3.0);
        double la;
        if (Ts >= 0) {
            la = 45068.7 - 42.8428 * Ts;
        }
        else {
            la = 51078.69 - 4.338 * Ts - 0.06367 * Ts * Ts;
        }
        double Da = satvapCpp2(Ta) * (1.0 - rh / 100.0);
        double De = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5);
        Ts = Ta + ((Rabs - Rema - ((la * ph) / (pk * rV)) * Da - G) /
            (Rer + ph * (cp / rHa + ((la * De) / (pk * rV)))));
    }
    return (Ts);
}
// Calculate rain transmission
// Calculates how much rainfall penetrates to each canopy depth.
// Wind tilts the rain trajectory, so the same canopy extinction geometry used for radiation is applied along the rain path.
static rainmodel rainintercept(const std::vector<double>& wcm, const std::vector<double>& pia,
    double uh, double rain, double wdir, double x, double sloper, double aspectr)
{
    // calculate rain downward velocity
    double vr = 3.78 * std::pow(rain, 0.067);
    // create output variables
    int n = static_cast<int>(wcm.size());
    std::vector<double> kd(n);
    std::vector<double> raintr(n);
    for (int i = 0; i < n; ++i) {
        // calculate average wind speed above
        double uzm = wcm[i] * uh;
        // calculate rain angle from vertical (in radians)
        double rainZ = std::atan(uzm / vr);
        // calculate rain extinction coefficient
        double si;
        if (sloper == 0.0) {
            si = std::cos(rainZ);
        }
        else {
            si = std::cos(rainZ) * std::cos(sloper) + std::sin(rainZ) *
                std::sin(sloper) * std::cos(wdir - aspectr);
        }
        kstruct kp = cankCpp(rainZ, x, si);
        kd[i] = kp.kd;
        // Calculate tranmission
        raintr[i] = std::exp(-kd[i] * pia[i]);
    }
    rainmodel out;
    out.tr = raintr;
    out.kd = kd;
    return out;
}
// Aerodynamic resistance to sensible/latent heat exchange between the
// ground and height z within the canopy, from the closed-form integral of
// the within-canopy eddy diffusivity profile (K-theory).
// Integrates turbulent resistance from the ground to a specified height inside the canopy.
// It represents the within-canopy leg of sensible/latent transport; the above-canopy leg to zref is added separately.
static double rhcanopy(double a2, double uf, double h, double z)
{
    const double mu = 1.0 / (a2 * h * uf);
    double inth;
    if (z == h) {
        inth = 4.293251 * h;
    }
    else {
        const double invh = 1.0 / h;
        const double x = pi * z * invh;
        const double s = std::sin(x);
        const double c = std::cos(x);
        const double c1 = c + 1.0;
        const double sqrt5 = 2.2360679774997896964;
        const double five32 = 11.180339887498948482; // 5^(3/2)
        const double t = (sqrt5 * s) / c1;
        const double atan_term = std::atan(t);
        const double s2 = s * s;
        const double c1_2 = c1 * c1;
        const double denom = c1 * ((25.0 * s2) / c1_2 + 5.0);
        const double inner = (48.0 * atan_term) / five32 + (32.0 * s) / denom;
        inth = (2.0 * h * inner) / pi;
    }
    double rHa = inth * mu;
    if (rHa < 0.001) rHa = 0.001;
    return rHa;
}
// Calculate rHa from h to zref
// Returns aerodynamic resistance between canopy top and the reference atmosphere.
// Subtracting two MOST resistances isolates the above-canopy path, which is then combined with within-canopy resistance where required.
static double rh_hzref(const windmodel& uzw, double h, double pai, double zref)
{
    // resistance from hgt to zref
    double rhgt_zref = 0.0;
    if (zref > h) {
        double d = zeroplanedisCpp2(h, pai);
        double zh = 0.2 * uzw.zm;
        double psih_z0 = dpsihCpp2(zh / uzw.LL) - dpsihCpp2((zref - d) / uzw.LL);
        double rz0 = (std::log((zref - d) / zh) + psih_z0) / (ka * uzw.uf); // resistance from zref to heat exchange surface
        double psih_h0 = dpsihCpp2(zh / uzw.LL) - dpsihCpp2((h - d) / uzw.LL);
        double rh0 = (std::log((h - d) / zh) + psih_h0) / (ka * uzw.uf); // resistance from hgt to heat exchange surface
        rhgt_zref = rz0 - rh0; // resistance from zref to hgt
    }
    return rhgt_zref;
}
// Calculate total sensible heat flux from canopy elements and ground and resulting heat exchange surface temperature
// Aggregates sensible heat exchange from all canopy elements and the ground.
// The summed flux is converted to the effective canopy heat-exchange surface temperature that feeds back into the above-canopy stability calculation.
static Hstruct sumHCpp(double tref, double tground, double pk, double zref,
    const std::vector<double>& z, const std::vector<double>& tleaf, 
    const std::vector<double>& rz_zref, // resistance from z to zref
    const std::vector<double>& rLB, // Leaf boundary layer resistance
    double rg_zref, // resistance from ground to zref
    const vegpstruct& vegp, const windmodel& uzw)
{
    size_t n = z.size();
    double Htot = 0.0;
    double ph = phairCpp(tref, pk);
    double cp = cpairCpp(tref);
    // Compute flux form individual canopy elements
    for (size_t i = 0; i < n; ++i) {
        double rHa = rz_zref[i] + rLB[i];
        Htot += ((ph * cp) / rHa) * (tleaf[i] - tref) * vegp.paii[i];
    }
    // Compute flux from ground;
    double Hground = ((ph * cp) / rg_zref) * (tground - tref);
    Htot += Hground;
    // Compute heat exchange surface temperature
    double d = zeroplanedisCpp2(vegp.hgt, vegp.pai);
    double zh = 0.2 * uzw.zm;
    double rHa_zref = (std::log((zref - d) / zh) + uzw.psi_h) / (ka * uzw.uf);
    Hstruct out;
    out.Htot = Htot;
    out.Tsurf = tref + (rHa_zref / (ph * cp)) * Htot;
    return out;
}
// Leaf energy balance for every canopy layer: Penman-Monteith temperature
// and evaporation for woody tissue and sunlit/shaded leaves separately
// (stomatal conductance from leafgs), area-weighted into a per-layer
// sensible/latent heat source (Hz/Lz, feeding LangrangianOne) and a
// canopy rain-interception water balance.
// Closes the energy and water balance of every canopy layer for the current atmospheric state.
// Woody, sunlit and shaded leaf fractions are solved separately, then area-weighted to update leaf temperature, stomatal conductance, evaporation/transpiration, sensible/latent source terms and intercepted water.
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
    for (int i = 0; i < n; ++i) {
        rLB[i] = leafrHa(onestepin.tair[i], dTs[i], onestepin.uz[i], vegp.len, vegp.wid, vegp.x);
        envdata.tair = onestepin.tair[i];
        envdata.tleaf = onestepin.tleaf[i];
        envdata.rh = onestepin.rh[i];
        const double ph = phairCpp(onestepin.tair[i], envdata.pk);
        const double Tk = onestepin.tair[i] + 273.15;
        // Reused below across the woody/sunlit/shaded branches (same argument each time).
        const double satvap_tair = satvapCpp2(onestepin.tair[i]);
        // Woody vegetation
        double rV = 9999.99;
        if (onestepin.swaterdepth[i] > 0.0) rV = rLB[i];
        double Rabs = swout.RswLav[i] + lwout.RlwLabs[i];
        const double twood = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        // Vapour pressure deficit: es(Twood) - es(Tair)*RH/100, matching
        // the sunlit/shaded DD formula below (RH applies only to the Tair term).
        double DD = (satvapCpp2(twood) - satvap_tair * (onestepin.rh[i] / 100.0)) * 1000.0;
        const double Evwood = (Mw / (RgasC * Tk)) * (DD / rV) * timestep; // surface water evaporation
        // Sunlit leaves
        envdata.PARabs = swout.RPARsun[i];
        const double gssun = leafgs(envdata, vegp, vegp.hgt / 2.0, C3);
        rV = leafrV(rLB[i], gssun, vegp.Lfrac[i], ph, onestepin.swaterdepth[i], envdata.precip);
        Rabs = swout.RswLsun[i] + lwout.RlwLabs[i];
        const double tsun = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        DD = (satvapCpp2(tsun) - satvap_tair * (onestepin.rh[i] / 100.0)) * 1000.0;
        double rVt = rLB[i] + ph / gssun;
        const double Evsun = (Mw / (RgasC * Tk)) * (DD / rV) * timestep; // surface water evaporation
        const double Etsun = (Mw / (RgasC * Tk)) * (DD / rVt) * timestep; // transpiration
        // Shaded leaves
        envdata.PARabs = swout.RPARshade[i];
        const double gsshade = leafgs(envdata, vegp, vegp.hgt / 2.0, C3);
        rV = leafrV(rLB[i], gsshade, vegp.Lfrac[i], ph, onestepin.swaterdepth[i], envdata.precip);
        Rabs = swout.RswLshade[i] + lwout.RlwLabs[i];
        const double tshade = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        DD = (satvapCpp2(tshade) - satvap_tair * (onestepin.rh[i] / 100.0)) * 1000.0;
        rVt = rLB[i] + ph / gsshade;
        const double Evshade = (Mw / (RgasC * Tk)) * (DD / rV) * timestep; // surface water evaporation
        const double Etshade = (Mw / (RgasC * Tk)) * (DD / rVt) * timestep; // transpiration
        // Perform averaging
        onestepin.gs[i] = swout.sunfrac[i] * gssun + (1.0 - swout.sunfrac[i]) * gsshade; // average stomatal conductance
        const double Evleaf = swout.sunfrac[i] * Evsun + (1.0 - swout.sunfrac[i]) * Evshade; // evaporation from leaves
        Ez[i] = vegp.Lfrac[i] * Evleaf + (1.0 - vegp.Lfrac[i]) * Evwood; // total evaporation
        const double tgreen = swout.sunfrac[i] * tsun + (1.0 - swout.sunfrac[i]) * tshade; // temperature of leaves
        tleafn[i] = vegp.Lfrac[i] * tgreen + (1.0 - vegp.Lfrac[i]) * twood; // temperature of foliage including woody
        Ezt[i] = vegp.Lfrac[i] * vegp.paii[i] * (swout.sunfrac[i] * Etsun + (1.0 - swout.sunfrac[i]) * Etshade);
        double la;
        if (tleafn[i] >= 0.0)
            la = 45068.7 - 42.8428 * tleafn[i];
        else
            la = 51078.69 - 4.338 * tleafn[i] - 0.06367 * tleafn[i] * tleafn[i];
        const double la_Jkg = la / Mw;
        const double cp = cpairCpp(onestepin.tair[i]);
        Lz[i] = la_Jkg * (Ezt[i] / timestep);
        Hz[i] = ((ph * cp) / rLB[i]) * (tleafn[i] - onestepin.tair[i]);
    }
    // Rain interception, top layer down: each layer's surface water depth
    // from throughfall (attenuated by rainvars.tr) plus drip carried over
    // from the layer above, capped at max water film thickness (mwft).
    double dripfrac = 0.0; // fraction of precipitation that drops to lower down
    for (int i = n - 1; i >= 0; --i) {
        const double truetrans = 1.0 - (1.0 - rainvars.tr[i]) * (1.0 - dripfrac);
        const double rainl = truetrans * envdata.precip; // precipitation reaching leaf surface
        if (i == 0) onestepin.precipground = rainl;  
        swaterdepthn[i] = onestepin.swaterdepth[i] + rainl - Ez[i]; // leaf surface water depth
        if (swaterdepthn[i] > vegp.mwft && envdata.precip > 0.0) {
            dripfrac = (swaterdepthn[i] - vegp.mwft) / envdata.precip;
            if (dripfrac > 1.0) dripfrac = 1.0;
        }
        if (swaterdepthn[i] > vegp.mwft) swaterdepthn[i] = vegp.mwft;
        if (swaterdepthn[i] < 0.0) swaterdepthn[i] = 0.0;
    }
    // Calculate total transpiration
    double Et = 0.0;
    for (int i = 0; i < n; ++i) if (Ezt[i] > 0.0) Et += Ezt[i];
    onestepin.tleaf = std::move(tleafn);
    onestepin.swaterdepth = std::move(swaterdepthn);
    onestepin.Et = Et;
    onestepin.Hz = std::move(Hz);
    onestepin.Lz = std::move(Lz);
    onestepin.rLB = std::move(rLB);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Soil model *********************************************************** //
// Soil heat and water are dynamic state variables carried from one timestep to the next.
// The two solvers share temperature/moisture-dependent material properties and meet the atmosphere at the ground energy/water boundary.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Relative humidity of soil pore air from water content, via the Kelvin
// equation applied to the Campbell retention curve's water potential.
// Converts surface-soil water status and temperature to equilibrium relative humidity in the soil pore space.
// This supplies the humidity boundary condition for evaporation from the ground surface.
static double soilrelhumCpp(const soilpstruct& soilp, double Tsoil, double theta)
{
    double psiw = soilp.psie[0] * std::pow(theta / soilp.thetaS[0], -soilp.b[0]);
    double Tk = Tsoil + 273.15;
    double hr = std::exp(Mw * psiw / (RgasC * Tk));
    return hr;
}
// Ground surface energy balance residual (net radiation minus sensible
// and latent heat) -- this residual is what the ground heat flux G is set
// to, elsewhere, to close the surface energy budget.
// Evaluates the residual of the soil-surface energy balance at a trial surface temperature.
// Short/longwave input, sensible heat, latent heat and conduction into the soil are combined so SoilHeatCpp can solve for the temperature that closes the balance.
static double soilsurfaceEB(const soilpstruct& soilp, double Rabs, double Tref,
    double Tsurface, double pk, double relhum, double rHa, double theta)
{
    double sb = 5.67e-8;
    double Rnet = Rabs - soilp.groundem * sb * radem(Tsurface);
    // Sensible heat
    double cp = cpairCpp(Tref);
    double ph = phairCpp(Tref, pk);
    double H = ((ph * cp) / rHa) * (Tsurface - Tref);
    // Latent heat
    double hr = soilrelhumCpp(soilp, Tsurface, theta);
    double es = satvapCpp2(Tsurface) * hr;
    double ea = satvapCpp2(Tref) * relhum / 100;
    double la;
    if (Tsurface >= 0) {
        la = 45068.7 - 42.8428 * Tsurface;
    }
    else {
        la = 51078.69 - 4.338 * Tsurface - 0.06367 * Tsurface * Tsurface;
    }
    double L = ((la * ph) / (rHa * pk)) * (es - ea);
    double Ba = Rnet - H - L;
    return Ba;
}
// Calculate Kmean
// Returns the requested mean of adjacent-layer conductivities for transport across a soil interface.
// Different averaging rules can be selected without changing the heat/water solvers.
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
// Calculates soil thermal conductivity (W/m/K)
// Calculates bulk soil thermal conductivity from mineral composition, porosity and water/ice content.
// This makes heat transport respond to both soil texture and the evolving moisture/freezing state.
static double thermalConductivityCpp(
    double Vq, // Volumetric quartz fraction (m^3/m^3
    double Vm, // Volumetric mineral fraction (m^3/m^3)
    double Vo, // Volumetric organic fraction (m^3/m^3)
    double Vw, // Volmetric water+ice fraction (m^3/m^3)
    double Vw_ice, // Volumetric ICE fraction, subset of Vw (m^3/m^3, 0 if unfrozen)
    double Mc, // Mass fraction of clay (kg/kg
    double Tc, // temperayure of air in air space
    double pk) // pressure (kPa in air space
{
    // Calculate thermal conductivity of solids
    double Vsolid = Vq + Vm + Vo;
    // Quartz conducts heat much better than other minerals (~8.8 vs ~2.5
    // W/m/K); Vq and Vm must each be weighted by their own conductivity,
    // not swapped with each other.
    double kSolid = std::pow(8.8, Vq / Vsolid) * std::pow(2.5, Vm / Vsolid) * std::pow(0.25, Vo / Vsolid);
    // Calculate shape factor
    double Ga = 0.088 * ((Vq + Vm) / Vsolid) + 0.5 * (Vo / Vsolid);
    double Q = 7.25 * Mc + 2.52;
    double xwo = 0.33 * Mc + 0.078;
    double porosity = 1.0 - Vsolid;
    double gasPorosity = porosity - Vw;
    if (gasPorosity < 0.0)  gasPorosity = 0.0;
    double Tk = Tc + 273.15;
    double Lv = 45144.0 - 48.0 * Tc;
    double svp = 0.611 * std::exp(17.502 * Tc / (Tc + 240.97));
    double slope = 17.502 * 240.97 * svp / std::pow(240.97 + Tc, 2);
    double Dv = 0.0000212 * (101.3 / pk) * std::pow(Tk / 273.15, 1.75);
    double rhoAir = 44.65 * (pk / 101.3) * (273.15 / Tk);
    double stcor = 1.0 - svp / pk;
    if (stcor < 0.3) stcor = 0.3;
    double kWater = 0.56 + 0.0018 * Tc;
    // Ice conducts heat ~4x better than liquid water (De Vries); blend the
    // "fluid" phase conductivity by how much of Vw is actually frozen. Pure
    // liquid (Vw_ice=0) reduces exactly to the original kWater.
    double kIce = 2.18 - 0.0068 * Tc;
    double kWaterPhase = (Vw > 1e-9) ? (kWater * (Vw - Vw_ice) + kIce * Vw_ice) / Vw : kWater;
    double wf;
    if (Vw < 0.01 * xwo) {
        wf = 0.0;
    }
    else {
        wf = 1.0 / (1.0 + std::pow(Vw / xwo, -Q));
    }
    double kGas = 0.0242 + 0.00007 * Tc + wf * Lv * rhoAir * Dv * slope / (pk * stcor);
    double Gc = 1.0 - 2.0 * Ga;
    double kFluid = kGas + (kWaterPhase - kGas) * pow(Vw / porosity, 2);
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
// Calculates the effective volumetric heat capacity of the soil mixture.
// Mineral, organic, air, liquid-water and ice fractions contribute sensible heat storage, while partial freezing adds latent heat through an apparent-heat-capacity term.
static double heatCapacityCpp(
    double Vq, // Volumetric quartz fraction (m^3/m^3
    double Vm, // Volumetric mineral fraction (m^3/m^3)
    double Vo, // Volumetric organic fraction (m^3/m^3)
    double Vw, // Volmetric water+ice fraction (m^3/m^3)
    double Vw_ice, // Volumetric ICE fraction, subset of Vw (m^3/m^3, 0 if unfrozen)
    double dtheta_ice_dTc, // |d(ice fraction)/dTc| (m^3 m^-3 K^-1), from the
                            // unfrozen-water curve -- apparent-heat-capacity
                            // latent term; 0 if unfrozen or unused
    double Tc, // temperature of air in air space
    double pk) // pressure (kPa in air space
{
    double Vsum = Vq + Vm + Vo + Vw;
    double Va = 0.0;
    double VwIceNorm = Vw_ice;
    double dThetaIceNorm = dtheta_ice_dTc;
    if (Vsum > 1) {
        Vq = Vq / Vsum;
        Vm = Vm / Vsum;
        Vo = Vo / Vsum;
        Vw = Vw / Vsum;
        VwIceNorm = VwIceNorm / Vsum;
        dThetaIceNorm = dThetaIceNorm / Vsum;
    }
    else {
        Va = 1.0 - Vsum;
    }
    double CH;
    double CHa = cpairCpp(Tc) * phairCpp(Tc, pk) / 1e6; // Specific heat of air in Mj/m^3/K
    if (Tc >= 0 || VwIceNorm <= 0.0) {
        CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Vw * 4.18 + Va * CHa;
    }
    else {
        double CHi = 1.93 + 0.0067 * Tc; // Specific heat of ice
        // Partial liquid/ice split from the unfrozen-water curve, not an
        // all-or-nothing switch at Tc=0.
        double VwLiq = Vw - VwIceNorm;
        if (VwLiq < 0.0) VwLiq = 0.0;
        CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + VwLiq * 4.18 + VwIceNorm * CHi + Va * CHa;
    }
    // Apparent heat capacity: latent heat released/absorbed as the unfrozen
    // fraction changes with Tc (rho_w * Lf * |d(theta_ice)/dTc|). CH is in
    // MJ/m^3/K here, so this J/m^3/K term is divided by 1e6 to match.
    if (dThetaIceNorm > 0.0) {
        constexpr double rho_w = 1000.0; // kg/m^3
        constexpr double Lf = 334000.0;  // J/kg, latent heat of fusion
        CH += (rho_w * Lf * dThetaIceNorm) / 1e6;
    }
    return CH * 1e6; // Convert to J/m^3/K
}
// Thomas algorithm: direct solve of a tridiagonal linear system (aa/bb/cc
// sub/main/super-diagonals, dd right-hand side) between layers first..last.
// Solves the tridiagonal linear system generated by the one-dimensional soil diffusion equations.
// The boundary rows are already assembled by the caller; this routine is only the numerical linear solver.
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
// Calculate degree of saturation from water potential
// Converts soil matric potential to effective saturation using the soil hydraulic retention curve.
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
// Calculate water potential from theta (in J/kg)
// Converts volumetric water content to matric water potential for one soil layer.
static double waterPotential(const soilpstruct& soilp, double theta, int i)
{
    double Se = 1.0;
    if (theta < soilp.thetaS[i]) Se = theta / soilp.thetaS[i];
    double psiw = soilp.psie[i] * std::pow(Se, -soilp.b[i]);
    return psiw;
}
// Calculate theta from water potential
// Inverse retention relation: converts matric water potential back to volumetric water content.
static double thetaFromPsi(const soilpstruct& soilp, double psiw, int i)
{
    double Se = degreeOfSaturation(soilp, psiw, i);
    double theta = Se * soilp.thetaS[i];
    return theta;
}
// Unfrozen (liquid) water content below 0C: freezing-point depression
// treated as extra suction through the same Campbell retention curve as
// thetaFromPsi() (dev-notes/soil_freezing_design_notes.md). Capped at
// theta; equilibrium curve, no hysteresis or persistent ice state.
// Partitions total soil water into the liquid fraction that can remain unfrozen at the current temperature.
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
// Calculate hydraulic conductivity from theta
// Returns unsaturated liquid-water conductivity from current soil water content.
static double hydraulicConductivityFromTheta(const soilpstruct& soilp,
    double theta, int i) {
    double n = 2.0 * soilp.b[i] + 3.0;
    double k = soilp.Ksat[i] * std::pow(theta / soilp.thetaS[i], n);
    return k;
}
// Calculate vapour from water potential
// Returns equilibrium pore-space water-vapour density from soil water potential and temperature.
static double vaporFromPsi(const soilpstruct& soilp, double psiw, double theta, double Tk, int i) {
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vapor = (soilp.thetaS[i] - theta) * satVapDensityCpp2(Tk) * humidity;
    return vapor;
}
// calculate change in theta with psi
// Derivative of the soil water-retention curve, used to linearise Richards-equation updates.
static double dTheta_dPsi(const soilpstruct& soilp, double psiw, int i)
{
    double psie = soilp.psie[i]; // assumed negative
    if (psiw >= psie) return 0.0;
    double Se = std::pow(psiw / psie, -1.0 / soilp.b[i]);
    double theta = soilp.thetaS[i] * Se;
    return -theta / (soilp.b[i] * psiw);
}
// Calculate vapour conductivity from psi and theta
// Effective vapour-phase conductivity through the soil pore space.
// This provides the diffusive vapour component of vertical soil-water transport.
static double vaporConductivityFromPsiTheta(const soilpstruct& soilp,
    double psiw, double theta, double Tk, int i)
{
    double dv = 0.000024;
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vp = satVapDensityCpp2(Tk);
    double k = 0.66 * (soilp.thetaS[i] - theta) * dv * vp * humidity * Mw / (RgasC * Tk);
    return k;
}
// Calculate change in vapour conductivity with psi
// Derivative of equilibrium soil vapour density with respect to water potential for the implicit water solve.
static double dvapor_dPsi(const soilpstruct& soilp, double psiw,
    double theta, double Tk, int i)
{
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vp = satVapDensityCpp2(Tk);
    double capacity_vapor = (soilp.thetaS[i] - theta) * vp * humidity *
        (Mw / (RgasC * Tk)) - dTheta_dPsi(soilp, psiw, i) * vp * humidity;
    return capacity_vapor;
}
// Calculates vapour loss from the soil surface to the atmosphere for the current surface state.
// The flux is controlled jointly by soil vapour availability and the aerodynamic resistance above the ground.
static double evaporation_flux(const soilpstruct& soilp, double theta, double Tsurface,
    double Tair, double relhum, double rHa, double dT)
{
    double psiw = waterPotential(soilp, theta, 0); // J/kg
    double Tk = Tair + 273.15; // Temperature deg C -> K
    double hs = std::exp(Mw * psiw / (RgasC * Tk));   // dimensionless soil effective relative humidity (0 to 1 range)
    double es = 1000.0 * satvapCpp2(Tsurface) * hs;            // kPa -> Pa
    double ea = 1000.0 * satvapCpp2(Tair) * (relhum / 100.0);  // kPa -> Pa
    double Ev = (Mw / (rHa * RgasC * Tk)) * (es - ea) * dT; // Bare soil evaporation (kg m^-2 s^-1 -> mm)
    return Ev;
}
// Calculate transpiration
// Wet-end soil-water stress factor used to reduce root uptake outside the favourable matric-potential range.
static double alpha_wet(double psiw, double psie)
{
    if (psiw >= 0) return(0.0); // saturated
    if (psiw <= psie) return(1.0); // beyond air entry to well aerated
    return std::abs(psiw) / std::abs(psie);
}
// Dry-end soil-water stress factor used to reduce root uptake as soil approaches wilting conditions.
static double alpha_dry(double psiw, double psi_dry, double psi_wilt) {
    if (psiw <= psi_wilt) return(0.0);
    if (psiw >= psi_dry)  return(1.0);
    return (psiw - psi_wilt) / (psi_dry - psi_wilt); // linear ramp 0 to 1
}
// Distribute root water uptake according to root fraction and water potential
// Partitions canopy transpiration demand among soil layers containing roots.
// Root abundance and layer water potential jointly determine where uptake occurs, while preserving the requested whole-profile demand as far as water stress allows.
static std::vector<double> transpiration_distribute(const soilpstruct& soilp, const std::vector<double>& rootfrac,
    double totalTransp_mm, double dT, const std::vector<double>& psiw, double p = 0.5)
{
    int n = static_cast<int>(psiw.size());
    std::vector<double> S(n);
    if (totalTransp_mm > 0.0) {
        std::vector<double> w(n);
        double Trate = totalTransp_mm / dT;
        for (int i = 0; i < n; ++i) {
            double theta_dry = soilp.thetaR[i] + p * (soilp.thetaS[i] - soilp.thetaR[i]);
            double psi_wilt = waterPotential(soilp, soilp.thetaR[i], i);
            double psi_dry = waterPotential(soilp, theta_dry, i);
            double aw = alpha_wet(psiw[i], soilp.psie[i]);
            double ad = alpha_dry(psiw[i], psi_dry, psi_wilt);
            double alpha = aw * ad;
            w[i] = rootfrac[i] * alpha;
        }
        double sumw = std::accumulate(w.begin(), w.end(), 0.0);
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
// Soil layer boundary depths: a thin top layer, then geometrically
// widening with depth (weighted by i^2) down to totalDepth.
// Creates the vertically stretched soil-layer geometry used by both heat and water solvers.
// Thin near-surface layers resolve rapid changes while deeper layers provide the lower boundary with fewer cells.
// [[Rcpp::export]]
std::vector<double> geometricCpp(int n, double totalDepth) {
    std::vector<double> z(n + 2);
    double weightSum = 0.0;
    for (int i = 1; i <= n; ++i) {
        weightSum += static_cast<double>(i) * static_cast<double>(i);
    }
    double dz_unit = totalDepth / weightSum;
    z[0] = 0.0;
    z[1] = dz_unit;  // thin top layer to avoid zero thickness
    for (int i = 2; i <= n; ++i) {
        z[i] = z[i - 1] + dz_unit * static_cast<double>(i) * static_cast<double>(i);
    }
    return z;
}
// Soil temperature profile for one time step: implicit (Thomas-solved)
// heat diffusion through the soil column, with the surface boundary
// condition iterated to convergence against the surface energy balance
// (soilsurfaceEB), and thermal conductivity/heat capacity (including the
// sub-zero freezing terms) evaluated per layer each pass.
// Advances the vertical soil-temperature profile through one timestep.
// Thermal properties respond to water/ice state, while the upper boundary temperature is solved implicitly from the surface energy balance and the lower boundary is prescribed by the soil setup.
static soilmod SoilHeatCpp(soilmod state, const soilpstruct& soilp, double Rabs, double Tref, double relhum, double atmPressure,
    double rHa, double dT = 3600.0, double Fact = 0.5, int maxNrIterations = 100, double tolerance = 1e-2)
{
    int n = state.n;
    double boundaryT = state.oldTe[n];
    std::vector<double> ff(n + 1), CT(n + 1), lambda(n + 1);
    std::vector<double> aa(n + 1), bb(n + 1), cc(n + 1), dd(n + 1);
    const double gg = 1.0 - Fact;
    // FIX: old state for THIS timestep must be constant through the nonlinear iterations
    const std::vector<double> oldTe_fixed = state.oldTe;
    // Start guess (use last solution if you have it)
    std::vector<double> Te_new = state.Te;
    std::vector<double> Te_prev = Te_new;
    std::vector<double> wc = state.wc;
    std::vector<double> dz = state.dz;
    // Extract from soilp
    std::vector<double> Vq = soilp.Vq;
    std::vector<double> Vm = soilp.Vm;
    std::vector<double> Vo = soilp.Vo;
    std::vector<double> Mc = soilp.Mc;
    double max_qsurface = 0.0;
    int nrIterations = 0;
    double maxdT = 1e99;
    double qsurface = 0.0;
    double qsurf1;
    double qsurf2;
    // Aitken-damps the surface-temperature iterate feeding Tav/qsurface
    // below -- this Tav<->qsurface<->tridiagonal-solve loop can diverge
    // undamped when seeded from a too-hot oldTe_fixed[0].
    Aitken1DState st_Tsurf;
    double Tsurf_iter = state.Te[0];
    while (maxdT > tolerance && nrIterations < maxNrIterations) {
        // Compute surface energy balance using midpoint temperature
        // Compute qsurface dynamically if <=10 iterations otherwise hold at average of iter 8 & 9
        if (nrIterations < 10) {
            double Tav = 0.5 * (oldTe_fixed[0] + Tsurf_iter);
            qsurface = soilsurfaceEB(soilp, Rabs, Tref, Tav, atmPressure, relhum, rHa, wc[0]);
            // Limit qsurface to value obtained in first or second iteration
            if (nrIterations < 2 && std::abs(qsurface) > std::abs(max_qsurface)) max_qsurface = std::abs(qsurface);
            if (nrIterations > 1) {
                if (std::abs(qsurface) > std::abs(max_qsurface)) qsurface = std::copysign(std::abs(max_qsurface), qsurface);
            }
            if (nrIterations == 8) qsurf1 = qsurface;
            if (nrIterations == 9) qsurf2 = qsurface;
        }
        else if (nrIterations == 10) qsurface = (qsurf1 + qsurf2) / 2.0;
        for (int i = 0; i <= n; ++i) {
            double Vw_ice = 0.0;
            double dtheta_ice_dTc = 0.0;
            if (Te_new[i] < 0.0) {
                double theta_u = unfrozenTheta(soilp, wc[i], Te_new[i], i);
                Vw_ice = wc[i] - theta_u;
                if (Vw_ice < 0.0) Vw_ice = 0.0;
                // Centred finite difference (not an analytic derivative)
                // to stay consistent with unfrozenTheta() itself.
                const double dTprobe = 0.01;
                double theta_u_lo = unfrozenTheta(soilp, wc[i], Te_new[i] - dTprobe, i);
                double theta_u_hi = unfrozenTheta(soilp, wc[i], Te_new[i] + dTprobe, i);
                dtheta_ice_dTc = std::abs((theta_u_hi - theta_u_lo) / (2.0 * dTprobe));
            }
            lambda[i] = thermalConductivityCpp(Vq[i], Vm[i], Vo[i], wc[i], Vw_ice, Mc[i], Te_new[i], atmPressure);
            CT[i] = heatCapacityCpp(Vq[i], Vm[i], Vo[i], wc[i], Vw_ice, dtheta_ice_dTc, Te_new[i], atmPressure) * state.vol[i];
        }
        // flux coefficients
        ff[0] = kMeanCpp("LOGARITHMIC", lambda[0], lambda[1]) / dz[0];
        for (int i = 1; i < n; ++i) ff[i] = lambda[i] / dz[i];
        // assemble tri-di system using oldTe_fixed ALWAYS
        for (int i = 0; i < n; ++i) {
            if (i == 0) {
                aa[i] = 0.0;
                bb[i] = CT[i] / dT + ff[i];
                cc[i] = -ff[i];
                dd[i] = CT[i] / dT * oldTe_fixed[i] + qsurface;
            }
            else if (i < (n - 1)) {
                aa[i] = -ff[i - 1] * Fact;
                bb[i] = CT[i] / dT + (ff[i - 1] + ff[i]) * Fact;
                cc[i] = -ff[i] * Fact;
                dd[i] = CT[i] / dT * oldTe_fixed[i] + gg * (
                    ff[i - 1] * oldTe_fixed[i - 1] +
                    ff[i] * oldTe_fixed[i + 1] - (ff[i - 1] + ff[i]) * oldTe_fixed[i]);
            }
            else {
                aa[i] = 0.0;
                bb[i] = 1.0;
                cc[i] = 0.0;
                dd[i] = boundaryT;
            }
        }
        Te_prev = Te_new;
        Thomas TBC = ThomasBoundaryCondition(aa, bb, cc, dd, Te_new, 0, n - 1);
        Te_new = TBC.x;
        // monotonic limiter
        for (int i = 1; i < n - 1; ++i) {
            double lo = std::min(Te_new[i - 1], Te_new[i + 1]);
            double hi = std::max(Te_new[i - 1], Te_new[i + 1]);
            if (Te_new[i] < lo) Te_new[i] = lo;
            else if (Te_new[i] > hi) Te_new[i] = hi;
        }
        Tsurf_iter = aitken1d(Tsurf_iter, Te_new[0], st_Tsurf);
        // convergence: max change in the iterate
        maxdT = 0.0;
        for (int i = 0; i <= n; ++i) {
            double d = std::abs(Te_new[i] - Te_prev[i]);
            if (d > maxdT) maxdT = d;
        }
        ++nrIterations;
    }
    // Gflux is not computed here -- callers compute it post-convergence as
    // the Aitken-damped surface energy balance residual (soilsurfaceEB).
    state.Te = Te_new;
    state.iters = nrIterations;
    return state;
}
// distribute roots across soil profile
// Creates the prescribed root-fraction profile across soil layers from total rooting depth and skew.
// The resulting fractions are used to allocate transpiration uptake in SoilWaterCpp.
static std::vector<double> root_distribute(const std::vector<double>& dz, double totalDepth, double skew)
{
    int n = static_cast<int>(dz.size());
    if (skew == 0.0) skew = 1e-12;
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) {
        double z = (i / static_cast<double>(n)) * totalDepth;
        if (i < (n - 1)) {
            v[i] = std::exp(-z * skew) * dz[i];
        }
        else {
            v[i] = std::exp(-z * skew) * dz[i - 1];
        }
    }
    double sumv = std::accumulate(v.begin(), v.end(), 0.0);
    for (int i = 0; i < n; ++i) v[i] = v[i] / sumv;
    return v;
}
// Soil water potential profile for one time step: Richards' equation
// (mixed liquid Darcy flow + vapour diffusion, Campbell retention curve),
// linearised and solved implicitly (Newton-Raphson via a Thomas solve)
// against surface evaporation, transpiration uptake (by root fraction)
// and drainage, with an ice-impedance term reducing liquid conductivity
// where the soil is partly frozen.
// Advances liquid water and water vapour through the soil profile for one timestep.
// An implicit Richards-type solve handles vertical redistribution and phase-coupled vapour transport, while rainfall/evaporation set the surface boundary and transpiration removes water through the root profile.
static soilwaterout SoilWaterCpp(soilwatermod soilmod, const soilpstruct& soilp,
    const climforwaterstruct& climdata, double dT = 3600.0, double pTAW = 0.5,
    int maxNrIterations = 100, double tolerance = 1e-4, bool useDamping = true)
{
    const double rho = 1000.0;
    const int n = soilp.nLayers;
    // Freeze previous-timestep state
    const std::vector<double> oldtheta = soilmod.oldtheta;
    const std::vector<double> oldvapor = soilmod.oldvapor;
    // Working variables
    std::vector<double> psiw = soilmod.psiw;
    std::vector<double> theta = soilmod.theta;
    std::vector<double> vapor = soilmod.vapor;
    std::vector<double> k(n + 1, 0.0);
    std::vector<double> aa(n, 0.0), bb(n, 0.0), cc(n, 0.0), dd(n, 0.0);
    std::vector<double> ff(n, 0.0), u(n, 0.0), du(n, 0.0), Ca(n, 0.0), dpsi(n, 0.0);
    // Bottom boundary initialisation if required
    if (soilp.FreeDrain == false) {
        psiw[n] = soilp.psie[n - 1];
        theta[n] = soilp.thetaS[n - 1];
        k[n] = soilp.Ksat[n - 1];
    }
    int iter = 1;
    double massBalance = 1.0;
    double Evapmmhr = 0.0;
    while (massBalance > tolerance && iter < maxNrIterations) {
        // Surface flux from current state
        Evapmmhr = evaporation_flux(soilp, theta[0], soilmod.Tc[0], climdata.Tair,
            climdata.relhum, climdata.rHa, dT);
        double surfaceFlux = (Evapmmhr - climdata.precip) / dT;
        // Transpiration outside Newton update, but based on current state
        std::vector<double> STr = transpiration_distribute(soilp, soilmod.rootfrac,
            climdata.Et, dT, psiw, pTAW);

        if (soilp.FreeDrain == false) {
            psiw[n] = soilp.psie[n - 1];
            theta[n] = soilp.thetaS[n - 1];
            k[n] = soilp.Ksat[n - 1];
        }
        else {
            psiw[n] = psiw[n - 1];
            theta[n] = theta[n - 1];
            k[n] = k[n - 1];
        }
        // Hydraulic properties
        for (int i = 0; i < n; ++i) {
            double Tkelvin = soilmod.Tc[i] + 273.15;
            // Ice is immobile: liquid (Darcy) conductivity uses only the
            // unfrozen water content, reduced by an impedance factor for
            // the remaining ice (Lundin 1990 / Zhao et al. form). Vapour
            // conductivity stays on total theta -- ice blocks vapour
            // diffusion like liquid water does.
            double theta_u = unfrozenTheta(soilp, theta[i], soilmod.Tc[i], i);
            double theta_ice = theta[i] - theta_u;
            if (theta_ice < 0.0) theta_ice = 0.0;
            double kh = hydraulicConductivityFromTheta(soilp, theta_u, i);
            if (theta_ice > 0.0) {
                constexpr double Omega = 7.0;
                kh *= std::pow(10.0, -Omega * theta_ice);
            }
            double kv = vaporConductivityFromPsiTheta(soilp, psiw[i], theta[i], Tkelvin, i);
            k[i] = kh + kv;
            u[i] = g * k[i];
            // Campbell conductivity exponent, matching
            // hydraulicConductivityFromTheta exactly (this must stay in
            // sync with that function's own local n).
            double nCampbell = 2.0 * soilp.b[i] + 3.0;
            du[i] = -u[i] * nCampbell / psiw[i];
            double Cw = dTheta_dPsi(soilp, psiw[i], i);
            double Cv = dvapor_dPsi(soilp, psiw[i], theta[i], Tkelvin, i);
            Ca[i] = soilmod.vol[i] * (rho * Cw + Cv) / dT;
        }
        if (!soilp.FreeDrain == false) {
            k[n] = k[n - 1];
        }
        // Flux term
        for (int i = 0; i < n; ++i) {
            double nCampbell = 2.0 * soilp.b[i] + 3.0;
            ff[i] = ((psiw[i + 1] * k[i + 1] - psiw[i] * k[i]) /
                (soilmod.dz[i] * (1.0 - nCampbell))) - u[i];
        }
        // Assemble system
        massBalance = 0.0;
        aa[0] = 0.0;
        cc[0] = -k[1] / soilmod.dz[0];
        bb[0] = k[0] / soilmod.dz[0] + Ca[0] + du[0];
        dd[0] = surfaceFlux + STr[0] - ff[0]
            + soilmod.vol[0] * (rho * (theta[0] - oldtheta[0]) + (vapor[0] - oldvapor[0])) / dT;
        massBalance += std::abs(dd[0]);

        for (int i = 1; i < n; ++i) {
            aa[i] = -k[i - 1] / soilmod.dz[i - 1] - du[i - 1];
            cc[i] = -k[i + 1] / soilmod.dz[i];
            bb[i] = k[i] / soilmod.dz[i - 1] + k[i] / soilmod.dz[i] + Ca[i] + du[i];
            dd[i] = ff[i - 1] + STr[i] - ff[i]
                + soilmod.vol[i] * (rho * (theta[i] - oldtheta[i]) + (vapor[i] - oldvapor[i])) / dT;
            massBalance += std::abs(dd[i]);
        }
        // Solve tridiagonal system
        Thomas TBC = ThomasBoundaryCondition(aa, bb, cc, dd, dpsi, 0, n - 1);
        dpsi = TBC.x;
        // Update state
        for (int i = 0; i < n; ++i) {
            double lambda = 1.0;

            if (useDamping) {
                double dryness = (psiw[i] - soilp.psie[i]) / (soilp.psi_min[i] - soilp.psie[i]);
                if (dryness < 0.0) dryness = 0.0;
                if (dryness > 1.0) dryness = 1.0;
                lambda = 0.2 + 0.8 * dryness;
            }

            psiw[i] = psiw[i] - lambda * dpsi[i];
            psiw[i] = std::min(psiw[i], soilp.psie[i] - 1e-8);
            psiw[i] = std::max(psiw[i], soilp.psi_min[i]);

            theta[i] = thetaFromPsi(soilp, psiw[i], i);
            vapor[i] = vaporFromPsi(soilp, psiw[i], theta[i], soilmod.Tc[i] + 273.15, i);
        }

        ++iter;
    }
    soilmod.psiw = psiw;
    soilmod.theta = theta;
    soilmod.vapor = vapor;
    soilmod.k = k;
    soilwaterout out;
    out.swo = soilmod;
    out.success = (massBalance < tolerance);
    out.iterations = iter;
    out.Evapmmhr = Evapmmhr;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ******************************************* Below-canopy Langrangian model ********************************************** //
// Distributed canopy sources do not directly prescribe air temperature/humidity.  This section
// uses a Lagrangian dispersion representation of within-canopy turbulence to convert heat and vapour sources into vertical profiles.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Below-canopy air temperature and humidity at each layer, from a
// localised near-field/far-field Lagrangian dispersion solution (Raupach
// 1989): each layer's air state is the sum of a far-field contribution
// (diffusive, from the cumulative source/sink strength of every layer
// below or above it) and a near-field correction (non-diffusive, from the
// small number of nearby source layers, via the kn() dispersion kernel).
// Source/sink strengths (Hz/Lz, sensible/latent heat) come from the plant
// and soil models; called repeatedly as OneStepBelow iterates the whole
// coupled system to convergence.
// Translates distributed canopy/ground heat and moisture sources into vertical air temperature and humidity profiles.
// The Lagrangian dispersion kernel represents turbulent transport within the canopy; the resulting profiles provide the atmospheric state seen by leaves and the soil at the next coupling iteration.
static void LangrangianOne(onestep& onestepin, double pk, double tground, double soilrelhum,
    double th, double eh, const vegpstruct& vegp,
    const std::vector<double>& z, const windmodel& windvars,
    double a0 = 0.25, double a1 = 1.25)
{
    auto& tair = onestepin.tair;
    const auto& tleaf = onestepin.tleaf;
    auto& rh = onestepin.rh;
    const size_t nn = rh.size();
    const double nnd = static_cast<double>(nn);
    // Floor uf for TL only (not elsewhere): near-calm wind sends TL, and
    // the near-field kernel's zeta, to a log singularity. See
    // dev-notes/multilayer_TL_blowup_handover.md.
    constexpr double UF_MIN_FOR_TL = 0.3;
    const double uf_for_TL = std::max(windvars.uf, UF_MIN_FOR_TL);
    const double TL = windvars.a2 * vegp.hgt / uf_for_TL; // Lagrangian time scale
    const double dz = vegp.hgt / nnd;
    // Physically plausible temperature/vapour-pressure range for this
    // step: between ground and canopy-top values, widened to also cover
    // every leaf's own temperature/saturation vapour pressure. Solved-for
    // values are clamped to this range below as a safety net.
    double tmn = std::min(tground, th) - 2.0, tmx = std::max(tground, th) + 2.0;
    for (size_t i = 0; i < nn; ++i) { tmn = std::min(tmn, tleaf[i]); tmx = std::max(tmx, tleaf[i]); }
    double esg = satvapCpp2(tground) * soilrelhum;  // ground vapour pressure
    double emn = std::min(eh, esg);
    double emx = std::max(eh, esg);
    for (size_t i = 0; i < nn; ++i) {
        double el = satvapCpp2(tleaf[i]);
        emn = std::min(emn, el);
        emx = std::max(emx, el);
    }
    // ow: standard deviation of vertical velocity (sigma_w) at each layer,
    // following Raupach's cosine profile between a0*uf near the ground and
    // a1*uf at canopy top. KH = sigma_w^2*TL is the equivalent far-field
    // (K-theory) eddy diffusivity; inowTL is 1/(sigma_w*TL), the length
    // scale used to non-dimensionalise distances in the near-field kernel
    // below. ST/SL are each layer's sensible/latent heat source strength
    // (per-layer flux Hz/Lz weighted by that layer's foliage area).
    std::vector<double> ow(nn), inowTL(nn), KH(nn), ST(nn), SL(nn), ST_over_ow(nn), SL_over_ow(nn);
    const double mu1 = (a1 + a0) * 0.5 * windvars.uf;
    const double mu2 = (a1 - a0) * 0.5 * windvars.uf;
    for (size_t i = 0; i < nn; ++i) {
        ow[i] = mu1 + mu2 * std::cos(pi * (1.0 - z[i] / vegp.hgt));
        KH[i] = TL * ow[i] * ow[i];
        inowTL[i] = 1.0 / (ow[i] * TL);
        ST[i] = vegp.paii[i] * onestepin.Hz[i];
        SL[i] = vegp.paii[i] * onestepin.Lz[i];
        ST_over_ow[i] = ST[i] / ow[i];
        SL_over_ow[i] = SL[i] / ow[i];
    }
    // Finite-canopy correction to the near-field kernel integral below,
    // larger when the canopy is resolved with fewer layers (Raupach's
    // localised near-field theory).
    const double mu = 1.0 + 0.894 * std::exp(-0.01386 * nnd) + 9.82 * std::exp(-0.15 * nnd);
    // Near-field concentration at canopy top, used as the upper boundary
    // condition below (subtracted back out of each layer's own near-field
    // term, since that term already includes the canopy-top contribution).
    double CnTh = 0.0;
    double CnLh = 0.0;
    for (size_t i = 0; i < (nn - 1); ++i) {
        double dz1 = (vegp.hgt - z[i]) * inowTL[i];
        double dz2 = (vegp.hgt + z[i]) * inowTL[i];
        double e = std::exp(-dz1);
        double kn = -0.39894 * std::log(1.0 - e) - 0.15623 * e; // near-field dispersion kernel
        double common = kn * (dz1 + dz2);
        CnTh += ST[i] / ow[i] * common;
        CnLh += SL[i] / ow[i] * common;
    }
    CnTh *= mu;
    CnLh *= mu;
    // Far-field concentration at canopy top: the above-canopy air
    // temperature/vapour pressure expressed in the same flux-like units
    // (enthalpy/latent-heat-flux equivalent) as the source terms below.
    double lah = (th < 0.0) ? (51078.69 - 4.338 * th - 0.06367 * th * th)
        : (45068.7 - 42.8428 * th);
    const double phh = phairCpp(th, pk);
    const double CfTh = phh * cpairCpp(th) * th;
    const double CfLh = eh * phh * lah / pk;
    // For each layer: combine far-field (diffusive, integrated over all
    // other layers' + ground's source strength) and near-field
    // (non-diffusive, from the localised kernel) contributions to get the
    // layer's total scalar concentration, then convert back to air
    // temperature and vapour pressure/RH.
    double sumRH = 0.0;
    for (size_t i = 0; i < nn; ++i) {
        // Integrated far-field resistance from the ground up to this
        // layer, floored to avoid an unrealistically small value very
        // close to the ground.
        double RH = 1.0 / KH[i];
        sumRH += RH;
        double rHa = sumRH * dz;
        if (rHa < 2.0) rHa = 2.0;
        // Ground-to-air sensible/latent heat exchange at this layer.
        double ph = phairCpp(tair[i], pk);
        double cp = cpairCpp(tair[i]);
        double GT = (ph * cp / rHa) * (tground - tair[i]) * dz;
        double ea = satvapCpp2(tair[i]) * (rh[i] / 100.0);
        double la = (tground < 0.0) ? (51078.69 - 4.338 * tground - 0.06367 * tground * tground)
            : (45068.7 - 42.8428 * tground);
        double GL = (la / (rHa * pk)) * (esg - ea) * dz;
        // Total sensible/latent source strength from the ground and every
        // foliage layer up to and including this one.
        double H = 0.0;
        double L = 0.0;
        for (size_t j = 0; j <= i; ++j) {
            H += ST[j];
            L += SL[j];
        }
        H += GT;
        L += GL;
        // Far-field (diffusive) contribution: integrate that source
        // strength over eddy diffusivity from this layer to canopy top.
        double CfT = 0.0;
        double CfL = 0.0;
        for (size_t j = i; j < nn; ++j) {
            CfT += (H / KH[j]) * dz;
            CfL += (L / KH[j]) * dz;
        }
        // Near-field (non-diffusive) contribution: the localised kernel's
        // response to every other layer's own source strength.
        double CnT = 0.0;
        double CnL = 0.0;
        for (size_t j = 0; j < nn; ++j) if (i != j) {
            double dz1 = (z[i] - z[j]) * inowTL[j];
            double dz2 = (z[i] + z[j]) * inowTL[j];
            double Zeta = std::abs(dz1);
            double e = std::exp(-Zeta);
            double kn = -0.39894 * std::log(1.0 - e) - 0.15623 * e;
            double common = kn * (dz1 + dz2);
            CnT += ST_over_ow[j] * common;
            CnL += SL_over_ow[j] * common;
        }
        // Raupach's near-field/far-field superposition: canopy-top
        // boundary value, replacing its own near-field term with this
        // layer's, plus this layer's far-field contribution.
        double CT = CfTh - CnTh + CfT + CnT * mu;
        double CL = CfLh - CnLh + CfL + CnL * mu;
        tair[i] = CT / (cp * ph);
        double la_i = (tleaf[i] < 0.0) ? (51078.69 - 4.338 * tleaf[i] - 0.06367 * tleaf[i] * tleaf[i])
            : (45068.7 - 42.8428 * tleaf[i]);
        double ean = (CL * pk) / (la_i * ph);
        if (tair[i] > tmx) tair[i] = tmx;
        if (tair[i] < tmn) tair[i] = tmn;
        if (ean > emx) ean = emx;
        if (ean < emn) ean = emn;
        rh[i] = (ean / satvapCpp2(tair[i])) * 100.0;
        if (rh[i] > 100.0) rh[i] = 100.0;
        if (rh[i] < 20.0)  rh[i] = 20.0;
    }
}
// Applies adaptive Aitken under-relaxation to a vertical profile during nonlinear coupling.
// This damps oscillatory updates while allowing well-behaved parts of the profile to converge rapidly.
static inline void aitkin_weightdif(
    const std::vector<double>& oldv,   // unchanged
    std::vector<double>& newv,         // updated in place
    const std::vector<double>& z,
    double hgt,
    WAitkenState& st,
    double omega_min = 0.02,
    double omega_max = 0.90,
    double w_bot = 0.05,
    double w_top = 0.80
)
{
    const size_t n = oldv.size();
    if (st.r_prev.size() != n) {
        st.r_prev.assign(n, 0.0);
        st.have_prev = false;
        st.omega = std::max(omega_min, std::min(st.omega, omega_max));
    }
    const double dw = (w_top - w_bot);
    // ---- First iteration ----
    if (!st.have_prev) {
        double omega = std::max(omega_min, std::min(st.omega, omega_max));
        for (size_t i = 0; i < n; ++i) {
            double r = newv[i] - oldv[i];
            st.r_prev[i] = r;
            double s = z[i] / hgt;
            double wz = w_bot + dw * (s * s);
            newv[i] = oldv[i] + (omega * wz) * r;
        }
        st.omega = omega;
        st.have_prev = true;
        return;
    }
    // ---- Learn omega (weighted Aitken emphasising bottom) ----
    constexpr double beta = 10.0;
    double num = 0.0;
    double den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double r = newv[i] - oldv[i];
        double dr = r - st.r_prev[i];
        double s = z[i] / hgt;
        double t = 1.0 - s;
        double g = 1.0 + beta * (t * t);  // bottom emphasis
        num += g * st.r_prev[i] * dr;
        den += g * dr * dr;
    }
    double omega = st.omega;
    if (den > 0.0)
        omega = -st.omega * (num / den);
    omega = std::max(omega_min, std::min(omega, omega_max));
    // ---- Apply update to newv ----
    for (size_t i = 0; i < n; ++i) {
        double r = newv[i] - oldv[i];
        st.r_prev[i] = r;
        double s = z[i] / hgt;
        double wz = w_bot + dw * (s * s);
        newv[i] = oldv[i] + (omega * wz) * r;
    }
    st.omega = omega;
}
// Canopy-top air temperature/vapour pressure at the reference height
// (zref), solved to convergence against the combined ground and canopy
// source strength (same flux-to-concentration approach as LangrangianOne,
// collapsed to a single layer since only the canopy-top value is needed).
// Closes the coupling between the multilayer canopy and the atmosphere above it.
// From current distributed heat/moisture sources it updates canopy-top conditions, stability, wind and transfer resistances so that the within- and above-canopy solutions share a consistent flux state.
static cantop canopytop(vegpstruct& vegpc, windmodel& wind, climstruct climdata,
    std::vector<double>& Hz, std::vector<double>& Lz, double zref,
    double Th, double eh, double tground, double soilrh,
    double rH_g, double rH_h_zref, int maxIter, double tolerance)
{
    size_t nb = vegpc.paii.size();
    const double d = zeroplanedisCpp2(vegpc.hgt, vegpc.pai);
    const double zm = roughlengthCpp2(vegpc.hgt, vegpc.pai, d);
    const double zh = 0.2 * zm;
    const double psih_h = dpsihCpp2(zm / wind.LL) - dpsihCpp2((vegpc.hgt - d) / wind.LL);
    const double rH_h = (std::log((vegpc.hgt - d) / zh) + psih_h) / (ka * wind.uf); // canopy top to h
    double FcH = 0.0;
    double FcL = 0.0;
    for (size_t i = 0; i < nb; ++i) {
        FcH += Hz[i] * vegpc.paii[i];
        FcL += Lz[i] * vegpc.paii[i];
    }
    const double eground = satvapCpp2(tground) * soilrh;
    const double eref = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    double err = 1e99;
    int nrIterations = 0;
    while (err > tolerance && nrIterations < maxIter) {
        double ph = phairCpp(Th, 101.3);
        double cp = cpairCpp(Th);
        double la;
        if (tground >= 0) {
            la = 45068.7 - 42.8428 * tground;
        }
        else {
            la = 51078.69 - 4.338 * tground - 0.06367 * tground * tground;
        }
        double FgH = ((ph * cp) / (rH_g + rH_h)) * (tground - Th); // ground sensible heat flux
        double FgL = ((la * ph) / ((rH_g + rH_h) * climdata.pk)) * (eground - eh); // ground latent heat flux
        double FzH = FcH + FgH;
        double FzL = FcL + FgL;
        double CfT = FzH * rH_h_zref;
        double CfL = FzL * rH_h_zref;
        double Th_new = CfT / (ph * cp) + climdata.tref;
        double eh_new = (CfL * climdata.pk) / (ph * la) + eref;
        err = std::abs(Th_new - Th);
        double err2 = std::abs(eh_new - eh);
        if (err2 > err) err = err2;
        Th = Th_new;
        eh = eh_new;
        ++nrIterations;
    }
    // ensure vapour pressure not > 100% relative humidity
    double rh = (eh / satvapCpp2(Th));
    if (rh > 1.0) eh = satvapCpp2(Th);
    cantop out;
    out.Th = Th;
    out.eh = eh;
    return out;
 }
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************** Run model for one step ********************************************** //
// These are the core timestep couplers.  A timestep is not a simple sequence: radiation, leaves,
// canopy air, ground, soil and above-canopy stability are iterated until their shared fluxes and boundary conditions agree.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// One vegetated timestep, fully coupled: radiation (fixed for the step) is
// computed once; wind, rain interception, the plant energy balance, the
// soil heat/water solves and the Lagrangian dispersion model are then
// iterated together until air temperature stops changing between passes
// (each depends on the others' previous-pass output -- e.g. wind depends
// on H, which depends on the plant model, which depends on wind's
// resistances). At least 3 passes always run, since dTs (leaf-air
// temperature difference, used to seed the plant model's convection
// terms) is only updated over the first 3 and held fixed after, for
// stability.
// Solves one complete timestep of the multilayer vegetated microclimate model.
// Radiation, canopy energy/water balance, turbulent profiles, soil heat/water and above-canopy stability are iterated because each supplies boundary conditions or fluxes to the others; convergence produces the timestep state passed forward in time.
static onestep OneStepBelow(onestep onestepin, const obsstruct& obsdata, const climstruct& climdata, vegpstruct& vegpc, const soilpstruct& soilpc,
    const std::vector<double>& z, const tsvegstruct& tspveg, const tsvegstruct& tspvegPAR, const tsdifstruct& tspdif,
    const tsdifstruct& tspdifPAR, const LWweights& wgts, const std::vector<double>& wc,
    double Ca, double latr, double lonr, double zref, int maxIter = 100, double  tolerance = 1e-3,
    double a0 = 0.25, double a1 = 1.25, bool C3 = true)
{
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double si = solarindexCpp2(soilpc.slope, soilpc.aspect, solp.zenr, solp.azir);
    kstruct kp = cankCpp(solp.zenr, vegpc.x, si);
    tsdirstruct tspdir = twostreamdirCpp(vegpc.pai, kp.kd, soilpc.gref, tspveg);
    tsdirstruct tspdirPAR = twostreamdirCpp(vegpc.pai, kp.kd, soilpc.grefPAR, tspvegPAR);
    radmodel swrad = shortwavemodelCpp(vegpc.pia, vegpc.pai, soilpc.gref, soilpc.grefPAR, vegpc.lref, vegpc.lrefp,
        climdata.Rsw, climdata.Rdif, si, solp, kp, tspveg, tspvegPAR, tspdif, tspdifPAR, tspdir, tspdirPAR);
    onestepin.Rdirdown = swrad.Rdirdown;
    onestepin.Rdifdown = swrad.Rdifdown;
    onestepin.Rswup = swrad.Rswup;
    double tdif = 1e99;
    int nrIterations = 0;
    size_t na = onestepin.tair.size();
    double tground = onestepin.soilheatvars.Te[0];
    double Th = onestepin.tair[na - 1];
    double rcanh = onestepin.rh[na - 1];
    double eh = satvapCpp2(Th) * (rcanh / 100.0);
    double psir = 0.0; // root-zone mean water potential
    size_t nb = onestepin.soilwatervars.rootfrac.size();
    for (size_t i = 0; i < nb; ++i) psir += onestepin.soilwatervars.rootfrac[i] * onestepin.soilwatervars.psiw[i];
    psir = psir / 1000.0; // J/kg -> MPa
    WAitkenState st_tair;
    WAitkenState st_rh;
    WAitkenState st_leaf;
    std::vector<double> dTs(na); // leaf-air temperature difference, held fixed after iteration 3 (see function doc comment)
    std::vector<double> rz_zref(na); // resistance from each layer to zref
    // Aitken-damps H feeding into windmodelCpp, breaking the self-referential
    // H -> wind/stability -> new resistances -> new H loop. The raw,
    // undamped H is still what's stored in onestepin.H and used elsewhere.
    Aitken1DState st_H;
    double H_iter = onestepin.H;
    // G (ground heat flux): Aitken-damped surface energy balance residual (soilsurfaceEB).
    Aitken1DState st_G;
    double G_iter = onestepin.soilheatvars.Gflux;
    while ((nrIterations < 3 || tdif > tolerance) && nrIterations < maxIter) {
        std::vector<double> oldTe_fixed = onestepin.soilheatvars.oldTe; // previous timestep's soil state, not touched by this pass's iteration
        radmodel2 lwrad = longwavemodelCpp(wgts, climdata.Rlw, tground, soilpc.groundem, vegpc.vegem, onestepin.tleaf);
        onestepin.Rlwdown = lwrad.Rlwdown;
        onestepin.Rlwup = lwrad.Rlwup;
        windmodel wind = windmodelCpp(wc, climdata.uref, vegpc.hgt, vegpc.pai, zref, H_iter, climdata.tref, climdata.pk, maxIter, a1,
            onestepin.psih, onestepin.psim, onestepin.phih);
        onestepin.uz = wind.uz;
        onestepin.psih = wind.psi_h;
        onestepin.psim = wind.psi_m;
        onestepin.phih = wind.phi_h;
        onestepin.LL = wind.LL;
        rainmodel rainvars = rainintercept(wc, vegpc.pia, wind.uh, climdata.precip, climdata.winddir, vegpc.x, soilpc.slope, soilpc.aspect);
        envstruct envdata{};
        envdata.pk = climdata.pk;
        envdata.psi_r = psir;
        envdata.Ca = Ca;
        envdata.precip = climdata.precip;
        if (nrIterations < 3) {
            for (size_t i = 0; i < na; ++i) dTs[i] = std::abs(onestepin.tleaf[i] - onestepin.tair[i]);
        }
        std::vector<double> tleaf = onestepin.tleaf;
        double rhg = rhcanopy(wind.a2, wind.uf, vegpc.hgt, vegpc.hgt); // rHa from ground to top of canopy
        double rhz = rh_hzref(wind, vegpc.hgt, vegpc.pai, zref); // rHa from top of canopy to zref
        double rHa = rhg + rhz; // resistance from ground to zref
        for (size_t i = 0; i < na; ++i) {
            rz_zref[i] = rhg - rhcanopy(wind.a2, wind.uf, vegpc.hgt, z[i]) + rhz;
        }
        plantmodelCpp(onestepin, envdata, vegpc, rainvars, swrad, lwrad, z, dTs, 3600.0, C3); // updates onestepin in place
        aitkin_weightdif(tleaf, onestepin.tleaf, z, vegpc.hgt, st_leaf);
        double Rabs = swrad.RswGabs + lwrad.RlwGabs;
        std::vector<double> stemp = onestepin.soilheatvars.Te;
        soilmod soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600, 0.5, maxIter);
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = onestepin.precipground; cfw.Et = onestepin.Et;
        onestepin.soilwatervars.oldTc = soilheat.oldTe;
        onestepin.soilwatervars.Tc = soilheat.Te;
        soilwaterout soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, 0.5, maxIter, tolerance);
        onestepin.witers = soilwater.iterations;
        tground = soilheat.Te[0];
        double soilrh = soilrelhumCpp(soilpc, tground, soilwater.swo.theta[0]);
        soilheat.wc = soilwater.swo.theta;
        double G_raw = soilsurfaceEB(soilpc, Rabs, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rHa, soilwater.swo.theta[0]);
        G_iter = aitken1d(G_iter, G_raw, st_G);
        soilheat.Gflux = G_iter;
        onestepin.soilheatvars = soilheat;
        onestepin.soilwatervars = soilwater.swo;
        onestepin.Ev = soilwater.Evapmmhr;
        onestepin.theta = soilwater.swo.theta[0];
        Hstruct HT = sumHCpp(climdata.tref, tground, climdata.pk, zref, z, onestepin.tleaf, rz_zref, onestepin.rLB, rHa, vegpc, wind);
        onestepin.H = HT.Htot;
        H_iter = aitken1d(H_iter, onestepin.H, st_H);
        onestepin.L = swrad.RswCabs - vegpc.vegem * sb * radem(HT.Tsurf) - onestepin.H - soilheat.Gflux;
        if (zref > vegpc.hgt) {
            cantop Theh = canopytop(vegpc, wind, climdata, onestepin.Hz, onestepin.Lz, zref, Th, eh, tground, soilrh,
                rhg, rhz, maxIter, tolerance);
            Th = Theh.Th;
            eh = Theh.eh;
        }
        else {
            Th = climdata.tref;
            eh = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
        }
        std::vector<double> tair = onestepin.tair;
        std::vector<double> rh = onestepin.rh;
        LangrangianOne(onestepin, climdata.pk, tground, soilrh, Th, eh, vegpc, z, wind, a0, a1); // updates onestepin.tair/rh in place
        aitkin_weightdif(tair, onestepin.tair, z, vegpc.hgt, st_tair);
        aitkin_weightdif(rh, onestepin.rh, z, vegpc.hgt, st_rh);
        tdif = 0.0; // max air-temperature change this pass, drives the convergence check above
        for (size_t i = 0; i < na; ++i) {
            double dif = std::abs(tair[i] - onestepin.tair[i]);
            if (dif > tdif) tdif = dif;
        }
        psir = 0.0;
        for (size_t i = 0; i < nb; ++i) psir += soilwater.swo.rootfrac[i] * soilwater.swo.psiw[i];
        psir = psir / 1000.0; // J/kg -> MPa
        onestepin.soilheatvars.oldTe = oldTe_fixed; // don't let SoilHeatCpp's own state carry the last pass's oldTe into the next pass
        ++nrIterations;
    }
    onestepin.iters = nrIterations;
    onestepin.error = tdif;
    return onestepin;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************************** Above canopy ************************************************************** //
// Above the canopy, temperature, humidity and wind are diagnosed with Monin-Obukhov similarity
// from the reference atmosphere and the canopy-top state produced by the coupled model.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Extrapolates air temperature above canopy from the Langrangian profile
// via a diabatically-corrected log profile between (hgt, th) and (zref,
// tref). za can exceed zref but must stay below h; LL defaults to 1e99 (neutral).
// Extrapolates air temperature from the reference level to another height above the canopy using the MOST heat profile.
// [[Rcpp::export]]
double Tabove(double za, double zref, double th, double tref, double hgt, double pai, double LL = 1e99)
{
    double d = 0.0;
    if (hgt > 0.0) d = zeroplanedisCpp2(hgt, pai);
    double num = std::log((za - d) / (hgt - d)) + dpsihCpp2((hgt - d) / LL) - dpsihCpp2((za - d) / LL);
    double den = std::log((zref - d) / (hgt - d)) + dpsihCpp2((hgt - d) / LL) - dpsihCpp2((zref - d) / LL);
    double Tz = th + (tref - th) * (num / den);
    return Tz;
}
// Derive humidity above canopy by extrapolating Langrangian profile. za can be > zref, but must be less than h
// See Tabove's doc comment immediately above -- same update, same rationale.
// Extrapolates humidity above the canopy consistently with the temperature and MOST heat-transfer profile.
// Vapour pressure rather than relative humidity is transported, then converted back at the target temperature.
// [[Rcpp::export]]
double RHabove(double za, double zref, double rh, double th, double tref, double tz, double relhum, double hgt, double pai, double LL = 1e99)
{
    double d = 0.0;
    if (hgt > 0.0) d = zeroplanedisCpp2(hgt, pai);
    const double eh = satvapCpp2(th) * (rh / 100.0);
    const double eref = satvapCpp2(tref) * (relhum / 100.0);
    double num = std::log((za - d) / (hgt - d)) + dpsihCpp2((hgt - d) / LL) - dpsihCpp2((za - d) / LL);
    double den = std::log((zref - d) / (hgt - d)) + dpsihCpp2((hgt - d) / LL) - dpsihCpp2((zref - d) / LL);
    const double ez = eh + (eref - eh) * (num / den);
    double rz = (ez / satvapCpp2(tz)) * 100.0;
    if (rz > 100.0) rz = 100.0;
    return rz;
}
// Uses dpsimCpp2 -- see Tabove's doc comment above; now smoothly tapered,
// so no separate shape form is needed here either.
// Extrapolates wind speed between heights above the canopy using the MOST momentum profile.
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
// One bare-ground timestep: the single-layer counterpart of OneStepBelow
// with no canopy, no plant model and no Lagrangian dispersion -- wind
// stability and the soil heat/water solves are iterated together until
// the ground surface temperature stops changing, then the converged
// state is used for a one-shot post-convergence profile at each height z.
// Solves one timestep for a bare surface using the same atmospheric and soil physics but no canopy layers.
// Ground radiation/energy balance, soil heat/water and atmospheric stability are iterated to a mutually consistent surface state.
static onestepbare OneStepBare(onestepbare onestepin, const obsstruct& obsdata, const climstruct& climdata, const soilpstruct& soilpc,
    const std::vector<double>& z, double latr, double lonr, double zref, double zm = 0.004, int maxIter = 100,
    double  tolerance = 1e-8)
{
    double Rswabs = 0.0;
    double Rb0 = 0.0;
    if (climdata.Rsw > 0.0) {
        solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
        if (solp.zenr < pi / 2.0) {
            Rb0 = (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr);
        }
        double si = solarindexCpp2(soilpc.slope, soilpc.aspect, solp.zenr, solp.azir);
        if (si < 0.0) si = 0.0;
        Rswabs = (1.0 - soilpc.gref) * (climdata.Rdif + si * Rb0);
    }
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
    double zmd = zm * std::exp(ka * psi_h);
    double zh = 0.2 * zmd;
    double uf = (ka * climdata.uref) / (std::log(zref / zmd) + psi_m);
    // Monin-Obukhov length; recomputed fresh from the loop's own live LL each pass below.
    double LL = (cpph * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
    double dif = 1e99;
    int nrIterations = 0;
    soilmod soilheat;
    soilwaterout soilwater;
    // Seeded from the pre-loop uf/zh/psi_h above, same formula the loop
    // itself uses (line below) -- keeps this well-defined even in the
    // degenerate maxIter=0 case, matching how uf/LL are already seeded
    // before the loop.
    double rHa = (std::log(zref / zh) + psi_h) / (ka * uf);
    // G: surface energy balance residual, Aitken-damped across outer passes.
    Aitken1DState st_G;
    double G_iter = onestepin.soilheatvars.Gflux;
    while (dif > tolerance && nrIterations < maxIter) {
        if (H != 0.0) {
            LL = (cpph * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
            double Lsafe = clipMOlength(LL, zref, 0.0, zmd);
            if (H > 0) {
                if (LL < Lsafe) LL = Lsafe;
            }
            else {
                if (LL > Lsafe) LL = Lsafe;
            }
            psi_m = dpsimCpp2(zmd / LL) - dpsimCpp2(zref / LL);
            psi_h = dpsihCpp2(zh / LL) - dpsihCpp2(zref / LL);
            zmd = zm * std::exp(ka * psi_h);
            
        }
        else {
            psi_m = 0.0;
            psi_h = 0.0;
            zmd = zm;
        }
        zh = 0.2 * zmd;
        double uf = (ka * climdata.uref) / (std::log(zref / zmd) + psi_m);
        rHa = (std::log(zref / zh) + psi_h) / (ka * uf);
        soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600, 0.5, maxIter);
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = climdata.precip; cfw.Et = 0.0;
        onestepin.soilwatervars.oldTc = soilheat.oldTe;
        onestepin.soilwatervars.Tc = soilheat.Te;
        soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, 0.5, maxIter, tolerance);
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
    double ea = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    double la = (Ts < 0.0) ? (51078.69 - 4.338 * Ts - 0.06367 * Ts * Ts) : (45068.7 - 42.8428 * Ts);
    onestepin.L = ((la * ph) / (climdata.pk * rHa)) * (eg - ea);
    // One-shot, post-convergence height sweep: turns the converged
    // LL/uf/H/L into profile values (tair, rh, uz) at each height z[i].
    size_t n = z.size();
    std::vector<double> uz(n); // wind speed
    for (size_t i = 0; i < n; ++i) {
        if (z[i] > zmd) {
            double psimz = dpsimCpp2(zmd / LL) - dpsimCpp2(z[i] / LL);
            uz[i] = (uf / ka) * (std::log(z[i] / zmd) + psimz);
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
            double rHz = (std::log(z[i] / zh) + psihz) / (ka * uf);
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
    onestepin.Rswup = (1.0 - soilpc.gref) * climdata.Rsw;
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
// The big-leaf model is a deliberately cheaper representation used to generate realistic soil
// initial conditions.  It retains surface/soil/stability feedbacks but avoids the multilayer canopy-air solve.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Canopy albedo (white-sky/diffuse and black-sky/direct-beam) and ground
// weighting factors (direct/diffuse transmission to the ground), from the
// two-stream coefficients -- the single-layer analogue of shortwavemodelCpp's
// per-layer albedo terms.
// Combines direct and diffuse two-stream reflectance into the effective canopy albedo for the current beam geometry.
albf albedo(double kd, double pai, const tsdifstruct& tspdif, const tsdirstruct& tspdir)
{
    albf out;
    out.albd = tspdif.p1 + tspdif.p2; // white sky albedo
    out.albb = tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7; // black sky albedo
    out.wgtg = std::exp(-kd * pai); // ground weighting(direct)
    out.wgti = std::exp(-pai); // ground weighting(diffuse)
    return out;
}
// Bulk (whole-canopy) surface resistance: leafgs' sunlit/shaded stomatal
// conductance, area-weighted by the Beer-Lambert-derived sunlit fraction,
// converted to a single resistance for the single-layer big-leaf model.
// Collapses the leaf stomatal model to a canopy-scale resistance for the big-leaf spin-up model.
// It integrates the same physiological controls approximately rather than resolving individual canopy layers.
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
    double gs_sun = leafgs(envdata, vegp, z, C3);
    envdata.PARabs = Rshade_abs;
    double gs_shade = leafgs(envdata, vegp, z, C3);
    double Gs = gs_sun * Lsun + gs_shade * Lshade;
    // Convert to bulk surface resistance
    double ph = phairCpp(tair, pk);
    return ph / Gs;
}
// Calculate dewpoint temperature
// Converts air temperature and relative humidity to dew-point temperature for humidity bookkeeping in the big-leaf model.
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
// Aitken-damped update for a self-referential fixed-point iteration (e.g.
// H feeding wind, which feeds resistances, which feeds H back): blends the
// raw new value with the old one by a weight (omega) that adapts each call
// from how the residual (r = newv-oldv) is shrinking, damping oscillation
// without needing a fixed damping constant tuned per use site.
// Scalar Aitken relaxation used to stabilise self-coupled iterative variables such as friction velocity and surface temperature.
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
// One timestep of the single-layer ("big leaf") model: canopy treated as
// one bulk surface (bulkstomatalresist) rather than resolved per-layer --
// used only to spin up soil heat/water state before a full multilayer
// run, not for the main model output. Iterates wind stability, canopy
// temperature and the surface energy balance together to convergence.
// Solves one vegetated big-leaf timestep used for inexpensive soil spin-up.
// Canopy radiation, bulk transpiration, ground/canopy energy balance, soil state and atmospheric stability are iterated together, but without the full multilayer canopy-air calculation.
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
    double sloper = soilp.slope * torad;
    double aspectr = soilp.aspect * torad;
    double si = solarindexCpp2(sloper, aspectr, solp.zenr, solp.azir);
    kstruct kk = cankCpp(solp.zenr, vegp.x, si);
    tsdirstruct tspdir = twostreamdirCpp(vegp.pai, kk.kd, soilp.gref, tsvegp);
    albf albs = albedo(kk.kd, vegp.pai, tspdif, tspdir);
    double Rabs_sw = 0.0;
    double RabsG_sw = 0.0;
    double amx = soilp.gref; if (amx < vegp.lref) amx = vegp.lref;
    if (climdata.Rsw > 0.0) {
        double Rb0 = (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr);
        if (Rb0 > 900.0) Rb0 = 900.0;
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
        double Rdifdowng = Rdddg * climdata.Rdif + Rdbdg * Rb0;
        RabsG_sw = Rdifdowng + Rdirdowng;
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
    double ea = satvapCpp2(climdata.tref);
    // Initialize values
    double psih = 0.0;
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
    while (error > 1e-2 && itr < maxiter) {
        double zm = roughlengthCpp2(vegp.hgt, vegp.pai, d);
        uf = (ka * climdata.uref) / (std::log((zref - d) / zm) + psim);
        double zh = 0.2 * zm;
        double rHa = (std::log((zref - d) / zh) + psih) / (ka * uf);
        // Bulk surface stomatal resistance. psiw here needs to be
        // the root-zone mean water potential (rootfrac-weighted across all
        // soil layers), matching OneStepBelow's treatment -- not just the
        // top layer's value, which would make transpiration blind to
        // moisture held deeper in the profile than the very surface.
        double psiw = 0.0;
        for (size_t i = 0; i < soilwater.swo.rootfrac.size(); ++i) {
            psiw += soilwater.swo.rootfrac[i] * soilwater.swo.psiw[i];
        }
        double rS = bulkstomatalresist(climdata.tref, tcanopy, climdata.relhum, climdata.pk, psiw,
            Ca, climdata.precip, climdata.Rsw, climdata.Rdif, kk.k, LAIfrac, vegp, solp, C3);
        rS = std::min(rS, 1e10);
        double rV = rS + rHa;
        double rhz = 0.0; // resistance from top of canopy to zref
        if (zref > vegp.hgt) {
            // Heat-exchange resistance (paired with rHa above): uses the
            // heat diabatic correction dpsihCpp2, matching rHa and
            // rh_hzref's equivalent for the multilayer canopy model.
            double psihh = dpsihCpp2(zh / LL) - dpsihCpp2((vegp.hgt - d) / LL);
            double rHh = (std::log((vegp.hgt - d) / zh) + psihh) / (ka * uf);
            rhz = rHa - rHh;
        }
        // Resistance from ground to top of canopy. phih reflects stability
        // at canopy top (a2 describes within-canopy exchange), not at zref.
        double phih = dphihCpp2((vegp.hgt - d) / LL);
        double a2 = (phih * 0.41 * (1 - d / vegp.hgt)) / 1.5625;
        double rhg = rhcanopy(a2, uf, vegp.hgt, vegp.hgt);
        double rGz = rhg + rhz; // resistance from ground to zref
        double RabsG_lw = (tr * climdata.Rlw + (1.0 - tr) * sb * radem(tcanopy)) * soilp.groundem;
        double RabsG = RabsG_sw + RabsG_lw;
        double Tkc = tcanopy + 273.15;
        double es = satvapCpp2(tcanopy);
        Tr = Mw / (RgasC * Tkc * rV) * (es - ea) * 1000.0 * 3600.0 * vegp.pai * LAIfrac; // transpiration
        if (Tr < 0.0) Tr = 0.0;
        soilheat.wc = soilwater.swo.theta;
        soilheat = SoilHeatCpp(soilheat, soilp, RabsG, climdata.tref, climdata.relhum, climdata.pk, rGz, 3600, 0.5, maxiter);
        climforwaterstruct cfw = {};
        cfw.Rabs = RabsG; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rGz;
        cfw.precip = climdata.precip; cfw.Et = Tr;
        soilwater.swo.oldTc = soilheat.oldTe;
        soilwater.swo.Tc = soilheat.Te;
        soilwater = SoilWaterCpp(soilwater.swo, soilp, cfw, 3600, 0.5, maxiter, 1e-4);
        double G_raw = soilsurfaceEB(soilp, RabsG, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rGz, soilwater.swo.theta[0]); // surface energy balance residual
        G = aitken1d(G, G_raw, st);
        soilheat.Gflux = G;
        double Te = (tcanopy + climdata.tref) / 2.0;
        double Rer = 4.0 * em * sb * std::pow(Te + 273.15, 3.0);
        double la = (tcanopy >= 0.0) ? (45068.7 - 42.8428 * tcanopy) : (51078.69 - 4.338 * tcanopy - 0.06367 * tcanopy * tcanopy);
        double De = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5);
        double oldtcanopy = tcanopy;
        tcanopy = climdata.tref + ((Rabs - Rema - ((la * ph) / (climdata.pk * rV)) * Da - G) /
            (Rer + ph * (cp / rHa + ((la * De) / (climdata.pk * rV)))));
        if (tcanopy < tdew) tcanopy = tdew;
        H = (ph * cp / rHa) * (tcanopy - climdata.tref);
        LL = (ph * cp * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
        double Lsafe = clipMOlength(LL, zref, d, zm);
        if (H > 0) {
            if (LL < Lsafe) LL = Lsafe;
        }
        else {
            if (LL > Lsafe) LL = Lsafe;
        }
        psim = dpsimCpp2(zm / LL) - dpsimCpp2((zref - d) / LL);
        psih = dpsihCpp2(zh / LL) - dpsihCpp2((zref - d) / LL);
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
// One timestep of the single-layer bare-ground model (no canopy at all) --
// the big-leaf counterpart of OneStepBare, used for soil spin-up.
// Bare-ground big-leaf timestep used to spin up soil state when vegetation is absent.
// It retains the coupled surface energy balance, soil heat/water and atmospheric stability calculations without canopy physiology or radiation.
bigleafone solveonestepbare(const obsstruct& obsdata, const climstruct& climdata, const soilpstruct soilp,
    soilmod soilheat, soilwaterout soilwater, double zmr, double latr, double lonr, double zref = 2.0,
    int maxiter = 20)
{
    solmodel solp = solpositionCpp2(latr, lonr, obsdata.year, obsdata.month, obsdata.day, obsdata.hour);
    double sloper = soilp.slope * torad;
    double aspectr = soilp.aspect * torad;
    double si = solarindexCpp2(sloper, aspectr, solp.zenr, solp.azir);
    double Rabs_sw = 0.0;
    if (climdata.Rsw > 0.0) {
        double Rb0 = (climdata.Rsw - climdata.Rdif) / std::cos(solp.zenr);
        if (Rb0 > 900.0) Rb0 = 900.0;
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
    // G: surface energy balance residual, Aitken-damped across outer passes.
    Aitken1DState st_G;
    double G = soilheat.Gflux;
    while (error > 1e-2 && itr < maxiter) {
        double zm = zmr * std::exp(-psih);
        uf = (ka * climdata.uref) / (std::log(zref / zm) + psim);
        double zh = 0.2 * zm;
        double rHa = (std::log(zref / zh) + psih) / (ka * uf);
        soilheat.wc = soilwater.swo.theta;
        soilheat = SoilHeatCpp(soilheat, soilp, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600, 0.5, maxiter);
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = climdata.precip; cfw.Et = 0.0;
        soilwater.swo.oldTc = soilheat.oldTe;
        soilwater.swo.Tc = soilheat.Te;
        soilwater = SoilWaterCpp(soilwater.swo, soilp, cfw, 3600, 0.5, maxiter, 1e-4);
        double G_raw = soilsurfaceEB(soilp, Rabs, climdata.tref, soilheat.Te[0], climdata.pk,
            climdata.relhum, rHa, soilwater.swo.theta[0]);
        G = aitken1d(G, G_raw, st_G);
        soilheat.Gflux = G;
        newtsoil = aitken1d(newtsoil, soilheat.Te[0], st);
        soilheat.Te[0] = newtsoil;
        H = (ph * cp / rHa) * (newtsoil - climdata.tref);
        LL = (ph * cp * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
        double Lsafe = clipMOlength(LL, zref, 0.0, zm);
        if (H > 0) {
            if (LL < Lsafe) LL = Lsafe;
        }
        else {
            if (LL > Lsafe) LL = Lsafe;
        }
        psim = dpsimCpp2(zm / LL) - dpsimCpp2(zref / LL);
        psih = dpsihCpp2(zh / LL) - dpsihCpp2(zref / LL);
        error = std::abs(newtsoil - oldtsoil);
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
// The ectotherm module takes modelled microclimate as forcing and solves the animal energy balance.
// Body geometry/orientation alter both radiative interception and convective exchange, allowing behavioural microclimate use to be represented.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Calculates the projected animal silhouette exposed to direct radiation for the specified body geometry and orientation.
// This converts solar angle and posture into the area that intercepts direct shortwave.
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
    // Calculate surface area using Knud Thomsens formula for ellipsoid
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
// Compute characteristic dimension of animal
// Returns the characteristic body dimension normal to the wind for convective heat transfer.
// Body orientation therefore affects aerodynamic resistance as well as radiative exposure.
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
    // Compute characteristic 
    double d = sa.V / sa.silA;
    return d;
}
// Animal boundary-layer resistance to heat transfer -- same forced/free
// convection Nusselt-number correlation as leafrHa above, but over the
// animal's characteristic dimension (d) rather than a leaf's.
// Boundary-layer resistance for an animal body, combining forced and free convection at the body scale.
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
// Calculate water diffusivity
// Molecular diffusivity of water vapour in air at the specified film temperature and pressure.
double Dw_waterVapour(double Tf, double pk)
{
    const double D0 = 2.26e-5;   // m^2 s^-1 at 273.15 K and 101325 Pa
    const double T0 = 273.15;    // reference temperature (K)
    const double P0 = 101325.0;  // reference pressure (Pa)
    double Tfk = Tf + 273.15;
    double Pa = pk * 1000.0;
    return D0 * std::pow(Tfk / T0, 1.75) * (P0 / Pa);
}
// Animal body temperature from a Penman-Monteith-style steady-state energy
// balance: radiative and evaporative exchange with the air (fraction cd,
// via rHa/rc) plus, where the animal is in contact with a substrate
// (fraction confrac), radiative, conductive (keff) and evaporative
// exchange through that contact interface, all solved simultaneously
// against absorbed radiation (Rabs) and metabolic heat (M).
// Solves animal surface temperature and evaporative heat exchange from the whole-body energy balance.
// Radiation, convection, conduction/contact, metabolism and evaporative resistance are combined for the specified morphology and microclimate.
// [[Rcpp::export]]
double PenmanMonteith_animal(double Rabs, double Ta, double Ts, double Te, double Tf, double pk, double rh,
    double rHa, double height, double rc, double confrac, double M, double em = 0.97, double k = 0.5,
    double surfrh = 1.0)
{
    double Tbody = Ta;
    double cp = cpairCpp(Ta);
    double ph = phairCpp(Ta, pk);
    double dT = Ts - Ta; // contact surface minus air temperature
    double cd = 1.0 - confrac;
    double ea = satvapCpp2(Ta);
    double Da = ea * (1.0 - rh / 100.0); // vapour pressure deficit, air
    double Dac = surfrh * satvapCpp2(Tf) - ea; // vapour pressure deficit, contact interface
    double keff = k / (0.5 * height); // conductive coefficient through the contact interface
    double mr = cd * em * sb; // radiative exchange coefficient, air-facing
    double mrc = confrac * em * sb; // radiative exchange coefficient, contact-facing
    double mc = keff * confrac; // conductive exchange coefficient, contact-facing
    double mv = ph * cp * cd / rHa; // sensible heat exchange coefficient, air-facing
    double la = (Tbody >= 0.0) ? 45068.7 - 42.8428 * Tbody : 51078.69 - 4.338 * Tbody - 0.06367 * Tbody * Tbody;
    double rv = rHa + rc; // total vapour resistance (s/m)
    double ml = (ph * la * cd) / (rv * pk); // latent heat exchange coefficient, air-facing
    double Dw = Dw_waterVapour(Tf, pk);
    double rv_contact = ((0.5 * height * (Ts + 273.15) * RgasC) / (Dw * 1000.0)) + rc;
    double mlc = (la * confrac) / rv_contact; // latent heat exchange coefficient, contact-facing
    double DeV = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5); // slope of the saturation vapour pressure curve, air-facing
    double DeC = satvapCpp2(Tf + 0.5) - satvapCpp2(Tf - 0.5); // slope of the saturation vapour pressure curve, contact-facing
    double top = (Rabs + M - mr * radem(Ta) - mrc * radem(Ta) - mc * dT - ml * Da - mlc * Dac);
    double btm = 4.0 * mr * std::pow(Te + 273.15, 3.0) + 4.0 * mrc * std::pow(Tf + 273.15, 3.0) +
        mc + mv + ml * DeV + surfrh * mlc * DeC;
    Tbody = Ta + top / btm;
    return (Tbody);
}
// Scalar version of aitkin_weightdif (Aitken-style adaptive damping),
// parameterised by backweight (= 1-omega) instead of omega directly.
// Adaptive scalar under-relaxation for the ectotherm body-temperature iteration.
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
// Animal metabolic rate (W): allometric scaling with body mass (m^b),
// scaled for body temperature via a Q10 temperature coefficient.
// Returns whole-animal metabolic heat production from body size and a Q10 temperature response.
double metabolic_rate(double volume, double rho, double Q10, double a0, double b, double Tref, double Tbody)
{
    double m = volume * rho;
    double M = a0 * std::pow(m, b) * std::pow(Q10, (Tbody - Tref) / 10.0);
    return M;
}
// Position details:
// one of fixed (constant adir and atilt assumed), 
// max (animal seeks to maximize radiation absorption)
// min (animal seeks to minimize radiation absorption)
// randomdir (animal maintains tilt, but adir is random)
// random (animal has random atilt and adir)
// Animal body temperature time series: at each timestep, absorbed
// shortwave (via the animal's solar-facing silhouette area) and longwave
// radiation feed PenmanMonteith_animal's energy balance, iterated to
// convergence (the animal's own temperature affects its emitted longwave
// and boundary-layer exchange, both inputs to the same balance).
// Runs the ectotherm energy-balance model through a time series for one supplied microclimate.
// At each timestep, body orientation controls radiation and convection, while metabolism and evaporative/conductive exchange are iterated to obtain body temperature.
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
                rc, confrac, M, em, k, surfrh[i]);
            aitkin_weightdif_scalar(Told, Tbody[i], st);
            double dTold = dT;
            dT = Tbody[i] - Ta[i];
            tdif = std::abs(dTold - dT);
            ++nrIterations;
        }
    }
    return Tbody;
}
// Ectotherm's own energy balance, over a time x vertical-layer grid of
// microclimate inputs instead of a single time series (e.g. body
// temperature at every canopy height simultaneously).
// Runs the ectotherm model across multiple candidate microclimates and selects the thermally relevant result at each timestep.
// This is the multi-location counterpart of Ectotherm for behavioural use of heterogeneous microclimate.
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
                    rc, confrac, M, em, k, 1.0);
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
// R wrappers translate flexible R lists/data frames into typed C++ state, run the sequential
// timestep model, and package either requested-height outputs or full diagnostic profiles back to R.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Plant area index above each canopy node (cumulative sum from the top down).
// Converts layer PAI increments to cumulative PAI above each layer, the vertical coordinate required by canopy radiative transfer.
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
// Creates c++ veg struct
// Translates the R vegetation parameter list into the typed C++ structure used throughout the model.
// Derived/cache fields and layer vectors are initialised here so downstream physics can work without repeated R lookups.
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
    out.rootskew = vegp["rootskew"]; // Skew towards top of root layer (unitless)
    int n = static_cast<int>(paii.size());
    out.pai = 0.0; // Total one-sided plant area index (m^2 / m^2)
    for (int i = 0; i < n; ++i) out.pai += paii[i];
    out.paii = paii; // one-sided plant area index at each canopy node
    out.pia = reverseCumsum(paii); // one-sided plant area index above each canopy node
    out.Lfrac = Lfrac; // Fraction of paii that is living plant material
    return out;
}
// Creates c++ soilc struct
// Translates the R soil parameter list into the typed hydraulic/thermal parameter structure used by the soil solvers.
static soilpstruct tosoilpstruct(List soilp) {
    soilpstruct out;
    // Scalars (force extraction to double)
    out.slope = Rcpp::as<double>(soilp["slope"]) * torad; // radians
    out.aspect = Rcpp::as<double>(soilp["aspect"]) * torad; // radians
    out.gref = Rcpp::as<double>(soilp["gref"]);
    out.grefPAR = Rcpp::as<double>(soilp["grefPAR"]);
    out.groundem = Rcpp::as<double>(soilp["groundem"]);
    out.nLayers = Rcpp::as<int>(soilp["nLayers"]);
    // Vectors (explicit conversion to std::vector<double>)
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
    // Ensure psi_e is negative (air entry potential)
    for (double& v : out.psie) v = -std::abs(v);
    // Calculate minimum psie
    out.psi_min = std::vector<double>(out.psie.size());
    for (size_t i = 0; i < out.psie.size(); ++i) {
        out.psi_min[i] = waterPotential(out, out.thetaR[i], i);
    }
    return out;
}
// Creates input file for soil heat model
// Builds the initial dynamic soil-heat state from supplied temperature and water-content profiles.
static soilmod toSoilheatmod(List soilp, std::vector<double> Te, std::vector<double> wc)
{
   
    int n = Rcpp::as<int>(soilp["nLayers"]);
    double totalDepth = Rcpp::as<double>(soilp["totalDepth"]);
    std::vector<double> z = geometricCpp(n, totalDepth);
    std::vector<double> dz(n + 1);
    std::vector<double> zCenter(n + 1);
    for (int i = 0; i <= n; ++i) {
        dz[i] = z[i + 1] - z[i];
        zCenter[i] = z[i] + dz[i] * 0.5;
    }
    soilmod state;
    state.n = n;
    state.z = z; // node depths
    state.dz = dz; // layer thickness
    state.zCenter = zCenter; // centre between nodes
    state.vol = dz; // volume of layer
    state.wc = wc; // volumetric water fraction of layer
    state.Te = Te; // temperature of layer
    state.oldTe = Te; // old temperature of layer
    state.Gflux = 0.0;
    state.iters = 0;
    return state;
}
// Creates input file for soil water model
// Builds the initial dynamic soil-water state, including matric potential, vapour and the root distribution derived from the soil/vegetation setup.
static soilwatermod toSoilwatermod(soilpstruct soilpc, soilmod state, double rootskew)
{
    soilwatermod out;
    out.n = state.n; // number of layers
    out.z = state.z; // node depths
    out.dz = state.dz; // layer thickness
    out.zCenter = state.zCenter; // centre between nodes
    out.Tc = state.Te; // temperature of soil profile
    out.oldTc = state.oldTe; // old temperature of soil profile
    out.theta = state.wc; // volumetric water content of soil profile
    // Initialize
    std::vector<double> vol(out.n + 1);
    std::vector<double> psiw(out.n + 1);
    std::vector<double> k(out.n + 1);
    std::vector<double> vapor(out.n + 1);
    // Node-centred control volume: half the distance to each neighbour.
    // The surface node (i=0, true physical boundary, no i-1 neighbour)
    // gets a one-sided half-cell instead -- leaving it at zero starves
    // the surface Newton row of capacitance, causing non-convergence as
    // k[0]->0 while drying.
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
    double totalDepth = 0.0; // deepest node depth, used to scale root_distribute below
    for (int i = 0; i <= out.n; ++i) {
        if (out.z[i] > totalDepth) totalDepth = out.z[i];
    }
    out.rootfrac = root_distribute(out.dz, totalDepth, rootskew);
    return out;
}
// Initial onestep state: every canopy layer set to the same starting
// temperature/humidity/wind (temp0/rh0/Rlw0), refined once the model
// actually starts iterating.
// Assembles the initial multilayer timestep state from spun-up soil conditions and reference atmospheric values.
// This provides the starting leaf/air/soil fields from which OneStepBelow begins its coupled iteration.
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
    onestepin.Hz = v0;
    onestepin.Lz = v0;
    onestepin.gs = v0;
    onestepin.precipground = 0.0;
    onestepin.H = 0.0;
    onestepin.L = 0.0;
    onestepin.Et = 0.0;
    onestepin.Ev = 0.0;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    onestepin.psim = 0.0;
    onestepin.psih = 0.0;
    onestepin.phih = 1.0;
    onestepin.LL = 999.99;
    onestepin.iters = 0;
    onestepin.witers = 0;
    onestepin.error = 0.0;
    return onestepin;
}
// Flattens a onestep struct into the R-facing named list returned by the
// profile/timeseries R functions.
// Packages the internal multilayer timestep state into an R list for diagnostic/profile output.
static Rcpp::List OneStepCpptoList(onestep onestepin, std::vector<double> z)
{
    Rcpp::List out;
    // Soil heat variables
    out["nb"] = Rcpp::wrap(onestepin.soilheatvars.n);
    out["zb"] = Rcpp::wrap(onestepin.soilheatvars.z);
    out["dz"] = Rcpp::wrap(onestepin.soilheatvars.dz);
    out["zCenter"] = Rcpp::wrap(onestepin.soilheatvars.zCenter);
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
    out["uz"] = Rcpp::wrap(onestepin.uz);
    out["rLB"] = Rcpp::wrap(onestepin.rLB);
    out["swaterdepth"] = Rcpp::wrap(onestepin.swaterdepth);
    out["Hz"] = Rcpp::wrap(onestepin.Hz);
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
    out["LL"] = onestepin.LL;
    out["iters"] = onestepin.iters;
    out["witers"] = Rcpp::wrap(onestepin.witers);
    out["error"] = onestepin.error;
    // Above ground z
    out["z"] = Rcpp::wrap(z);
    return out;
}
// Bare-ground vertical profile at a single requested hour: first spins the
// soil state forward hour-by-hour up to hourtoplot using a coarse
// two-point height grid (z2), then re-runs the final hour at the full
// resolution of the requested z to return that hour's profile.
// R-facing diagnostic runner for a bare-ground profile at a selected hour.
// It advances the model to that timestep and returns the detailed internal state rather than only standard time-series outputs.
// [[Rcpp::export]]
List profilebareR(size_t hourtoplot, DataFrame obstime, DataFrame climdata, List soilc, std::vector<double> z,
    double zref, double lat, double lon, std::vector<double> SoilTempIni,
    std::vector<double> SoilThetaIni, double zm = 0.004, int maxNrIterations = 100, double tolerance = 1e-3)
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
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 10.0);
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
    out["zCenter"] = onestepin.soilheatvars.zCenter;
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
// Wrapper function for returning profile
// R-facing diagnostic runner for the vegetated multilayer model at a selected hour.
// It is intended for inspecting the full vertical state and flux structure of one timestep.
// [[Rcpp::export]]
List profileR(size_t hourtoplot, DataFrame obstime, DataFrame climdata, List soilc, List vegp,
    std::vector<double> paii20, std::vector<double> paii, std::vector<double> Lfrac20, std::vector<double> Lfrac, 
    double zref, double Ca, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double tolerance = 1e-3, double a0 = 0.25, double a1 = 1.25, bool C3 = true)
{
    // ** ------------------ Run model up to hour with 20 layers ------------------------- ** //
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
    std::vector<double> wc = windprofileCpp(vegpc);
    size_t na = paii20.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
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
    // ** --------------------------------- Run full profile for hour  --------------------------------- ** //
    vegpc = tovegpstruct(vegp, paii, Lfrac);
    wgts = lwradweights(vegpc.paii);
    wc = windprofileCpp(vegpc);
    size_t hr = hourtoplot;
    obsstruct obsdata;
    climstruct climin;
    obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
    climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
    climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
    climin.winddir = wdir[hr]; climin.precip = precip[hr];
    na = paii.size();
    std::vector<double> z2(na);
    for (size_t i = 0; i < na; ++i) z2[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[hr], relhum[hr], Rlw[hr], na);
    onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z2, tspveg, tspvegPAR, tspdif,
        tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
    return OneStepCpptoList(onestepin, z2);
}
// Bare-ground time series at a single requested height (reqhgt): runs
// OneStepBare hour-by-hour across the whole input time series, returning
// one row per hour.
// Main R-facing time-series driver for the bare-ground microclimate model.
// It advances soil and surface state sequentially through the forcing data and extracts conditions at the requested height for each timestep.
// [[Rcpp::export]]
DataFrame RunBareR(double reqhgt, DataFrame obstime, DataFrame climdata, List soilc, double zref,
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    double zm = 0.004, int maxNrIterations = 100, double  tolerance = 1e-3)
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
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 10.0);
    onestepin.H = 0.0; onestepin.L = 0.0; onestepin.psim = 0.0; onestepin.psih = 0.0; 
    onestepin.Rswup = 0.0; onestepin.Rb0 = 0.0; onestepin.Rlwup = 0.0;
    onestepin.LL = 999.99; onestepin.iters = 1; onestepin.error = 1e99;
    onestepin.theta = onestepin.soilwatervars.theta[0];
    // identifiy which value to select if below ground
    size_t sel = 0;
    if (reqhgt < 0.0) {
        std::vector<double> zz = onestepin.soilheatvars.z;
        double mn = std::abs(reqhgt + zz[0]);
        for (size_t i = 1, n = zz.size(); i < n; ++i) {
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
    std::vector<int> iters(n);
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
        // Error variables
        Named("iters") = Rcpp::wrap(iters),
        Named("error") = Rcpp::wrap(error)
    );
}
// Vegetated time series at a single requested height (reqhgt): runs
// OneStepBelow hour-by-hour across the whole input time series, returning
// one row per hour -- the multilayer-canopy counterpart of RunBareR.
// Main R-facing time-series driver for the multilayer vegetated microclimate model.
// It prepares static canopy/radiation geometry once, advances the fully coupled OneStepBelow state hour by hour, and returns requested-height weather plus soil/canopy diagnostics.
// [[Rcpp::export]]
Rcpp::List RunModelR(double reqhgt, Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List soilc,
    Rcpp::List vegp, std::vector<double> paii, std::vector<double> Lfrac, double zref, double Ca,
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni,
    int maxNrIterations = 100, double tolerance = 1e-3, double a0 = 0.25,
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
    std::vector<double> wc = windprofileCpp(vegpc);
    size_t n = temp.size();
    size_t na = paii.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i)
        z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    double latr = lat * torad;
    double lonr = lon * torad;
    std::vector<double> zz = z;
    if (reqhgt < 0.0) {
        zz = soilheatvars.z;
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
    std::vector<double> tair(n), tground(n), tleaf(n), rz(n), uz(n);
    std::vector<double> H(n), L(n), G(n), Ev(n), Et(n), error(n);
    std::vector<int> iters(n);
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
        Et[hr] = onestepin.Et;
        iters[hr] = onestepin.iters;
        error[hr] = onestepin.error;
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
            if (reqhgt > vegpc.hgt) {
                Rdifdown[hr] = Rdif[hr];
                Rlwdown[hr] = Rlw[hr];
                double Th = onestepin.tair[na - 1];
                double rh = onestepin.rh[na - 1];
                double uh = onestepin.uz[na - 1];
                tair[hr] = Tabove(reqhgt, zref, Th, temp[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
                tleaf[hr] = -999.99;
                rz[hr] = RHabove(reqhgt, zref, rh, Th, temp[hr], tair[hr], relhum[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
                uz[hr] = Uabove(reqhgt, zref, uh, wspeed[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
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
                    rz[hr] = onestepin.rh[sel];
                    uz[hr] = onestepin.uz[sel];
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
    out["relhum"] = rz;
    out["windspeed"] = uz;
    out["H"] = H;
    out["L"] = L;
    out["G"] = G;
    out["Evaporation"] = Ev;
    out["Transpiration"] = Et;
    out["iters"] = iters;
    out["error"] = error;
    out.attr("class") = "data.frame";
    out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(n));
    return out;
}
// Re-expresses a weather time series measured at height zin as the
// equivalent series at height zout, above canopy: runs OneStepBelow at
// zin each hour, then extrapolates temperature/humidity/wind/pressure to
// zout via the above-canopy profile functions (Tabove/RHabove, a
// hypsometric correction for pressure).
// Runs the multilayer model and converts weather from one height to another, using the modelled canopy/above-canopy profile rather than a fixed lapse or wind adjustment.
// It is the height-adjustment interface built on the same timestep physics as RunModelR.
// [[Rcpp::export]]
List WeatherhgtCpp2(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zin, double zout, double lat, double lon, 
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, double CO2ppm = 430.0)
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
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc);
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
            tspdifPAR, wgts, wc, CO2ppm, latr, lonr, zin, 100, 1e-2);
        double th = onestepin.tair[na - 1];
        double rh = onestepin.rh[na - 1];
        temp_new[hr] = Tabove(zout, zin, th, temp[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
        relhum_new[hr] = RHabove(zout, zin, rh, th, temp[hr], temp_new[hr], relhum[hr], vegpc.hgt, vegpc.pai, onestepin.LL);
        // One-shot post-convergence evaluation of onestepin.LL (already
        // finished above, not part of OneStepBelow's own iteration).
        double zm = roughlengthCpp2(vegpc.hgt, vegpc.pai, d);
        double uf = (ka * wspeed[hr]) / (std::log((zin - d) / zm) + onestepin.psim);
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

// Writes vector v into row hr of mat (one hour's profile into a
// time x height output matrix).
// Small output helper that copies a std::vector profile into one row of an R numeric matrix.
void FillMyMatrix(NumericMatrix& mat, std::vector<double>& v, size_t hr)
{
    for (size_t i = 0; i < v.size(); ++i) {
        mat(hr, i) = v[i];
    }
}
// Bare-ground time x height matrices for the whole input time series and
// requested z: runs OneStepBare hour-by-hour, writing each hour's profile
// into a row via FillMyMatrix.
// Runs the bare-ground model and retains the complete vertical/soil state at every timestep.
// This is the full-field diagnostic counterpart of the reduced requested-height output in RunBareR.
// [[Rcpp::export]]
List RunBelowFullBare(DataFrame obstime, DataFrame climdata, List soilc, std::vector<double> z,
    double zref, double zm, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double  tolerance = 1e-2)
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
    onestepin.soilwatervars = toSoilwatermod(soilpc, onestepin.soilheatvars, 10.0);
    onestepin.H = 0.0; onestepin.L = 0.0; onestepin.psim = 0.0; onestepin.psih = 0.0;
    onestepin.Rswup = 0.0; onestepin.Rb0 = 0.0; onestepin.Rlwup = 0.0;
    onestepin.LL = 999.99; onestepin.iters = 1; onestepin.error = 1e99;
    onestepin.theta = onestepin.soilwatervars.theta[0];
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
        Named("theta") = Rcpp::wrap(theta)
    );
}
// Vegetated time x height matrices for the whole input time series and
// requested z: the multilayer-canopy counterpart of RunBelowFullBare.
// Runs the vegetated multilayer model and retains full vertical fields through time.
// It exposes the internal canopy-air state for diagnostics or downstream calculations instead of collapsing to a single requested height.
// [[Rcpp::export]]
List RunBelowFull(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zref, double Ca, double lat, double lon,
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, int maxNrIterations = 100,
    double  tolerance = 1e-2, double a0 = 0.25, double a1 = 1.25, bool C3 = true)
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
    std::vector<double> wc = windprofileCpp(vegpc);
    size_t n = temp.size();
    size_t na = paii.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    // Create Numeric Matrices for above ground
    NumericMatrix Rdirdown(n, na);
    NumericMatrix Rdifdown(n, na);
    NumericMatrix Rswup(n, na);
    NumericMatrix Rlwdown(n, na);
    NumericMatrix Rlwup(n, na);
    NumericMatrix tair(n, na);
    NumericMatrix tleaf(n, na);
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
        onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR, tspdif,
            tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
        // Extract values
        FillMyMatrix(Rdirdown, onestepin.Rdirdown, hr);
        FillMyMatrix(Rdifdown, onestepin.Rdifdown, hr);
        FillMyMatrix(Rswup, onestepin.Rswup, hr);
        FillMyMatrix(Rlwdown, onestepin.Rlwdown, hr);
        FillMyMatrix(Rlwup, onestepin.Rlwup, hr);
        FillMyMatrix(tair, onestepin.tair, hr);
        FillMyMatrix(tleaf, onestepin.tleaf, hr);
        FillMyMatrix(rh, onestepin.rh, hr);
        FillMyMatrix(uz, onestepin.uz, hr);
        FillMyMatrix(tsoil, onestepin.soilheatvars.Te, hr);
        FillMyMatrix(theta, onestepin.soilwatervars.theta, hr);
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    return List::create(
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
        // Below ground temp and theta
        Named("tsoil") = Rcpp::wrap(tsoil),
        Named("theta") = Rcpp::wrap(theta)
    );
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Other useful functions *********************************************** //
// Convenience utilities derived from the same solar/radiation machinery, plus the public big-leaf
// spin-up entry points used to prepare soil initial conditions for the main model.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Clear-sky direct-beam shortwave irradiance, from optical depth terms
// for Rayleigh scattering + permanent gases, water vapour absorption
// (via dewpoint) and aerosols.
// Calculates clear-sky shortwave radiation through the supplied date/time series.
// The result provides a solar-geometry benchmark independent of the observed/modelled cloud attenuation.
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
// Diffuse fraction of shortwave radiation: a clearness-index (k =
// swrad/top-of-atmosphere beam) separation model, with an additional
// correction (sigma3/delta terms) for sky variability.
// Estimates the diffuse fraction of incoming shortwave from total radiation and solar geometry.
// This supplies the direct/diffuse partition required by the two-stream canopy radiation model.
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
// Solar altitude (degrees above horizon).
// Returns solar altitude through a date/time series for the supplied location.
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
// Linearly upsamples a matrix from nin to nout rows per column (e.g. for
// smoother plotting of hourly output), interpolating column-wise.
// Linearly resamples time-by-variable output to a requested number of rows.
// This is a presentation/interpolation utility and does not alter the physical model state.
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
// Spins up soil heat/water state ahead of a full multilayer run: runs the
// single-layer big-leaf model (solveonestep) hour-by-hour across the
// whole input time series from a linearly-interpolated initial soil
// profile, and returns the converged final soil state as a warm start.
// Spins up soil temperature and moisture beneath vegetation before the full multilayer run.
// The cheaper big-leaf model is marched through the forcing series so the returned soil profiles are dynamically consistent rather than arbitrary initial conditions.
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
    // Soil profile initial-condition seed: temperature interpolates from
    // temp[0] to boundaryT; water content interpolates from 60% to 100%
    // of thetaS[n-1] (matching weatherhgt_adjust()'s R-level initial guess).
    double end1 = boundaryT;
    double end2 = soilpc.thetaS[n - 1];
    double start1 = temp[0];
    double start2 = 0.6 * end2;
    for (size_t i = 0; i < n; ++i) {
        SoilTempIni[i] = start1 + (end1 - start1) * static_cast<double>(i) / static_cast<double>(n - 1);
        SoilThetaIni[i] = start2 + (end2 - start2) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    soilmod soilheat = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwater = toSoilwatermod(soilpc, soilheat, vegpc.rootskew);
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
// Bare-ground counterpart of BigLeafCpp2: spins up soil heat/water state
// via solveonestepbare instead of solveonestep, for runs with no canopy.
// Bare-ground soil spin-up counterpart of BigLeafCpp2.
// It equilibrates soil heat and water with the forcing series without paying the cost of the full multilayer canopy calculation.
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
    // Soil profile initial-condition seed: temperature interpolates from
    // temp[0] to boundaryT; water content interpolates from 60% to 100%
    // of thetaS[n-1] (matching weatherhgt_adjust()'s R-level initial guess).
    double end1 = boundaryT;
    double end2 = soilpc.thetaS[n - 1];
    double start1 = temp[0];
    double start2 = 0.6 * end2;
    for (size_t i = 0; i < n; ++i) {
        SoilTempIni[i] = start1 + (end1 - start1) * static_cast<double>(i) / static_cast<double>(n - 1);
        SoilThetaIni[i] = start2 + (end2 - start2) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    soilmod soilheat = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwater = toSoilwatermod(soilpc, soilheat, 10.0);
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
