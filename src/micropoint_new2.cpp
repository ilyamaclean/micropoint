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
constexpr double vp = 0.017;
constexpr double torad = 3.14159265358979323846 / 180.0;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************** Solar model ***************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculates Astronomical Julian day ** //
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
static double soltimeCpp(int jd, double lt, double lond)
{
    double m = 6.24004077 + 0.01720197 * (static_cast<double>(jd) - 2451545.0);
    double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2.0 * m + 3.5932);
    double st = lt + (4.0 * lond + eot) / 60.0;
    return st;
}
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculate canopy extinction coefficient for sloped ground surfaces** //
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
        k = std::tan(zenr);
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
static tsvegstruct twostreamvegCpp(double pai, double x, double lref, double ltra, double gref)
{
    // Calculate base parameters
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
    // Calculate extended parameters
    params.S1 = std::exp(-params.h * pai);
    params.u1 = params.a + params.gma * (1.0 - 1.0 / gref);
    params.u2 = params.a + params.gma * (1.0 - gref);
    params.D1 = (params.a + params.gma + params.h) * (params.u1 - params.h) *
        1.0 / params.S1 - (params.a + params.gma - params.h) * (params.u1 + params.h) * params.S1;
    params.D2 = (params.u2 + params.h) * 1.0 / params.S1 - (params.u2 - params.h) * params.S1;
    return params;
}
// ** Calculates parameters for diffuse radiation using two-stream model ** //
static tsdifstruct twostreamdifCpp(tsvegstruct& tsvegp)
{
    tsdifstruct params{};
    // Calculate parameters
    params.p1 = (tsvegp.gma / (tsvegp.D1 * tsvegp.S1)) * (tsvegp.u1 - tsvegp.h);
    params.p2 = (-tsvegp.gma * tsvegp.S1 / tsvegp.D1) * (tsvegp.u1 + tsvegp.h);
    params.p3 = (1.0 / (tsvegp.D2 * tsvegp.S1)) * (tsvegp.u2 + tsvegp.h);
    params.p4 = (-tsvegp.S1 / tsvegp.D2) * (tsvegp.u2 - tsvegp.h);
    return params;
}
// ** Calculates parameters for direct radiation using two - stream model * * //
static tsdirstruct twostreamdirCpp(double pai, double kd, double gref, tsvegstruct tsvegp)
{
    tsdirstruct params{};
    // Calculate base parameters
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
static radmodel shortwavemodelCpp(const std::vector<double>& pia, double pai, double gref, double grefPAR, double lref, double lrefp, double Rswdown, double Rdif, double si, const solmodel& solp, const kstruct& kp, const tsvegstruct& tspveg, const tsvegstruct& tspvegPAR, const tsdifstruct& tspdif, const tsdifstruct& tspdifPAR, const tsdirstruct& tspdir, const tsdirstruct& tspdirPAR) {
    const int n = static_cast<int>(pia.size());
    radmodel out;
    out.RswLsun.assign(n, 0.0); out.RswLshade.assign(n, 0.0); out.RswLav.assign(n, 0.0);
    out.RPARsun.assign(n, 0.0); out.RPARshade.assign(n, 0.0); out.sunfrac.assign(n, 0.0);
    out.Rdirdown.assign(n, 0.0); out.Rdifdown.assign(n, 0.0); out.Rswup.assign(n, 0.0);
    out.RswGabs = 0.0; out.RswCabs = 0.0;
    if (Rswdown <= 0.0) return out;
    const double cosz = std::cos(solp.zenr);
    double Rbeam = (Rswdown - Rdif) / cosz;
    if (Rbeam > 1352.0) Rbeam = 1352.0;
    const double Rb = Rbeam * cosz;
    double amx = gref; if (amx < lref)amx = lref;
    double amxp = grefPAR; if (amxp < lrefp)amxp = lrefp;
    const double kd = kp.kd, h = tspveg.h, hp = tspvegPAR.h;
    const double sil = kp.k * cosz;
    const double exp_kd_pai = std::exp(-kd * pai);
    const double exp_mh_pai = std::exp(-h * pai);
    const double exp_ph_pai = std::exp(h * pai);
    double Rdddg = tspdif.p3 * exp_mh_pai + tspdif.p4 * exp_ph_pai;
    if (Rdddg > 1.0) Rdddg = 1.0;
    if (Rdddg < 0.0) Rdddg = 1.0;
    double Rdbdg = 0.0;
    if (Rb > 0.0) {
        Rdbdg = (tspdir.p8 / tspdir.sig) * exp_kd_pai + tspdir.p9 * exp_mh_pai + tspdir.p10 * exp_ph_pai;
        if (Rdbdg > amx) Rdbdg = amx;
        if (Rdbdg < 0.0) Rdbdg = 0.0;
    }
    const double Rdirdowng = Rbeam * exp_kd_pai;
    const double Rdifdowng = Rdddg * Rdif + Rdbdg * Rbeam;
    out.RswGabs = (1.0 - gref) * (Rdifdowng + si * Rdirdowng);
    const double albd = tspdif.p1 + tspdif.p2;
    const double albb = tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7;
    const double trg = exp_kd_pai;
    if (Rb > 0.0) {
        const double Rbc = (trg * si + (1.0 - trg) * cosz) * Rbeam;
        out.RswCabs = (1.0 - albd) * Rdif + (1.0 - albb) * Rbc;
    }
    else out.RswCabs = (1.0 - albd) * Rdif;
    for (int i = 0; i < n; ++i) {
        const double p = pia[i];
        const double exp_kd = std::exp(-kd * p);
        const double exp_mh = std::exp(-h * p);
        const double exp_ph = std::exp(h * p);
        double Rddu = tspdif.p1 * exp_mh + tspdif.p2 * exp_ph;
        if (Rddu > 1.0) Rddu = 1.0;
        if (Rddu < 0.0) Rddu = 0.0;
        double Rddd = tspdif.p3 * exp_mh + tspdif.p4 * exp_ph;
        if (Rddd > 1.0) Rddd = 1.0;
        if (Rddd < 0.0) Rddd = 1.0;
        double Rdbu = 0.0, Rdbd = 0.0;
        if (Rb > 0.0) {
            Rdbu = (tspdir.p5 / -tspdir.sig) * exp_kd + tspdir.p6 * exp_mh + tspdir.p7 * exp_ph;
            if (Rdbu > amx) Rdbu = amx;
            if (Rdbu < 0.0) Rdbu = 0.0;
            Rdbd = (tspdir.p8 / tspdir.sig) * exp_kd + tspdir.p9 * exp_mh + tspdir.p10 * exp_ph;
            if (Rdbd > amx) Rdbd = amx;
            if (Rdbd < 0.0) Rdbd = 0.0;
        }
        out.Rdirdown[i] = Rbeam * exp_kd;
        out.Rdifdown[i] = Rddd * Rdif + Rdbd * Rb;
        out.Rswup[i] = Rddu * Rdif + Rdbu * Rb;
        out.RswLsun[i] = tspveg.a * 0.5 * (sil * Rbeam + out.Rdifdown[i] + out.Rswup[i]);
        out.RswLshade[i] = tspveg.a * 0.5 * (out.Rdifdown[i] + out.Rswup[i]);
        out.RswLav[i] = tspveg.a * 0.5 * (sil * out.Rdirdown[i] + out.Rdifdown[i] + out.Rswup[i]);
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
        out.RPARshade[i] = tspvegPAR.a * (out.Rdifdown[i] + out.Rswup[i]);
        out.RPARsun[i] = tspvegPAR.a * (sil * Rbeam + out.Rdifdown[i] + out.Rswup[i]);
        out.sunfrac[i] = exp_kd;
    }
    return out;
}
// ** Longwave radiation model
// ** Calculate longwave radiation weights ** //
static LWweights lwradweights(const std::vector<double>& paii) {
    // Calculate total pai
    int n = static_cast<int>(paii.size());
    double pait = 0.0;
    for (int i = 0; i < n; ++i) pait += paii[i];
    // Calculate pai above and below
    std::vector<double> paia(n);
    std::vector<double> paib(n);
    paia[0] = pait;
    paib[0] = 0.0;
    for (int i = 1; i < n; ++i) {
        paib[i] = paib[i - 1] + paii[i];
        paia[i] = pait - paib[i];
    }
    // For each element i, calculate weight to every other element
    // weight is paii of j * transmission between i & j
    NumericMatrix wgts(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Calculate leaf area between i and j (used for scaling canopy gaps)
                double pij = std::abs(paia[i] - paia[j]); // pai between i and j
                double tr = std::exp(-pij); // transmission 
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
    // ** calculate total canopy weight
    std::vector<double> wgt(n);
    for (int i = 0; i < n; ++i) {
        wgt[i] = 0.0;
        for (int j = 0; j < n; ++j) wgt[i] += wgts(i, j);
    }
    // ** Check how close result is to two
    std::vector<double> wsum(n);
    for (int i = 0; i < n; ++i) wsum[i] = wgt[i] + trh[i] + trg[i];
    // ** calculate adjustment factor
    std::vector<double> mu(n);
    for (int i = 0; i < n; ++i) {
        double mu1 = wsum[i] / 2.0;
        mu[i] = 1.0 + (2.0 - 2.0 * mu1) / wgt[i];
    }
    // ** apply adjustment factor
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) wgts(i, j) = wgts(i, j) * mu[i];
    }
    // Calculate weight of canopy elements for ground absorbed
    std::vector<double> wgtg(n);
    for (int i = 0; i < n; ++i) {
        wgtg[i] = trg[i] * paii[i];
    }
    // Check how close summed weight plus transmission to sky is to 1
    // ** summed weight
    double gwsum = 0.0;
    for (int i = 0; i < n; ++i) gwsum += wgtg[i];
    // ** transmission to sky
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
static double radem(double tc)
{
    return std::pow(tc + 273.15, 4.0);
}
// Calculate longwave radiation below canopy
static radmodel2 longwavemodelCpp(const LWweights& wgts, double lwdown, double tground, double groundem, double vegem,
    const std::vector<double>& tleaf)
{
    // Initialize output variables
    int n = static_cast<int>(wgts.trh.size());
    std::vector<double> RlwLabs(n);
    std::vector<double> Rlwdown(n);
    std::vector<double> Rlwup(n);
    // Compute leaf and ground longwave (speeds up function)
    std::vector<double> leafLW(n);
    for (int i = 0; i < n; ++i) leafLW[i] = vegem * sb * radem(tleaf[i]);
    double groundLW = groundem * sb * radem(tground);
    for (int i = 0; i < n; ++i) {
        // longwave radiation downward from sky
        double lwsky = lwdown * wgts.trh[i];
        // longwave radiation upward from ground
        double lwgro = groundLW * wgts.trg[i];
        // longwave radiation from canopy elements
        double lwcand = 0.0; // down from canopy
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
    // Ground longwave radiation
    double lwgrfromcan = 0.0;
    for (int i = 0; i < n; ++i) {
        lwgrfromcan += wgts.wgtg[i] * leafLW[i];
    }
    // longwave radiation downward from sky
    double lwgrfromsky = wgts.trsky * lwdown;
    // Ground absorbed
    double RlwGabs = groundem * (lwgrfromsky + lwgrfromcan);
    // Canopy absorbed
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Worker functions **
static double zeroplanedisCpp2(double h, double pai)
{
    if (pai < 0.001) pai = 0.001;
    double d = (1.0 - (1.0 - exp(-std::sqrt(7.5 * pai))) / std::sqrt(7.5 * pai)) * h;
    return d;
}
// ** Calculate roughness length ** //
static double roughlengthCpp2(double h, double pai, double d, double psi_h)
{
    double Be = std::sqrt(0.003 + (0.2 * pai) / 2.0);
    double zm = (h - d) * std::exp(-ka / Be) * std::exp(ka * psi_h);
    if (zm < 0.0005) zm = 0.0005;
    // safety check to stop diabatic coefficient reversing profile
    if (zm > (0.9 * (h - d))) zm = 0.9 * (h - d);
    return zm;
}
// ** Calculate molar density of air ** //
static double phairCpp(double tc, double pk)
{
    double tk = tc + 273.15;
    double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
    return ph;
}
// ** Calculate specific heat of air at constant pressure ** //
static double cpairCpp(double tc)
{
    double cp = 2e-05 * tc * tc + 0.0002 * tc + 29.119;
    return cp;
}
// **  Calculate integrated diabatic correction coefficient for momentum ** //
static double dpsimCpp2(double ze)
{
    double psim;
    // unstable
    if (ze < 0.0) {
        double x = std::pow((1.0 - 15.0 * ze), 0.25);
        psim = std::log(std::pow((1.0 + x) / 2.0, 2.0) * (1.0 + x * x) / 2.0) - 2.0 * std::atan(x) + pi / 2.0;
    }
    // stable
    else {
        psim = -4.7 * ze;
    }
    if (psim < -4.0) psim = -4.0;
    if (psim > 3.0) psim = 3.0;
    return psim;
}
// **  Calculate integrated diabatic correction coefficient for heat ** //
static double dpsihCpp2(double ze)
{
    double psih;
    // unstable
    if (ze < 0.0) {
        double y = std::sqrt(1.0 - 9.0 * ze);
        psih = std::log(std::pow((1.0 + y) / 2.0, 2.0));
    }
    // stable
    else {
        psih = -(4.7 * ze) / 0.74;
    }
    if (psih < -4.0) psih = -4.0;
    if (psih > 3.0) psih = 3.0;
    return psih;
}
// **  Calculate diabatic influencing factor for heat ** //  
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
// Clip Monin Obukhov length to keep it within limits
static double clipMOlength(double L, double zref, double d, double zm, double beta = 0.9)
{
    const double ln_z = std::log((zref - d) / zm);
    const double psim_min = -beta * ln_z;
    const double tol = 1e-4;
    // Compute current psim
    double psim = dpsimCpp2(zm / L) - dpsimCpp2((zref - d) / L);
    if (L < 0.0 && psim < psim_min) {
        // Need to find new L such that psim_new = psim_min
        // Bracket range for L: search between -0.001 and original L
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
    return L;  // already safe
}
// ** Calculate scaled wind profile ** //
// Calculate canopy wind profile
static std::vector<double> windprofileCpp(const vegpstruct& vegp) {
    int n = static_cast<int>(vegp.paii.size());
    if (n < 10) Rcpp::stop("Wind profile requires at least 10 layers");
    // Calculate whole canopy attenuation coefficient
    double Be = std::sqrt(0.003 + 0.1 * vegp.pai);
    double a = vegp.pai / vegp.hgt;
    double Lc = std::pow(0.25 * a, -1.0);
    double Lm = 2.0 * std::pow(Be, 3.0) * Lc;
    double at = Be * vegp.hgt / Lm;
    // Calculate attenuation coefficient for canopy elements
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
    // Adjust attenuation coefficient
    for (int i = 0; i < n; ++i) ati[i] = (ati[i] / sati) * at;
    // Calculate canopy wind shelter coefficient
    int n2 = std::trunc(n / 10);
    std::vector<double> ui(n, 1.0);
    for (int i = n - 1; i >= n2; --i) {
        ui[i - 1] = ui[i] * (1.0 - ati[i - 1]);
    }
    // Calculate bottom 10 percent
    double zm = vegp.hgt / (20.0 * n2);
    for (int i = 0; i < n2; ++i) {
        double z2 = (i + 1) * vegp.hgt / (10 * n2);
        double uf = (ka * ui[n2 - 1]) / std::log(vegp.hgt / (10.0 * zm));
        ui[i] = (uf / ka) * std::log(z2 / zm);
    }
    return ui;
}
static windmodel windmodelCpp(const std::vector<double>& wc, double uref, double hgt, double pai,
    double zref, double H = 0.0, double tc = 15, double pk = 101.3, int maxiter = 100,
    double a1 = 1.25, double psi_h = 0.0, double psi_m = 0.0, double phi_h = 1.0)
{
    if (zref < hgt) {
        Rcpp::stop("Height of wind speed measurement must be above canopy");
    }
    double Tk = tc + 273.15;
    // get base parameters
    double cp = cpairCpp(tc);
    double ph = phairCpp(tc, pk);
    double d = zeroplanedisCpp2(hgt, pai);
    double zm = roughlengthCpp2(hgt, pai, d, psi_h);
    double zh = 0.2 * zm;
    double uf = (ka * uref) / (std::log((zref - d) / zm) + psi_m);
    // Calculate safe limits for Monin Obukhov length 
    double LL = 1e99;
    if (H > 0.0) LL = (ph * cp * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
    // L negative when H positive etc
    double Lsafe = clipMOlength(LL, zref, d, zm);
    // Derive diabatic coefficients iteratively
    double dif = 10.0;
    double olduf = -100.00;
    int iter = 1;
    if (H != 0.0) {
        while (dif > 0.00000001) {
            zm = roughlengthCpp2(hgt, pai, d, psi_h);
            zh = 0.2 * zm;
            LL = (ph * cp * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
            if (H > 0) {
                if (LL < Lsafe) LL = Lsafe;
            }
            else {
                if (LL > Lsafe) LL = Lsafe;
            }
            psi_m = dpsimCpp2(zm / LL) - dpsimCpp2((zref - d) / LL);
            psi_h = dpsihCpp2(zh / LL) - dpsihCpp2((zref - d) / LL);
            uf = (ka * uref) / (std::log((zref - d) / zm) + psi_m);
            dif = std::abs(olduf - uf);
            if (iter > maxiter) dif = 0.0;
            olduf = uf;
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
    phi_h = dphihCpp2((zref - d) / LL);
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculate saturated vapour pressure * * //  
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
// Compute leaf boundary layer resistance
static double leafrHa(double tair, double dT, double uz, double len, double wid,
    double x, double rHmax = 300.0)
{
    double Tk = tair + 273.15;
    // Compute thermal diffusivity
    double Kh = (1.6667e-10 * Tk * Tk + 2.9935e-8 * Tk - 1.7128e-6);
    // Compute kinematic viscosity
    double v = 1.326 * std::pow(10.0, -5.0) * std::pow(Tk / 273.15, 1.5) * (393.55 / (Tk + 120.0));
    // Compute projected area in direction of wind flow (horizontal)
    double y = 1.0 / x;
    double Gy = std::sqrt(y * y + (1.0 - y * y)) / (y + 1.774 * std::pow((y + 1.182), -0.733));
    // Compute characteristic dimension
    double d = std::sqrt(len * wid) * Gy;
    // Compute Reynolds number
    double Re = (uz * d) / v;
    // Compute Prandlt number
    double Pr = v / Kh;
    // Compute Nusselt number for forced convection
    double Nuf;
    if (Re > 2e5) {
        Nuf = 0.37 * std::pow(Re, 0.6) * std::pow(Pr, 1.0 / 3.0) + 9.08;
    }
    else if (Re > 1000) {
        Nuf = 2.0 + (0.48 * std::pow(Re, 0.6) - 11.31) * std::pow(Pr, 1.0 / 3.0);
    }
    else {
        Nuf = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0 / 3.0);
    }
    // Compute Grashof number 
    double Gr = (g * std::pow(d, 3.0) * dT) / (Tk * v * v);
    // Compute Nusselt number for free convection
    double Nun = std::pow(0.825 + (0.387 * std::pow(Gr * Pr, 1.0 / 6.0)) /
        std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0), 2.0);
    // Compute Nusselt number (mixed forced and free)
    double Nu = std::pow(std::pow(Nuf, 3.0) + std::pow(Nun, 3.0), 1.0 / 3.0);
    // Compute boundary layer resistance
    double rHa = d / (Kh * Nu);
    if (rHa > rHmax) rHa = rHmax;
    return rHa;
}
// Compute minimum resistances for stomatal model
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
// Compute leaf stomatal conductance
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
        double A;
        double Acol;
        double cicol;
        if (C3) { // C3 photosynthetic pathway
            // Calculate Kc and Ko
            double Q10Kc = 2.1;
            double Kc = 30.0 * std::pow(Q10Kc, 0.1 * (envdata.tair - 25.0));
            double Q10Ko = 1.2;
            double Ko = 30000.0 * std::pow(Q10Ko, 0.1 * (envdata.tair - 25.0));
            // Calculate gross assimilation
            double Wc = Vcmax * ((ci - photocomp) /(ci + Kc * (1 + Oa / Ko))); // Limitating rate due to carbon
            double Wl = vegp.alpha * IPAR * ((ci - photocomp) / (ci + 2.0 * photocomp)); // Limitating rate due to light
            double We = 0.5 * Vcmax; // Limiting rate due to transport
            if (Wc < 0.0) Wc = 0.0;
            if (Wl < 0.0) Wl = 0.0;
            if (We < 0.0) We = 0.0;
            // Rate limited gross assimilation
            double W = Wc;
            if (Wl < W) W = Wl;
            if (We < W) W = We;
            A = W - Rd; // net assimilation (can be negative)
            // Co - limiting assimilation
            double Wcol = ((We + Wl) - std::sqrt(std::pow(We + Wl, 2.0) - 4.0 * 0.93 * (We * Wl)))
                / (2.0 * 0.93);  // Point at which no longer limited by ci
            Acol = Wcol - Rd;
            // Co-limiting CO2 concentration
            cicol = (-Vcmax * photocomp - Kc * (1.0 + Oa / Ko) * Wcol)
                / (Wcol - Vcmax);
        }
        else { // C4 pathway
            double k = 2.0e-4;
            double Wc = Vcmax;
            double Wl = vegp.alpha * IPAR;
            double We = k * Vcmax * (ci / (envdata.pk * 1000.0));
            double W = Wc;
            if (Wl < W) W = Wl;
            if (We < W) W = We;
            A = W - Rd; // net assimilation (can be negative)
            // Co - limiting assimilation
            double Wcol = ((Wc + Wl) - std::sqrt(std::pow(Wc + Wl, 2.0) - 4.0 * 0.83 * (We * Wl)))
                / (2.0 * 0.83);
            Acol = Wcol - Rd;
            // Co-limiting CO2 concentration
            cicol = (Wcol * envdata.pk * 1000.0) / (k * Vcmax);
        }
        // Compute change in assimulation per ci gradient
        double dadc = (A - Acol) / (ci - cicol);
        // Compute change in conductivity per psi gradient
        if (vegp.apsi < 0.0) {
            double stem_slope = 65.15 * pow(-vegp.psi50, -1.25);
            vegp.apsi = -4.0 * stem_slope / 100.0 * vegp.psi50;
        }
        double rhow = 1000.0 * (1 - (envdata.tleaf + 288.9414) * std::pow(envdata.tleaf - 3.9863, 2.0) /
            (508929.2 * (envdata.tleaf + 68.12963)));
        double psi_pd = envdata.psi_r - z * g * rhow * 1e-6;
        double K_psi_pd = 1.0 / (1.0 + std::pow(psi_pd / vegp.psi50, vegp.apsi));
        double K_50f = 0.5 / (1.0 + std::pow((psi_pd + vegp.psi50) / vegp.psi50, vegp.apsi));
        double psi_50f = (psi_pd + vegp.psi50) / 2.0;
        double dKdpKi = ((K_psi_pd - K_50f) / (psi_pd - psi_50f)) * (1.0 / K_psi_pd);
        //  # Compute rp
        double rpmin = rpmin_calc(vegp.hgt, vegp.hv, vegp.Kxmx);
        double rp = rpmin / K_psi_pd;
        // Compute zeta
        double DDm = (es - ea) / envdata.pk; // divide by pk to convert to mol / mol
        if (DDm < 0.05) DDm = 0.05;
        double zeta = 2.0 / (dKdpKi * rp * 1.6 * DDm);
        // compute limits
        if (dadc < 1e-99) dadc = 1e-99;
        double mu = 1.0 + (4.0 * zeta) / dadc;
        if (mu < 1.0) mu = 1.0;
        gs = 0.5 * dadc * (std::sqrt(mu) - 1.0);
    }
    return gs;
}
// Compute leaf vapour resistance
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
// Calculates rHa from ground to z
static double rhcanopy(double a2, double uf, double h, double z)
{
    // mu part is cheap; keep as-is
    const double mu = 1.0 / (a2 * h * uf);
    double inth;
    if (z == h) {
        // constant fallback
        inth = 4.293251 * h;
    }
    else {
        const double invh = 1.0 / h;
        const double x = pi * z * invh;
        const double s = std::sin(x);
        const double c = std::cos(x);
        const double c1 = c + 1.0;
        // constants
        const double sqrt5 = 2.2360679774997896964;          // sqrt(5)
        const double five32 = 11.180339887498948482;         // 5*sqrt(5) = 5^(3/2)
        // common subexpressions
        const double t = (sqrt5 * s) / c1;                   // argument of atan
        const double atan_term = std::atan(t);
        const double s2 = s * s;
        const double c1_2 = c1 * c1;
        const double denom = c1 * ((25.0 * s2) / c1_2 + 5.0);
        const double inner = (48.0 * atan_term) / five32 + (32.0 * s) / denom;
        inth = (2.0 * h * inner) / pi;
    }
    double rHa = inth * mu;
    // clamp
    if (rHa < 0.001) rHa = 0.001;
    return rHa;
}
// Calculate rHa from h to zref
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
// Runs plant model - one iteration (in-place update)
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
        // Woody vegetation
        double rV = 9999.99;
        if (onestepin.swaterdepth[i] > 0.0) rV = rLB[i];
        double Rabs = swout.RswLav[i] + lwout.RlwLabs[i];
        const double twood = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        double DD = (satvapCpp2(twood) - satvapCpp2(onestepin.tair[i])) * (onestepin.rh[i] / 100.0) * 1000.0;
        const double Evwood = (Mw / (RgasC * Tk)) * (DD / rV) * timestep; // surface water evaporation
        // Sunlit leaves
        envdata.PARabs = swout.RPARsun[i];
        const double gssun = leafgs(envdata, vegp, vegp.hgt / 2.0, C3);
        rV = leafrV(rLB[i], gssun, vegp.Lfrac[i], ph, onestepin.swaterdepth[i], envdata.precip);
        Rabs = swout.RswLsun[i] + lwout.RlwLabs[i];
        const double tsun = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        DD = (satvapCpp2(tsun) - satvapCpp2(onestepin.tair[i]) * (onestepin.rh[i] / 100.0)) * 1000.0;
        double rVt = rLB[i] + ph / gssun;
        const double Evsun = (Mw / (RgasC * Tk)) * (DD / rV) * timestep; // surface water evaporation
        const double Etsun = (Mw / (RgasC * Tk)) * (DD / rVt) * timestep; // transpiration
        // Shaded leaves
        envdata.PARabs = swout.RPARshade[i];
        const double gsshade = leafgs(envdata, vegp, vegp.hgt / 2.0, C3);
        rV = leafrV(rLB[i], gsshade, vegp.Lfrac[i], ph, onestepin.swaterdepth[i], envdata.precip);
        Rabs = swout.RswLshade[i] + lwout.RlwLabs[i];
        const double tshade = PenmanMonteithCpp2(Rabs, onestepin.tair[i], envdata.pk, onestepin.rh[i], vegp.vegem, rLB[i], rV, 0.0, 4);
        DD = (satvapCpp2(tshade) - satvapCpp2(onestepin.tair[i]) * (onestepin.rh[i] / 100.0)) * 1000.0;
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
    // Rain intercept model model
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
    // Write results back into onestepin
    onestepin.tleaf = std::move(tleafn);
    onestepin.swaterdepth = std::move(swaterdepthn);
    onestepin.Et = Et;
    onestepin.Hz = std::move(Hz);
    onestepin.Lz = std::move(Lz);
    onestepin.rLB = std::move(rLB);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************** Soil model *********************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
static double soilrelhumCpp(const soilpstruct& soilp, double Tsoil, double theta)
{
    double psiw = soilp.psie[0] * std::pow(theta / soilp.thetaS[0], -soilp.b[0]);
    double Tk = Tsoil + 273.15;
    double hr = std::exp(Mw * psiw / (RgasC * Tk));
    return hr;
}
static double soilsurfaceEB(const soilpstruct& soilp, double Rabs, double Tref,
    double Tsurface, double pk, double relhum, double rHa, double theta)
{
    // Net radiation
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
static double thermalConductivityCpp(
    double Vq, // Volumetric quartz fraction (m^3/m^3
    double Vm, // Volumetric mineral fraction (m^3/m^3)
    double Vo, // Volumetric organic fraction (m^3/m^3)
    double Vw, // Volmetric water fraction (m^3/m^3)
    double Mc, // Mass fraction of clay (kg/kg
    double Tc, // temperayure of air in air space
    double pk) // pressure (kPa in air space
{
    // Calculate thermal conductivity of solids
    double Vsolid = Vq + Vm + Vo;
    double kSolid = std::pow(2.5, Vq / Vsolid) * std::pow(8.8, Vm / Vsolid) * std::pow(0.25, Vo / Vsolid);
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
    double wf;
    if (Vw < 0.01 * xwo) {
        wf = 0.0;
    }
    else {
        wf = 1.0 / (1.0 + std::pow(Vw / xwo, -Q));
    }
    double kGas = 0.0242 + 0.00007 * Tc + wf * Lv * rhoAir * Dv * slope / (pk * stcor);
    double Gc = 1.0 - 2.0 * Ga;
    double kFluid = kGas + (kWater - kGas) * pow(Vw / porosity, 2);
    double ka = (2.0 / (1.0 + (kGas / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kGas / kFluid - 1.0) * Gc)) / 3.0;
    double kw = (2.0 / (1.0 + (kWater / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kWater / kFluid - 1.0) * Gc)) / 3.0;
    double ks = (2.0 / (1.0 + (kSolid / kFluid - 1.0) * Ga) +
        1.0 / (1.0 + (kSolid / kFluid - 1.0) * Gc)) / 3.0;
    double out = (kw * kWater * Vw + ka * kGas * gasPorosity + ks * kSolid * Vsolid)
        / (kw * Vw + ka * gasPorosity + ks * Vsolid);
    return out;
}
static double heatCapacityCpp(
    double Vq, // Volumetric quartz fraction (m^3/m^3
    double Vm, // Volumetric mineral fraction (m^3/m^3)
    double Vo, // Volumetric organic fraction (m^3/m^3)
    double Vw, // Volmetric water fraction (m^3/m^3)
    double Tc, // temperature of air in air space
    double pk) // pressure (kPa in air space
{
    double Vsum = Vq + Vm + Vo + Vw;
    double Va = 0.0;
    if (Vsum > 1) {
        Vq = Vq / Vsum;
        Vm = Vm / Vsum;
        Vo = Vo / Vsum;
        Vw = Vw / Vsum;
    }
    else {
        Va = 1.0 - Vsum;
    }
    double CH;
    double CHa = cpairCpp(Tc) * phairCpp(Tc, pk) / 1e6; // Specific heat of air in Mj/m^3/K
    if (Tc >= 0) {
        CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Vw * 4.18 + Va * CHa;
    }
    else {
        double CHi = 1.93 + 0.0067 * Tc; // Specific heat of ice
        CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Vw * CHi + Va * CHa;
    }
    return CH * 1e6; // Convert to J/m^3/K
}
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
// Calculate effective saturation from theta
static double SeFromTheta(const soilpstruct& soilp, double theta, int i)
{
    if (theta > soilp.thetaS[i]) return 1.0;
    return theta / soilp.thetaS[i];
}
// Calculate theta from effective saturation
static double thetaFromSe(const soilpstruct& soilp, double Se, int i)
{
    return Se * soilp.thetaS[i];
}
// Calculate degree of saturation from water potential
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
static double waterPotential(const soilpstruct& soilp, double theta, int i)
{
    double Se = SeFromTheta(soilp, theta, i);
    double psiw = soilp.psie[i] * std::pow(Se, -soilp.b[i]);
    return psiw;
}
static double thetaFromPsi(const soilpstruct& soilp, double psiw, int i)
{
    double Se = degreeOfSaturation(soilp, psiw, i);
    double theta = thetaFromSe(soilp, Se, i);
    return theta;
}
// Calculate hydraulic conductivity from theta
static double hydraulicConductivityFromTheta(const soilpstruct& soilp,
    double theta, int i) {
    double psiw = waterPotential(soilp, theta, i);
    double k = soilp.Ksat[i] * std::pow(soilp.psie[i] / psiw, soilp.n[i]);
    return k;
}
// Calculate vapour from water potential
static double vaporFromPsi(const soilpstruct& soilp, double psiw, double theta, double Tk, int i) {
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double vapor = (soilp.thetaS[i] - theta) * vp * humidity;
    return vapor;
}
// calculate change in theta with psi
static double dTheta_dPsi(const soilpstruct& soilp, double psiw, int i)
{
    double theta = soilp.thetaS[i] * degreeOfSaturation(soilp, psiw, i);
    return -theta / (soilp.b[i] * psiw);
}
static double vaporConductivityFromPsiTheta(const soilpstruct& soilp,
    double psiw, double theta, double Tk, int i)
{
    double dv = 0.000024;
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double k = 0.66 * (soilp.thetaS[i] - theta) * dv * vp * humidity * Mw / (RgasC * Tk);
    return k;
}
static double dvapor_dPsi(const soilpstruct& soilp, double psiw,
    double theta, double Tk, int i)
{
    double humidity = std::exp(Mw * psiw / (RgasC * Tk));
    double capacity_vapor = (soilp.thetaS[i] - theta) * vp * humidity *
        (Mw / (RgasC * Tk)) - dTheta_dPsi(soilp, psiw, i) * vp * humidity;
    return capacity_vapor;
}
static double evaporation_flux(const soilpstruct& soilp, double theta, double Tsurface,
    double Tair, double relhum, double rHa, double dT)
{
    double psiw = waterPotential(soilp, theta, 0) * 1000.0; // water potential (J/kg -> Pa)
    double Tk = Tair + 273.15; // Temperature deg C -> K
    double hs = std::exp(Mw * psiw / (RgasC * Tk));   // dimensionless soil effective relative humidity (0 to 1 range)
    double es = 1000.0 * satvapCpp2(Tsurface) * hs;            // kPa -> Pa
    double ea = 1000.0 * satvapCpp2(Tair) * (relhum / 100.0);  // kPa -> Pa
    double Ev = (Mw / (rHa * RgasC * Tk)) * (es - ea) * dT; // Bare soil evaporation (kg m^-2 s^-1 -> mm)
    return Ev;
}
// Calculate transpiration
static double alpha_wet(double psiw, double psie)
{
    if (psiw >= 0) return(0.0); // saturated
    if (psiw <= psie) return(1.0); // beyond air entry to well aerated
    return std::abs(psiw) / std::abs(psie);
}
static double alpha_dry(double psiw, double psi_dry, double psi_wilt) {
    if (psiw <= psi_wilt) return(0.0);
    if (psiw >= psi_dry)  return(1.0);
    return (psiw - psi_wilt) / (psi_dry - psi_wilt); // linear ramp 0 to 1
}
// Distribute root water uptake according to root fraction and water potontial
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
    while (maxdT > tolerance && nrIterations < maxNrIterations) {
        // Compute surface energy balance using midpoint temperature
        // Compute qsurface dynamically if <=10 iterations otherwise hold at average of iter 8 & 9
        if (nrIterations < 10) {
            double Tav = 0.5 * (oldTe_fixed[0] + Te_new[0]);
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
            lambda[i] = thermalConductivityCpp(Vq[i], Vm[i], Vo[i], wc[i], Mc[i], Te_new[i], atmPressure);
            CT[i] = heatCapacityCpp(Vq[i], Vm[i], Vo[i], wc[i], Te_new[i], atmPressure) * state.vol[i];
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
        // convergence: max change in the iterate
        maxdT = 0.0;
        for (int i = 0; i <= n; ++i) {
            double d = std::abs(Te_new[i] - Te_prev[i]);
            if (d > maxdT) maxdT = d;
        }
        ++nrIterations;
    }
    // compute G using fixed old state
    double G = 0.0;
    for (int i = 1; i < n; ++i) {
        double heatflux = CT[i] * (Te_new[i] - oldTe_fixed[i]) / dT;
        G += heatflux;
    }
    state.Te = Te_new;
    state.Gflux = G;
    state.iters = nrIterations;
    return state;
}
// distribute roots across soil profile
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
static soilwaterout SoilWaterCpp(soilwatermod soilmod, const soilpstruct& soilp, const climforwaterstruct& climdata,
    double dT = 3600.0, double pTAW = 0.5, int maxNrIterations = 100, double tolerance = 1e-9)
{
    // constants
    const double rho = 1000.0;
    const int n = soilp.nLayers;
    // --- FIX: freeze previous-timestep state ---
    const std::vector<double> oldtheta = soilmod.oldtheta;
    const std::vector<double> oldvapor = soilmod.oldvapor;
    // working variables (iterates)
    std::vector<double> psiw = soilmod.psiw;
    std::vector<double> theta = soilmod.theta;
    std::vector<double> vapor = soilmod.vapor;
    std::vector<double> k(n + 1);
    std::vector<double> aa(n), bb(n), cc(n), dd(n);
    std::vector<double> ff(n), u(n), du(n), Ca(n), dpsi(n);
    // bottom boundary
    if (soilp.deepSaturated) {
        psiw[n] = soilp.psie[n - 1];
        k[n] = soilp.Ksat[n - 1];
        theta[n] = soilp.thetaS[n - 1];
    }
    int iter = 0;
    double maxStep = 1e99;
    double Evapmmhr = 0.0;
    double max_Evap = 0.0;
    double Evap1;
    double Evap2;
    while (maxStep > tolerance && iter < maxNrIterations) {
        // hydraulic properties
        for (int i = 0; i < n; ++i) {
            k[i] = hydraulicConductivityFromTheta(soilp, theta[i], i) + 
                vaporConductivityFromPsiTheta(soilp, psiw[i], theta[i], soilmod.Tc[i] + 273.15, i);
            u[i] = g * k[i];
            du[i] = -u[i] * soilp.n[i] / psiw[i];
            double Cw = dTheta_dPsi(soilp, psiw[i], i);
            double Cv = dvapor_dPsi(soilp, psiw[i], theta[i], soilmod.Tc[i] + 273.15, i);
            Ca[i] = soilmod.vol[i] * (rho * Cw + Cv) / dT;
        }
        // Calculate surface flux
        if (iter < 10) {
            double theta_av = (theta[0] + soilmod.oldtheta[0]) / 2.0;
            Evapmmhr = evaporation_flux(soilp, theta_av, soilmod.Tc[0], climdata.Tair, climdata.relhum, climdata.rHa, dT);
            // Limit qsurface to value obtained in first or second iteration
            if (iter < 2 && std::abs(Evapmmhr) > std::abs(max_Evap)) max_Evap = std::abs(Evapmmhr);
            if (iter > 1) {
                if (std::abs(Evapmmhr) > std::abs(max_Evap)) Evapmmhr = std::copysign(std::abs(max_Evap), Evapmmhr);
            }
            if (iter == 8) Evap1 = Evapmmhr;
            if (iter == 9) Evap2 = Evapmmhr;
        }
        else if (iter == 10) Evapmmhr = (Evap1 + Evap2) / 2.0;
        double surfaceFlux = (Evapmmhr - climdata.precip) / dT;
        std::vector<double> STr = transpiration_distribute(soilp, soilmod.rootfrac, climdata.Et, dT, psiw, pTAW);
        // assemble system
        for (int i = 0; i < n; ++i) {
            ff[i] = ((psiw[i + 1] * k[i + 1] - psiw[i] * k[i]) / (soilmod.dz[i] * (1.0 - soilp.n[i]))) - u[i];
            if (i == 0) {
                aa[i] = 0.0;
                cc[i] = -k[i + 1] / soilmod.dz[i];
                bb[i] = k[i] / soilmod.dz[i] + Ca[i] + du[i];
                dd[i] = surfaceFlux + STr[i] - ff[i] + soilmod.vol[i] * (rho * (theta[i] - oldtheta[i]) + (vapor[i] - oldvapor[i])) / dT;
            }
            else {
                aa[i] = -k[i - 1] / soilmod.dz[i - 1] - du[i - 1];
                cc[i] = -k[i + 1] / soilmod.dz[i];
                bb[i] = k[i] / soilmod.dz[i - 1] + k[i] / soilmod.dz[i] + Ca[i] + du[i];
                dd[i] =  ff[i - 1] + STr[i] - ff[i]  + rho * soilmod.vol[i] * (theta[i] - oldtheta[i]) / dT;
            }
        }
        // solve
        Thomas TBC = ThomasBoundaryCondition(aa, bb, cc, dd, dpsi, 0, n - 1);
        dpsi = TBC.x;
        // update iterate
        maxStep = 0.0;
        for (int i = 0; i < n; ++i) {
            double psi_new = psiw[i] - dpsi[i];
            psi_new = std::min(psi_new, soilp.psie[i]);
            psi_new = std::max(psi_new, soilp.psi_min[i]);
            maxStep = std::max(maxStep, std::abs(psi_new - psiw[i]));
            psiw[i] = psi_new;
            theta[i] = thetaFromPsi(soilp, psiw[i], i);
            vapor[i] = vaporFromPsi(soilp, psiw[i], theta[i], soilmod.Tc[i], i);
        }
        ++iter;
    }
    soilmod.psiw = psiw;
    soilmod.theta = theta;
    soilmod.vapor = vapor;
    soilmod.k = k;
    soilwaterout out;
    out.swo = soilmod;
    out.success = (maxStep < tolerance);
    out.iterations = iter;
    out.Evapmmhr = Evapmmhr;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ******************************************* Below-canopy Langrangian model ********************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Run Langrangian model for one iteration
static void LangrangianOne(onestep& onestepin, double pk, double tground, double soilrelhum,
    double th, double eh, const vegpstruct& vegp,
    const std::vector<double>& z, const windmodel& windvars,
    double a0 = 0.25, double a1 = 1.25)
{
    // ** Create base ariables
    auto& tair = onestepin.tair;
    const auto& tleaf = onestepin.tleaf;
    auto& rh = onestepin.rh;
    const size_t nn = rh.size();
    const double nnd = static_cast<double>(nn);
    const double TL = windvars.a2 * vegp.hgt / windvars.uf;
    const double dz = vegp.hgt / nnd;
    // ** Calculate limits to avoid model going out of range
    // Calculate limits: temperature
    double tmn = std::min(tground, th) - 2.0, tmx = std::max(tground, th) + 2.0;
    for (size_t i = 0; i < nn; ++i) { tmn = std::min(tmn, tleaf[i]); tmx = std::max(tmx, tleaf[i]); }
    // Calculate limits: vapour pressure
    double esg = satvapCpp2(tground) * soilrelhum;  // ground vapour pressure
    double emn = std::min(eh, esg);
    double emx = std::max(eh, esg);
    for (size_t i = 0; i < nn; ++i) {
        double el = satvapCpp2(tleaf[i]);
        emn = std::min(emn, el);
        emx = std::max(emx, el);
    }
    // ** Compute vectors that are re-used
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
    // ** Compute near-field correction factor for small sample size
    const double mu = 1.0 + 0.894 * std::exp(-0.01386 * nnd) + 9.82 * std::exp(-0.15 * nnd);
    // ** Compute near-field concentrations at the top of the canopy
    double CnTh = 0.0;
    double CnLh = 0.0;
    for (size_t i = 0; i < (nn - 1); ++i) {
        // Compute kernal function
        double dz1 = (vegp.hgt - z[i]) * inowTL[i];
        double dz2 = (vegp.hgt + z[i]) * inowTL[i];
        double e = std::exp(-dz1);
        double kn = -0.39894 * std::log(1.0 - e) - 0.15623 * e;
        // Compute concentrations
        double common = kn * (dz1 + dz2);
        CnTh += ST[i] / ow[i] * common;
        CnLh += SL[i] / ow[i] * common;
    }
    CnTh *= mu;
    CnLh *= mu;
    // ** Compute far-field concentrations at the top of the canopy
    double lah = (th < 0.0) ? (51078.69 - 4.338 * th - 0.06367 * th * th)
        : (45068.7 - 42.8428 * th);
    const double phh = phairCpp(th, pk);
    const double CfTh = phh * cpairCpp(th) * th;
    const double CfLh = eh * phh * lah / pk;
    // ** Compute near and far-field concentrations for each canopy element
    double sumRH = 0.0;
    for (size_t i = 0; i < nn; ++i) {
        // Compute resistance from ground to z
        double RH = 1.0 / KH[i];
        sumRH += RH;
        double rHa = sumRH * dz;
        if (rHa < 2.0) rHa = 2.0;
        // Compute the ground sensible and latent heat flux
        double ph = phairCpp(tair[i], pk);
        double cp = cpairCpp(tair[i]);
        double GT = (ph * cp / rHa) * (tground - tair[i]) * dz;
        double ea = satvapCpp2(tair[i]) * (rh[i] / 100.0);
        double la = (tground < 0.0) ? (51078.69 - 4.338 * tground - 0.06367 * tground * tground)
            : (45068.7 - 42.8428 * tground);
        double GL = (la / (rHa * pk)) * (esg - ea) * dz;
        // Compute Hz and Lz
        double H = 0.0;
        double L = 0.0;
        for (size_t j = 0; j <= i; ++j) {
            H += ST[j];
            L += SL[j];
        }
        H += GT;
        L += GL;
        // Compute far-field concentration
        double CfT = 0.0;
        double CfL = 0.0;
        for (size_t j = i; j < nn; ++j) {
            CfT += (H / KH[j]) * dz;
            CfL += (L / KH[j]) * dz;
        }
        // Compute near-field concentration
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
        // Compute total source Concetration
        double CT = CfTh - CnTh + CfT + CnT * mu;
        double CL = CfLh - CnLh + CfL + CnL * mu;
        // Compute temperature and vapour pressure
        tair[i] = CT / (cp * ph);
        double la_i = (tleaf[i] < 0.0) ? (51078.69 - 4.338 * tleaf[i] - 0.06367 * tleaf[i] * tleaf[i])
            : (45068.7 - 42.8428 * tleaf[i]);
        double ean = (CL * pk) / (la_i * ph);
        // Impose limits and temperature and vapour pressure
        if (tair[i] > tmx) tair[i] = tmx;
        if (tair[i] < tmn) tair[i] = tmn;
        if (ean > emx) ean = emx;
        if (ean < emn) ean = emn;
        // Convert to relative humidity and impose limits
        rh[i] = (ean / satvapCpp2(tair[i])) * 100.0;
        if (rh[i] > 100.0) rh[i] = 100.0;
        if (rh[i] < 20.0)  rh[i] = 20.0;
    }
}
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
static cantop canopytop(vegpstruct& vegpc, windmodel& wind, climstruct climdata,
    std::vector<double>& Hz, std::vector<double>& Lz, double zref, 
    double Th, double eh, double tground, double soilrh, 
    double rH_g, double rH_h_zref, int maxIter, double tolerance)
{
    size_t nb = vegpc.paii.size();
    // ** Compute resistances
    const double d = zeroplanedisCpp2(vegpc.hgt, vegpc.pai);
    const double zm = roughlengthCpp2(vegpc.hgt, vegpc.pai, d, wind.psi_h);
    const double zh = 0.2 * zm;
    // ** resistance from canopy top to zref
    // from canopy hes to h
    const double psih_h = dpsihCpp2(zm / wind.LL) - dpsihCpp2((vegpc.hgt - d) / wind.LL);
    const double rH_h = (std::log((vegpc.hgt - d) / zh) + psih_h) / (ka * wind.uf);
    // ** Compute total canopy source concentrations
    double FcH = 0.0;
    double FcL = 0.0;
    for (size_t i = 0; i < nb; ++i) {
        FcH += Hz[i] * vegpc.paii[i];
        FcL += Lz[i] * vegpc.paii[i];
    }
    // ** Compute ground effective vapour pressure at reference height
    const double eground = satvapCpp2(tground) * soilrh;
    const double eref = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    // ** Compute vapour pressure at 
    // ** Set up convergence loop
    double err = 1e99;
    int nrIterations = 0;
    while (err > tolerance && nrIterations < maxIter) {
        double ph = phairCpp(Th, 101.3);
        double cp = cpairCpp(Th);
        // ** Compute ground fluxes
        double la;
        if (tground >= 0) {
            la = 45068.7 - 42.8428 * tground;
        }
        else {
            la = 51078.69 - 4.338 * tground - 0.06367 * tground * tground;
        }
        double FgH = ((ph * cp) / (rH_g + rH_h)) * (tground - Th);
        double FgL = ((la * ph) / ((rH_g + rH_h) * climdata.pk)) * (eground - eh);
        // ** Compute total flux
        double FzH = FcH + FgH;
        double FzL = FcL + FgL;
        // Compute concentrations
        double CfT = FzH * rH_h_zref;
        double CfL = FzL * rH_h_zref;
        // Convert back to temperature and vapour pressure
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
    // Set up iterations
    double tdif = 1e99; // error used to check convergence
    int nrIterations = 0;
    // Initialize variables for first run
    size_t na = onestepin.tair.size();
    double tground = onestepin.soilheatvars.Te[0];
    double Th = onestepin.tair[na - 1];
    double rcanh = onestepin.rh[na - 1];
    double eh = satvapCpp2(Th) * (rcanh / 100.0);
    double psir = 0.0;
    size_t nb = onestepin.soilwatervars.rootfrac.size();
    for (size_t i = 0; i < nb; ++i) psir += onestepin.soilwatervars.rootfrac[i] * onestepin.soilwatervars.psiw[i];
    psir = psir / 1000.0; // conversion from J/kg to MPa (1 J/kg = 1 kPa) 
    // set differences to guide convergence
    // used to check convergence
    WAitkenState st_tair;
    WAitkenState st_rh;
    WAitkenState st_leaf;
    // Leaf air temperature difference (held fixed after 3rd iteration to ensure convergence
    std::vector<double> dTs(na);
    // Resistance from z to zref
    std::vector<double> rz_zref(na);
    while ((nrIterations < 3 || tdif > tolerance) && nrIterations < maxIter) {
        // Ensure values in previous timestep stay fixed
        std::vector<double> oldTe_fixed = onestepin.soilheatvars.oldTe;
        // ** run longwave radiation model
        radmodel2 lwrad = longwavemodelCpp(wgts, climdata.Rlw, tground, soilpc.groundem, vegpc.vegem, onestepin.tleaf);
        onestepin.Rlwdown = lwrad.Rlwdown;
        onestepin.Rlwup = lwrad.Rlwup;
        // run wind model
        windmodel wind = windmodelCpp(wc, climdata.uref, vegpc.hgt, vegpc.pai, zref, onestepin.H, climdata.tref, climdata.pk, maxIter, a1,
            onestepin.psih, onestepin.psim, onestepin.phih);
        onestepin.uz = wind.uz;
        onestepin.psih = wind.psi_h;
        onestepin.psim = wind.psi_m;
        onestepin.phih = wind.phi_h;
        onestepin.LL = wind.LL;
        // ** run rain model
        rainmodel rainvars = rainintercept(wc, vegpc.pia, wind.uh, climdata.precip, climdata.winddir, vegpc.x, soilpc.slope, soilpc.aspect);
        // ** run plant model
        // create env data
        envstruct envdata{};
        envdata.pk = climdata.pk; // atmospheric pressure(kPa)
        envdata.psi_r = psir; // mean water potential in root zone(MPa)
        envdata.Ca = Ca; // CO2 ppm
        envdata.precip = climdata.precip; // precipitation (mm) 
        // Leaf air temperature difference (heald fixed after third iteration
        if (nrIterations < 3) {
            for (size_t i = 0; i < na; ++i) dTs[i] = std::abs(onestepin.tleaf[i] - onestepin.tair[i]);
        }
        std::vector<double> tleaf = onestepin.tleaf;
        // Compute resistances
        double rhg = rhcanopy(wind.a2, wind.uf, vegpc.hgt, vegpc.hgt); // rHa from ground to top of canopy
        double rhz = rh_hzref(wind, vegpc.hgt, vegpc.pai, zref); // rHa from top of canopy to zref
        double rHa = rhg + rhz; // resistance from ground to zref
        // Compute resistance from z to zref
        for (size_t i = 0; i < na; ++i) {
            rz_zref[i] = rhg - rhcanopy(wind.a2, wind.uf, vegpc.hgt, z[i]) + rhz;
        }
        // Modified inside function as passed by reference
        plantmodelCpp(onestepin, envdata, vegpc, rainvars, swrad, lwrad, z, dTs, 3600.0, C3);
        // Apply dynamic weighting between old and new to return onestep.tleaf partially weighted by old
        aitkin_weightdif(tleaf, onestepin.tleaf, z, vegpc.hgt, st_leaf);
        // Calculate additional inputs for soil model
        double Rabs = swrad.RswGabs + lwrad.RlwGabs;
        // ** Run soil heat and water model
        std::vector<double> stemp = onestepin.soilheatvars.Te;
        soilmod soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600, 0.5, maxIter);
        soilheat.iters = 0;
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = onestepin.precipground; cfw.Et = onestepin.Et;
        // Slot in soil temperature
        onestepin.soilwatervars.Tc = soilheat.Te;
        soilwaterout soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, 0.5, maxIter, tolerance);
        onestepin.witers = soilwater.iterations;
        tground = soilheat.Te[0];
        double soilrh = soilrelhumCpp(soilpc, tground, soilwater.swo.theta[0]);
        soilheat.wc = soilwater.swo.theta;
        onestepin.soilheatvars = soilheat;
        onestepin.soilwatervars = soilwater.swo;
        onestepin.Ev = soilwater.Evapmmhr;
        onestepin.theta = soilwater.swo.theta[0];
        // ** compute Sensible and Latent heat flux
        Hstruct HT = sumHCpp(climdata.tref, tground, climdata.pk, zref, z, onestepin.tleaf, rz_zref, onestepin.rLB, rHa, vegpc, wind);
        onestepin.H = HT.Htot;
        onestepin.L = swrad.RswCabs - vegpc.vegem * sb * radem(HT.Tsurf) - onestepin.H - soilheat.Gflux;
        // Compute temperature and vapour pressure at top of canopy
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
        // ** Compute air temperature and relative humidity using Langrangian model
        std::vector<double> tair = onestepin.tair;
        std::vector<double> rh = onestepin.rh;
        // Does infunction modification of original
        LangrangianOne(onestepin, climdata.pk, tground, soilrh, Th, eh, vegpc, z, wind, a0, a1);
        // Apply dynamic weighting between old and new to return onestep.tair and rh partially weighted by old
        aitkin_weightdif(tair, onestepin.tair, z, vegpc.hgt, st_tair);
        aitkin_weightdif(rh, onestepin.rh, z, vegpc.hgt, st_rh);
        // Check max air temperature difference - used to stop while loop
        tdif = 0.0;
        for (size_t i = 0; i < na; ++i) {
            double dif = std::abs(tair[i] - onestepin.tair[i]);
            if (dif > tdif) tdif = dif;
        }
        // compute psir
        psir = 0.0;
        for (size_t i = 0; i < nb; ++i) psir += soilwater.swo.rootfrac[i] * soilwater.swo.psiw[i];
        psir = psir / 1000.0; // conversion from J/kg to MPa (1 J/kg = 1 kPa) 
        // Ensure oldTe not updated
        onestepin.soilheatvars.oldTe = oldTe_fixed;
        ++nrIterations;
    }
    // Update variables
    onestepin.iters = nrIterations;;
    onestepin.error = tdif;
    return onestepin;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ***************************************** Above canopy ************************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Derive temperature above canopy by extrapolating Langrangian profile. za can be > zref, but must be less than h
// [[Rcpp::export]]
double Tabove(double za, double zref, double th, double tref, double hgt, double pai)
{
    double d = 0.0;
    if (hgt > 0.0) d = zeroplanedisCpp2(hgt, pai);
    double Tz = tref + (th - tref) * std::log((za - d) / (zref - d)) / std::log((hgt - d) / (zref - d));
    return Tz;
}
// Derive humidity above canopy by extrapolating Langrangian profile. za can be > zref, but must be less than h
// [[Rcpp::export]]
double RHabove(double za, double zref, double rh, double th, double tref, double tz, double relhum, double hgt, double pai)
{
    double d = 0.0;
    if (hgt > 0.0) d = zeroplanedisCpp2(hgt, pai);
    const double eh = satvapCpp2(th) * (rh / 100.0);
    const double eref = satvapCpp2(tref) * (relhum / 100.0);
    const double ez = eref + (eh - eref) * std::log((za - d) / (zref - d)) / std::log((hgt - d) / (zref - d));
    double rz = (ez / satvapCpp2(tz)) * 100.0;
    if (rz > 100.0) rz = 100.0;
    return rz;
}
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
// Run the model when there is no vegetation
static onestepbare OneStepBare(onestepbare onestepin, const obsstruct& obsdata, const climstruct& climdata, const soilpstruct& soilpc,
    const std::vector<double>& z, double latr, double lonr, double zref, double zm = 0.004, int maxIter = 100, 
    double  tolerance = 1e-8) 
{
    // ** Calculate Rabs
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
    // Calculate base parameters
    double Tk = climdata.tref + 273.15;
    double cp = cpairCpp(climdata.tref);
    double ph = phairCpp(climdata.tref, climdata.pk);
    double cpph = cp * ph;
    // Get values that require update on each iteration
    double Ts = onestepin.soilheatvars.Te[0];
    double psi_m = onestepin.psim;
    double psi_h = onestepin.psih;
    double H = onestepin.H;
    double zmd = zm * std::exp(ka * psi_h);
    double zh = 0.2 * zmd;
    double uf = (ka * climdata.uref) / (std::log(zref / zmd) + psi_m);
    // Calculate safe limits for Monin Obukhov length 
    double LL = (cpph * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
    // L negative when H positive etc
    double Lsafe = clipMOlength(LL, zref, 0.0, zmd);
    double dif = 1e99;
    int nrIterations = 0;
    // Initialize soil heat and water model
    soilmod soilheat;
    soilwaterout soilwater;
    double rHa;
    while (dif > tolerance && nrIterations < maxIter) {
        // ** Calculate rHa
        if (H != 0.0) {
            LL = (cpph * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
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
        // ** Run soil heat model
        soilheat = SoilHeatCpp(onestepin.soilheatvars, soilpc, Rabs, climdata.tref, climdata.relhum, climdata.pk, rHa, 3600, 0.5, maxIter);
        // Construct climate input to water model
        climforwaterstruct cfw = {};
        cfw.Rabs = Rabs; cfw.Tair = climdata.tref; cfw.relhum = climdata.relhum; cfw.pk = climdata.pk; cfw.rHa = rHa;
        cfw.precip = climdata.precip; cfw.Et = 0.0;
        // Slot in soil temperature
        onestepin.soilwatervars.Tc = soilheat.Te;
        // Run water model
        soilwater = SoilWaterCpp(onestepin.soilwatervars, soilpc, cfw, 3600, 0.5, maxIter, tolerance);
        onestepin.soilheatvars.wc = soilwater.swo.theta;
        // Recalculate H
        dif = std::abs(Ts - soilheat.Te[0]);
        Ts = soilheat.Te[0];
        H = (cpph / rHa) * (Ts - climdata.tref);
        ++nrIterations;
    }
    // Compute L
    double soilrh = soilrelhumCpp(soilpc, Ts, soilwater.swo.theta[0]);
    double eg = satvapCpp2(Ts) * soilrh;
    double ea = satvapCpp2(climdata.tref) * (climdata.relhum / 100.0);
    double la = (Ts < 0.0) ? (51078.69 - 4.338 * Ts - 0.06367 * Ts * Ts) : (45068.7 - 42.8428 * Ts);
    onestepin.L = ((la * ph) / (climdata.pk * rHa)) * (eg - ea);
    // Compute tair, relhum and wind profiles
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
    // ** Calculate radiation streams
    onestepin.Rb0 = Rb0;
    onestepin.Rswup = (1.0 - soilpc.gref) * climdata.Rsw;
    onestepin.Rlwup = soilpc.groundem * sb * radem(Ts);
    // ** Update onstepin
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ****************************************** Ecotherm model *********************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
// Compute animal boundary layer resistance
static double animalrHa(double tair, double dT, double uz, double d, double rHmax = 300.0)
{
    double Tk = tair + 273.15;
    // Compute thermal diffusivity
    double Kh = (1.6667e-10 * Tk * Tk + 2.9935e-8 * Tk - 1.7128e-6);
    // Compute kinematic viscosity
    double v = 1.326 * std::pow(10.0, -5.0) * std::pow(Tk / 273.15, 1.5) * (393.55 / (Tk + 120.0));
    // Compute Reynolds number
    double Re = (uz * d) / v;
    // Compute Prandlt number
    double Pr = v / Kh;
    // Compute Nusselt number for forced convection
    double Nuf;
    if (Re > 2e5) {
        Nuf = 0.37 * std::pow(Re, 0.6) * std::pow(Pr, 1.0 / 3.0) + 9.08;
    }
    else if (Re > 1000) {
        Nuf = 2.0 + (0.48 * std::pow(Re, 0.6) - 11.31) * std::pow(Pr, 1.0 / 3.0);
    }
    else {
        Nuf = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0 / 3.0);
    }
    // Compute Grashof number 
    double Gr = (g * std::pow(d, 3.0) * dT) / (Tk * v * v);
    // Compute Nusselt number for free convection
    double Nun = std::pow(0.825 + (0.387 * std::pow(Gr * Pr, 1.0 / 6.0)) /
        std::pow(1.0 + std::pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0), 2.0);
    // Compute Nusselt number (mixed forced and free)
    double Nu = std::pow(std::pow(Nuf, 3.0) + std::pow(Nun, 3.0), 1.0 / 3.0);
    // Compute boundary layer resistance
    double rHa = d / (Kh * Nu);
    if (rHa > rHmax) rHa = rHmax;
    return rHa;
}
// Calculate water diffusivity
double Dw_waterVapour(double Tf, double pk)
{
    const double D0 = 2.26e-5;   // m^2 s^-1 at 273.15 K and 101325 Pa
    const double T0 = 273.15;    // reference temperature (K)
    const double P0 = 101325.0;  // reference pressure (Pa)
    double Tfk = Tf + 273.15;
    double Pa = pk * 1000.0;
    return D0 * std::pow(Tfk / T0, 1.75) * (P0 / Pa);
}
// Compute body temperature of animal
// [[Rcpp::export]]
double PenmanMonteith_animal(double Rabs, double Ta, double Ts, double Te, double Tf, double pk, double rh,
    double rHa, double height, double wetfrac, double confrac, double M, double em = 0.97, double k = 0.5,
    double surfrh = 1.0)
{
    double Tbody = Ta;
    double cp = cpairCpp(Ta);
    double ph = phairCpp(Ta, pk);
    double dT = Ts - Ta; // difference contact surface and air temperature
    double cd = 1.0 - confrac;
    double ea = satvapCpp2(Ta);
    double Da = ea * (1.0 - rh / 100.0); // Vapour pressure deficit
    double Dac = surfrh * satvapCpp2(Tf) - ea; // contact Vapour pressure deficit
    // latent heat exchange through a contact interface
    double keff = k / (0.5 * height);
    // Calculate mus
    double mr = cd * em * sb; // emmited ratiation mu
    double mrc = confrac * em * sb; // emmited ratiation mu
    double mc = keff * confrac; // conductance my
    double mv = ph * cp * cd / rHa;
    double la = (Tbody >= 0.0) ? 45068.7 - 42.8428 * Tbody : 51078.69 - 4.338 * Tbody - 0.06367 * Tbody * Tbody;
    double ml = (ph * la * wetfrac * cd) / (rHa * pk);
    double Dw = Dw_waterVapour(Tf, pk);
    double mlc = (la * Dw * wetfrac * confrac * 1000.0) / (0.5 * height * (Ts + 273.15) * RgasC);
    double DeV = satvapCpp2(Te + 0.5) - satvapCpp2(Te - 0.5); // Slope of the saturated vapour pressure curve (Lv)
    double DeC = satvapCpp2(Tf + 0.5) - satvapCpp2(Tf - 0.5); // Slope of the saturated vapour pressure curve (Lc)
        // Calculate various contributions
    double top = (Rabs + M - mr * radem(Ta) - mrc * radem(Ta) - mc * dT - ml * Da - mlc * Dac);
    double btm = 4.0 * mr * std::pow(Te + 273.15, 3.0) + 4.0 * mrc * std::pow(Tf + 273.15, 3.0) +
        mc + mv + ml * DeV + surfrh * mlc * DeC;
    Tbody = Ta + top / btm;
    return (Tbody);
}
static inline void aitkin_weightdif_scalar(
    double oldv,
    double& newv,   // raw new value on input, updated in place
    WAitkenStateScalar& st,
    double backweight_min = 0.10,
    double backweight_max = 0.98
)
{
    // Residual for current iteration
    double r = newv - oldv;
    // Convert stored omega limits to backweight limits
    // backweight = 1 - omega
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
// Computes animal metabolic rate in W
double metabolic_rate(double volume, double rho, double Q10, double a0, double b, double Tref, double Tbody)
{
    // compute mass
    double m = volume * rho;
    // Compute metabolic rate
    double M = a0 * std::pow(m, b) * std::pow(Q10, (Tbody - Tref) / 10.0);
    return M;
}
// Position details:
// one of fixed (constant adir and atilt assumed), 
// max (animal seeks to maximize radiation absorption)
// min (animal seeks to minimize radiation absorption)
// randomdir (animal maintains tilt, but adir is random)
// random (animal has random atilt and adir)
// [[Rcpp::export]]
Rcpp::NumericVector Ectotherm(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List animal,
    double lat, // Latitude (decimal degrees)
    double lon, // longitude decimal degrees
    double leaft, // leaf transmittance
    int maxIter = 100, // max number of iterations
    double tolerance = 1e-2) // error tolerence for convergence
{
    // Extract from obstime data.frame
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
    double skinwetfrac = animal["skinwetfrac"]; // fraction of skin surface acting like a freely evaporating surface (0 for most animals, 1 for amphibians)
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
    // Varibles needed
    size_t n = Ta.size();
    double latr = lat * torad;
    double lonr = lon * torad;
    Rcpp::NumericVector Tbody(n);
    for (size_t i = 0; i < n; ++i) {
        Tbody[i] = Ta[i] + 3.0;
        double dT = Tbody[i] - Ta[i]; // initial animal air temperature difference
        // Calculate solar position
        solmodel solp = solpositionCpp2(latr, lonr, year[i], month[i], day[i], hour[i]);
        // Calculate shortwave radiation absorption
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
        // Calculate total radiation absorption
        double Rabs = (1.0 - refl) * Rsw_flux + 0.5 * em * (Rlwdown[i] + Rlwup[i]);
        // Calculate characteristic dimension
        double danim = chardim(wdir[i], solp.zenr, solp.azir, height, width, length, adir, atilt, position);
        double tdif = 1e99;
        WAitkenStateScalar st;
        int nrIterations = 0;
        while (tdif > tolerance && nrIterations < maxIter) {
            // Calculate average temperatures
            double Te = (Tbody[i] + Ta[i]) / 2.0;
            double Tf = (Tbody[i] + Ts[i]) / 2.0;
            // Calculate rHa
            double rHa = animalrHa(Ta[i], std::abs(dT), uz[i], danim, 300.0);
            // Calculate Metabolic rate
            double M = metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody[i]) / area; 
            // Calculate body temperature
            double Told = Tbody[i];
            Tbody[i] = PenmanMonteith_animal(Rabs, Ta[i], Ts[i], Te, Tf, pk[i], rh[i], rHa, height,
                skinwetfrac, confrac, M, em, k, surfrh[i]);
            aitkin_weightdif_scalar(Told, Tbody[i], st);
            double dTold = dT;
            dT = Tbody[i] - Ta[i];
            tdif = std::abs(dTold - dT);
            ++nrIterations;
        }
    }
    return Tbody;
}
// Run Ectotherm model using Matrix Inputs
// [[Rcpp::export]]
Rcpp::NumericVector EctothermM(Rcpp::DataFrame obstime, Rcpp::List climdata, Rcpp::List animal,
    double lat, // Latitude (decimal degrees)
    double lon, // longitude decimal degrees
    double leaft, // leaf transmittance
    int maxIter = 100, // max number of iterations
    double tolerance = 1e-2) // error tolerence for convergence
{
    // Extract from obstime data.frame
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Extract from climdata data.frame
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
    // Extract animal parameters
    double height = animal["height"]; // height of animal (cm)
    double width = animal["width"]; // width of animal (cm)
    double length = animal["len"]; // length of animal (cm)
    double refl = animal["refl"]; // reflectance of animal (0-1)
    double confrac = animal["confrac"]; // fraction of animal in direct contact with surface
    double skinwetfrac = animal["skinwetfrac"]; // fraction of skin surface acting like a freely evaporating surface (0 for most animals, 1 for amphibians)
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
    // Get number of time steps and vertical segments
    int tsteps = Ta.nrow();
    int nlayers = Ta.ncol();
    // Varibles needed
    Rcpp::NumericMatrix Tbody(tsteps, nlayers);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (int i = 0; i < tsteps; ++i) { 
        // Calculate solar position
        solmodel solp = solpositionCpp2(latr, lonr, year[i], month[i], day[i], hour[i]);
        // Calculate solar coefficient
        double sc = 1.0;
        if (solp.zenr < pi / 2.0) {
            silstruct sa = silhouette(solp.zenr, solp.azir, height, width, length, adir, atilt, position);
            sc = sa.silA / sa.A; // solar coefficient
        }
        // Calculate characteristic dimension
        double danim = chardim(wdir[i], solp.zenr, solp.azir, height, width, length, adir, atilt, position);
        for (int j = 0; j < nlayers; ++j) {
            // Calculate shortwave radiation absorption
            double Rsw_flux = 0.0;
            if (solp.zenr < pi / 2.0) {
                if (atilt < 0.0) { // animal suspended below leaf
                    Rsw_flux = leaft * (sc * Rdirdown(i,j) + 0.5 * Rdifdown(i, j)) + 0.5 * Rswup(i, j);
                }
                else {
                    Rsw_flux = sc * Rdirdown(i, j) + 0.5 * Rdifdown(i, j) + 0.5 * Rswup(i, j);
                }
            }
            // Calculate total radiation absorption
            double Rabs = (1.0 - refl) * Rsw_flux + 0.5 * em * (Rlwdown(i, j) + Rlwup(i, j));
            // Assume initial values
            Tbody(i, j) = Ta(i, j) + 3.0;
            double dT = Tbody(i, j) - Ta(i, j); // initial animal air temperature difference
            double tdif = 1e99;
            WAitkenStateScalar st;
            int nrIterations = 0;
            // Iteratively calculate body temperature
            while (tdif > tolerance && nrIterations < maxIter) {
                // Calculate average temperatures
                double Te = (Tbody(i, j) + Ta(i, j)) / 2.0;
                double Tf = (Tbody(i, j) + Ts(i, j)) / 2.0;
                // Calculate rHa
                double rHa = animalrHa(Ta(i, j), std::abs(dT), uz(i, j), danim, 300.0);
                // Calculate Metabolic rate
                double M = metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody(i, j)) / area;
                // Calculate body temperature
                double Told = Tbody(i, j);
                Tbody(i, j) = PenmanMonteith_animal(Rabs, Ta(i, j), Ts(i, j), Te, Tf, pk[i], rh(i, j), rHa, height,
                    skinwetfrac, confrac, M, em, k, 1.0);
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Calculate cumulative paii
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
    out.n = Rcpp::as<std::vector<double>>(soilp["n"]);
    out.deepSaturated = Rcpp::as<bool>(soilp["deepSaturated"]);
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
static soilwatermod toSoilwatermod(soilpstruct soilpc, soilmod state, double rootskew)
{
    soilwatermod out;
    out.n = state.n; // number of layers
    out.z = state.z; // node depths
    out.dz = state.dz; // layer thickness
    out.zCenter = state.zCenter; // centre between nodes
    out.Tc = state.Te; // temperature of soil profile
    out.theta = state.wc; // volumetric water content of soil profile
    // Initialize
    std::vector<double> vol(out.n + 1);
    std::vector<double> psiw(out.n + 1);
    std::vector<double> k(out.n + 1);
    std::vector<double> vapor(out.n + 1);
    vol[0] = 0.0;
    for (int i = 0; i <= out.n; ++i) {
        if (i > 0) vol[i] = (out.z[i + 1] - out.z[i - 1]) / 2.0;
        psiw[i] = waterPotential(soilpc, out.theta[i], i);
        k[i] = hydraulicConductivityFromTheta(soilpc, out.theta[i], i);
        vapor[i] = vaporFromPsi(soilpc, psiw[i], out.theta[i], out.Tc[i] + 273.15, i);
    }
    // assign to out
    out.vol = vol;
    out.psiw = psiw;
    out.vapor = vapor;
    out.k = k;
    out.oldvapor = vapor;
    out.oldtheta = out.theta;
    double totalDepth = 0.0;
    for (int i = 0; i <= out.n; ++i) {
        if (out.z[i] > totalDepth) totalDepth = out.z[i];
    }
    out.rootfrac = root_distribute(out.dz, totalDepth, rootskew);
    return out;
}
// Function for initalising one step
static onestep CreateOneStep(soilmod soilheatvars, soilwatermod soilwatervars, double temp0, double rh0, double Rlw0, int na) {
    // ** Create input vectors for onestep initialization
    std::vector<double> v0(na, 0.0);
    std::vector<double> v1(na, 1.0);
    std::vector<double> Rlwv(na, Rlw0);
    std::vector<double> tv(na, temp0);
    std::vector<double> rv(na, rh0);
    // ** initialize onestep
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
// Wrapper function for returning profile with bare ground
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
// [[Rcpp::export]]
List profileR(size_t hourtoplot, DataFrame obstime, DataFrame climdata, List soilc, List vegp,
    std::vector<double> paii20, std::vector<double> paii, std::vector<double> Lfrac20, std::vector<double> Lfrac, 
    double zref, double Ca, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double tolerance = 1e-3, double a0 = 0.25, double a1 = 1.25, bool C3 = true)
{
    // ** ------------------ Run model up to hour with 20 layers ------------------------- ** //
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    // ** Convert vegp and soilc
    vegpstruct vegpc = tovegpstruct(vegp, paii20, Lfrac20);
    soilpstruct soilpc = tosoilpstruct(soilc);
    // ** Create additional inputs
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc);
    // ** Compute z
    size_t na = paii20.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    // ** Initialize soil heat and water model
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    // ** Onestep initialization
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    double latr = lat * torad;
    double lonr = lon * torad;
    if (hourtoplot > 0) {
        for (size_t hr = 0; hr < hourtoplot; ++hr) {
            // Create input data
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
    // ** Create new vegp and other time-invarient variables
    vegpc = tovegpstruct(vegp, paii, Lfrac);
    wgts = lwradweights(vegpc.paii);
    wc = windprofileCpp(vegpc);
    // ** Run model
    // Create input data
    size_t hr = hourtoplot;
    obsstruct obsdata;
    climstruct climin;
    obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
    climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
    climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
    climin.winddir = wdir[hr]; climin.precip = precip[hr];
    // ** Recompute z
    na = paii.size();
    std::vector<double> z2(na);
    for (size_t i = 0; i < na; ++i) z2[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    // Initialize onestepin
    onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[hr], relhum[hr], Rlw[hr], na);
    onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z2, tspveg, tspvegPAR, tspdif,
        tspdifPAR, wgts, wc, Ca, latr, lonr, zref, maxNrIterations, tolerance, a0, a1, C3);
    // Convert to list and return
    return OneStepCpptoList(onestepin, z2);
}
// Wrapper function for running the full model with bare ground
// [[Rcpp::export]]
DataFrame RunBareR(double reqhgt, DataFrame obstime, DataFrame climdata, List soilc, double zref, 
    double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    double zm = 0.004, int maxNrIterations = 100, double  tolerance = 1e-3)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    // ** Convert soilc
    soilpstruct soilpc = tosoilpstruct(soilc);
    // ** Derive additional variables
    size_t n = temp.size();
    // ** Compute z
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
        // Create input data
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
        // Update inputs
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
// Wrapper function for running the full model
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
                tair[hr] = Tabove(reqhgt, zref, Th, temp[hr], vegpc.hgt, vegpc.pai);
                tleaf[hr] = -999.99;
                rz[hr] = RHabove(reqhgt, zref, rh, Th, temp[hr], tair[hr], relhum[hr], vegpc.hgt, vegpc.pai);
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
// Wrapper function for performing height adjustment to weather data
// [[Rcpp::export]]
List WeatherhgtCpp2(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zin, double zout, double lat, double lon, 
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, double CO2ppm = 430.0)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    // ** Convert vegp and soilc
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    size_t n = temp.size();
    // ** Compute z
    size_t na = vegpc.paii.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    // ** Initialize soil heat and water models 
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    // ** Create additional inputs
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc);
    // ** Onestep initialization
    onestep onestepin = CreateOneStep(soilheatvars, soilwatervars, temp[0], relhum[0], Rlw[0], na);
    // Create output variables for storing
    std::vector<double> temp_new(n);
    std::vector<double> relhum_new(n);
    std::vector<double> windspeed_new(n);
    std::vector<double> pres_new(n);
    // Run model for every hr
    double d = zeroplanedisCpp2(vegpc.hgt, vegpc.pai);
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        // Create input data
        obsstruct obsdata;
        climstruct climin;
        obsdata.year = year[hr]; obsdata.month = month[hr]; obsdata.day = day[hr]; obsdata.hour = hour[hr];
        climin.tref = temp[hr]; climin.relhum = relhum[hr]; climin.pk = pres[hr];
        climin.Rsw = Rsw[hr]; climin.Rdif = Rdif[hr]; climin.Rlw = Rlw[hr]; climin.uref = wspeed[hr];
        climin.winddir = wdir[hr]; climin.precip = precip[hr];
        onestepin = OneStepBelow(onestepin, obsdata, climin, vegpc, soilpc, z, tspveg, tspvegPAR, tspdif,
            tspdifPAR, wgts, wc, CO2ppm, latr, lonr, zin, 100, 1e-2);
        // height adjust temperature and relative humidity
        double th = onestepin.tair[na - 1];
        double rh = onestepin.rh[na - 1];
        temp_new[hr] = Tabove(zout, zin, th, temp[hr], vegpc.hgt, vegpc.pai);
        relhum_new[hr] = RHabove(zout, zin, rh, th, temp[hr], temp_new[hr], relhum[hr], vegpc.hgt, vegpc.pai);
        // height adjust windspeed;
        double zm = roughlengthCpp2(vegpc.hgt, vegpc.pai, d, onestepin.psih);
        double uf = (ka * wspeed[hr]) / (std::log((zin - d) / zm) + onestepin.psim);
        double psi_m = dpsimCpp2(zm / onestepin.LL) - dpsimCpp2((zout - d) / onestepin.LL);
        windspeed_new[hr] = (uf / ka) * (std::log((zout - d) / zm) + psi_m);
        // height adjust pressure
        double Tv = ((temp_new[hr] + climin.tref) / 2.0) + 273.15;
        pres_new[hr] = pres[hr] * std::exp(-(g * (zout - zin)) / (287.05 * Tv));
        // Update inputs
        onestepin.soilheatvars.oldTe = onestepin.soilheatvars.Te;
        onestepin.soilwatervars.oldtheta = onestepin.soilwatervars.theta;
        onestepin.soilwatervars.oldvapor = onestepin.soilwatervars.vapor;
    }
    return List::create(
        // Obstime
        Named("new_temp") = Rcpp::wrap(temp_new),
        Named("new_relhum") = Rcpp::wrap(relhum_new),
        Named("new_windspeed") = Rcpp::wrap(windspeed_new),
        Named("new_pressure") = Rcpp::wrap(pres_new)
    );
}

void FillMyMatrix(NumericMatrix& mat, std::vector<double>& v, size_t hr)
{
    for (size_t i = 0; i < v.size(); ++i) {
        mat(hr, i) = v[i];
    }
}
// Wrapper function for running the full model for bare ground
// [[Rcpp::export]]
List RunBelowFullBare(DataFrame obstime, DataFrame climdata, List soilc, std::vector<double> z, 
    double zref, double zm, double lat, double lon, std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, 
    int maxNrIterations = 100, double  tolerance = 1e-2)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    // ** Convert soilc
    soilpstruct soilpc = tosoilpstruct(soilc);
    // ** Derive additional variables
    size_t n = temp.size();
    size_t na = z.size();
    // ** Initialize onstepin
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
        // Create input data
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
        // Update inputs
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
// Wrapper function for running the full model below canopy
// [[Rcpp::export]]
List RunBelowFull(DataFrame obstime, DataFrame climdata, List soilc, List vegp, std::vector<double> paii,
    std::vector<double> Lfrac, double zref, double Ca, double lat, double lon,
    std::vector<double> SoilTempIni, std::vector<double> SoilThetaIni, int maxNrIterations = 100,
    double  tolerance = 1e-2, double a0 = 0.25, double a1 = 1.25, bool C3 = true)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> precip = climdata["precip"];
    // ** Convert vegp and soilc
    vegpstruct vegpc = tovegpstruct(vegp, paii, Lfrac);
    soilpstruct soilpc = tosoilpstruct(soilc);
    // ** Create additional inputs
    tsvegstruct  tspveg = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lref, vegpc.ltra, soilpc.gref);
    tsvegstruct  tspvegPAR = twostreamvegCpp(vegpc.pai, vegpc.x, vegpc.lrefp, vegpc.ltrap, soilpc.gref);
    tsdifstruct tspdif = twostreamdifCpp(tspveg);
    tsdifstruct tspdifPAR = twostreamdifCpp(tspvegPAR);
    LWweights wgts = lwradweights(vegpc.paii);
    std::vector<double> wc = windprofileCpp(vegpc);
    size_t n = temp.size();
    // ** Compute z
    size_t na = paii.size();
    std::vector<double> z(na);
    for (size_t i = 0; i < na; ++i) z[i] = (static_cast<double>(i + 1) / static_cast<double>(na)) * vegpc.hgt;
    // ** Initialize soil heat and water model
    soilmod soilheatvars = toSoilheatmod(soilc, SoilTempIni, SoilThetaIni);
    soilwatermod soilwatervars = toSoilwatermod(soilpc, soilheatvars, vegpc.rootskew);
    // ** Onestep initialization
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
        // Create input data
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
        // Update inputs
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
std::vector<double> clearskyradCpp2(DataFrame obstime, DataFrame climdata, double lat, double lon)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // ** Access columns of climdata
    std::vector<double> temp = climdata["temp"];
    std::vector<double> relhum = climdata["relhum"];
    std::vector<double> pres = climdata["pres"];
    // ** compute solar 
    double latr = lat * torad;
    double lonr = lon * torad;
    size_t n = hour.size();
    std::vector<double> csr(n, 0.0);
    for (size_t hr = 0; hr < n; ++hr) {
        solmodel solp = solpositionCpp2(latr, lonr, year[hr], month[hr], day[hr], hour[hr]);
        if (solp.zenr <= pi / 2.0) {
            double m = 35 * std::cos(solp.zenr) * std::pow(1224 * std::cos(solp.zenr) * std::cos(solp.zenr) + 1, -0.5);
            double TrTpg = 1.021 - 0.084 * std::sqrt(m * 0.00949 * pres[hr] + 0.051);
            double xx = std::log(relhum[hr] / 100.0) + ((17.27 * temp[hr]) / (237.3 + temp[hr]));
            double Td = (237.3 * xx) / (17.27 - xx);
            double u = std::exp(0.1133 - std::log(3.78) + 0.0393 * Td);
            double Tw = 1.0 - 0.077 * std::pow(u * m, 0.3);
            double Ta = 0.935 * m;
            double od = TrTpg * Tw * Ta;
            csr[hr] = 1352.778 * std::cos(solp.zenr) * od;
        }
    }
    return csr;
}
// ** Calculates diffuse fraction ** //
// [[Rcpp::export]]
std::vector<double> difpropCpp(DataFrame obstime, std::vector<double> swrad, double lat, double lon)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    size_t n = hour.size();
    std::vector<double> dp(n, 1.0);
    // ** compute solar 
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
// [[Rcpp::export]]
std::vector<double> solaltCpp(DataFrame obstime, double lat, double lon)
{
    // ** Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    size_t n = hour.size();
    std::vector<double> sa(n, 1.0);
    // ** compute solar 
    double latr = lat * torad;
    double lonr = lon * torad;
    for (size_t hr = 0; hr < n; ++hr) {
        solmodel solp = solpositionCpp2(latr, lonr, year[hr], month[hr], day[hr], hour[hr]);
        double zd = solp.zenr * 180.0 / pi;
        sa[hr] = 90.0 - zd;
    }
    return sa;
}
// Code for spline interpolating a Numeric Matrix for better visual representation
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
