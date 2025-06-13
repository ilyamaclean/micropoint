#include <Rcpp.h>
#include "micropointheaders.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
using namespace Rcpp;
// ** Calculates Astronomical Julian day ** //
int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// ** Calculates solar time ** //
double soltimeCpp(int jd, double lt, double lon)
{

    double m = 6.24004077 + 0.01720197 * (jd - 2451545);
    double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
    double st = lt + (4 * lon + eot) / 60;
    return st;
}
// ** Calculates solar position ** //
// [[Rcpp::export]]
std::vector<double> solpositionCpp(double lat, double lon, int year, int month, int day, double lt)
{
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lon);
    // Calculate solar zenith (degrees)
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh) * (180 / M_PI);
    // Calculate solar azimuth (degrees)
    double sh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double hh = atan(sh / sqrt(1 - sh * sh));
    double sazi = cos(dec) * sin(tt) / cos(hh);
    double cazi = (sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec)) /
        sqrt(pow(cos(dec) * sin(tt), 2) + pow(sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec), 2));
    double sqt = 1 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    double azi = 180 + (180 * atan(sazi / sqrt(sqt))) / M_PI;
    if (cazi < 0) {
        if (sazi < 0) {
            azi = 180 - azi;
        }
        else {
            azi = 540 - azi;
        }
    }
    // Define and return output variable
    std::vector<double> solpos(2, 0.0);
    solpos[0] = z;
    solpos[1] = azi;
    return solpos;
}
// ** Calculate solar index ** //
// [[Rcpp::export]]
double solarindexCpp(double slope, double aspect, double zen, double azi, bool shadowmask = false)
{
    double si;
    if (zen > 90.0 && !shadowmask) {
        si = 0;
    }
    else {
        if (slope == 0.0) {
            si = cos(zen * M_PI / 180);
        }
        else {
            si = cos(zen * M_PI / 180) * cos(slope * M_PI / 180) + sin(zen * M_PI / 180) *
                sin(slope * M_PI / 180) * cos((azi - aspect) * M_PI / 180);
        }
    }
    if (si < 0.0) si = 0.0;
    return si;
}
// ** Calculate clear sky radiation ** //
// [[Rcpp::export]]
std::vector<double> clearskyradCpp(std::vector<int> year, std::vector<int> month, std::vector<int> day,
    std::vector<double> lt, double lat, double lon, std::vector<double> tc, std::vector<double> rh,
    std::vector<double> pk)
{
    std::vector<double> Ic(year.size());
    for (size_t i = 0; i < Ic.size(); ++i) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
        double zen = sp[0];
        double z = zen * M_PI / 180.0;
        if (zen <= 90.0) {
            double m = 35 * cos(z) * pow(1224 * cos(z) * cos(z) + 1, -0.5);
            double TrTpg = 1.021 - 0.084 * sqrt(m * 0.00949 * pk[i] + 0.051);
            double xx = log(rh[i] / 100) + ((17.27 * tc[i]) / (237.3 + tc[i]));
            double Td = (237.3 * xx) / (17.27 - xx);
            double u = exp(0.1133 - log(3.78) + 0.0393 * Td);
            double Tw = 1 - 0.077 * pow(u * m, 0.3);
            double Ta = 0.935 * m;
            double od = TrTpg * Tw * Ta;
            Ic[i] = 1352.778 * cos(z) * od;
        }
    }
    return Ic;
}
// ** Calculate canopy extinction coefficient for sloped ground surfaces ** //
std::vector<double> cankCpp(double zen, double x, double si) {
    double k;
    if (zen > 90.0) zen = 90.0;
    if (si < 0.0) si = 0.0;
    double Z = zen * M_PI / 180.0;
    // Calculate normal canopy extinction coefficient
    if (x == 1.0) {
        k = 1 / (2 * cos(Z));
    }
    else if (std::isinf(x)) {
        k = 1.0;
    }
    else if (x == 0.0) {
        k = tan(Z);
    }
    else {
        k = sqrt(x * x + (tan(Z) * tan(Z))) / (x + 1.774 * pow((x + 1.182), -0.733));
    }
    if (k > 6000.0) k = 6000.0;
    // Calculate adjusted k
    double kd = k * cos(Z) / si;
    if (si == 0) kd = 1.0;
    double Kc = 1.0 / si;
    if (si == 0.0) Kc = 600.0;
    std::vector<double> kparams(3, 0.0);
    kparams[0] = k;
    kparams[1] = kd;
    kparams[2] = Kc;
    return kparams;
}
// ** Calculates parameters for diffuse radiation using two-stream model ** //
std::vector<double> twostreamdifCpp(double pait, double om, double a, double gma, double h, double gref)
{
    // Adjust leaf area for clumping factor
    // Calculate base parameters
    double S1 = exp(-h * pait);
    double u1 = a + gma * (1 - 1 / gref);
    double u2 = a + gma * (1 - gref);
    double D1 = (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h) * (u1 + h) * S1;
    double D2 = (u2 + h) * 1 / S1 - (u2 - h) * S1;
    // Calculate parameters
    double p1 = (gma / (D1 * S1)) * (u1 - h);
    double p2 = (-gma * S1 / D1) * (u1 + h);
    double p3 = (1 / (D2 * S1)) * (u2 + h);
    double p4 = (-S1 / D2) * (u2 - h);
    // Define and return output variable
    std::vector<double> params(8, 0.0);
    params[0] = p1;
    params[1] = p2;
    params[2] = p3;
    params[3] = p4;
    params[4] = u1;
    params[5] = S1;
    params[6] = D1;
    params[7] = D2;
    return params;
}
// ** Calculates parameters for direct radiation using two-stream model ** //
std::vector<double> twostreamdirCpp(double pait, double om, double a, double gma, double J, double del, double h, double gref,
    double k, double sig, double u1, double S1, double D1, double D2)
{
    // Calculate base parameters
    double ss = 0.5 * (om + J * del / k) * k;
    double sstr = om * k - ss;
    double S2 = exp(-k * pait);
    double u2 = a + gma * (1 - gref);
    double p5 = -ss * (a + gma - k) - gma * sstr;
    double v1 = ss - (p5 * (a + gma + k)) / sig;
    double v2 = ss - gma - (p5 / sig) * (u1 + k);
    double p6 = (1 / D1) * ((v1 / S1) * (u1 - h) - (a + gma - h) * S2 * v2);
    double p7 = (-1 / D1) * ((v1 * S1) * (u1 + h) - (a + gma + h) * S2 * v2);
    sig = -sig;
    double p8 = sstr * (a + gma + k) - gma * ss;
    double v3 = (sstr + gma * gref - (p8 / sig) * (u2 - k)) * S2;
    double p9 = (-1 / D2) * ((p8 / (sig * S1)) * (u2 + h) + v3);
    double p10 = (1 / D2) * (((p8 * S1) / sig) * (u2 - h) + v3);
    // Define and return output variable
    std::vector<double> params(6, 0.0);
    params[0] = p5;
    params[1] = p6;
    params[2] = p7;
    params[3] = p8;
    params[4] = p9;
    params[5] = p10;
    return params;
}
// ** Calculate absorbed shortwave radiation ** //
radmodel RadswabsCpp(double pai, double x, double lref, double ltra, double clump, double gref,
    double slope, double aspect, double lat, double lon,
    std::vector<int> year, std::vector<int> month, std::vector<int> day, std::vector<double> lt,
    std::vector<double> Rsw, std::vector<double> Rdif)
{
    // ** Define variables
    std::vector<double> radGsw(Rsw.size());
    std::vector<double> radCsw(Rsw.size());
    std::vector<double> albedo(Rsw.size());
    if (pai > 0.0) {
        // Calculate time-invariant variables
        double pait = pai;
        if (clump > 0.0) pait = pai / (1 - clump);
        double om = lref + ltra;
        double a = 1 - om;
        double del = lref - ltra;
        double J = 1.0 / 3.0;
        if (x != 1.0) {
            double mla = 9.65 * pow((3 + x), -1.65);
            if (mla > M_PI / 2) mla = M_PI / 2;
            J = cos(mla) * cos(mla);
        }
        double gma = 0.5 * (om + J * del);
        double h = sqrt(a * a + 2 * a * gma);
        // Calculate two-stream parameters (diffuse)
        std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref);
        double p1 = tspdif[0];
        double p2 = tspdif[1];
        double p3 = tspdif[2];
        double p4 = tspdif[3];
        double u1 = tspdif[4];
        double S1 = tspdif[5];
        double D1 = tspdif[6];
        double D2 = tspdif[7];
        double trd = clump * clump;
        // Calculate albedo and ground flux
        double amx = gref;
        if (amx < lref) amx = lref;
        double albd = gref * trd + (1.0 - trd) * (p1 + p2);
        if (albd > amx) albd = amx;
        if (albd < 0.01) albd = 0.01;
        double groundRdd = trd + (1.0 - trd) * (p3 * exp(-h * pait) + p4 * exp(h * pait));
        // Calculate time-variant two-stream parameters (direct)
        for (size_t i = 0; i < Rsw.size(); ++i) {
            if (Rsw[i] > 0.0) {
                // Calculate solar variables
                std::vector<double> solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double zen = solp[0];
                double azi = solp[1];
                double si = solarindexCpp(slope, aspect, zen, azi);
                if (zen > 90.0) zen = 90.0;
                // Calculate canopy extinction coefficient
                double cosz = cos(zen * M_PI / 180);
                std::vector<double> kp = cankCpp(zen, x, si);
                double kd = kp[1];
                double Kc = kp[2];
                // Calculate two-stream parameters (direct)      
                double sig = (kd * kd + gma * gma - pow((a + gma), 2.0));
                std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref, kd, sig, u1, S1, D1, D2);
                double p5 = tspdir[0];
                double p6 = tspdir[1];
                double p7 = tspdir[2];
                double p8 = tspdir[3];
                double p9 = tspdir[4];
                double p10 = tspdir[5];
                // Calculate beam normalisations
                double Rbeam = (Rsw[i] - Rdif[i]) / cosz;
                if (Rbeam > 1352.0) Rbeam = 1352.0;
                double trb = pow(clump, Kc);
                if (trb > 0.999) trb = 0.999;
                if (trb < 0.0) trb = 0.0;
                double Rb = Rbeam * cosz;
                double trg = trb + exp(-kd * pait); // transmission to ground though gaps and leaves
                double Rbc = (trg * si + (1 - trg) * cosz) * Rbeam;
                // Calculate albedo and ground flux
                double albb = trd * gref + (1.0 - trd) * (p5 / sig + p6 + p7);
                if (albb > amx) albb = amx;
                if (albb < 0.01) albb = 0.01;
                double groundRbdd = trb + (1.0 - trb) * ((p8 / -sig) * exp(-kd * pait) +
                    p9 * exp(-h * pait) + p10 * exp(h * pait));
                if (groundRbdd > amx) groundRbdd = amx;
                if (groundRbdd < 0.0) groundRbdd = 0.0;
                // Calculate canopy absorbed
                radCsw[i] = (1.0 - albd) * Rdif[i] + (1.0 - albb) * Rbc;
                // Calculate ground absorbed
                double Rgdif = groundRdd * Rdif[i] + groundRbdd * Rb;
                radGsw[i] = (1.0 - gref) * (Rgdif + exp(-kd * pait) * Rbeam * si);
                // Calculate albedo
                albedo[i] = 1.0 - (radCsw[i] / (Rdif[i] + Rb));
                if (albedo[i] > amx) albedo[i] = amx;
                if (albedo[i] < 0.01) albedo[i] = 0.01;
            }
            else {
                radGsw[i] = 0;
                radCsw[i] = 0;
                albedo[i] = lref;
            }
        }
    }
    else {
        for (size_t i = 0; i < Rsw.size(); ++i) {
            albedo[i] = gref;
            if (Rsw[i] > 0) {
                std::vector<double> solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double zen = solp[0];
                double azi = solp[1];
                double si = solarindexCpp(slope, aspect, zen, azi);
                if (zen > 90.0) zen = 90.0;
                double dirr = (Rsw[i] - Rdif[i]) / cos(zen * M_PI / 180.0);
                radGsw[i] = (1 - gref) * (Rdif[i] + si * dirr);
                radCsw[i] = radGsw[i];
            }
            else {
                radGsw[i] = 0;
                radCsw[i] = 0;
            }
        }
    }
    radmodel out;
    out.ground = radGsw;
    out.canopy = radCsw;
    out.albedo = albedo;
    return out;
}
// ** Calculate molar density of air ** //
double phairCpp(double tc, double pk)
{
    double tk = tc + 273.15;
    double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
    return ph;
}
// ** Calculate specific heat of air at constant pressure ** //
double cpairCpp(double tc)
{
    double cp = 2e-05 * pow(tc, 2) + 0.0002 * tc + 29.119;
    return cp;
}
// ** Calculate zero-plane displacement ** //
// [[Rcpp::export]]
double zeroplanedisCpp(double h, double pai)
{
    if (pai < 0.001) pai = 0.001;
    double d = (1 - (1 - exp(-sqrt(7.5 * pai))) / sqrt(7.5 * pai)) * h;
    return d;
}
// ** Calculate roughness length ** //
// [[Rcpp::export]]
double roughlengthCpp(double h, double pai, double d, double psi_h)
{
    double Be = sqrt(0.003 + (0.2 * pai) / 2);
    double zm = (h - d) * exp(-0.4 / Be) * exp(psi_h);
    if (zm < 0.0005) zm = 0.0005;
    return zm;
}
// **  Calculate integrated diabatic correction coefficient for momentum ** //
// [[Rcpp::export]]
double dpsimCpp(double ze)
{
    double psim;
    // unstable
    if (ze < 0) {
        double x = pow((1 - 15 * ze), 0.25);
        psim = log(pow((1 + x) / 2, 2) * (1 + pow(x, 2)) / 2) - 2 * atan(x) + M_PI / 2;
    }
    // stable
    else {
        psim = -4.7 * ze;
    }
    if (psim < -4) psim = -4;
    if (psim > 3) psim = 3;
    return psim;
}
// **  Calculate integrated diabatic correction coefficient for heat ** //
// [[Rcpp::export]]
double dpsihCpp(double ze)
{
    double psih;
    // unstable
    if (ze < 0) {
        double y = sqrt(1 - 9 * ze);
        psih = log(pow((1 + y) / 2, 2));
    }
    // stable
    else {
        psih = -(4.7 * ze) / 0.74;
    }
    if (psih < -4) psih = -4;
    if (psih > 3) psih = 3;
    return psih;
}
// **  Calculate diabatic influencing factor for heat ** //  
// [[Rcpp::export]]
double dphihCpp(double ze)
{
    double phih;
    // unstable
    if (ze < 0) {
        double phim = 1 / pow((1.0 - 16.0 * ze), 0.25);
        phih = pow(phim, 2.0);
    }
    // stable
    else {
        phih = 1 + ((6.0 * ze) / (1.0 + ze));
    }
    if (phih > 1.5) phih = 1.5;
    if (phih < 0.5) phih = 0.5;
    return phih;
}
// **  Calculate free convection ** //  
double gfreeCpp(double leafd, double H)
{
    double d = 0.71 * leafd;
    double dT = 0.7045388 * pow((d * pow(H, 4)), 0.2);
    double gha = 0.0375 * pow(dT / d, 0.25);
    if (gha < 0.1) gha = 0.1;
    return gha;
}
// **  Calculate molar conductance above canopy ** //  
double gturbCpp(double uf, double d, double zm, double zref, double ph, double psi_h, double gmin)
{
    double z0 = 0.2 * zm + d; // heat exchange surface height
    double ln = log((zref - d) / (z0 - d));
    double g = (0.4 * ph * uf) / (ln + psi_h);
    if (g < gmin) g = gmin;
    return g;
}
// **  Stomatal conductance ** //  
double stomcondCpp(double Rsw, double gsmax, double q50)
{
    double rpar = Rsw * 4.6;
    double gs = (gsmax * rpar) / (rpar + q50);
    return gs;
}
// **  Saturated vapour pressure ** //
// [[Rcpp::export]]  
double satvapCpp(double tc)
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
// **  Dewpoint temperature ** // 
// [[Rcpp::export]]   
double dewpointCpp(double tc, double ea)
{
    double e0;
    double L;
    double it;
    // Dew point
    if (tc >= 0) {
        e0 = 611.2 / 1000;
        L = (2.501 * pow(10, 6)) - (2340 * tc);
        it = 1 / 273.15 - (461.5 / L) * log(ea / e0);
    }
    // Frost point
    else {
        e0 = 610.78 / 1000;
        L = 2.834 * pow(10, 6);
        it = 1 / 273.16 - (461.5 / L) * log(ea / e0);
    }
    double Tdew = 1 / it - 273.15;
    return Tdew;
}
// **  Penman-Monteith equation ** //
// [[Rcpp::export]]  
double PenmanMonteithCpp(double Rabs, double gHa, double gV, double tc, double te, double pk, double ea, double em, double G, double erh)
{
    double sb = 5.67 * pow(10, -8);
    double Rema = em * sb * pow(tc + 273.15, 4);
    double la = 0;
    if (te >= 0) {
        la = 45068.7 - 42.8428 * te;
    }
    else {
        la = 51078.69 - 4.338 * te - 0.06367 * te * te;
    }
    double cp = cpairCpp(te);
    double Da = satvapCpp(tc) - ea;
    double gR = (4 * em * sb * pow(te + 273.15, 3)) / cp;
    double De = satvapCpp(te + 0.5) - satvapCpp(te - 0.5);
    double Ts = tc + ((Rabs - Rema - la * (gV / pk) * Da * erh - G) / (cp * (gHa + gR) + la * (gV / pk) * De * erh));
    return Ts;
}
// **  Function to compute daily from hourly** // 
// [[Rcpp::export]]   
std::vector<double> hourtodayCpp(std::vector<double> hourly, std::string stat) {
    int numDays = hourly.size() / 24;
    std::vector<double> daily(numDays * 24, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double dailyStat = hourly[i * 24];
        if (stat == "max") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::max(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "min") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::min(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "mean") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
            dailyStat /= 24;
        }
        else {
            throw std::invalid_argument("Invalid statistic. Please use 'max', 'min', or 'mean'.");
        }
        // Fill daily with replicated values for each day
        for (int j = 0; j < 24; ++j) {
            daily[i * 24 + j] = dailyStat;
        }
    }
    return daily;
}
// **  Function to compute daily from hourly without replicating each value 24 times** //  
std::vector<double> hourtodayCpp2(std::vector<double> hourly, std::string stat) {
    int numDays = hourly.size() / 24;
    std::vector<double> daily(numDays, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double dailyStat = hourly[i * 24];
        if (stat == "max") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::max(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "min") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::min(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "mean") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
            dailyStat /= 24;
        }
        else if (stat == "sum") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
        }
        else {
            throw std::invalid_argument("Invalid statistic. Please use 'max', 'min', 'mean' or 'sum'.");
        }
        // Fill daily with replicated values for each day
        daily[i] = dailyStat;
    }
    return daily;
}
// **  Function to compute rolling mean temp ** // 
// [[Rcpp::export]] 
std::vector<double> maCpp(std::vector<double> x, int n) {
    std::vector<double> y(x.size());
    int m = x.size();
    for (int i = 0; i < m; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += x[(i - j + m) % m];
        }
        y[i] = sum / n;
    }
    return y;
}
// **  Function to compute rolling mean yearly ** // 
// [[Rcpp::export]]  
std::vector<double> mayCpp(std::vector<double> x) {
    // Calculate daily mean
    int numDays = x.size() / 24;
    std::vector<double> d(numDays, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 24; ++j) {
            sum += x[i * 24 + j];
        }
        d[i] = sum / 24.0;
    }
    // Calculate rolling mean
    std::vector<double> y = maCpp(d, 91);
    std::vector<double> z;
    for (double val : y) {
        for (int i = 0; i < 24; ++i) {
            z.push_back(val);
        }
    }
    return z;
}
// **  Function to compute ground heat flux ** //  
Gmodel GFluxCpp(std::vector<double> Tg, std::vector<double> soilm, double rho, double Vm, double Vq, double Mc,
    std::vector<double> Gmax, std::vector<double> Gmin, int iter, bool yearG = true) {
    // Time invariant variables
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1 - frs);
    double c3 = 1 + 2.6 * pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * frs * frs;
    double mu1 = 2400 * rho / 2.64;
    double mu2 = 1.06 * rho;
    // Calculate daily mean soil surface temperature
    std::vector<double> Td = hourtodayCpp(Tg, "mean");
    // Initalise variables that need retaining
    std::vector<double> Gmu(Tg.size());
    std::vector<double> dT(Tg.size());
    std::vector<double> k(Tg.size());
    std::vector<double> ka(Tg.size());
    for (size_t i = 0; i < Tg.size(); ++i) {
        // Find soil diffusivity
        double cs = mu1 + 4180 * soilm[i];
        double ph = (rho * (1 - soilm[i]) + soilm[i]) * 1000;
        double c2 = mu2 * soilm[i];
        k[i] = c1 + c2 * soilm[i] - (c1 - c4) * exp(-pow(c3 * soilm[i], 4));
        ka[i] = k[i] / (cs * ph);
        double omdy = (2 * M_PI) / (24 * 3600);
        double DD = sqrt(2 * ka[i] / omdy);
        Gmu[i] = sqrt(2) * (k[i] / DD) * 0.5;
        // Calculate T fluctuation from daily mean
        dT[i] = Tg[i] - Td[i];
    }
    // Calculate 6 hour back rolling mean of Gmu and dT to get 3 hour lag
    std::vector<double> Gmud = maCpp(Gmu, 6);
    std::vector<double> G = maCpp(dT, 6);
    // remove this line **
    for (size_t i = 0; i < G.size(); ++i) G[i] = G[i] * Gmud[i] * 1.1171;
    if (iter == 0) {
        Gmin = hourtodayCpp(G, "min");
        Gmax = hourtodayCpp(G, "max");
    }
    for (size_t i = 0; i < G.size(); ++i) {
        if (G[i] < Gmin[i]) G[i] = Gmin[i];
        if (G[i] > Gmax[i]) G[i] = Gmax[i];
    }
    // Calculate yearly ground flux
    if (yearG) {
        // Calculate moving average of k
        std::vector<double> kma = mayCpp(k);
        std::vector<double> kama = mayCpp(ka);
        std::vector<double> Gmuy(Tg.size());
        for (size_t i = 0; i < G.size(); ++i) {
            double omyr = (2 * M_PI) / (G.size() * 3600);
            Gmuy[i] = sqrt(2) * kma[i] / sqrt(2 * kama[i] / omyr);
        }
        // Calculate mean Td
        double sumTd = 0.0;
        for (double num : Td) sumTd += num;
        std::vector<double> dTy(Tg.size());
        for (size_t i = 0; i < dTy.size(); ++i) dTy[i] = Td[i] - sumTd / Td.size();
        std::vector<double> madTy = mayCpp(dTy);
        for (size_t i = 0; i < G.size(); ++i) {
            G[i] = G[i] + madTy[i] * Gmuy[i] * 1.1171;
        }
    }
    Gmodel out;
    out.G = G;
    out.Gmin = Gmin;
    out.Gmax = Gmax;
    return out;
}
// Run point model in progress
// [[Rcpp::export]]
Rcpp::List BigLeafCpp(DataFrame obstime, DataFrame climdata, std::vector<double> vegp, std::vector<double> groundp, std::vector<double> soilm,
    double lat, double lon, double dTmx = 25.0, double zref = 2.0, int maxiter = 100, double bwgt = 0.5, double tol = 0.5, double gmn = 0.1, bool yearG = true)
{
    // Access items of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> prec = climdata["precip"];
    // Access items of vegp
    double h = vegp[0];
    double pai = vegp[1];
    double vegx = vegp[2];
    double clump = vegp[3];
    double lref = vegp[4];
    double ltra = vegp[5];
    double leafd = vegp[6];
    double em = vegp[7];
    double gsmax = vegp[8];
    double q50 = vegp[9];
    // Access items of groundp
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    double groundem = groundp[3];
    double rho = groundp[4];
    double Vm = groundp[5];
    double Vq = groundp[6];
    double Mc = groundp[7];
    double Smin = groundp[11];
    double Smax = groundp[10];
    // Calculate SW radiation
    radmodel swabs = RadswabsCpp(pai, vegx, lref, ltra, clump, gref, slope, aspect, lat, lon, year, month, day, hour, Rsw, Rdif);
    // Calculate time-invarient variables
    double pait = pai / (1 - clump);
    double trd = (1 - clump * clump) * exp(-pait) + clump * clump;
    double sb = 5.67 * pow(10, -8);
    double d = zeroplanedisCpp(h, pai);
    // Used to avoid (h-d)/zm being less than one, meaning log((h-d)/zm) becomes negative
    double Belim = 0.4 / sqrt(0.003 + (0.2 * pai) / 2);
    // Initialise variables set to zero on first run
    std::vector<double> Tg = tc;
    std::vector<double> Tc = tc;
    std::vector<double> tcc = tc;
    std::vector<double> tcg = tc;
    std::vector<double> psim(Rsw.size(), 0.0);
    std::vector<double> psih(Rsw.size(), 0.0);
    std::vector<double> phih(Rsw.size(), 0.0);
    std::vector<double> LL(Rsw.size(), 0.0);
    std::vector<double> G(Rsw.size(), 0.0);
    std::vector<double> Gmin(Rsw.size(), -999.0);
    std::vector<double> Gmax(Rsw.size(), 999.0);
    std::vector<double> uf(Rsw.size(), 999.0);
    std::vector<double> RabsG(Rsw.size(), 999.0);
    // New variables for storing
    // Initalise H
    std::vector<double> H(Rsw.size());
    for (size_t i = 0; i < H.size(); ++i) H[i] = 0.5 * Rsw[i] - em * sb * pow(tc[i] + 273.15, 4);
    double tstf = tol * 2;
    int iter = 0;
    double tst = 0;
    while (tstf > tol) {
        tst = 0;
        for (size_t i = 0; i < Rsw.size(); ++i) {
            // Calculate longwave radiation
            double RemC = em * sb * pow(Tc[i] + 273.15, 4);
            double radClw = em * Rlw[i]; // Longwave radiation down from sky
            double radGlw = groundem * (trd * radClw + (1 - trd) * RemC);
            // Calculate absorbed radiation
            RabsG[i] = swabs.ground[i] + radGlw;
            double RabsC = swabs.canopy[i] + radClw;
            // Calculate canopy temperature
            double zm = roughlengthCpp(h, pai, d, psih[i]);
            uf[i] = (0.4 * wspeed[i]) / (log((zref - d) / zm) + psim[i]);
            if (uf[i] < 0.0002) uf[i] = 0.0002;
            double gmin = gfreeCpp(leafd, abs(H[i])) * 2 * pai;
            double ph = phairCpp(tcc[i], pk[i]);
            double gHa = gturbCpp(uf[i], d, zm, zref, ph, psih[i], gmin);
            double gC = stomcondCpp(Rsw[i], gsmax * 3, q50 * 3);
            double gV = 1 / (1 / gHa + 1 / gC);
            if (gC == 0) gV = 0;
            double ea = satvapCpp(tc[i]) * rh[i] / 100;
            double Tcn = PenmanMonteithCpp(RabsC, gHa, gV, tc[i], tcc[i], pk[i], ea, em, G[i], 1);
            double tdew = dewpointCpp(tc[i], ea);
            if (Tcn < tdew) Tcn = tdew;
            // Calculate ground surface temperature
            double srh = (soilm[i] - Smin) / (Smax - Smin);
            double Tgn = PenmanMonteithCpp(RabsG[i], gHa, gHa, tcg[i], tcc[i], pk[i], ea, em, G[i], srh);
            if (Tgn < tdew) Tgn = tdew;
            // Cap values
            double dTc = Tcn - tc[i];
            double dTg = Tgn - tc[i];
            if (dTc > dTmx) dTc = dTmx;
            if (dTg > dTmx) dTg = dTmx;
            Tcn = tc[i] + dTc;
            Tgn = tc[i] + dTg;
            // Run tests for convergence
            double tst2 = abs(Tcn - Tc[i]);
            double tst3 = abs(Tgn - Tg[i]);
            if (tst2 > tst) tst = tst2;
            if (tst3 > tst) tst = tst3;
            // Reassign Tc and Tg using bwgt
            Tc[i] = bwgt * Tc[i] + (1 - bwgt) * Tcn;
            Tg[i] = bwgt * Tg[i] + (1 - bwgt) * Tgn;
            // Recalculate variables
            tcc[i] = (Tc[i] + tc[i]) / 2;
            tcg[i] = (Tg[i] + tc[i]) / 2;
            double Tk = 273.15 + tcc[i];
            ph = phairCpp(tcc[i], pk[i]);
            double cp = cpairCpp(tcc[i]);
            // Calculate H
            H[i] = bwgt * H[i] + (1 - bwgt) * (gHa * cp * (Tcn - tc[i]));
            // Set limits to H
            double Rnet = RabsC - sb * em * pow(Tc[i] + 273.15, 4);
            if (Rnet > 0 && H[i] > Rnet) H[i] = Rnet;
            // Recalculate stablity variables
            // Stability
            LL[i] = (ph * cp * pow(uf[i], 3) * Tk) / (-0.4 * 9.81 * H[i]);
            //if (LL[i] > 10000.0) LL[i] = 10000.0;
            //if (LL[i] < -10000.0) LL[i] = -10000.0;
            psim[i] = dpsimCpp(zm / LL[i]) - dpsimCpp((zref - d) / LL[i]);
            psih[i] = dpsihCpp((0.2 * zm) / LL[i]) - dpsihCpp((zref - d) / LL[i]);
            phih[i] = dphihCpp((zref - d) / LL[i]);
            // Set limits to diabatic coefficients
            double ln1 = log((zref - d) / zm);
            double ln2 = log((zref - d) / (0.2 * zm));
            if (psim[i] < -0.9 * ln1) psim[i] = -0.9 * ln1;
            if (psih[i] < -0.9 * ln2) psih[i] = -0.9 * ln2;
            if (psim[i] > 0.9 * ln1) psim[i] = 0.9 * ln1;
            if (psih[i] > 0.9 * ln2) psih[i] = 0.9 * ln2;
            if (psih[i] > 0.9 * Belim) psih[i] = 0.9 * Belim;
        }
        // Recalculate Ground heat flux
        Gmodel GG = GFluxCpp(Tg, soilm, rho, Vm, Vq, Mc, Gmax, Gmin, iter, yearG);
        G = GG.G;
        Gmin = GG.Gmin;
        Gmax = GG.Gmax;
        tstf = tst;
        ++iter;
        if (iter >= maxiter) tstf = 0;
    }
    // Return outputs
    Rcpp::List out;
    out["Tc"] = Rcpp::wrap(Tc);
    out["Tg"] = Rcpp::wrap(Tg);
    out["H"] = Rcpp::wrap(H);
    out["G"] = Rcpp::wrap(G);
    out["psih"] = Rcpp::wrap(psih);
    out["psim"] = Rcpp::wrap(psim);
    out["phih"] = Rcpp::wrap(phih);
    out["OL"] = Rcpp::wrap(LL);
    out["uf"] = Rcpp::wrap(uf);
    out["RabsG"] = Rcpp::wrap(RabsG);
    out["err"] = Rcpp::wrap(tst);
    out["albedo"] = Rcpp::wrap(swabs.albedo);
    return out;
}
// Perform weather height adjustment
// [[Rcpp::export]]
DataFrame weatherhgtCpp(DataFrame obstime, DataFrame climdata, double zin, double uzin, double zout, double lat, double lon, bool yearG = true)
{
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> ws = climdata["windspeed"];
    std::vector<double> vegp({ 0.12,1,1,0.1,0.4,0.2,0.05,0.97,0.33,100.0 });
    std::vector<double> groundp({ 0.15,0.0, 180.0,0.97,1.529643,0.509,0.06,0.5422,5.2,-5.6,0.419,0.074 });
    std::vector<double> soilm(tc.size(), 0.2); // Initialize soilm with size and value 0.2
    int hrs = tc.size();
    if (hrs < 8760) yearG = false;
    // Run point model
    Rcpp::List bigleafp = BigLeafCpp(obstime, climdata, vegp, groundp, soilm, lat, lon, 25, 2, 20, 0.5, 0.5, 0.1, yearG);
    // Extract things needed from list
    std::vector<double> psih = Rcpp::as<std::vector<double>>(bigleafp["psih"]);
    std::vector<double> Tc = Rcpp::as<std::vector<double>>(bigleafp["Tc"]);
    // Define variables
    std::vector<double> Tz(tc.size());
    std::vector<double> Rh(tc.size());
    std::vector<double> Uz(tc.size());
    double d = zeroplanedisCpp(0.12, 1);
    for (size_t i = 0; i < tc.size(); ++i) {
        double zm = roughlengthCpp(0.12, 1, d, psih[i]);
        double zh = 0.2 * zm;
        double lnr = log((zout - d) / zh) / log((zin - d) / zh);
        // Temperature
        Tz[i] = (Tc[i] - tc[i]) * (1 - lnr) + tc[i];
        // Humidity
        double ea = satvapCpp(tc[i]) * rh[i] / 100;
        double es = satvapCpp(Tc[i]) * sqrt(rh[i] / 100);
        double ez = ea + (es - ea) * (1 - lnr);
        es = satvapCpp(Tz[i]);
        Rh[i] = (ez / es) * 100;
        // Cap Rh
        if (Rh[i] < 0.25 * rh[i]) Rh[i] = 0.25 * rh[i];
        if (Rh[i] > 100.0) Rh[i] = 100.0;
        // Wind speed
        double lnru = log((zout - d) / zm) / log((uzin - d) / zm);
        Uz[i] = ws[i] * lnru;
    }
    // Clone climdata to stop original input being over-written in R
    DataFrame climdata_copy = clone(climdata); // Make a copy of climdata
    climdata_copy["temp"] = Tz;
    climdata_copy["relhum"] = Rh;
    climdata_copy["windspeed"] = Uz;
    return climdata_copy;
}
// Calculate mean dtr
// [[Rcpp::export]]
double meandtrCpp(std::vector<double> temp) {
    int ndays = temp.size() / 24;
    // Calculate daily diurnal temperature range
    std::vector<double> dtr(ndays);
    int idx = 0;
    for (int d = 0; d < ndays; ++d) {
        double tmx = -999.9;
        double tmn = 999.9;
        for (int h = 0; h < 24; ++h) {
            if (temp[idx] > tmx) tmx = temp[idx];
            if (temp[idx] < tmn) tmn = temp[idx];
            ++idx;
        }
        dtr[d] = tmx - tmn;
    }
    // Calculate mean daily diurnal temperature range
    double dtrm = 0.0;
    for (int d = 0; d < ndays; ++d) dtrm += dtr[d];
    dtrm = dtrm / ndays;
    return dtrm;
}
// Adjust dtr by a fixed amount
// [[Rcpp::export]]
std::vector<double> adjustdtrCpp(std::vector<double> temp, double dtrc) {
    int ndays = temp.size() / 24;
    // Create new vector
    std::vector<double> tempadj(temp.size());
    int idx1 = 0;
    int idx2 = 0;
    for (int d = 0; d < ndays; ++d) {
        // Calculate daily max and min
        double tmx = -999.9;
        double tmn = 999.9;
        for (int h = 0; h < 24; ++h) {
            if (temp[idx1] > tmx) tmx = temp[idx1];
            if (temp[idx1] < tmn) tmn = temp[idx1];
            ++idx1;
        }
        // Perform adjustment
        double tmean = (tmx + tmn) / 2.0;
        for (int h = 0; h < 24; ++h) {
            double dif = temp[idx2] - tmean;
            tempadj[idx2] = (dif * dtrc) + tmean;
            ++idx2;
        }
    }
    return tempadj;
}
// Adjust relative humidity by new temperature
std::vector<double> adjustrelhum(std::vector<double> tc,
    std::vector<double> tcn, std::vector<double> rh)
{
    // Initialise new relative humidity variable
    std::vector<double> rhn(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) {
        // Calculate vapour pressure
        double ea = satvapCpp(tc[i]) * (rh[i] / 100.0);
        // Calculate new relative humidity
        double esn = satvapCpp(tcn[i]);
        rhn[i] = (ea / esn) * 100.0;
        // Set limits
        if (rhn[i] > 100.0) rhn[i] = 100.0;
        if (rhn[i] < 25.0) rhn[i] = 25.0;
    }
    return rhn;
}
// Correct dtr when applying weather height adjustment
// [[Rcpp::export]]
DataFrame dtr_correctCpp(DataFrame obstime, DataFrame climdata,
    double zin, double uzin, double zout, double lat, double lon,
    bool yearG = true)
{
    int iter = 0;
    int tst = 1;
    double dtrc = 1.0;
    DataFrame climdatan = clone(climdata);
    std::vector<double> tc = climdatan["temp"];
    std::vector<double> tcn = tc;
    std::vector<double> rh = climdatan["relhum"];
    while (tst > 0) {
        DataFrame climdata2 = weatherhgtCpp(obstime, climdatan, zin, uzin, zout,
            lat, lon, yearG);
        std::vector<double> tc2 = climdata2["temp"];
        double dtrr = meandtrCpp(tc2) / meandtrCpp(tcn);
        if (dtrr > 1.0) {
            dtrc += 0.1;
            tcn = adjustdtrCpp(tc, dtrc);
            climdatan["temp"] = tcn;
        }
        else {
            tst = 0;
        }
        ++iter;
        if (iter > 20) tst = 0;
    }
    rh = adjustrelhum(tc, tcn, rh);
    climdatan["relhum"] = rh;
    return climdatan;
}
// Soil moisture model
// [[Rcpp::export]]
std::vector<double> soilmCpp(DataFrame climdata, double rmu, double mult, double pwr, double Smax, double Smin, double Ksat, double a)
{
    // Extract data from columns
    std::vector<double> temp = climdata["temp"];
    std::vector<double> swdown = climdata["swdown"];
    std::vector<double> lwdown = climdata["lwdown"];
    std::vector<double> rainh = climdata["precip"];
    // Calculate net radiation
    std::vector<double> rnet(temp.size());
    for (size_t i = 0; i < temp.size(); ++i) {
        double swrad = (1 - 0.15) * swdown[i];
        double lwout = 5.67 * pow(10, -8) * 0.95 * pow(temp[i] + 273.15, 4);
        double lwnet = lwout - lwdown[i];
        rnet[i] = swrad - lwnet;
        if (rnet[i] < 0) rnet[i] = 0;
    }
    // Convert to daily
    std::vector<double> rnetd = hourtodayCpp2(rnet, "mean");
    std::vector<double> rain = hourtodayCpp2(rainh, "sum");
    // Run soil model
    std::vector<double> soilm(rain.size());
    double s1 = Smax;
    double s2 = Smax;
    soilm[0] = Smax;
    for (size_t i = 1; i < rain.size(); ++i) {
        double sav = (s1 + s2) / 2;
        double dif = s2 - s1;
        s1 = s1 + rmu * rain[i] - mult * rnetd[i];
        double k = Ksat * pow(sav / Smax, pwr);
        s1 = s1 + a * k * dif;
        s2 = s2 - ((a * k * dif) / 10);
        if (s1 > Smax) s1 = Smax;
        if (s2 > Smax) s2 = Smax;
        if (s1 < Smin) s1 = Smin;
        if (s2 < Smin) s2 = Smin;
        soilm[i] = (s1 + s2) / 2;
    }
    return soilm;
}
// ========================================================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Small leaf model from here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
// ========================================================================================================================= #
// Calculate cumulative paii
std::vector<double> reverseCumsum(std::vector<double> paii) {
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
// Calculate shortwave radiation below canopy
radmodel2 RadiationSmallLeafSWCpp(double lat, double lon, int year, int month, int day, double lt,
    double reqhgt, double hgt, std::vector<double> paii, double x, double lref, double lrefp, double ltra,
    double ltrap, double clump, double vegem, double gref, double groundem, double slope, double aspect,
    double Rsw, double Rdif)
{
    // Create output variables
    std::vector<double> Rswabs(paii.size());
    std::vector<double> PAR(paii.size());
    std::vector<double> Rdirdown(paii.size());
    std::vector<double> Rdifdown(paii.size());
    std::vector<double> Rswup(paii.size());
    if (Rsw > 0) {
        // Calculate heights and paia, setting one z to reqhgt
        double mindif = hgt;
        int whichz = 0;
        std::vector<double> z(paii.size());
        std::vector<double> paic = reverseCumsum(paii);
        double nn = paii.size();
        double pai = paii[0];
        z[0] = (1 / nn) * hgt;
        for (size_t i = 1; i < paii.size(); ++i) {
            z[i] = ((i + 1) / nn) * hgt;
            pai += paii[i];
            double newdif = abs(z[i] - reqhgt);
            if (newdif < mindif) {
                whichz = i;
                mindif = newdif;
            }
        }
        z[whichz] = reqhgt;
        // Calculate and adjust pait
        double pait = pai / (1 - clump);
        // Calculate additional variables  
        double om = lref + ltra;
        double del = lref - ltra;
        double a = 1 - om;
        double J = 1.0 / 3.0;
        if (x != 1.0) {
            double mla = 9.65 * pow((3 + x), -1.65);
            if (mla > M_PI / 2) mla = M_PI / 2;
            J = cos(mla) * cos(mla);
        }
        double gma = 0.5 * (om + J * del);
        double h = sqrt(a * a + 2 * a * gma);
        // Calculate additional variables (PAR) 
        double omp = lrefp + ltrap;
        double delp = lrefp - ltrap;
        double ap = 1 - om;
        double gmap = 0.5 * (omp + J * delp);
        double hp = sqrt(ap * ap + 2 * ap * gmap);
        // Calculate solar variables
        std::vector<double> solp = solpositionCpp(lat, lon, year, month, day, lt);
        double zen = solp[0];
        double azi = solp[1];
        double si = solarindexCpp(slope, aspect, zen, azi);
        if (zen > 90.0) zen = 90.0;
        // Calculate two-stream parameters (diffuse)
        std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref);
        std::vector<double> tspdifPAR = twostreamdifCpp(pait, omp, ap, gmap, hp, gref);
        // Extract two-stream diffuse parameters
        double p1 = tspdif[0];
        double p2 = tspdif[1];
        double p3 = tspdif[2];
        double p4 = tspdif[3];
        double u1 = tspdif[4];
        double S1 = tspdif[5];
        double D1 = tspdif[6];
        double D2 = tspdif[7];
        // Extract two-stream diffuse PAR parameters
        double p1R = tspdifPAR[0];
        double p2R = tspdifPAR[1];
        double p3R = tspdifPAR[2];
        double p4R = tspdifPAR[3];
        double u1R = tspdifPAR[4];
        double S1R = tspdifPAR[5];
        double D1R = tspdifPAR[6];
        double D2R = tspdifPAR[7];
        // Calculate canopy extinction coefficient
        double cosz = cos(zen * M_PI / 180);
        std::vector<double> kp = cankCpp(zen, x, si);
        double k = kp[0];
        double kd = kp[1];
        double Kc = kp[2];
        // Calculate two-stream parameters (direct)      
        double sig = kd * kd + gma * gma - pow((a + gma), 2);
        double sigp = kd * kd + gmap * gmap - pow(ap + gmap, 2);
        std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref, kd, sig, u1R, S1R, D1R, D2R);
        std::vector<double> tspdirPAR = twostreamdirCpp(pait, omp, ap, gmap, J, delp, hp, gref, kd, sigp, u1, S1, D1, D2);
        // Extract two-stream direct parameters
        double p5 = tspdir[0];
        double p6 = tspdir[1];
        double p7 = tspdir[2];
        double p8 = tspdir[3];
        double p9 = tspdir[4];
        double p10 = tspdir[5];
        // Extract two-stream parameters direct (PAR)
        double p5R = tspdirPAR[0];
        double p6R = tspdirPAR[1];
        double p7R = tspdirPAR[2];
        double p8R = tspdirPAR[3];
        double p9R = tspdirPAR[4];
        double p10R = tspdirPAR[5];
        // Calculate transmissions ot bottom of the canopy
        double trdn = clump * clump;
        // ~~ Direct beam above canopy
        double Rbeam = (Rsw - Rdif) / cosz;
        if (Rbeam > 1352.0) Rbeam = 1352.0;
        double Rb = Rbeam * cosz;
        // Iterate through to calaculate radiation for each canopy element
        double amx = gref;
        double amxp = gref;
        if (amx < lref) amx = lref;
        if (amxp < lrefp) amxp = lref;
        for (size_t i = 0; i < paii.size(); ++i) {
            double Pi = paic[i]; // leaf area above
            // Calculate gap fractions and transmissions looking up and down
            double gi = pow(clump, Pi / pai);
            double giu = pow(clump, (pai - Pi) / pai);
            if (gi > 0.99) gi = 0.99;
            if (giu > 0.99) giu = 0.99;
            double trdd = gi * gi;
            double trdu = giu * giu;
            double trbd = pow(gi, Kc);
            if (trbd > 0.999) trbd = 0.999;
            if (trbd < 0.0) trbd = 0.0;
            double Pia = Pi / (1.0 - gi); // leaf area adjusted for gap fraction
            // **************** Calculate Shortwave ************************
            // ~~ Calculate normalised diffuse only upward
            double Rddu = trdu * gref + (1.0 - trdn) * (p1 * exp(-h * Pia) +
                p2 * exp(h * Pia));
            if (Rddu > 1.0) Rddu = 1.0;
            if (Rddu < 0.0) Rddu = 0.0;
            // ~~ Calculate normalised diffuse only downward
            double Rddd = trdd + (1.0 - trdd) * (p3 * exp(-h * Pia) +
                p4 * exp(h * Pia));
            if (Rddd > 1.0) Rddd = 1.0;
            if (Rddd < 0.0) Rddd = 1.0;
            // ~~ Calculate normalised contribution of direct to upward diffuse
            double Rdbu = trdu * gref + (1.0 - trdu) * ((p5 / sig) * exp(-kd * Pia) +
                p6 * exp(-h * Pia) + p7 * exp(h * Pia));
            if (Rdbu > amx) Rdbu = amx;
            if (Rdbu < 0.0) Rdbu = 0.0;
            // ~~ Calculate normalised contribution of direct to downward diffuse
            double Rdbd = (1.0 - trbd) * ((p8 / -sig) * exp(-kd * Pia) +
                p9 * exp(-h * Pia) + p10 * exp(h * Pia));
            if (Rdbd > amx) Rdbd = amx;
            if (Rdbd < 0.0) Rdbd = 0.0;
            // ~~ Calculate actual upward and downward fluxes
            Rdirdown[i] = Rbeam * trbd + ((1.0 - trbd) * exp(-kd * Pia));
            Rdifdown[i] = Rddd * Rdif + Rdbd * Rb;
            Rswup[i] = Rddu * Rdif + Rdbu * Rb;
            // ~~ Calculate shortwave radiation absorbed by leaf
            double sil = k * cosz;
            double swupper = (1 - lref) * (Rdifdown[i] + sil * Rdirdown[i]);
            double swunder = (1 - lref) * Rswup[i];
            Rswabs[i] = 0.5 * (swupper + swunder);
            // **************** Calculate PAR ************************
            // ~~ Calculate normalised diffuse only upward
            Rddu = trdu * gref + (1.0 - trdn) * (p1R * exp(-hp * Pia) +
                p2R * exp(hp * Pia));
            if (Rddu > 1.0) Rddu = 1.0;
            if (Rddu < 0.0) Rddu = 0.0;
            // ~~ Calculate normalised diffuse only downward
            Rddd = trdd + (1.0 - trdd) * (p3R * exp(-hp * Pia) +
                p4R * exp(hp * Pia));
            if (Rddd > 1.0) Rddd = 1.0;
            if (Rddd < 0.0) Rddd = 0.0;
            // ~~ Calculate normalised contribution of direct to upward diffuse
            Rdbu = trdu * gref + (1.0 - trdu) * ((p5R / sigp) * exp(-kd * Pia) +
                p6R * exp(-hp * Pia) + p7R * exp(hp * Pia));
            if (Rdbu > amxp) Rdbu = amxp;
            if (Rdbu < 0.0) Rdbu = 0.0;
            // ~~ Calculate normalised contribution of direct to downward diffuse
            Rdbd = (1.0 - trbd) * ((p8R / -sigp) * exp(-kd * Pia) +
                p9R * exp(-hp * Pia) + p10R * exp(h * Pia));
            if (Rdbd > amxp) Rdbd = amxp;
            if (Rdbd < 0.0) Rdbd = 0.0;
            // ~~ Calculate actual upward and downward fluxes
            double PARdn = Rddd * Rdif + Rdbd * Rb;
            double PARup = Rddu * Rdif + Rdbu * Rb;
            // ~~ Calculate PAR absorbed by leaf
            swupper = (1 - lrefp) * (PARdn + sil * Rdirdown[i]);
            swunder = (1 - lrefp) * PARup;
            PAR[i] = 0.5 * (swupper + swunder);
        }
    }
    else {
        for (size_t i = 0; i < paii.size(); ++i) {
            Rswabs[i] = 0.0;
            PAR[i] = 0.0;
            Rdirdown[i] = 0.0;
            Rdifdown[i] = 0.0;
            Rswup[i] = 0.0;
        }
    }
    radmodel2 out;
    out.Rswabs = Rswabs;
    out.PAR = PAR;
    out.Rdirdown = Rdirdown;
    out.Rdifdown = Rdifdown;
    out.Rswup = Rswup;
    return out;
}
// ** Calculate longwave radiation weights ** //
LWweights lwradweights(std::vector<double> paii) {
    // Calculate pai above and below
    std::vector<double> paib(paii.size());
    paib[0] = paii[0];
    for (size_t i = 1; i < paii.size(); ++i) paib[i] = paib[i - 1] + paii[i];
    double pait = paib[paib.size() - 1];
    std::vector<double> paia(paii.size());
    for (size_t i = 0; i < paii.size(); ++i) paia[i] = pait - paib[i];
    // Initialise wgt matrix
    int n = paii.size();
    // Calculate relative weighting from each foliage element to every other
    NumericMatrix wgt(n, n);
    for (size_t i = 0; i < paii.size(); ++i) {
        for (size_t j = 0; j < paii.size(); ++j) {
            double pai = 0;
            if (i > j) pai = paib[i] - paib[j];
            if (i < j) pai = paia[i] - paia[j];
            double tr = exp(-pai);
            wgt(i, j) = tr * paii[j];
        }
    }
    // Calculate transmission from ground and canopy
    std::vector<double> trg(paii.size());
    std::vector<double> trh(paii.size());
    for (size_t i = 0; i < paii.size(); ++i) {
        trg[i] = exp(-paib[i]);
        trh[i] = exp(-paia[i]);
    }
    LWweights out;
    out.trg = trg;
    out.trh = trh;
    out.wgt = wgt;
    return out;
}
// Calculate longwave radiation below canopy
radmodel3 RadiationSmallLeafLWCpp(std::vector<double> paii, double lwdown, double tground, double groundem, double vegem,
    std::vector<double> tleaf)
{
    // Calculate longwave radiation weights
    LWweights wgts = lwradweights(paii);
    // ** Longwave 
    double sb = 5.67 * pow(10, -8);
    std::vector<double> lwupper(paii.size());
    std::vector<double> lwunder(paii.size());
    for (size_t i = 0; i < paii.size(); ++i)
    {
        // longwave radiation downward from sky
        double lwsky = lwdown * wgts.trh[i];
        // longwave radiation upward from ground
        double lwgro = groundem * sb * pow(tground + 273.15, 4) * wgts.trg[i];
        // Calculate multipliers for foliage radiation 
        // ** Calculate sums of weights
        double smd = 0;
        double smu = 0;
        for (size_t j = i; j < paii.size(); ++j) smd = smd + wgts.wgt(i, j);
        for (size_t j = 0; j <= i; ++j) smu = smu + wgts.wgt(i, j);
        double mua = 1.0;
        double mub = 1.0;
        if (smd > 0) mua = (1 / smd) * (1 - wgts.trh[i]);
        if (smu > 0) mub = (1 / smu) * (1 - wgts.trg[i]);
        // Calculate summed longwave radiation for foliage
        double lwfd = 0;
        double lwfu = 0;
        for (size_t j = i; j < paii.size(); ++j) lwfd = lwfd + wgts.wgt(i, j) * sb * pow(tleaf[j] + 273.15, 4);
        for (size_t j = 0; j <= i; ++j) lwfu = lwfu + wgts.wgt(i, j) * sb * pow(tleaf[j] + 273.15, 4);
        lwupper[i] = lwsky + lwfd * mua * vegem;
        lwunder[i] = lwgro + lwfu * mub * vegem;
    }
    radmodel3 out;
    out.Rlwdown = lwupper;
    out.Rlwup = lwunder;
    return out;
}
// Calculate canopy wind profile
// [[Rcpp::export]]
std::vector<double> CanopyWindCpp(double hgt, std::vector<double> paii) {
    // Calculate whole canopy attenuation coefficient
    double pai = 0;
    for (size_t i = 0; i < paii.size(); ++i) pai += paii[i];
    double Be = 0.205 * pow(pai, 0.445) + 0.1;
    double a = pai / hgt;
    double Lc = pow(0.25 * a, -1);
    double Lm = 2 * pow(Be, 3) * Lc;
    double at = Be * hgt / Lm;
    // Calculate attenuation coefficient for canopy elements
    std::vector<double> ati(paii.size());
    double sati = 0;
    for (size_t i = 0; i < paii.size(); ++i) {
        double Bei = 0.205 * pow(paii[i], 0.445) + 0.1;
        double ai = paii[i] / hgt;
        double Lci = pow(0.25 * ai, -1);
        double Lmi = 2 * pow(Bei, 3) * Lci;
        ati[i] = Bei * hgt / Lmi;
        sati += ati[i];
    }
    // Adjust attenuation coefficient
    for (size_t i = 0; i < paii.size(); ++i) ati[i] = (ati[i] / sati) * at;
    // Calculate canopy wind shelter coefficient
    int n = paii.size();
    int n2 = std::trunc(n / 10);
    std::vector<double> ui(n, 1.0);
    for (int i = n - 1; i >= n2; --i) {
        ui[i - 1] = ui[i] * (1 - ati[i - 1]);
    }
    // Calculate bottom 10 percent
    double zm = hgt / (20.0 * n2);
    for (int i = 0; i < n2; ++i) {
        double z2 = (i + 1) * hgt / (10 * n2);
        double uf = (0.4 * ui[n2 - 1]) / log(hgt / (10 * zm));
        ui[i] = (uf / 0.4) * log(z2 / zm);
    }
    return ui;
}
canHL CanopyHL(double uh, double pk, radmodel2 swrad, radmodel3 lwrad, std::vector<double> wc,
    double leafd, double gsmax, double q50, double vegem, std::vector<double> tleaf, std::vector<double> tair, std::vector<double> ea,
    double surfwet)
{
    // Calculate sensible heat flux
    std::vector<double> H(tair.size());
    std::vector<double> L(tair.size());
    std::vector<double> tleafn(tair.size());
    std::vector<double> uz(tair.size());
    for (size_t i = 0; i < tair.size(); ++i) {
        // Conductance for heat
        uz[i] = wc[i] * uh;
        double gHa = 1.4 * 0.135 * std::sqrt(uz[i] / (0.71 * leafd));
        double gFo = 1.4 * 0.05 * pow(abs(tleaf[i] - tair[i]) / (0.71 * leafd), 0.25);
        if (gFo > gHa) gHa = gFo;
        if (gHa < 0.25) gHa = 0.25;
        // Conductance for vapour
        double gS = stomcondCpp(swrad.PAR[i], gsmax, q50);
        double gV = 0.0;
        if (gS > 0.0) gV = 1 / (1 / gHa + 1 / gS);
        // PenmanMonteith
        double lwabs = 0.5 * vegem * (lwrad.Rlwdown[i] + lwrad.Rlwup[i]);
        double Rabs = swrad.Rswabs[i] + lwabs;
        double te = (tair[i] + tleaf[i]) / 2;
        tleafn[i] = PenmanMonteithCpp(Rabs, gHa, gV, tair[i], te, pk, ea[i], vegem, 0.0, surfwet);
        double tdew = dewpointCpp(tair[i], ea[i]);
        if (tleafn[i] < tdew) tleafn[i] = tdew;
        te = (tair[i] + tleafn[i]) / 2;
        // Sensible heat flux
        double cp = cpairCpp(te);
        H[i] = cp * gHa * (tleafn[i] - tair[i]);
        // Latent heat flux
        double la = 45068.7 - 42.8428 * te;
        if (te < 0) la = 51078.69 - 4.338 * te - 0.06367 * te * te;
        double es = satvapCpp(te);
        L[i] = la * gV / pk * (es - ea[i]) * surfwet;
    }
    canHL out;
    out.H = H;
    out.L = L;
    out.tleaf = tleafn;
    out.uz = uz;
    return out;
}
// Run Langrangian model for one iteration
Lang LangrangianOne(double reqhgt, double uh, double th, double tlh, double eh, double pk, double lwdown, radmodel2 swrad, 
    std::vector<double> wc, std::vector<double> vegp, std::vector<double> paii, double groundem, std::vector<double> tleaf, 
    std::vector<double> tair, std::vector<double> ea, double tground, double surfwet, double theta, double psim, double psih, 
    double phih, std::vector<double> z, double a0 = 0.25, double a1 = 1.25)
{
    // Extract vegp
    double hgt = vegp[0];
    double pai = vegp[1];
    //double x = vegp[2];
    //double clump = vegp[3];
    //double lref = vegp[4];
    //double ltra = vegp[5];
    double leafd = vegp[6];
    double vegem = vegp[7];
    double gsmax = vegp[8];
    double q50 = vegp[9];
    // Calculate longwave rad
    radmodel3 lwrad = RadiationSmallLeafLWCpp(paii, lwdown, tground, groundem, vegem, tleaf);
    // Calculate canopy sensible and latent heat
    canHL HL = CanopyHL(uh, pk, swrad, lwrad, wc, leafd, gsmax, q50, vegem, tleaf, tair, ea, surfwet);
    std::vector<double> tleafn = HL.tleaf;
    // Calculate Langrangian timescale
    double d = zeroplanedisCpp(hgt, pai);
    double zm = roughlengthCpp(hgt, pai, d, psih);
    double uf = (0.4 * uh) / (log((hgt - d) / zm) + psim);
    double a2 = 0.4 * (1 - d / hgt) / (phih * pow(a1, 2.0));
    double TL = a2 * hgt / uf;
    std::vector<double> ow(tair.size());
    std::vector<double> KH(tair.size());
    std::vector<double> H(tair.size());
    std::vector<double> L(tair.size());
    std::vector<double> ST(tair.size());
    std::vector<double> SL(tair.size());
    double sumRH = 0.0;
    double sbtm = 0;
    double CnzrT = 0;
    double CnzrL = 0;
    double nn = paii.size();
    for (size_t i = 0; i < tair.size(); ++i) {
        // Calculate thermal diffusivity
        double x = M_PI * (1 - z[i] / hgt);
        ow[i] = uf * (0.5 * (a1 + a0) + 0.5 * (a1 - a0) * cos(x));
        KH[i] = TL * ow[i] * ow[i];
        double RH = 1 / KH[i];
        sumRH = sumRH + RH;
        double rHa = sumRH * (hgt / nn);
        if (rHa < 2.0) rHa = 2.0;
        // Calculate ground heat flux
        double ph = phairCpp(tair[i], pk);
        double cp = cpairCpp(tair[i]);
        double GT = (ph * cp / rHa) * (tground - tair[i]);
        // Calculate ground evapotranspiration
        double es = satvapCpp(tleafn[i]);
        double te = (tleafn[i] - tair[i]) / 2;
        double la = 45068.7 - 42.8428 * te;
        if (te < 0) la = 51078.69 - 4.338 * te - 0.06367 * te * te;
        double GL = (la / (rHa * pk)) * (es - ea[i]) * theta;
        // Add to total flux
        H[i] = HL.H[i] + GT;
        L[i] = HL.L[i] + GL;
        // Calculate source concentration
        ST[i] = (paii[i] / hgt) * HL.H[i];
        SL[i] = (paii[i] / hgt) * HL.L[i];
        // For near-field correction factor
        double btm = -0.399 * hgt * log(1 - exp(-abs(hgt / 2.0 - z[i]))) / nn;
        if (hgt / 2.0 != z[i]) sbtm = sbtm + btm;
        // Compute reference height near - field and far - field
        double Zeta = abs((hgt - z[i]) / (ow[i] * TL));
        double kn = -0.39894 * log(1.0 - exp(-Zeta)) - 0.15623 * exp(-Zeta);
        if (z[i] < hgt) {
            CnzrT = CnzrT + (ST[i] / ow[i]) * (kn * ((hgt - z[i]) / (ow[i] * TL)) + kn * ((hgt + z[i]) / (ow[i] * TL)));
            CnzrL = CnzrL + (SL[i] / ow[i]) * (kn * ((hgt - z[i]) / (ow[i] * TL)) + kn * ((hgt + z[i]) / (ow[i] * TL)));
        }
    }
    // Calculate near-field correction factor for small sample size
    double alpha = 1.31;
    if (hgt < 10) alpha = 0.756 + 0.3012 * hgt;
    double mu = alpha / sbtm;
    CnzrT = CnzrT * mu;
    CnzrL = CnzrL * mu;
    // Far-field concentration at top of canopy
    double te = (tlh + th) / 2;
    double cph = cpairCpp(th);
    double phh = phairCpp(th, pk);
    double lah = 45068.7 - 42.8428 * te;
    if (te < 0) lah = 51078.69 - 4.338 * te - 0.06367 * te * te;
    double CfTh = (th * cph * phh);
    double CfLh = (eh * phh * lah / pk);
    // Calculate near and far-field concentrations for all layers
    std::vector<double> tairn(tair.size());
    std::vector<double> ean(tair.size());
    double mxdif = 0.0;
    for (size_t i = 0; i < tair.size(); ++i) {
        // Near field
        double CnT = 0.0;
        double CnL = 0.0;
        for (size_t j = 0; j < tair.size(); ++j) {
            if (i != j) {
                double Zeta = abs((z[i] - z[j]) / (ow[j] * TL));
                double kn = -0.39894 * log(1.0 - exp(-Zeta)) - 0.15623 * exp(-Zeta);
                CnT = CnT + (ST[j] / ow[j]) * (kn * ((z[i] - z[j]) / (ow[j] * TL)) + kn * ((z[i] + z[j]) / (ow[j] * TL)));
                CnL = CnL + (SL[j] / ow[j]) * (kn * ((z[i] - z[j]) / (ow[j] * TL)) + kn * ((z[i] + z[j]) / (ow[j] * TL)));
            }
        }
        CnT = CnT * mu;
        CnL = CnL * mu;
        // Far-field
        double CfsT = 0.0;
        double CfsL = 0.0;
        for (size_t j = i; j < tair.size(); ++j) {
            CfsT = CfsT + H[j] / KH[j];
            CfsL = CfsL + L[j] / KH[j];
        }
        CfsT = CfsT * (hgt / nn);
        CfsL = CfsL * (hgt / nn);
        double CfT = CfTh - CnzrT + CfsT;
        double CfL = CfLh - CnzrL + CfsL;
        double CT = CfT + CnT;
        double CL = CfL + CnL;
        double cp = cpairCpp(tair[i]);
        double ph = phairCpp(tair[i], pk);
        double te = (tair[i] + tleafn[i]) / 2;
        double la = 45068.7 - 42.8428 * te;
        if (te < 0) la = 51078.69 - 4.338 * te - 0.06367 * te * te;
        tairn[i] = CT / (cp * ph);
        ean[i] = (CL * pk) / (la * ph);
        // Calculate tmx and tmn
        double tmx = tground;
        double tmn = tground;
        if (th > tmx) tmx = th;
        if (th < tmn) tmn = th;
        if (tleafn[i] > tmx) tmx = tleafn[i];
        if (tleafn[i] < tmn) tmn = tleafn[i];
        if (tairn[i] > tmx) tairn[i] = tmx;
        if (tairn[i] < tmn) tairn[i] = tmn;
        // Cap at dewpoint
        //double tdew = dewpointCpp(tairn[i], ean[i]);
        //if (tairn[i] < tdew) tairn[i] = tdew;
        // Calculate emx
        double emx = satvapCpp(tairn[i]);
        double emn = satvapCpp(tairn[i]) * 0.1;
        if (ean[i] > emx) ean[i] = emx;
        if (ean[i] < emn) ean[i] = emn;
        // Calculate mxdif
        if (abs(tairn[i] - tair[i]) > mxdif) mxdif = abs(tairn[i] - tair[i]);
        if (abs(tleafn[i] - tleaf[i]) > mxdif) mxdif = abs(tleafn[i] - tleaf[i]);
        if (abs(ean[i] - ea[i]) > mxdif) mxdif = abs(ean[i] - ea[i]);
    }
    Lang out;
    out.tleaf = tleafn;
    out.tair = tairn;
    out.ea = ean;
    out.uz = HL.uz;
    out.Rlwdown = lwrad.Rlwdown;
    out.Rlwup = lwrad.Rlwup;
    out.mxdif = mxdif;
    return out;
}
// Run below canopy model for one step;
// [[Rcpp::export]]
Rcpp::List SmallLeafOne(double reqhgt, double zref, double lat, double lon, std::vector<double> obsvars,
    std::vector<double> climvars, std::vector<double> bigleafvars, int maxiters, std::vector<double> wc, std::vector<double> vegp, 
    std::vector<double> paii, std::vector<double> groundp, std::vector<double> tleaf, std::vector<double> tair, std::vector<double> ea, 
    double surfwet, std::vector<double> z, double a0 = 0.25, double a1 = 1.25, double bwgt = 0.5)
{
    // Extract obsvars
    int year = obsvars[0];
    int month = obsvars[1];
    int day = obsvars[2];
    double hour = obsvars[3];
    // Extract climvars
    double tc = climvars[0]; 
    double rh = climvars[1]; 
    double pk = climvars[2]; 
    double Rsw = climvars[3];
    double Rdif = climvars[4];
    double lwdown = climvars[5];
    //double u2 = climvars[6];
    // Extract vegp
    double hgt = vegp[0];
    double pai = vegp[1];
    double x = vegp[2];
    double clump = vegp[3];
    double lref = vegp[4];
    double ltra = vegp[5];
    //double leafd = vegp[6];
    double vegem = vegp[7];
    //double gsmax = vegp[8];
    //double q50 = vegp[9];
    double lrefp = vegp[10];
    double ltrap = vegp[11];
    // Extract groundp
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    double groundem = groundp[3];
    //double rho = groundp[4];
    //double Vm = groundp[5];
    //double Vq = groundp[6];
    //double Mc = groundp[7];
    //double b = groundp[8];
    //double Psie = groundp[9];
    double Smax = groundp[10];
    double Smin = groundp[11];
    // Extract bigleafvars
    double Tc = bigleafvars[0]; // Canopy heat exchange surface temperature (deg C)
    double Tg = bigleafvars[1]; // Ground surface temperature (deg C)
    double psih = bigleafvars[2]; // Diabatic coefficient for heat
    double psim = bigleafvars[3]; // Diabatic coefficient for momentum
    double phih = bigleafvars[4]; // Diabatic influencing factor for heat
    double uf = bigleafvars[5]; // Friction velocity (m/s)
    double sm = bigleafvars[6]; // Fractional volumetric soil moisture
    // Calculate theta
    double theta = (sm - Smin) / (Smax - Smin);
    if (theta > 1.0) theta = 1.0;
    // Calculate variables at top of canopy
    if (hgt < 0.001) hgt = 0.001;
    double d = zeroplanedisCpp(hgt, pai);
    double zm = roughlengthCpp(hgt, pai, d, psih);
    double zh = 0.2 * zm;
    double lnm = log((hgt - d) / zm);
    double lnh = log((hgt - d) / zh);
    double lnrm = lnm / log((zref - d) / zm);
    double lnrh = lnh / log((zref - d) / zh);
    psim = lnrm * psim;
    psih = lnrh * psih;
    if (psim > 0.9 * lnm) psim = 0.9 * lnm;
    if (psih > 0.9 * lnh) psih = 0.9 * lnh;
    phih = pow(phih, lnrh);
    double uh = (uf / 0.4) * (log((hgt - d) / zm) + psim);
    double th = tc + (Tc - tc) * (1 - lnrh);
    double es = satvapCpp(Tc);
    double es2 = satvapCpp(tc);
    double e2 = es2 * (rh / 100);
    double eh = e2 + (es - e2) * (1 - lnrh) * surfwet;
    // Calculate shortwave
    radmodel2 swrad = RadiationSmallLeafSWCpp(lat, lon, year, month, day, hour, reqhgt, hgt, paii, x, lref, lrefp, ltra,
        ltrap, clump, vegem, gref, groundem, slope, aspect, Rsw, Rdif);
    double tlh = tleaf[tleaf.size() - 1];
    Lang mout = LangrangianOne(reqhgt, uh, th, tlh, eh, pk, lwdown, swrad, wc, vegp, paii, groundem, tleaf, tair, ea, Tg, surfwet,
        theta, psim, psih, phih, z, a0, a1);
    tleaf = mout.tleaf;
    tair = mout.tair;
    ea = mout.ea;
    double mxdif = mout.mxdif;
    int iters = 0;
    int tst = 1;
    if (maxiters > 1) {
        while (tst > 0) {
            mout = LangrangianOne(reqhgt, uh, th, tlh, eh, pk, lwdown, swrad, wc, vegp, paii, groundem, tleaf, tair, ea, Tg, surfwet,
                theta, psim, psih, phih, z, a0, a1);
            if (mout.mxdif > mxdif) bwgt = 1-(1-bwgt) * 0.8;
            for (size_t j = 0; j < tair.size(); ++j) {
                tleaf[j] = tleaf[j] * bwgt + (1 - bwgt) * mout.tleaf[j];
                tair[j] = tair[j] * bwgt + (1 - bwgt) * mout.tair[j];
                ea[j] = ea[j] * bwgt + (1 - bwgt) * mout.ea[j];
            }
            mxdif = mout.mxdif;
            ++iters;
            if (iters >= maxiters) tst = 0;
            if (mxdif < 0.1) tst = 0;
        }
    }
    Rcpp::List out;
    out["tleaf"] = Rcpp::wrap(tleaf);
    out["tair"] = Rcpp::wrap(tair);
    out["ea"] = Rcpp::wrap(ea);
    out["uz"] = Rcpp::wrap(mout.uz);
    out["Rdirdown"] = Rcpp::wrap(swrad.Rdirdown);
    out["Rdifdown"] = Rcpp::wrap(swrad.Rdifdown);
    out["Rlwdown"] = Rcpp::wrap(mout.Rlwdown);
    out["Rswup"] = Rcpp::wrap(swrad.Rswup);
    out["Rlwup"] = Rcpp::wrap(mout.Rlwup);
    out["mxdif"] = Rcpp::wrap(mxdif);
    return out;
}
// Run below canopy model for all time steps;
// [[Rcpp::export]]
DataFrame BelowCanopy(double reqhgt, double zref, double lat, double lon, DataFrame obstime, DataFrame climdata, DataFrame bigleafvars,
    int iters, std::vector<double> vegp, std::vector<double> paii, std::vector<double> groundp, double a0 = 0.25,
    double a1 = 1.25, double bwgt = 0.5)
{
    // Access items of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    // Access columns of bigleafvars
    std::vector<double> Tc = bigleafvars["Tc"]; // Canopy heat exchange surface temperature (deg C)
    std::vector<double> Tg = bigleafvars["Tg"]; // Ground surface temperature (deg C)
    std::vector<double> psih = bigleafvars["psih"]; // Diabatic coefficient for heat
    std::vector<double> psim = bigleafvars["psim"]; // Diabatic coefficient for momentum
    std::vector<double> phih = bigleafvars["phih"]; // Diabatic influencing factor for heat
    std::vector<double> uf = bigleafvars["uf"]; // Friction velocity (m/s)
    std::vector<double> sm = bigleafvars["sm"]; // Fractional volumetric soil moisture
    std::vector<double> surfwet = bigleafvars["surfwet"]; // surface wet
    // Calculate canopy wind shelter
    std::vector<double> wc = CanopyWindCpp(vegp[0], paii);
    // Calculate whichz
    std::vector<double> z(paii.size());
    double mindif = vegp[0];
    int whichz = 0;
    double nn = paii.size();
    for (size_t i = 0; i < paii.size(); ++i) {
        z[i] = ((i + 1) / nn) * vegp[0];
        double newdif = abs(z[i] - reqhgt);
        if (newdif < mindif) {
            whichz = i;
            mindif = newdif;
        }
    }
    // Create output data variables
    std::vector<double> tleafo(tc.size());
    std::vector<double> tairo(tc.size());
    std::vector<double> eao(tc.size());
    std::vector<double> uzo(tc.size());
    std::vector<double> Rdirdowno(tc.size());
    std::vector<double> Rdifdowno(tc.size());
    std::vector<double> Rlwdowno(tc.size());
    std::vector<double> Rswupo(tc.size());
    std::vector<double> Rlwupo(tc.size());
    // Initialise variables
    // Initialise variables
    double ea1 = satvapCpp(tc[0]) * rh[0] / 100.0;
    std::vector<double> tleaf(paii.size(), Tc[0]);
    std::vector<double> tair(paii.size(), tc[0]);
    std::vector<double> ea(paii.size(), ea1);
    // Run model in hourly time-steps
    for (size_t hr = 0; hr < tc.size(); ++hr) {
        // create inputs
        std::vector<double> obsvarsone{
            static_cast<double>(year[hr]),
            static_cast<double>(month[hr]),
            static_cast<double>(day[hr]),hour[hr] };
        std::vector<double> climvarsone{ tc[hr],rh[hr],pk[hr],Rsw[hr],Rdif[hr],Rlw[hr] };
        std::vector<double> blvarsone{ Tc[hr],Tg[hr],psih[hr],psim[hr],phih[hr],uf[hr],sm[hr] };
        Rcpp::List onerun = SmallLeafOne(reqhgt, zref, lat, lon, obsvarsone, climvarsone, blvarsone, iters, wc, vegp, paii, groundp,
            tleaf, tair, ea, surfwet[hr], z, a0, a1, bwgt);
        // Extract outputs
        tleaf = Rcpp::as<std::vector<double>>(onerun["tleaf"]);
        tair = Rcpp::as<std::vector<double>>(onerun["tair"]);
        ea = Rcpp::as<std::vector<double>>(onerun["ea"]);
        std::vector<double> uz = Rcpp::as<std::vector<double>>(onerun["uz"]);
        std::vector<double> Rdirdown = Rcpp::as<std::vector<double>>(onerun["Rdirdown"]);
        std::vector<double> Rdifdown = Rcpp::as<std::vector<double>>(onerun["Rdifdown"]);
        std::vector<double> Rlwdown = Rcpp::as<std::vector<double>>(onerun["Rlwdown"]);
        std::vector<double> Rswup = Rcpp::as<std::vector<double>>(onerun["Rswup"]);
        std::vector<double> Rlwup = Rcpp::as<std::vector<double>>(onerun["Rlwup"]);
        // store outputs
        tleafo[hr] = tleaf[whichz];
        tairo[hr] = tair[whichz];
        eao[hr] = ea[whichz];
        uzo[hr] = uz[whichz];
        Rdirdowno[hr] = Rdirdown[whichz];
        Rdifdowno[hr] = Rdifdown[whichz];
        Rlwdowno[hr] = Rlwdown[whichz];
        Rswupo[hr] = Rswup[whichz];
        Rlwupo[hr] = Rlwup[whichz];
    }
    std::vector<double> relhum(tc.size());
    for (size_t hr = 0; hr < tc.size(); ++hr) {
        relhum[hr] = (eao[hr] / satvapCpp(tairo[hr])) * 100.0;
        if (relhum[hr] > 100.0) relhum[hr] = 100.0;
    }
    Rcpp::DataFrame out;
    out["tair"] = Rcpp::wrap(tairo);
    out["tleaf"] = Rcpp::wrap(tleafo);
    out["relhum"] = Rcpp::wrap(relhum);
    out["windspeed"] = Rcpp::wrap(uzo);
    out["Rdirdown"] = Rcpp::wrap(Rdirdowno);
    out["Rdifdown"] = Rcpp::wrap(Rdifdowno);
    out["Rlwdown"] = Rcpp::wrap(Rlwdowno);
    out["Rswup"] = Rcpp::wrap(Rswupo);
    out["Rlwup"] = Rcpp::wrap(Rlwupo);
    return out;
}
// Run above canopy model for all time steps;
// [[Rcpp::export]]
Rcpp::DataFrame AboveCanopy(double reqhgt, double zref, double lat, double lon, DataFrame obstime, DataFrame climdata, DataFrame bigleafvars, std::vector<double> vegp)
{
    // Access items of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    // Access columns of bigleafvars
    std::vector<double> Tc = bigleafvars["Tc"]; // Canopy heat exchange surface temperature (deg C)
    std::vector<double> psih = bigleafvars["psih"]; // Diabatic coefficient for heat
    std::vector<double> psim = bigleafvars["psim"]; // Diabatic coefficient for momentum
    std::vector<double> uf = bigleafvars["uf"]; // Friction velocity (m/s)
    std::vector<double> albedo = bigleafvars["albedo"]; // albedo
    std::vector<double> surfwet = bigleafvars["surfwet"]; // albedo
    // Extract vegp
    double hgt = vegp[0];
    double pai = vegp[1];
    double vegem = vegp[7];
    // Calculate tair and ea
    std::vector<double> tair(tc.size());
    std::vector<double> ez(tc.size());
    std::vector<double> uz(tc.size());
    std::vector<double> Rswup(tc.size());
    std::vector<double> Rlwup(tc.size());
    std::vector<double> Rdir(tc.size());
    if (hgt < 0.001) hgt = 0.001;
    double d = zeroplanedisCpp(hgt, pai);
    for (size_t hr = 0; hr < tc.size(); ++hr) {
        double zm = roughlengthCpp(hgt, pai, d, psih[hr]);
        double zh = 0.2 * zm;
        double lnrh = log((reqhgt - d) / zh) / log((zref - d) / zh);
        double lnrm = log((reqhgt - d) / zm) / log((zref - d) / zm);
        // Calculate wind speed
        double psimz = lnrm * psim[hr];
        uz[hr] = (uf[hr] / 0.4) * (log((reqhgt - d) / zm) + psimz);
        // Calculate temperature
        tair[hr] = tc[hr] + (Tc[hr] - tc[hr]) * (1 - lnrh);
        // Calculate ez
        double ea = satvapCpp(tc[hr]) * rh[hr] / 100;
        double estl = satvapCpp(Tc[hr]);
        ez[hr] = ea + (estl - ea) * surfwet[hr] * (1 - lnrh);
        // Rswup and Rlwup
        Rswup[hr] = (1 - albedo[hr]) * Rsw[hr];
        Rlwup[hr] = vegem * 5.67 * pow(10, -8) * pow(Tc[hr] + 273.15, 4);
        // solar position
        std::vector<double> solp = solpositionCpp(lat, lon, year[hr], month[hr], day[hr], hour[hr]);
        double zen = solp[0];
        Rdir[hr] = (Rsw[hr] - Rdif[hr]) / cos(zen * M_PI / 180);
    }
    std::vector<double> relhum(tc.size());
    for (size_t hr = 0; hr < tc.size(); ++hr) relhum[hr] = (ez[hr] / satvapCpp(tair[hr])) * 100.0;
    Rcpp::DataFrame out;
    out["tair"] = Rcpp::wrap(tair);
    out["tcanopy"] = Rcpp::wrap(Tc);
    out["relhum"] = Rcpp::wrap(ez);
    out["windspeed"] = Rcpp::wrap(uz);
    out["Rdirdown"] = Rcpp::wrap(Rdir);
    out["Rdifdown"] = Rcpp::wrap(Rdif);
    out["Rlwdown"] = Rcpp::wrap(Rlw);
    out["Rswup"] = Rcpp::wrap(Rswup);
    out["Rlwup"] = Rcpp::wrap(Rlwup);
    return out;
}
// **  Function to compute rolling mean yearly ** // 
// [[Rcpp::export]]  
std::vector<double> manCpp(std::vector<double> x, int n) {
    std::vector<double> z(x.size());
    // Calculate rolling mean if n < 48
    if (n <= 48) {
        z = maCpp(x, n);
    }
    // Average to daily first if n > 48
    else {
        // Calculate daily mean
        int numDays = x.size() / 24;
        std::vector<double> d(numDays, 0.0);
        for (int i = 0; i < numDays; ++i) {
            double sum = 0.0;
            for (int j = 0; j < 24; ++j) {
                sum += x[i * 24 + j];
            }
            d[i] = sum / 24.0;
        }
        int n2 = n / 24;
        // Calculate rolling mean of daily
        std::vector<double> y = maCpp(d, n2);
        // Expand back out to hourly
        for (int i = 0; i < numDays; ++i) {
            for (int j = 0; j < 24; ++j) {
                z[i * 24 + j] = y[i];
            }
        }
        z = maCpp(z, 24);
    }
    return z;
}
// Run below ground model for all time steps;
// [[Rcpp::export]]
Rcpp::DataFrame Belowground(double reqhgt, DataFrame bigleafvars, std::vector<double> groundp)
{
    // extract ground params
    double rho = groundp[4];
    double Vm = groundp[5];
    double Vq = groundp[6];
    double Mc = groundp[7];
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1 - frs);
    double c3 = 1 + 2.6 * pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * frs * frs;
    // Extract bigleafparams
    std::vector<double> Tg = bigleafvars["Tg"]; // Ground surface temperature (deg C)
    std::vector<double> sm = bigleafvars["sm"]; // Fractional volumetric soil moisture
    // Calculate soil variables
    double sumDD = 0.0;
    for (size_t hr = 0; hr < sm.size(); ++hr) {
        double cs = 2400 * rho / 2.64 + 4180.0 * sm[hr]; // specific heat of soil in J / kg / K
        double ph = (rho * (1 - sm[hr]) + sm[hr]) * 1000.0; // bulk density in kg / m3
        double c2 = 1.06 * rho * sm[hr];
        double k = c1 + c2 * sm[hr] - (c1 - c4) * exp(-pow(c3 * sm[hr], 4.0)); // Thermal conductivity W / m / K
        double ka = k / (cs * ph);
        // Calculate damping depth
        double omdy = (2.0 * M_PI) / (24.0 * 3600.0);
        sumDD = sumDD + sqrt(2 * ka / omdy);
    }
    double meanDD = sumDD / sm.size();
    std::vector<double> Tz(Tg.size());
    // Calculate n (smoothing parameter)
    double nb = -118.35 * reqhgt / meanDD;
    int n = static_cast<int>(std::round(nb));
    if (n < static_cast<int>(Tg.size())) {
        Tz = manCpp(Tg, n);
    }
    else {
        double sumT = 0.0;
        for (size_t i = 0; i < Tg.size(); ++i) sumT = sumT + Tg[i];
        double meanT = sumT / Tg.size();
        for (size_t i = 0; i < Tg.size(); ++i) Tz[i] = meanT;
    }
    Rcpp::DataFrame out;
    out["tground"] = Rcpp::wrap(Tz);
    out["soilm"] = Rcpp::wrap(sm);
    return out;
}
// Run model for ground surface all time steps;
// [[Rcpp::export]]
Rcpp::DataFrame Atground(double lat, double lon, DataFrame obstime, DataFrame climdata, DataFrame bigleafvars,
    std::vector<double> vegp, std::vector<double> groundp)
{
    // Access items of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    //std::vector<double> wspeed = climdata["windspeed"];
    // Access columns of bigleafvars
    std::vector<double> Tc = bigleafvars["Tc"]; // Canopy heat exchange surface temperature (deg C)
    std::vector<double> Tg = bigleafvars["Tg"]; // Ground surface temperature (deg C)
    std::vector<double> sm = bigleafvars["sm"]; // Fractional volumetric soil moisture
    // Extract vegetation parameters
    double pai = vegp[1];
    double x = vegp[2];
    double clump = vegp[3];
    double lref = vegp[4];
    double ltra = vegp[5];
    double vegem = vegp[7];
    // extract ground params
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    double groundem = groundp[3];
    // Run shortwave radiation model
    // ** Calculate time-invariant variables
    double pait = pai / (1 - clump);
    double om = lref + ltra;
    double a = 1 - om;
    double del = lref - ltra;
    double J = 1.0 / 3.0;
    if (x != 1.0) {
        double mla = 9.65 * pow((3 + x), -1.65);
        if (mla > M_PI / 2) mla = M_PI / 2;
        J = cos(mla) * cos(mla);
    }
    double gma = 0.5 * (om + J * del);
    double h = sqrt(a * a + 2 * a * gma);
    // Calculate two-stream parameters (diffuse)
    std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref);
    double p1 = tspdif[0];
    double p2 = tspdif[1];
    double p3 = tspdif[2];
    double p4 = tspdif[3];
    double u1 = tspdif[4];
    double S1 = tspdif[5];
    double D1 = tspdif[6];
    double D2 = tspdif[7];
    // ** Upward diffuse stream (exluding contributions from direct and gaps
    double trd = clump * clump;
    double Rudm = (1 - trd) * (p1 * exp(-h * pait) + p2 * exp(h * pait)) + trd * gref;
    // ** Downward diffuse stream (exluding contributions from direct)
    double Rddm = (1 - trd) * (p3 * exp(-h * pait) + p4 * exp(h * pait)) + trd;
    // ** Define variables
    std::vector<double> Rdirdown(Rsw.size());
    std::vector<double> Rdifdown(Rsw.size());
    std::vector<double> Rlwdown(Rsw.size());
    std::vector<double> Rswup(Rsw.size());
    std::vector<double> Rlwup(Rsw.size());
    // ** Calculate limits on contribution of direct to diffuse
    double amx = gref;
    if (amx < lref) amx = lref;
    for (size_t i = 0; i < Rsw.size(); ++i) {
        if (Rsw[i] > 0) {
            // Calculate solar variables
            std::vector<double> solp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
            double zen = solp[0];
            double azi = solp[1];
            double si = solarindexCpp(slope, aspect, zen, azi);
            // Calculate canopy extinction coefficient
            double cosz = cos(zen * M_PI / 180);
            std::vector<double> kp = cankCpp(zen, x, si);
            double kd = kp[1];
            double Kc = kp[2];
            // Calculate two-stream parameters (direct)      
            double sig = kd * kd + gma * gma - pow((a + gma), 2);
            std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref, kd, sig, u1, S1, D1, D2);
            double p5 = tspdir[0];
            double p6 = tspdir[1];
            double p7 = tspdir[2];
            double p8 = tspdir[3];
            double p9 = tspdir[4];
            double p10 = tspdir[5];
            // ~~ Direct beam above canopy
            double Rbeam = (Rsw[i] - Rdif[i]) / cosz;
            if (Rbeam > 1352.0) Rbeam = 1352.0;
            double trb = pow(clump, Kc);
            double Rb = Rbeam * cosz;
            // Contribution of direct to upward diffuse stream
            double Rudbm = (1 - trd) * ((p5 / sig) * exp(-kd * pait) +
                p6 * exp(-h * pait) + p7 * exp(h * pait)) + trd * gref;
            if (Rudbm > amx) Rudbm = amx;
            if (Rudbm < 0.0) Rudbm = 0.0;
            // Contribution of direct to downward diffuse stream
            double Rddbm = (1 - trb) * ((p8 / -sig) * exp(-kd * pait) +
                p9 * exp(-h * pait) + p10 * exp(h * pait));
            if (Rddbm > amx) Rddbm = amx;
            if (Rddbm < 0.0) Rddbm = 0.0;
            // Radiation streams
            Rdirdown[i] = ((1.0 - trb) * exp(-kd * pait) + trb) * Rbeam;
            Rdifdown[i] = Rddm * Rdif[1] + Rb * Rddbm;
            Rswup[i] = Rudm * Rdif[1] + Rb * Rudbm;
        }
        else {
            Rdirdown[i] = 0.0;
            Rdifdown[i] = 0.0;
            Rswup[i] = 0.0;
        }
        double lwcan = vegem * 5.67 * pow(10, -8.0) * pow(Tc[i] + 273.15, 4.0);
        Rlwdown[i] = trd * Rlw[i] + (1 - trd) * lwcan;
        Rlwup[i] = groundem * 5.67 * pow(10, -8.0) * pow(Tg[i] + 273.15, 4.0);

    }
    Rcpp::DataFrame out;
    out["tground"] = Rcpp::wrap(Tg);
    out["soilm"] = Rcpp::wrap(sm);
    out["Rdirdown"] = Rcpp::wrap(Rdirdown);
    out["Rdifdown"] = Rcpp::wrap(Rdifdown);
    out["Rlwdown"] = Rcpp::wrap(Rlwdown);
    out["Rswup"] = Rcpp::wrap(Rswup);
    out["Rlwup"] = Rcpp::wrap(Rlwup);
    return out;
}
// Run below or above canopy model for all time steps;
// [[Rcpp::export]]
Rcpp::DataFrame runmodel(double reqhgt, double zref, double lat, double lon, DataFrame obstime, DataFrame climdata, DataFrame bigleafvars,
    int iters, std::vector<double> vegp, std::vector<double> paii, std::vector<double> groundp, double a0 = 0.25,
    double a1 = 1.25, double bwgt = 0.5)
{
    double hgt = vegp[0];
    Rcpp::DataFrame out;
    if (reqhgt >= hgt) {
        out = AboveCanopy(reqhgt, zref, lat, lon, obstime, climdata, bigleafvars, vegp);
    }
    else {
        if (reqhgt > 0) {
            out = BelowCanopy(reqhgt, zref, lat, lon, obstime, climdata, bigleafvars, iters, vegp, paii, groundp, a0, a1, bwgt);
        }
        else {
            if (reqhgt == 0) {
                out = Atground(lat, lon, obstime, climdata, bigleafvars, vegp, groundp);
            }
            else {
                out = Belowground(reqhgt, bigleafvars, groundp);
            }
        } // end else reqhgt <= 0
    } // end else reqhgt < hgt
    return out;
}