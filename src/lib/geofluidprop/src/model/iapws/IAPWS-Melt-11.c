
#include "IAPWS-Melt-11.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "util/utility.h"

static const double pt = 611.6547711; // [Pa]
static const double Tt = 273.16;  // [K]

// TODO: multiple solution can exist at T<= Tt
// double iapws11_melt_presssure(double T)
// {
//     assert(T >= 251.165);
//     if (251.165 <= T <= Tt)
//         return iapws11_melt_presssure_ice_Ih(T);
//     if (T < 256.164)
//         return iapws11_melt_presssure_ice_III(T);
//     if (T < 273.31)
//         return iapws11_melt_presssure_ice_V(T);
//     if (T < 355)
//         return iapws11_melt_presssure_ice_VI(T);
//     if (T <= 715)
//         return iapws11_melt_presssure_ice_VII(T);
//     assert(T > 715);
//     return -1;
// }

// 273.16 ~ 251.165
double iapws11_melt_presssure_ice_Ih(double T)
{
    assert(251.165 <= T && T <= Tt);
    double theta = T/Tt;
    double pi = 1.
                + 0.119539337e7*(1.-pow(theta,0.3e1))
                + 0.808183159e5*(1.-pow(theta,0.2575e2))
                + 0.333826860e4*(1.-pow(theta,0.103750e3));
    return pi*pt;
}

// 251.165 ~ 256.164
double iapws11_melt_presssure_ice_III(double T)
{
    assert(251.165 <= T && T <= 256.164);
    double theta = T/251.165;
    double pi = 1.- 0.299948*(1-pow(theta, 60));
    return pi*208.566e6;
}

// 256.164~ 273.31
double iapws11_melt_presssure_ice_V(double T)
{
    assert(256.164 <= T && T <= 273.31);
    double theta = T/256.164;
    double pi = 1.- 1.18721*(1-pow(theta, 8));
    return pi*350.1e6;
}

// 273.31~355
double iapws11_melt_presssure_ice_VI(double T)
{
    assert(273.31 <= T && T <= 355);
    double theta = T/ 273.31;
    double pi = 1.- 1.07476*(1-pow(theta, 4.6));
    return pi*632.4e6;
}

// 355~715
double iapws11_melt_presssure_ice_VII(double T)
{
    assert(355 <= T && T <= 715);
    double theta = T/355;
    double pi = 0.173683e1*(1-pow(theta, -1))
                -0.544606e-1*(1-pow(theta, 5))
                +0.806106e-7*(1-pow(theta, 22));
    return exp(pi)*2216.e6;
}

// 50 ~ 273.16
double iapws11_sublimation_pressure(double T)
{
    //printf("T=%e\n",T);
    assert(50 <= T && T <= Tt);
    double theta = T/Tt;
    double pi = 1./theta *(
                    -0.212144006e2*pow(theta,0.333333333e-2)
                    +0.273203819e2*pow(theta,0.120666667e1)
                    -0.610598130e1*pow(theta,0.170333333e1)
                );
    return exp(pi)*pt;
}

double iapws11_invert(double p, double (*f)(double), double Tmin, double Tmax)
{
    double T = 0.5*(Tmin + Tmax);
    for (int i=0; i<50; i++)
    {
        double r = (*f)(T) - p;
        //printf("%d: T=%.4e, r=%.3e\n", i, T, r);
        if (fabs(r)<p*1e-8){
            //printf("converged with %d iterations: r=%.3e, err=%.3e\n", i, r, fabs(r)/p*1e6);
            return T;
        }
        double h = T*1e-8;
        double df = ((*f)(T+h)-(*f)(T-h))/(2*h);
        T -= r/df;
        //printf("%d: T0=%.4e\n", i, T);
        if (T<Tmin) T = Tmin*(1.+1e-6);
        if (Tmax<T) T = Tmax*(1.-1e-6);
    }
    printf("**ERROR: invert() did not converge\n");
    return -1;
}

double iapws11_melt_temperature(double p)
{
    if (p < pt) { //6.116550000e+02) {
        double Tm = iapws11_invert(p, iapws11_sublimation_pressure, 50, Tt);
        return min(Tm, Tt);
    } else if (p < 2.085658841e+08) {
        // Ih
        double Tm = iapws11_invert(p, iapws11_melt_presssure_ice_Ih, 251.165, Tt);
        //return Tm;
        return max(min(Tm, Tt), 251.165);
    } else if (p < 3.501000157e+08) {
        // III
        return iapws11_invert(p, iapws11_melt_presssure_ice_III, 251.165, 256.164);
    } else if (p < 6.323993474e+08) {
        // V
        return iapws11_invert(p, iapws11_melt_presssure_ice_V, 256.164, 273.31);
    } else if (p < 2.216002257e+09) {
        // VI
        return iapws11_invert(p, iapws11_melt_presssure_ice_VI, 273.31, 355);
    } else {
        // VII
        return iapws11_invert(p, iapws11_melt_presssure_ice_VII, 355, 715);
    }

    return 0;
}
