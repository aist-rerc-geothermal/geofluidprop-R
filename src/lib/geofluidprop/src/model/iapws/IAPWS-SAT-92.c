
#include "IAPWS-SAT-92.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/utility.h"

#define Tc 647.096
#define pc 22.064e6
#define rhoc 322.0

#define alpha0 1000 // [J/kg]
#define phi0 1000./647.096

// static const double Tc = 647.096;
// static const double pc = 22.064e6;
// static const double rhoc = 322.0;
// // //
// static const double alpha0 = 1000; // [J/kg]
// static const double phi0 = 1000./647.096;
// //
#define a SAT92_a
#define b SAT92_b
#define c SAT92_c
#define d SAT92_d

const double SAT92_a[6] = {
    -7.85951783,
     1.84408259,
    -11.7866497,
     22.6807411,
    -15.9618719,
     1.80122502
};

const double SAT92_b[6] = {
     1.99274064,
     1.09965342,
    -0.510839303,
    -1.75493479,
    -45.5170352,
    -6.74694450e5
};

const double SAT92_c[6] = {
    -2.03150240,
    -2.68302940,
    -5.38626492,
    -17.2991605,
    -44.7586581,
    -63.9201063
};

const double SAT92_d[5] = {
    -5.65134998e-8,
     2690.66631,
     127.287297,
    -135.003439,
     0.981825814
};

#define d_alpha -1135.905627715
#define d_phi 2319.5246
// static const double d_alpha = -1135.905627715;
// static const double d_phi = 2319.5246;
// //
double iapws92_sat_p_T(double T)
{
    double theta = T/Tc;
    double tau = 1.-theta;

    double p = a[0]*tau+a[1]*pow(tau,1.5)+a[2]*pow(tau,3)
                +a[3]*pow(tau,3.5)+a[4]*pow(tau,4)+a[5]*pow(tau,7.5);
    p = pc*exp(Tc/T*p);
    return p;
}

double iapws92_sat_T_p(double p)
{
    double T = 273.15+100;
    for (int i=0; i<25; i++)
    {
        double p0 = iapws92_sat_p_T(T);
        double r = p0 - p;
//         printf("%d: T=%.10e, p=%.10e, r=%g\n", i, T, p0, r);
        if (fabs(r)<1.e-8*p)
            return T;
        double drdT = iapws92_sat_dp_dT_T(T);
        // double d = 100.*1e-8;
        // double drdT = (T+d<Tc) ? (vapor_pressure(T+d) - p0)/d : (p0-vapor_pressure(T-d))/d;
        T -= r/drdT;
        T = max(T, 273.16);
        T = min(T, Tc);
    }
    printf("%s: diverged for p=%.10e\n", __FUNCTION__, p);
    abort();
    return -1;
}

double iapws92_sat_dp_dT_T(double T)
{
    double theta = T/Tc;
    double tau = 1.-theta;
    double dtau = -1./Tc;

    double tmp = a[0]*tau+a[1]*pow(tau,1.5)+a[2]*pow(tau,3)
                +a[3]*pow(tau,3.5)+a[4]*pow(tau,4)+a[5]*pow(tau,7.5);
    double dtmp = a[0]*dtau+1.5*a[1]*pow(tau,0.5)*dtau+3*a[2]*pow(tau,2)*dtau
                +3.5*a[3]*pow(tau,2.5)*dtau+4*a[4]*pow(tau,3)*dtau
                +7.5*a[5]*pow(tau,6.5)*dtau;
    double p = pc*exp(Tc/T*tmp)*(-Tc/(T*T)*tmp + Tc/T*dtmp);
    return p;
}

double iapws92_sat_rhol_T(double T)
{
    double theta = T/Tc;
    double tau = 1.-theta;
    return rhoc*(1+b[0]*pow(tau,1./3)+b[1]*pow(tau,2./3)+b[2]*pow(tau,5./3)
            +b[3]*pow(tau,16./3)+b[4]*pow(tau,43./3)+b[5]*pow(tau,110./3));
}

double iapws92_sat_drhol_dT_T(double T)
{
    double theta = T/Tc;
    double tau = 1.-theta;
    double h = T*1e-6;
    return (iapws92_sat_rhol_T(T+h)-iapws92_sat_rhol_T(T-h))/(2*h);
    // return rhoc*(1+b[0]*pow(tau,1./3)+b[1]*pow(tau,2./3)+b[2]*pow(tau,5./3)
    //         +b[3]*pow(tau,16./3)+b[4]*pow(tau,43./3)+b[5]*pow(tau,110./3));
}

double iapws92_sat_rhov_T(double T)
{
    double theta = T/Tc;
    double tau = 1.-theta;
    return rhoc*exp(c[0]*pow(tau,2./6)+c[1]*pow(tau,4./6)+c[2]*pow(tau,8./6)
            +c[3]*pow(tau,18./6)+c[4]*pow(tau,37./6)+c[5]*pow(tau,71./6));
}

double iapws92_sat_alpha(double T)
{
    double theta = T/Tc;
    return alpha0*(d_alpha+d[0]*pow(theta, -19)+d[1]*theta+d[2]*pow(theta, 4.5)
            +d[3]*pow(theta, 5)+d[4]*pow(theta, 54.5));
}

double iapws92_sat_phi(double T)
{
    double theta = T/Tc;
    return phi0*(d_phi+19./20*d[0]*pow(theta, -20)+d[1]*log(theta)+9./7*d[2]*pow(theta, 3.5)
            +5./4*d[3]*pow(theta, 4)+109./107*d[4]*pow(theta, 53.5));
}

double iapws92_sat_hl_T(double T)
{
    return (iapws92_sat_alpha(T)+T/iapws92_sat_rhol_T(T)*iapws92_sat_dp_dT_T(T));
}

double iapws92_sat_hv_T(double T)
{
    return (iapws92_sat_alpha(T)+T/iapws92_sat_rhov_T(T)*iapws92_sat_dp_dT_T(T));
}

double iapws92_sat_sl_T(double T)
{
    return (iapws92_sat_phi(T)+1/iapws92_sat_rhol_T(T)*iapws92_sat_dp_dT_T(T));
}

double iapws92_sat_sg_T(double T)
{
    return (iapws92_sat_phi(T)+1/iapws92_sat_rhov_T(T)*iapws92_sat_dp_dT_T(T));
}
