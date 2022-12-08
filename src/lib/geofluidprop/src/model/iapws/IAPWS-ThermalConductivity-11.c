
#include "IAPWS-ThermalConductivity-11.h"

#include <math.h>

#include "util/utility.h"

#include "IAPWS-95.h"
#include "IAPWS-Viscosity-08.h"


// Reference parameter values
static const double Tr = 647.096;
static const double pr = 22.064e6;
static const double rhor = 322.0;
static const double lambdar = 1.e-3;
static const double mur = 1.e-6;

// Constants
static const double R = 0.46151805e3;

// L_ij
static const double L[5*6] = {
    1.60397357,    -0.646013523,    0.111443906,    0.102997357,    -0.0504123634,    0.00609859258,
    2.33771842,    -2.78843778,    1.53616167,        -0.463045512,    0.0832827019,    -0.00719201245,
    2.19650529,    -4.54580785,    3.55777244,        -1.40944978,    0.275418278,    -0.0205938816,
    -1.21051378,    1.60812989,    -0.621178141,    0.0716373224,    0,                0,
    -2.7203370,    4.57586331,        -3.18369245,    1.1168348,        -0.19268305,    0.012913842
};

// Critical region constants
static const double Lambda = 177.8514;
static const double qd = 1./0.4;
static const double v = 0.630;
static const double gamma_ = 1.239;
static const double xi0 = 0.13;
static const double Gamma0 = 0.06;
static const double TbarR = 1.5;


double iapws11_lambda0bar(double Tbar)
{
    double a =   2.443221e-3  / pow(Tbar,0)
               + 1.323095e-2  / pow(Tbar,1)
               + 6.770357e-3  / pow(Tbar,2)
               - 3.454586e-3  / pow(Tbar,3)
               + 4.096266e-4  / pow(Tbar,4);
    return sqrt(Tbar)/a;
}

double iapws11_lambda1bar(double Tbar, double rhobar)
{
    double a = 0;
    for (int i=0; i<5; i++)
    {
        double b = 0;
        for (int j=0; j<6; j++)
            b += L[i*6+j] * pow(rhobar-1, j);
        a += pow(1/Tbar-1,i)*b;
    }
    a *= rhobar;
    return exp(a);
}


static double sigma(double Tbar, double rhobar)
{
    // (drhobar/dpbar)_T = pr/rhor*(dp/drho)^-1
    return pr/rhor/iapws95_dp_drho_rhoT(rhobar*rhor, Tbar*Tr);
}

double iapws11_lambda2bar(double Tbar, double rhobar)
{
    if (rhobar==0.0)
        return 0.0;
    double dchi = rhobar*(sigma(Tbar, rhobar) - sigma(TbarR,rhobar)*TbarR/Tbar);
    dchi = max(dchi, 0);
    double xi = xi0*pow(dchi/Gamma0, v/gamma_);
    double y = qd*xi;
    if (y<1.2e-7)
        return 0.0;

    double cp1 = iapws95_cp_rhoT(rhobar*rhor, Tbar*Tr);
    double cv1 = iapws95_cv_rhoT(rhobar*rhor, Tbar*Tr);
    double cpbar = cp1/R;
    double mubar = iapws08_viscosity_rhoT(rhobar*rhor,Tbar*Tr) / mur;

    double Z = 0.0;
    double invk = cv1/cp1;
    Z = (1-invk)*atan(y) + invk*y - (1-exp(-1/(1/y+y*y/(3*rhobar*rhobar))));
    Z *= 2./(M_PI*y);

    return Lambda*rhobar*cpbar*Tbar/mubar*Z;
}

double iapws11_thermal_conductivity_rhoT(double rho, double T)
{
    double Tbar = T/Tr;
    double rhobar = rho/rhor;
    return lambdar*(iapws11_lambda0bar(Tbar)*iapws11_lambda1bar(Tbar,rhobar)+iapws11_lambda2bar(Tbar,rhobar));
}
