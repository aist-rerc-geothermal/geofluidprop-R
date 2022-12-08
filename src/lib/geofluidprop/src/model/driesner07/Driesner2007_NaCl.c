
#include "Driesner2007_NaCl.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "util/utility.h"

extern double driesner07_H2O_NaCl_h_pTx(double p_Pa, double T_K, double x);

//-------------------------------------------------------------------
// Pure NaCl constants
//-------------------------------------------------------------------
static double const Tt_NaCl = 800.7 + 273.15; // [K]
static double const Tt_NaCl_C = 800.7;        // [C]
static double const pt_NaCl = 50.0;           // [Pa]
static double const pt_NaCl_bar = 50.0e-5;    // [bar]

//-------------------------------------------------------------------
// Phase boundary
//-------------------------------------------------------------------
double driesner07_NaCl_VLH_p() {return pt_NaCl; }
double driesner07_NaCl_VLH_T() {return Tt_NaCl; }
double driesner07_NaCl_VLH_TC() {return Tt_NaCl_C; }
double driesner07_NaCl_VLH_h_halite()
{
    double h_moltenNaCl = driesner07_H2O_NaCl_h_pTx(driesner07_NaCl_VLH_p(), driesner07_NaCl_VLH_T(), 1.);
    double h_fusion = 29.8e3;
    return h_moltenNaCl + h_fusion;
}

double driesner07_NaCl_LH_T_p(double p)
{
    double Thm = Tt_NaCl_C + 2.47260e-2*(to_bar(p - pt_NaCl));
    return to_K(Thm);
}

double driesner07_NaCl_VH_p_T(double T)
{
    double logp = log10(pt_NaCl_bar) + 1.18061e4*(1/Tt_NaCl-1/T);
    return to_Pa(pow(10, logp));
}

double driesner07_NaCl_VL_p_T(double T)
{
    double logp = log10(pt_NaCl_bar) + 0.941812e4*(1/Tt_NaCl-1/T);
    return to_Pa(pow(10, logp));
}

double driesner07_NaCl_sat_p_T(double T)
{
    if (T<Tt_NaCl)
        return driesner07_NaCl_VH_p_T(T);
    else
        return driesner07_NaCl_VL_p_T(T);
}

//-------------------------------------------------------------------
// Halite property
//-------------------------------------------------------------------
double driesner07_NaCl_H_rho_pT(double p_Pa, double T_K)
{
    double p = to_bar(p_Pa);
    double T = to_C(T_K);
    double l0 = 2.1704e3;
    double l1 = -2.4599e-1;
    double l2 = -9.5797e-5;
    double l3 = 5.727e-3;
    double l4 = 2.715e-3;
    double l5 = 733.4;
    double rho0 = l0 + l1*T + l2*T*T;
    double l = l3 + l4*exp(T/l5);
    double rho = rho0 + l*p;
    return rho;
}

double driesner07_NaCl_H_h_pT(double p_Pa, double T_K)
{
    double h_tr = driesner07_NaCl_VLH_h_halite();
    double cp = driesner07_NaCl_H_cp_pT(p_Pa, T_K);

    double h = h_tr;
    int n = 10;
    double dT = (T_K - driesner07_NaCl_VLH_T())/n;
    for (int i = 0; i < n; i++)
    {
        double T2 = driesner07_NaCl_VLH_T() + dT*(i+1);
        double cp = driesner07_NaCl_H_cp_pT(p_Pa, T2);
        h += cp*dT;
    }

    return h;
}

double driesner07_NaCl_H_cp_pT(double p_Pa, double T_K)
{
    double p = to_bar(p_Pa);
    double T = to_C(T_K);
    double Ttr = to_C(driesner07_NaCl_VLH_T());

    double r0 = 1148.81;
    double r1 = 0.275774;
    double r2 = 8.8103e-5;
    double r3 = -1.7099e-3 - 3.82734e-6*T - 8.65455e-9*T*T;
    double r4 = 5.29063e-8 - 9.63084e-11 + 6.50745e-13;

    double cp = r0 + 2*r1*(T-Ttr) + 3*r2*(T-Ttr)*(T-Ttr) + r3*p + r4*p*p;
    return cp;
}

//-------------------------------------------------------------------
// Liquid NaCl property
//-------------------------------------------------------------------
double driesner07_NaCl_L_rho_pT(double p_Pa, double T_K)
{
    double p = to_bar(p_Pa);
    double T = to_C(T_K);

    double m0 = 58443;
    double m1 = 23.772;
    double m2 = 0.018639;
    double m3 = -1.9687e-6;
    double rho0 = m0/(m1 + m2*T + m3*T*T);
    double K = driesner07_NaCl_L_compressibility_pT(T_K);
    double rho = rho0/(1. - 0.1*log(1 + 10*p*K));
    return rho;
}

double driesner07_NaCl_L_compressibility_pT(double T_K)
{
    double T = to_C(T_K);
    double K = -1.5259e-5 + 5.5058e-8*T;
    return K;
}

