
#include "Driesner2007_H2ONaCl_prost.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <steam4.h>

#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-IF97.h"

#include "util/utility.h"
#include "Driesner2007_NaCl.h"


#pragma region water
static double h2o_sat_p(double T)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(T, propl, propg);
    double p = propl->p;
    freeProp(propl);
    freeProp(propg);
    return p;
}

static double h2o_sat_T(double p)
{
    Prop* propl = newProp('p', 'p', 0);
    Prop* propg = newProp('p', 'p', 0);
    sat_p(p, propl, propg);
    double T = propl->T;
    freeProp(propl);
    freeProp(propg);
    return T;
}

static void h2o_sat_Trholv_p(double p, double* Tsat, double* rhol, double* rhov)
{
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_p(p, pliq, pvap);
    *Tsat = pliq->T;
    *rhol = pliq->d;
    *rhov = pvap->d;
    freeProp(pliq);
    freeProp(pvap);
}

static double h2o_rho_pT(double p, double T)
{
    Prop* prop = newProp('p', 'T', 0);
    water_tp(T, p, 1.e3, 1.e-8, prop);
    double rho = prop->d;
    freeProp(prop);
    return rho;
}

static double h2o_h_pT(double p, double T)
{
    Prop* prop0 = newProp('p', 'T', 0);
    water_tp(T, p, 0, 1.e-8, prop0);
    double h = prop0->h;
    freeProp(prop0);
    return h;
}

static double h2o_cp_pT(double p, double T)
{
    Prop* prop0 = newProp('p', 'T', 0);
    water_tp(T, p, 0, 1.e-8, prop0);
    //dumpProp(stdout, prop0);
    //double rho = prop0->d;
    double cp = prop0->cp;
    freeProp(prop0);
    return cp;
}

#pragma endregion


//-------------------------------------------------------------------
// Thermodynamic constants
//-------------------------------------------------------------------
#pragma region constants
static const double H2O_Tc = 647.126; //373.976+273.15;
static const double H2O_pc = 22.054915e6;
static const double H2O_molar_mass  = 18.015e-3; // kg/mol
static const double NaCl_molar_mass = 58.443e-3; // kg/mol

// critical curve
static const double c[14] = {
    -2.36,
     1.28534e-1,
    -2.3707e-2,
     3.20089e-3,
    -1.38917e-4,
     1.02789e-7,
    -4.8376e-11,
     2.36,
    -1.31417e-2,
     2.98491e-3,
    -1.30114e-4,
     5.810101498620e+02, // pc at 500C
     3.405488304672,     // dpc/dT at 500C
    -4.88336e-4
};

static const double cA[11] = {
    1.,
    1.5,
    2.,
    2.5,
    3.,
    4.,
    5.,
    1.,
    2.,
    2.5,
    3.
};


static const double d[11] = {
    8.00000e-05,
    1.00000e-05,
    -1.37125e-07,
    9.46822e-10,
    -3.50549e-12,
    6.57369e-15,
    -4.89423e-18,
    7.77761e-2,
    2.7042e-4,
    -4.244821e-7,
    2.580872e-10
};
#pragma endregion

//-------------------------------------------------------------------
// convinient functions
//-------------------------------------------------------------------
#pragma region util
// convert kg/m3 to m3/mol
static inline double to_molar_volume(double rho, double molar_mass)
{
    // mol/m3 = rho/Mi
    return 1./(rho/molar_mass);
}

// drho to dV
static inline double to_dmolar_volume(double rho, double drho, double molar_mass)
{
    return -1./(rho*rho/molar_mass)*drho;
}

static double to_H2O_molar_volume(double rho)
{
    return to_molar_volume(rho, H2O_molar_mass);
}

static double to_dH2O_molar_volume(double rho, double drho)
{
    return to_dmolar_volume(rho, drho, H2O_molar_mass);
}

static double to_NaCl_molar_volume(double rho)
{
    return to_molar_volume(rho, NaCl_molar_mass);
}

// convert mass volume to molar volume
static double H2O_to_vm(double v_m3_kg)
{
    double H2O_M = 18.015e-3; // kg/mol
    return v_m3_kg * H2O_M;
}

#pragma endregion



#pragma region critical
double driesner07_prost_H2O_NaCl_pc_T(double T)
{
    return driesner07_prost_H2O_NaCl_pc_T2(T, false);
}

// Eqs. (5a)~(5c)
double driesner07_prost_H2O_NaCl_pc_T2(double T, bool use_satp)
{
    double const Tc_H2O = H2O_Tc;
    double const pc_H2O = H2O_pc;
    double pc = 0;
    if (T < Tc_H2O) {
        if (use_satp) {
            pc = h2o_sat_p(T);
        } else {
            pc = 0;
            for (int i = 0; i < 7; i++) {
                pc += c[i] * pow(to_C(Tc_H2O) - to_C(T), cA[i]);
            }
            pc = to_Pa(pc) + pc_H2O;
        }
    } else if (T<=to_K(500)) {
        for (int i=7; i<11; i++)
            pc += c[i]*pow(T-Tc_H2O, cA[i]);
        pc = pc*1e5 + pc_H2O;
        // // calculate c12, c13
        // double dpc = 0;
        // for (int i=7; i<11; i++)
        //     dpc += c[i]*cA[i]*pow(T-Tc_H2O, cA[i]-1.);
        // dpc = dpc*1e5;
        // printf("T=%g[C]: pc=%15.12e, dpc/dT=%15.12e\n", T-273.15, pc*1e-5, dpc*1e-5);
    } else {
        double const dT = to_C(T)-500;
        for (int i=11; i<14; i++)
            pc += c[i]*pow(dT, i+1-12);
        pc = pc*1e5;
    }
    return pc;
}


double driesner07_prost_H2O_NaCl_xc_T(double T)
{
    double const Tc_H2O = H2O_Tc;
    double xc = 0;
    if (Tc_H2O <= T && T <= to_K(600)) {
        for (int i=0; i<7; i++)
            xc += d[i]*pow(T-Tc_H2O, i+1);
    } else if (to_K(600) < T && T <= to_K(1000)) {
        for (int i=7; i<11; i++)
            xc += d[i]*pow(T-to_K(600), i+1-8);
    } else {
        return 0;
        assert(false);
    }
    return xc;
}

double driesner07_prost_H2O_NaCl_Tc_x(double x)
{
    if (x==0.) return H2O_Tc;
    double T = H2O_Tc + 100.;
    const double d = 1e-8;
    for (int i=0; i<100; i++)
    {
        double x1 = driesner07_prost_H2O_NaCl_xc_T(T);
        double r = x1 - x;
        //printf("%d: T=%g, x=%g, r=%g\n", i, T, x1, r);
        if (fabs(r) < 1e-15)
            return T;
        double drdT = (driesner07_prost_H2O_NaCl_xc_T(T+d) - x1)/d;
        T += -r/drdT;
        T = max(T, H2O_Tc+1.e-12);
    }
    printf("%s: diverged for x=%g\n", __FUNCTION__, x);
    return -1;
}
#pragma endregion

//-------------------------------------------------------------------
// vapor-liquid
//-------------------------------------------------------------------
#pragma region VL

// Normalized K fom Eqs. (17)
double driesner07_prost_H2O_NaCl_VL_log10_Kn(double pn, double T_C)
{
    if (pn > 1.0 && pn < (1.0+1e-6)) pn = 1.0;
    //assert(pn <= 1.0);
    double k0 = -0.235694;
    double k1 = -0.188838;
    double k2 = 0.004;
    double k3 = 0.0552466;
    double k4 = 0.66918;
    double k5 = 396.848;
    double k6 = 45.0;
    double k7 = -3.2719e-7;
    double k8 = 141.699;
    double k9 = -0.292631;
    double k10 = -0.00139991;
    double k11 = 1.95965e-6;
    double k12 = -7.3653e-10;
    double k13 = 0.904411;
    double k14 = 0.000769766;
    double k15 = -1.18658e-6;

    double j0 = k0 + k1*exp(-k2*T_C);
    double j1 = k4 + (k3-k4)/(1.+exp((T_C-k5)/k6))
                + k7 * pow(T_C + k8, 2.);
    double j2 = k9 + k10*T_C + k11*T_C*T_C + k12*T_C*T_C*T_C;
    double j3 = k13 + k14*T_C + k15*T_C*T_C;

    double v = 1. + j0 * pow(1.-pn, j1)
                  + j2 * (1.-pn)
                  + j3 * pow(1.-pn, 2.)
                  - (1.+j0+j2+j3) * pow(1.-pn, 3.);



//    assert(0<=v && v<=1.);

    return v;
}

// Modified K from Eqs. (15)
double driesner07_prost_H2O_NaCl_VL_log10_Km(double p_bar, double T_K)
{
    double T_C = to_C(T_K);
    double p_NaCl = to_bar(driesner07_NaCl_VH_p_T(T_K));
    double pc = to_bar(driesner07_prost_H2O_NaCl_pc_T(T_K));
    double pn = (p_bar - p_NaCl)/(pc - p_NaCl);
    double log10Kn = driesner07_prost_H2O_NaCl_VL_log10_Kn(pn, T_C);
    double x_sat_l = driesner07_prost_H2O_NaCl_LH_xl_pT(to_Pa(p_NaCl), T_K);

    double v = log10Kn*(log10(p_NaCl/pc)-log10(x_sat_l)) + log10(x_sat_l);
    // printf("p_NaCl=%g, pc=%g, pn=%g\n", p_NaCl, pc, pn);
    // printf("log10Kn=%g, x_sat_l=%g, log10Km=%g\n", log10Kn, x_sat_l, v);
    return v;
}

// K from Eqs. (14)
double driesner07_prost_H2O_NaCl_VL_K(double p, double T)
{
    double p_bar = to_bar(p);
    double p_NaCl = to_bar(driesner07_NaCl_VH_p_T(T));
    double log_Km = driesner07_prost_H2O_NaCl_VL_log10_Km(p_bar, T);
    double log_K = log_Km - log10(p_NaCl/p_bar);
    return pow(10., log_K);
}


double driesner07_prost_H2O_vm_pT(double p, double T)
{
    double rho = h2o_rho_pT(p, T);
    double v = 1./rho; // m3/kg
    return H2O_to_vm(v);
}

double driesner07_prost_H2O_NaCl_VL_xl_pT(double p_Pa, double T_K)
{
    double const p = to_bar(p_Pa);
    double const T = to_C(T_K);

    double const h[11] = {
        1.68486e-3, 2.19379e-4, 4.3858e2, 1.84508e1, -5.6765e-10,
        6.73704e-6, 1.44951e-7, 3.84904e2, 7.07477e0, 6.06896e-5, 7.62859e-3};

    double const g1 = h[1] + (h[0]-h[1])/(1.+exp((T-h[2])/h[3])) + h[4]*T*T;
    double const g2 = h[6] + (h[5]-h[6])/(1.+exp((T-h[7])/h[8])) + h[9]*exp(-h[10]*T);
    double const pc = to_bar(driesner07_prost_H2O_NaCl_pc_T(T_K));
    double xc = driesner07_prost_H2O_NaCl_xc_T(T_K);

    // get g0 from xVLH = f(pVLH, T)
    double g0;
    if (T_K < driesner07_NaCl_VLH_T())
    {
        // below NaCl triple temperature
        double p_VLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
        double x_VLH = driesner07_prost_H2O_NaCl_LH_xl_pT(p_VLH, T_K);
        double dp_VLH = pc - to_bar(p_VLH);

        if (T_K < H2O_Tc) {
            // another constrain: x = 0 at p=pH2O
            double pH2O = h2o_sat_p(T_K);
            double dp_H2O = pc - to_bar(pH2O);
            g0 = x_VLH - g1*(dp_VLH - dp_H2O) - g2*(dp_VLH*dp_VLH - dp_H2O*dp_H2O);
            g0 /= (sqrt(dp_VLH) - sqrt(dp_H2O));
            xc = - g0*sqrt(dp_H2O) - g1*dp_H2O - g2*dp_H2O*dp_H2O;
            //printf("g0=%g\n", g0);
        } else {
            g0 = (x_VLH - xc - g1*dp_VLH - g2*dp_VLH*dp_VLH)/sqrt(dp_VLH);
        }
    } else {
        // at and above NaCl triple temperature
        double p_NaCl_boil = driesner07_NaCl_sat_p_T(T_K);
        double x_NaCl_boil = 1.0; //driesner07_prost_H2O_NaCl_LH_xl_pT(p_NaCl_boil, T_K);
        //printf("x_LH,sat=%e\n", driesner07_prost_H2O_NaCl_LH_xl_pT(p_NaCl_boil, T_K));
        double dp_NaCl_boil = pc - to_bar(p_NaCl_boil);
        g0 = (x_NaCl_boil - xc - g1*dp_NaCl_boil - g2*dp_NaCl_boil*dp_NaCl_boil)/sqrt(dp_NaCl_boil);
    }


    // Eqs (11)
    double dp = pc - p;
    //assert(pc >= p);
    double x = xc + g0*sqrt(dp) + g1*dp + g2*dp*dp;


    //printf("g0=%12.9e\n", g0);
    //printf("g1=%12.9e\n", g1);
    //printf("g2=%12.9e\n", g2);
    //printf("dp=%12.9e\n", dp);

    return x;
}


double driesner07_prost_H2O_NaCl_VL_xv_pT(double p, double T)
{
    double x = driesner07_prost_H2O_NaCl_VL_xl_pT(p,T)/driesner07_prost_H2O_NaCl_VL_K(p,T);
    return x;
}

#pragma endregion

#pragma region LH
// H2O(l)-NaCl(s)
double driesner07_prost_H2O_NaCl_LH_xl_pT(double p_Pa, double T_K)
{
    double p = to_bar(p_Pa);
    double T_C = to_C(T_K);
    double e[6];
    e[0] = 0.0989944  + 3.30796e-6*p - 4.71759e-10*p*p;
    e[1] = 0.00947257 - 8.66460e-6*p + 1.69417e-9 *p*p;
    e[2] = 0.610863   - 1.51716e-5*p + 1.19290e-8 *p*p;
    e[3] = -1.64994   + 2.03441e-4*p - 6.46015e-8 *p*p;
    e[4] = 3.36474    - 1.54023e-4*p + 8.17048e-8 *p*p;
    e[5] = 1.0 - e[0] - e[1] - e[2] - e[3] - e[4];
    double const Thm = to_C(driesner07_NaCl_LH_T_p(p_Pa));
    double x = 0;
    for (int i=0; i<6; i++)
        x += e[i]*pow(T_C/Thm, i);

    //assert(0<=x && x<=1.);
    return x;
}
#pragma endregion

#pragma region VH
// H2O(v)-NaCl(s)
double driesner07_prost_H2O_NaCl_VH_xv_pT(double p, double T)
{
    double x_v_sat = driesner07_prost_H2O_NaCl_LH_xl_pT(p, T)
                    / driesner07_prost_H2O_NaCl_VL_K(p,T);
    return x_v_sat;
}
#pragma endregion

#pragma region VLH

// H2O(v)-H2O(l)-NaCl(s)
double driesner07_prost_H2O_NaCl_VLH_p_T(double T_K)
{
    double const T_C = to_C(T_K);
    assert(T_C <= driesner07_NaCl_VLH_TC());
    double const f[11] = {
        4.64e-3,
        5.0e-7,
        1.69078e1,
        -2.69148e2,
        7.63204e3,
        -4.95636e4,
        2.33119e5,
        -5.13556e5,
        5.49708e5,
        -2.84628e5,
        5.754079606e4
    };
    //f[10] = pt_NaCl_bar - (f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9]);
    //printf("f10=%12.9e\n",f[10]);

    double p = 0;
    for (int i=0; i<11; i++)
        p += f[i] * pow(T_C/driesner07_NaCl_VLH_TC(), i);
    assert(p>0);
    return to_Pa(p);
}

double driesner07_prost_H2O_NaCl_VLH_xl_T(double T)
{
    double p = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    double x = driesner07_prost_H2O_NaCl_LH_xl_pT(p ,T);
    return x;
}

double driesner07_prost_H2O_NaCl_VLH_xv_T(double T)
{
    double p = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    double x = driesner07_prost_H2O_NaCl_VH_xv_pT(p ,T);
    return x;
}

#pragma endregion

#pragma region phase

static bool isSuperCritical(double p, double T, double x)
{
    if (p < H2O_pc || T < H2O_Tc)
        return false;
    if (x==0)
        return (p>=H2O_pc && T>=H2O_Tc);
    const double xc = driesner07_prost_H2O_NaCl_xc_T(T);
    const double pc = driesner07_prost_H2O_NaCl_pc_T(T);
    return (p >= pc); // && x >= xc);
}

static bool isVapor(double p, double T, double x)
{
    // check pmax of V region
    if (T >= H2O_Tc && p >= driesner07_prost_H2O_NaCl_pc_T(T))
        return false;
    if (T < H2O_Tc && p > h2o_sat_p(T))
        return false;

    // check xmax of V region
    double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    double xv_max;
    if (p < pVLH)
        return true;
//        xv_max  = driesner07_prost_H2O_NaCl_VH_xv_pT(p, T);
    // p in VL
    xv_max = driesner07_prost_H2O_NaCl_VL_xv_pT(p, T);
    return (x <= xv_max);
}


static bool isLquid(double p, double T, double x)
{
    // check pmax
    if (T >= H2O_Tc && p >= driesner07_prost_H2O_NaCl_pc_T(T))
        return false;
    if (T < H2O_Tc && p > h2o_sat_p(T))
        return true;

    // check pmin
    if (T < driesner07_NaCl_VLH_T())
    {
        double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
        if (p < pVLH)
            return false;
    } else if (p < driesner07_NaCl_VLH_p()) {
        return false;
    }

    // check xmin
    double x_min = driesner07_prost_H2O_NaCl_VL_xl_pT(p, T);
    return (x_min <= x);
}



static bool isTwoPhase(double p, double T, double x)
{
    // check max pressure of VL region
    const double pc = driesner07_prost_H2O_NaCl_pc_T(T);
    if (p >= pc)
        return false;

    if (T < H2O_Tc) {
        const double H2O_ps = h2o_sat_p(T);
        if (p > H2O_ps)
            return false;
    }

    // check min pressure of VL region
    double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    if (p < pVLH)
        return false;

    // checl min xv of VL region
    double VL_x_min = driesner07_prost_H2O_NaCl_VL_xv_pT(p, T);
    double VL_x_max = driesner07_prost_H2O_NaCl_VL_xl_pT(p, T);
    if (VL_x_min < x && x < VL_x_max)
        return true;

    return (p == pVLH && x >= VL_x_max);
}


H2O_NaCl_PhaseType driesner07_prost_H2O_NaCl_phase_type(double p, double T, double x)
{
    if (isLquid(p, T, x))
    {
        double xsat = driesner07_prost_H2O_NaCl_NaCl_solubility_xl_pT(p,T);
        if (x < xsat || (xsat==1.0 && x==1.0))
            return H2O_NaCl_PhaseType_L;
        return H2O_NaCl_PhaseType_LH;
    }

    if (isTwoPhase(p, T, x))
    {
        if (p > driesner07_prost_H2O_NaCl_VLH_p_T(T))
            return H2O_NaCl_PhaseType_VL;
        return H2O_NaCl_PhaseType_VLH;
    }

    if (isVapor(p, T, x)) {
        if (x < driesner07_prost_H2O_NaCl_NaCl_solubility_xv_pT(p,T))
            return H2O_NaCl_PhaseType_V;
        return H2O_NaCl_PhaseType_VH;
    }
    if (isSuperCritical(p, T, x))
    {
        if (x < driesner07_prost_H2O_NaCl_NaCl_solubility_xl_pT(p,T))
            return H2O_NaCl_PhaseType_SC;
        return H2O_NaCl_PhaseType_SCH;
    }
    assert(false);
    return H2O_NaCl_PhaseType_INVALID;
}

char* driesner07_prost_H2O_NaCl_phase_type_to_str(H2O_NaCl_PhaseType phase, char*name)
{
    if (phase == H2O_NaCl_PhaseType_V)
        strcpy(name, "V");
    else if (phase == H2O_NaCl_PhaseType_VH)
        strcpy(name, "VH");
    else if (phase == H2O_NaCl_PhaseType_L)
        strcpy(name, "L");
    else if (phase == H2O_NaCl_PhaseType_LH)
        strcpy(name, "LH");
    else if (phase == H2O_NaCl_PhaseType_VL)
        strcpy(name, "VL");
    else if (phase == H2O_NaCl_PhaseType_VLH)
        strcpy(name, "VLH");
    else if (phase == H2O_NaCl_PhaseType_SC)
        strcpy(name, "F");
    else if (phase == H2O_NaCl_PhaseType_SCH)
        strcpy(name, "FH");
    else
        strcpy(name, "Invalid");

    return name;
}

#pragma endregion

#pragma region vol_prop
// scaled temperature for molar volume
double driesner07_prost_H2O_NaCl_Tv_pTx(double p /*bar*/, double T /*C*/, double x)
{
    double n11 = -54.2958 - 45.7623*exp(-9.44785e-4*p);
    double n1_x1 = 330.47 + 0.942876*sqrt(p) + 0.0817193*p
                   - 2.47556e-8*p*p + 3.45052e-10*p*p*p;
    double n10 = n1_x1;
    double n12 = - n10 - n11;

    double n21 = -2.6142 - 0.000239092*p;
    double n22 = 0.0356828 + 4.37235e-6*p + 2.0566e-9*p*p;
    double n2_x1 = -0.0370751 + 0.00237723*sqrt(p) + 5.42049e-5*p
                   + 5.84709e-9*p*p - 5.99373e-13*p*p*p;
    double n20 = 1. - n21 * sqrt(n22);
    double n23 = n2_x1 - n20 - n21 * sqrt(1.+n22);

    double n1 = n10 + n11*(1.-x) + n12*(1.-x)*(1.-x);
    double n2 = n20 + n21*sqrt(x + n22) + n23*x;

    double Tv = n1 + n2 * T;

    // for low T, low p region
    double n300 = 7.60664e6/pow(p + 472.051, 2);
    double n301 = -50 - 86.1446*exp(-6.21128e-4*p);
    double n302 = 294.318*exp(-5.66735e-3*p);
    double n310 = -0.0732761*exp(-2.3772e-3*p) - 5.2948e-5*p;
    double n311 = -47.2747 + 24.3653*exp(-1.25533e-3*p);
    double n312 = -0.278529 - 0.00081381*p;

    double n30 = n300*(exp(n301*x)-1) + n302*x;
    double n31 = n310*exp(n311*x) + n312*x;
    double D = n30*exp(n31*T);
    Tv += D;

    return Tv; // deg C
}


// Molar volume function for most cases
double driesner07_prost_H2O_NaCl_vm_pTx0(double p_Pa, double T_K, double x)
{
    const double p = to_bar(p_Pa);
    const double T = to_C(T_K);

    double Tv = driesner07_prost_H2O_NaCl_Tv_pTx(p, T, x);
    double Tv_K = to_K(Tv);
    //printf("Tv=%g\n",Tv);
    assert(Tv>0);

    double rho = h2o_rho_pT(p_Pa, to_K(Tv));
    //printf("rho=%g\n", rho);
    double v = to_H2O_molar_volume(rho); // m3/mol
    return v;
}

// 1st workaround for the molar volume at low p (< ca 1.5MPa), low T
// ad highest xNaCl close to halite saturation. see Equation (17)
double driesner07_prost_H2O_NaCl_vm_pTx1(double p_Pa, double T_K, double x)
{
    const double p_bar = to_bar(p_Pa);
    const double T_C = to_C(T_K);
    // get Tsat, rho_l,sat for H2O
    double Tsat, rhol, rhov;
    h2o_sat_Trholv_p(p_Pa, &Tsat, &rhol, &rhov);
    double Tsat_C = to_C(Tsat);
    double d_sat = rhol;
    double Vm_sat = to_H2O_molar_volume(d_sat)*1e6; // cm3/mol

    // compute drho/dT for H2O
    double dT = Tsat*1e-6;
    //printf("p=%e, T=%g\n", p_Pa, Tsat - dT);
    double d1 = h2o_rho_pT(p_Pa, Tsat - dT);
    double dd_dT = (d_sat - d1)/dT;
    double dVm_dT_sat = to_dH2O_molar_volume(d_sat, dd_dT)*1e6;

    double o2 = 2.0125e-7 + 3.29977e-9*exp(-4.31279*log(p_bar))
                -1.17748e-7*log(p_bar) + 7.58009e-8*pow(log(p_bar),2);
    // compute o0, o1 from values at Tsat
    double o1 = dVm_dT_sat - 3*o2*Tsat_C*Tsat_C;
    double o0 = Vm_sat - o1*Tsat_C - o2*Tsat_C*Tsat_C*Tsat_C;
    // compute V at Tv
    double Tv = driesner07_prost_H2O_NaCl_Tv_pTx(p_bar, T_C, x);
    double v = o0 + o1*Tv + o2*Tv*Tv*Tv;
    //double v = o0 + o1*T_C + o2*T_C*T_C*T_C;
    return v* 1e-6; // cm3/mol to m3/mol
}

// 2nd workaround for the molar volume at high T (> ca600C), low p (< ca35 MPa),
// ad high xNaCl (> xVL,liquid). see Equation (18)
double driesner07_prost_H2O_NaCl_vm_pTx2(double p_Pa, double T_K, double x)
{
    const double p = to_bar(p_Pa);
    const double T = to_C(T_K);

    double Tv = driesner07_prost_H2O_NaCl_Tv_pTx(p, T, x);
    double Tv_K = to_K(Tv);
    //printf("Tv=%g\n",Tv);
    assert(Tv>0);

    // >ca600C, <ca350bar, >xVL,liquid
    double p1 = 390.147; //bar
    double p2 = 1000.0;  //bar
    //double H2O_molar_volume1; //390bar
    double dp = p1*1e-6;

    double Vm_p1 = driesner07_prost_H2O_NaCl_vm_pTx0(to_Pa(p1), T_K, x);
    double dVm_p1 = driesner07_prost_H2O_NaCl_vm_pTx0(to_Pa(p1+dp), T_K, x);
    dVm_p1 -= driesner07_prost_H2O_NaCl_vm_pTx0(to_Pa(p1-dp), T_K, x);
    dVm_p1 /= 2.*dp;
    double dVm_p2 = driesner07_prost_H2O_NaCl_vm_pTx0(to_Pa(p2+dp), T_K, x);
    dVm_p2 -= driesner07_prost_H2O_NaCl_vm_pTx0(to_Pa(p2-dp), T_K, x);
    dVm_p2 /= 2.*dp;

    double o4 = (dVm_p1 - dVm_p2);
    o4 /= (1./(p1+1000.)-1./(p2+1000.));
    double o5 = dVm_p1 - o4/(p1+1000.);
    double o3 = Vm_p1 - o4*log(p1+1000.) - o5*p1;
    double v = o3 + o4*log(p+1.e3) + o5*p;

    return v;
}

double driesner07_prost_H2O_NaCl_vm_singlephase_pTx(double p_Pa, double T_K, double x)
{
    const double p_bar = to_bar(p_Pa);
    const double T_C = to_C(T_K);

    const double Tv_C = driesner07_prost_H2O_NaCl_Tv_pTx(p_bar, T_C, x);
    assert(Tv_C>0);
    const double Tv_K = to_K(Tv_C);
    //printf("Tv=%g\n",Tv);

    H2O_NaCl_PhaseType phase_id = driesner07_prost_H2O_NaCl_phase_type(p_Pa, T_K,  x);
    double v = -1;
    if (phase_id/10 == H2O_NaCl_PhaseType_L/10 && p_Pa < 20e5)
    {
        if (Tv_K > h2o_sat_T(p_Pa)) // in case Tv exceeds boiling temperature
            v = driesner07_prost_H2O_NaCl_vm_pTx1(p_Pa, T_K, x);
    } else if (phase_id/10 == H2O_NaCl_PhaseType_L/10 && T_C > 600 && p_Pa < 40e6) {
        // at high T (> ca600C), low p (< ca35 MPa), high xNaCl (> xVL,liquid)
        double VL_xl = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
        if (x >= VL_xl)
            v = driesner07_prost_H2O_NaCl_vm_pTx2(p_Pa, T_K, x);
    }

    // default case
    if (v < 0)
        v = driesner07_prost_H2O_NaCl_vm_pTx0(p_Pa, T_K, x);

    return v;
}

double driesner07_prost_H2O_NaCl_rho_singlephase_pTx(double p_Pa, double T_K, double x)
{
    double mv = driesner07_prost_H2O_NaCl_vm_singlephase_pTx(p_Pa, T_K, x); //m3/mol
    //double total_mass = x*NaCl_molar_mass + (1 - x)*H2O_molar_mass; // kg/total mol;
    //double rho = total_mass * 1. / mv;
    double total_mol = 1./(x*NaCl_molar_mass+(1-x)*H2O_molar_mass); // mol/kg
    double rho = 1./total_mol * 1./mv;
    return rho;
}

double driesner07_prost_H2O_NaCl_VL_rhov_pT(double p_Pa, double T_K)
{
    double xv = driesner07_prost_H2O_NaCl_VL_xv_pT(p_Pa, T_K);
    double mv = driesner07_prost_H2O_NaCl_vm_singlephase_pTx(p_Pa, T_K, xv); //m3/mol
    double total_mol = 1./(xv*NaCl_molar_mass+(1.-xv)*H2O_molar_mass); // mol/kg
    double rho = 1./total_mol * 1./mv;
    return rho;
}

double driesner07_prost_H2O_NaCl_VL_rhol_pT(double p_Pa, double T_K)
{
    double x = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
    double molar_v = driesner07_prost_H2O_NaCl_vm_singlephase_pTx(p_Pa, T_K, x); //m3/mol
    double total_mol = 1./(x*NaCl_molar_mass+(1.-x)*H2O_molar_mass); // mol/kg
    double rho = 1./total_mol * 1./molar_v;
    return rho;
}

double driesner07_prost_H2O_NaCl_rhol_pTx1(double p_Pa, double T_K, double x)
{
    double mv = driesner07_prost_H2O_NaCl_vm_pTx1(p_Pa, T_K, x); //m3/mol
    double total_mol = 1./(x*NaCl_molar_mass+(1-x)*H2O_molar_mass); // mol/kg
    double rho = 1./total_mol * 1./mv;
    return rho;
}

double driesner07_prost_H2O_NaCl_rhol_pTx2(double p_Pa, double T_K, double x)
{
    double mv = driesner07_prost_H2O_NaCl_vm_pTx2(p_Pa, T_K, x); //m3/mol
    double total_mol = 1./(x*NaCl_molar_mass+(1-x)*H2O_molar_mass); // mol/kg
    double rho = 1./total_mol * 1./mv;
    return rho;
}

#pragma endregion


#pragma region solubility
double driesner07_prost_H2O_NaCl_NaCl_solubility_xv_pT(double p, double T)
{
    double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    if (p<=pVLH)
        return driesner07_prost_H2O_NaCl_VH_xv_pT(p, T);
    // two phase (TODO: check if pT is in liquid phase)
    return 1.0; // no saturation limit
}

double driesner07_prost_H2O_NaCl_NaCl_solubility_xl_pT(double p, double T)
{
    if (T > driesner07_NaCl_VLH_T())
        return 1.0; // NaCl melt
    if (T >= driesner07_NaCl_LH_T_p(p))
        return 1.0; // NaCl melt
    double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    if (p < pVLH)
        return -1.; // vapor phase TODO assert?
    return driesner07_prost_H2O_NaCl_LH_xl_pT(p, T);
}

double driesner07_prost_H2O_NaCl_NaCl_solubility_x_pT(double p, double T)
{
    if (T > driesner07_NaCl_VLH_T())
        return 1.;

    double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T);
    if (p >= pVLH)
        return driesner07_prost_H2O_NaCl_LH_xl_pT(p, T);
    return driesner07_prost_H2O_NaCl_VH_xv_pT(p, T);
}

#pragma endregion


#pragma region enthalpy

static double H2O_NaCl_q2(double p /*bar*/, double T /*C*/, double x)
{
    double q21 = -1.69513 - 4.52781e-4*p - 6.04279e-8*p*p;
    double q22 = 0.0612567 + 1.88082e-5*p;
    double q2_x1 = 0.241022 + 3.45087e-5*p - 4.28356e-9*p*p;

    double q20 = 1. - q21 * sqrt(q22);
    double q23 = q2_x1 - q20 - q21 * sqrt(1+q22);

    double q2 = q20 + q21*sqrt(x+q22)+q23*x;
    return q2;
}

// scaled temperature for enthalpy
static double H2O_NaCl_Th(double p /*bar*/, double TC /*C*/, double x)
{
    double q11 = -32.1724 + 0.0621255*p;
    double q1_x1 = 47.9048 - 9.36994e-3*p + 6.51059e-6*p*p;
    double q10 = q1_x1;
    double q12 = - q10 - q11;
    double q1 = q10 + q11*(1-x)+q12*(1-x)*(1-x);

    double q2 = H2O_NaCl_q2(p, TC, x);

    double Th = q1 + q2 * TC;
    return Th;
}

// J/kg
double driesner07_prost_H2O_NaCl_h_singlephase_pTx(double p_Pa, double T_K, double x)
{
    const double p = to_bar(p_Pa);
    const double T = to_C(T_K);
    double Th = H2O_NaCl_Th(p, T, x);
    return h2o_h_pT(p_Pa, to_K(Th));
}

double driesner07_prost_H2O_NaCl_cp_singlephase_pTx(double p_Pa, double T_K, double x)
{
    const double p = to_bar(p_Pa);
    const double T = to_C(T_K);
    double q2 = H2O_NaCl_q2(p, T, x);
    double Th = H2O_NaCl_Th(p, T, x);
    double cp1 = h2o_cp_pT(p_Pa, to_K(Th));
    return q2*cp1;
}

#pragma endregion

#pragma region saturation

double driesner07_prost_H2O_NaCl_mass_fraction_liquid(double p_Pa, double T_K, double x)
{
    if (isTwoPhase(p_Pa, T_K, x))
    {
        double xl = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
        double xv = driesner07_prost_H2O_NaCl_VL_xv_pT(p_Pa, T_K);
        double ml = (x-xv)/(xl-xv);
        return ml;
    } else if (isLquid(p_Pa, T_K, x)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double driesner07_prost_H2O_NaCl_mass_fraction_vapor(double p_Pa, double T_K, double x)
{
    double ml = driesner07_prost_H2O_NaCl_mass_fraction_liquid(p_Pa, T_K, x);
    return (1.-ml);
}

double driesner07_prost_H2O_NaCl_volume_fraction_liquid(double p_Pa, double T_K, double x)
{
    if (isTwoPhase(p_Pa, T_K, x))
    {
        double xl = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
        double xv = driesner07_prost_H2O_NaCl_VL_xv_pT(p_Pa, T_K);
        double rhol = driesner07_prost_H2O_NaCl_VL_rhol_pT(p_Pa, T_K);
        double rhov = driesner07_prost_H2O_NaCl_VL_rhov_pT(p_Pa, T_K);
        double satl = rhov*(xv-x)/(rhov*(xv-x)+rhol*(x-xl));
        return satl;
    } else if (isLquid(p_Pa, T_K, x)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double driesner07_prost_H2O_NaCl_volume_fraction_vapor(double p_Pa, double T_K, double x)
{
    double satl = driesner07_prost_H2O_NaCl_volume_fraction_liquid(p_Pa, T_K, x);
    return (1.-satl);
}

double driesner07_prost_H2O_NaCl_volume_fraction_halite(double p_Pa, double T_K, double x)
{
    //assert is V+H or L+H
    H2O_NaCl_PhaseType phase_type = driesner07_prost_H2O_NaCl_phase_type(p_Pa, T_K, x);
    assert(phase_type==H2O_NaCl_PhaseType_VH || phase_type==H2O_NaCl_PhaseType_LH);

    // V+H
    if (phase_type==H2O_NaCl_PhaseType_VH)
    {
        double xv = driesner07_prost_H2O_NaCl_VL_xv_pT(p_Pa, T_K);
        double xh = 1.;
        double rhov = driesner07_prost_H2O_NaCl_VL_rhol_pT(p_Pa, T_K);
        double rhoh = driesner07_NaCl_H_rho_pT(p_Pa, T_K);
        double sath = rhov*(xv-x)/(rhoh*(x-xh)+rhov*(xv-x));
        return sath;
    }

    // L+H
    double xl = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
    double xh = 1.;
    double rhol = driesner07_prost_H2O_NaCl_VL_rhov_pT(p_Pa, T_K);
    double rhoh = driesner07_NaCl_H_rho_pT(p_Pa, T_K);
    double sath = rhol*(xl-x)/(rhoh*(x-xh)+rhol*(xl-x));
    return sath;
}

double driesner07_prost_H2O_NaCl_T_rhopx(double rho, double p_Pa, double x, double T0)
{
    double T_K = T0;
    const double rtol = 1e-6;
    const int nitrmax = 50;
    for (int i=0; i<nitrmax; i++)
    {
        double rho_tmp = driesner07_prost_H2O_NaCl_rho_pTx(p_Pa, T_K, x);
        double r_rho = rho_tmp - rho;
#ifdef DEBUG_INVERSE
        printf("%d: r_rho=%.3e, rho=%g, T=%g\n", i, fabs(r_rho/rho), rho_tmp, T_K-273.15);
#endif
        if (fabs(r_rho)<rtol*rho)
        {
            return T_K;
        }

        double dT = 1e-3;
        double drhodT = (driesner07_prost_H2O_NaCl_rho_pTx(p_Pa, T_K+dT, x) - rho_tmp)/dT;
        double deltaT = - r_rho / drhodT;
        T_K += deltaT;
    }
#ifdef DEBUG_INVERSE
    printf("%s: not convered for rho=%g, pMPa=%g, x=%g\n", __FUNCTION__, rho, p_Pa*1e-6, x);
#endif
    return -1;
}

double driesner07_prost_H2O_NaCl_volume_fraction_liquid_rhopx(double rho, double p_Pa, double x, double T0)
{
    double T_K = driesner07_prost_H2O_NaCl_T_rhopx(rho, p_Pa, x, T0);
    double satl = driesner07_prost_H2O_NaCl_volume_fraction_liquid(p_Pa, T_K, x);
    return satl;
}

double driesner07_prost_H2O_NaCl_T_phx(double p_Pa, double h, double x, double T0)
{
    double T_K = T0;
    const double rtol = 1e-6;
    const int nitrmax = 50;
    for (int i=0; i<nitrmax; i++)
    {
        double h_tmp = driesner07_prost_H2O_NaCl_h_pTx(p_Pa, T_K, x);
        double r_h = h_tmp - h;
#ifdef DEBUG_INVERSE
        printf("%d: r=%.3e, h=%g, T=%g\n", i, fabs(r_h/h), h_tmp*1e-3, T_K-273.15);
#endif
        if (fabs(r_h)<rtol*h)
        {
            return T_K;
        }

        double dT = 1e-3;
        double drdT = (driesner07_prost_H2O_NaCl_h_pTx(p_Pa, T_K+dT, x) - h_tmp)/dT;
        double deltaT = - r_h / drdT;
        T_K += deltaT;
    }
#ifdef DEBUG_INVERSE
    printf("%s: not convered for h=%g, p=%g, x=%g\n", __FUNCTION__, h*1e-3, p_Pa*1e-6, x);
#endif
    return -1;
}

double driesner07_prost_H2O_NaCl_rho_pTx(double p_Pa, double T_K, double x)
{
    const H2O_NaCl_PhaseType phase_type = driesner07_prost_H2O_NaCl_phase_type(p_Pa, T_K, x);

    if (phase_type == H2O_NaCl_PhaseType_V
        || phase_type == H2O_NaCl_PhaseType_L
        || phase_type == H2O_NaCl_PhaseType_SC)
    {
        return driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p_Pa, T_K, x);
    }
    else if (phase_type == H2O_NaCl_PhaseType_VL
             || phase_type == H2O_NaCl_PhaseType_VLH)
    {
        double rhol = driesner07_prost_H2O_NaCl_VL_rhol_pT(p_Pa, T_K);
        double rhov = driesner07_prost_H2O_NaCl_VL_rhov_pT(p_Pa, T_K);
        double satl = driesner07_prost_H2O_NaCl_volume_fraction_liquid(p_Pa, T_K, x);
        double rho = satl*rhol + (1.-satl)*rhov;
        return rho;
    } else if (phase_type == H2O_NaCl_PhaseType_VH) {
        double xv = driesner07_prost_H2O_NaCl_VH_xv_pT(p_Pa, T_K);
        double rhov = driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p_Pa, T_K, xv);
        return rhov; //TODO halite precipitation
    } else if (phase_type == H2O_NaCl_PhaseType_LH
               || phase_type == H2O_NaCl_PhaseType_SCH) {
        double xl = driesner07_prost_H2O_NaCl_LH_xl_pT(p_Pa, T_K);
        double rhol = driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p_Pa, T_K, xl);
        return rhol; //TODO halite precipitation
    }
    assert(false);
    return -1;
}

double driesner07_prost_H2O_NaCl_h_pTx(double p_Pa, double T_K, double x)
{
    const H2O_NaCl_PhaseType phase_type = driesner07_prost_H2O_NaCl_phase_type(p_Pa, T_K, x);
    if (phase_type == H2O_NaCl_PhaseType_V || phase_type == H2O_NaCl_PhaseType_L)
    {
        return driesner07_prost_H2O_NaCl_h_singlephase_pTx(p_Pa, T_K, x);
    }
    else if (phase_type == H2O_NaCl_PhaseType_VL || phase_type == H2O_NaCl_PhaseType_VLH)
    {
        double xl = driesner07_prost_H2O_NaCl_VL_xl_pT(p_Pa, T_K);
        double rhol = driesner07_prost_H2O_NaCl_VL_rhol_pT(p_Pa, T_K);
        double hl = driesner07_prost_H2O_NaCl_h_pTx(p_Pa, T_K, xl);
        double xv = driesner07_prost_H2O_NaCl_VL_xv_pT(p_Pa, T_K);
        double rhov = driesner07_prost_H2O_NaCl_VL_rhov_pT(p_Pa, T_K);
        double hv = driesner07_prost_H2O_NaCl_h_pTx(p_Pa, T_K, xv);
        double satl = driesner07_prost_H2O_NaCl_volume_fraction_liquid(p_Pa, T_K, x);
        double rho = satl*rhol + (1.-satl)*rhov;
        double h = (satl*rhol*hl + (1.-satl)*rhov*hv)/rho;
        return h;
    } else if (phase_type == H2O_NaCl_PhaseType_VH) {
        double xv = driesner07_prost_H2O_NaCl_VH_xv_pT(p_Pa, T_K);
        double rhov = driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p_Pa, T_K, xv);
        double hv = driesner07_prost_H2O_NaCl_h_singlephase_pTx(p_Pa, T_K, xv);
        return hv; //TODO halite precipitation
    } else if (phase_type == H2O_NaCl_PhaseType_LH) {
        double xl = driesner07_prost_H2O_NaCl_LH_xl_pT(p_Pa, T_K);
        double rhol = driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p_Pa, T_K, xl);
        double hl = driesner07_prost_H2O_NaCl_h_singlephase_pTx(p_Pa, T_K, xl);
        return hl; //TODO halite precipitation
    }
    assert(false);
    return -1;

}

double driesner07_prost_H2O_NaCl_VL_p_Txv(double T_K, double x, double p0)
{
    double p = p0 > 0 ? p0 : 1e7;
    double pVLH = 0.1;
    if (T_K < driesner07_NaCl_VLH_T())
        pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
    double pc = driesner07_prost_H2O_NaCl_pc_T(T_K);
    p = max(pVLH, min(p, pc));
    const double rtol = 1e-6;
    for (int i=0; i<100; i++)
    {
        double x_tmp = driesner07_prost_H2O_NaCl_VL_xv_pT(p, T_K);
        double r = x_tmp - x;
        double error = fabs( r / x );
        //printf("%d: pMPa=%g, x=%g, r=%.3e\n", i, p*1e-6, x_tmp, error);
        if (error < rtol)
            return p;
        double dp = p * 1e-5;
        double drdp = 0;
        if (p+dp > pc)
            drdp = (x_tmp - driesner07_prost_H2O_NaCl_VL_xv_pT(p-dp, T_K))/dp;
        else
            drdp = (driesner07_prost_H2O_NaCl_VL_xv_pT(p+dp, T_K) - x_tmp)/dp;
        double deltap = -r / drdp;
        double factor = 1.;
        if (deltap + p > pc)
            factor = 0.2;
        else if (deltap + p < pVLH)
            factor = 0.2;
        p += deltap * factor;
        if (p < pVLH) p = pVLH;
        if (p > pc) p = pc*0.999;
    }
    //printf("***error in %s: diverged at T=%g, x=%g\n", __FUNCTION__, T_K, x);
    return -1;
}

int driesner07_prost_H2O_NaCl_VL_p_all_Txv(double T_K, double x, double pmin, double pmax, double* p)
{
    double p1 = driesner07_prost_H2O_NaCl_VL_p_Txv(T_K, x, pmin*1.01);
    double p2 = driesner07_prost_H2O_NaCl_VL_p_Txv(T_K, x, pmax*0.99);
    if (p1 < 0 && p2 < 0) {
        p1 = driesner07_prost_H2O_NaCl_VL_p_Txv(T_K, x, (pmin+pmax)*0.5);
        if (p1 < 0) return 0;
    }
    if (p2 < 0) p2 = p1;
    if (p1 < 0) p1 = p2;

    if (fabs(p1-p2)<p1*1e-6)
    {
        p[0] = p1;
        return 1;
    }
    p[0] = p1;
    p[1] = p2;
    return 2;
}

double driesner07_prost_H2O_NaCl_VL_p_Txl(double T_K, double x, double p0)
{
    double p = p0 > 0 ? p0 : 1e7;
    double pVLH = -1;
    if (T_K < driesner07_NaCl_VLH_T())
        pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
    double pc = driesner07_prost_H2O_NaCl_pc_T(T_K);
    p = max(pVLH, min(p, pc));
    const double rtol = 1e-6;
    for (int i=0; i<50; i++)
    {
        double x_tmp = driesner07_prost_H2O_NaCl_VL_xl_pT(p, T_K);
        double r = x_tmp - x;
        //printf("%d: pMPa=%g, x=%g, r=%.3e, error=%.3e\n", i, p*1e-6, x_tmp, r, fabs(r/x));
        if (fabs(r) < rtol*x)
            return p;
        double dp = p * 1e-5;
        double drdp = 0;
        if (p+dp > pc)
            drdp = (x_tmp - driesner07_prost_H2O_NaCl_VL_xl_pT(p-dp, T_K))/dp;
        else
            drdp = (driesner07_prost_H2O_NaCl_VL_xl_pT(p+dp, T_K) - x_tmp)/dp;
        double deltap = -r / drdp;
        p += deltap;
        if (p < pVLH) p = pVLH;
        if (p > pc) p = pc;
    }
    //printf("***error in %s: diverged at T=%g, x=%g\n", __FUNCTION__, T_K, x);
    return -1;
}

int driesner07_prost_H2O_NaCl_sat_p_Tx(double T_K, double x, double p0, double* p)
{
    int n = 0; // the number of saturation pressure at the given condition
    // pure water
    if (x==0)
    {
        if (T_K > H2O_Tc)
            return 0;
        p[n++] = h2o_sat_p(T_K);
        return n;
    }
    // above melting point of NaCl
    if (T_K >= driesner07_NaCl_VLH_T())
    {
        double pv[2];
        int nv = driesner07_prost_H2O_NaCl_VL_p_all_Txv(T_K, x, 0.1, 250e6, pv);
        for (int i=0; i<nv; i++)
            p[n++] = pv[i];
        double pl = driesner07_prost_H2O_NaCl_VL_p_Txl(T_K, x, p0);
        if (pl > 0)
            p[n++] = pl;
        return n;
    }
    // below melting point of NaCl
    const double xlVLH = driesner07_prost_H2O_NaCl_VLH_xl_T(T_K);
    const double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
    if (x >= xlVLH)
    {
        p[n++] = pVLH;
        return n;
    }

    const double xc = driesner07_prost_H2O_NaCl_xc_T(T_K);
    const double xvVLH = driesner07_prost_H2O_NaCl_VLH_xv_T(T_K);
    const double pc = driesner07_prost_H2O_NaCl_pc_T(T_K);
    if (x < xvVLH)
    {
        double pv[2];
        int nv = driesner07_prost_H2O_NaCl_VL_p_all_Txv(T_K, x, pVLH, pc, pv);
        if (xc < x)
        {
            assert(nv==1);
            p[n++] = pv[0];
            p[n++] = driesner07_prost_H2O_NaCl_VL_p_Txl(T_K, x, p0);
            return n;
        }
        for (int i=0; i<nv; i++)
            p[n++] = pv[i];
        return n;
    }

    if (x > xvVLH && xc < xvVLH)
    {
        double pl = driesner07_prost_H2O_NaCl_VL_p_Txl(T_K, x, p0);
        double pv[2] = {0, 0};
        int nv = driesner07_prost_H2O_NaCl_VL_p_all_Txv(T_K, x, pVLH, pc, pv);

        if (nv==0)
        {
            p[n++] = pVLH;
            p[n++] = pl;
        } else if (nv == 2) {
            p[n++] = pVLH;
            p[n++] = pv[0];
            p[n++] = pv[1];
            p[n++] = pl;
        } else {
            // printf("nv=%d\n", nv);
            // assert(nv==0 || nv==2);
            return 0;
        }
        return n;
    }

    if (x > xvVLH && xc > xvVLH)
    {
        double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
        double pl = driesner07_prost_H2O_NaCl_VL_p_Txl(T_K, x, p0);
        double pv[2] = {0, 0};
        int nv = driesner07_prost_H2O_NaCl_VL_p_all_Txv(T_K, x, pVLH, pc, pv);
        if (nv==0)
        {
            p[n++] = pVLH;
            p[n++] = pl;
        } else if (nv == 1) {
            p[n++] = pVLH;
            p[n++] = pv[0];
        } else if (nv == 2) {
            p[n++] = pVLH;
            p[n++] = pv[0];
            p[n++] = pv[1];
            p[n++] = pl;
        } else {
            assert(nv==0 || nv==2);
        }
        return n;
    }
    assert(false);
    return -1;
}

// // return pressure at V-VL boundary or L-VL boundary
// double driesner07_prost_H2O_NaCl_sat_p_upper_Tx(double T_K, double x, double p0)
// {
//     if (T_K > driesner07_NaCl_VLH_T())
//         return -1;
//     if (x==0)
//         return h2o_sat_p(T_K);

//     // p at V-VL boundary and L-VL boundary
//     double pv = driesner07_prost_H2O_NaCl_VL_p_Txv(T_K, x, p0);
//     double pl = driesner07_prost_H2O_NaCl_VL_p_Txl(T_K, x, p0);
//     // if (pv > 0 && pl > 0)
//     //     return min(pv, pl);
//     return max(pv, pl); // return non-negative
// }

// // return pressure at VH-VL boundary
// double driesner07_prost_H2O_NaCl_sat_p_lower_Tx(double T_K, double x, double p0)
// {
//     if (T_K > driesner07_NaCl_VLH_T())
//         return -1;
//     if (x==0)
//         return h2o_sat_p(T_K);

//     double xvVLH = driesner07_prost_H2O_NaCl_VLH_xv_T(T_K);

//     if (x >= xvVLH)
//     {
//         double pVLH = driesner07_prost_H2O_NaCl_VLH_p_T(T_K);
//         return pVLH;
//     }
//     double xc = driesner07_prost_H2O_NaCl_xc_T(T_K);
//     if (x >= xc)
//     {
//         double pv = driesner07_prost_H2O_NaCl_VL_p_Txv(T_K, x, p0);
//         return pv;
//     }
//     return -1; // NaCl concentration is too low at the given temperature
// }

#pragma endregion
