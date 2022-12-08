
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utest.h"

#include "model/iapws/IAPWS-IF97.h"
#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-SAT-92.h"
#include "model/iapws/IAPWS-Melt-11.h"
#include "model/iapws/IAPWS-Viscosity-08.h"
#include "model/iapws/IAPWS-ThermalConductivity-11.h"
#include "model/iapws/IAPS-84.h"
#include "model/driesner07/Driesner2007_NaCl.h"
#include "model/driesner07/Driesner2007_H2ONaCl.h"

#include "util/utility.h"


double to_mole_fraction(double mass_frac)
{
    double mNaCl = mass_frac / 58.443e-3;
    double mH2O = (1 - mass_frac) / 18.015e-3;
    double x = mNaCl / (mNaCl + mH2O);
    return x;
}

double molality_to_mole_fraction(double molality)
{
    double mNaCl = molality;
    double mH2O = 1./18.015e-3;
    double x = mNaCl / (mNaCl + mH2O);
    return x;
}

double to_mass_frac(double mole_frac)
{
    double W_NaCl = mole_frac * 58.443e-3;
    double W_H2O = (1 - mole_frac) * 18.015e-3;
    double w = W_NaCl / (W_NaCl + W_H2O);
    return w;
}

UTEST(driesner07, phase_v)
{
    double p, T, x;
    T = to_K(300); p = 30.e5; x = 0;
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));
    x = to_mole_fraction(1e-10*1e-2);
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));
    p = 70.e5;
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(500); p = 200.e5; x = 0;
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));
    x = to_mole_fraction(1e-3*1e-2);
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));
    p = 500.e5;
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(800); p = 1.e5; x = to_mole_fraction(1e-1*1e-2);
    EXPECT_EQ(H2O_NaCl_PhaseType_V, driesner07_H2O_NaCl_phase_type(p, T, x));
}

UTEST(driesner07, phase_l)
{
    double p, T, x;
    T = to_K(300); p = 100.e5; x = 0;
    EXPECT_EQ(H2O_NaCl_PhaseType_L, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(500); p = 500.e5; x = to_mole_fraction(0.5);
    EXPECT_EQ(H2O_NaCl_PhaseType_L, driesner07_H2O_NaCl_phase_type(p, T, x));
}

UTEST(driesner07, phase_lh)
{
    double p, T, x;
    T = to_K(300); p = 100.e5; x = to_mole_fraction(0.8);
    EXPECT_EQ(H2O_NaCl_PhaseType_LH, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(500); p = 500.e5; x = to_mole_fraction(0.8);
    EXPECT_EQ(H2O_NaCl_PhaseType_LH, driesner07_H2O_NaCl_phase_type(p, T, x));
}

UTEST(driesner07, phase_vl)
{
    double p, T, x;
    T = to_K(300); p = 70.e5; x = to_mole_fraction(0.1);
    EXPECT_EQ(H2O_NaCl_PhaseType_VL, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(500); p = 500.e5; x = to_mole_fraction(0.2);
    EXPECT_EQ(H2O_NaCl_PhaseType_VL, driesner07_H2O_NaCl_phase_type(p, T, x));
}

UTEST(driesner07, phase_sc)
{
    double p, T, x;
    T = to_K(375); p = 230.e5; x = to_mole_fraction(0.);
    EXPECT_EQ(H2O_NaCl_PhaseType_SC, driesner07_H2O_NaCl_phase_type(p, T, x));

    T = to_K(500); p = 600.e5; x = to_mole_fraction(0.1);
    EXPECT_EQ(H2O_NaCl_PhaseType_SC, driesner07_H2O_NaCl_phase_type(p, T, x));
    //TODO above pc,Tc but not in p-x view
    x = to_mole_fraction(0.01);
    EXPECT_EQ(H2O_NaCl_PhaseType_SC, driesner07_H2O_NaCl_phase_type(p, T, x));
}

UTEST(driesner07, vol1)
{
    FILE* fp = fopen("vol_sat.csv", "w");
    fprintf(fp, "T, v, v0, v2, vwater, Tv, xsat\n");
    // double p = 200e5;
    // int nT = 100;
    // for (int i=0; i<nT; i++)
    // {
    //     double T = to_K(100. + (1000.-100.)/nT*i);
    //     double xsat = driesner07_H2O_NaCl_NaCl_solubility_x_pT(p, T);
    //     double v = driesner07_H2O_NaCl_vm_singlephase_pTx(p, T, xsat);
    //     double v0 = driesner07_H2O_NaCl_vm_pTx0(p, T, xsat);
    //     double v2 = driesner07_H2O_NaCl_vm_pTx2(p, T, xsat);
    //     double Tv = driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T), xsat);
    //     double vwater = driesner07_H2O_vm_pT(p, T);
    //     fprintf(fp, "%g, %g, %g, %g, %g, %g, %g\n", to_C(T), v*1e6, v0*1e6, v2*1e6, vwater*1e6, Tv, to_mass_frac(xsat));
    // }
    double p = 5e5;
    int nT = 100;
    for (int i=0; i<nT; i++)
    {
        double T = to_K(130. + (170.-130.)/nT*i);
        double xsat = driesner07_H2O_NaCl_NaCl_solubility_x_pT(p, T);
        double v = driesner07_H2O_NaCl_vm_singlephase_pTx(p, T, xsat);
        double v0 = driesner07_H2O_NaCl_vm_pTx0(p, T, xsat);
        double v1 = driesner07_H2O_NaCl_vm_pTx1(p, T, xsat);
        double Tv = driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T), xsat);
        double vwater = driesner07_H2O_vm_pT(p, T);
        fprintf(fp, "%g, %g, %g, %g, %g, %g, %g\n", to_C(T), v*1e6, v0*1e6, v1*1e6, vwater*1e6, Tv, to_mass_frac(xsat));
    }
    fclose(fp);
}

#if 0
UTEST(driesner07, write_phase)
{
    double Tmin = 20, Tmax = 1000, dT  = 20; //(Tmax-Tmin)/nT;
    int nT = (Tmax-Tmin)/dT;
    // int nx = 40;
    // double xmin = -11, xmax = 2.0, dx  = (xmax-xmin)/nx;
    int nx = 20;
    double xmin = 0, xmax = 1.0, dx  = (xmax-xmin)/nx;
    double xscale = 50., xs = 1000;
    int np = 400;
    double pmin = 1, dp  = 2000./np, pmax=3e3;

    FILE* fp = fopen("boundary.csv", "w");
    fprintf(fp, "pbar, TC, x, id\n");
    int id = 0;
    // VLH
    {
        for (int i=0; i<=nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            if (T > driesner07_NaCl_VLH_T())
                continue;
            double VLH_p = driesner07_H2O_NaCl_VLH_p_T(T);
            double VLH_xv = driesner07_H2O_NaCl_VH_xv_pT(VLH_p, T);
            double p = VLH_p, x = VLH_xv;
            fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            for (int j=0; j<=nx; j++)
            {
                double x = to_mole_fraction(xmin + dx*j);
                if (x < VLH_xv) continue;
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            }
            // for (int j=0; j<np; j++)
            // {
            //     double p = (pmin + dp*j)*1e5;
            //     if (p > VLH_p) continue;
            //     double x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            //     fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            // }
        }
    }
    id++;
    // Melt
    {
        int np = 10; double dp = 2000/np;
        for (int i=0; i<np; i++)
        {
            double p = (pmin + dp*i)*1e5;
            double T = driesner07_NaCl_LH_T_p(p);
            for (int j=0; j<nx; j++)
            {
                double x = to_mole_fraction(xmin + dx*j);
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            }

        }
    }
    id++;
    // C.P.
    {
        for (int i=0; i<nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            if (T<iapws95_get_Tc())
                continue;
            double p = driesner07_H2O_NaCl_pc_T(T);
            if (p*1e-5>pmax) continue;
            double x = driesner07_H2O_NaCl_xc_T(T);
            fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
        }
    }
    id++;
    // LH
    {
        int np = 20; double dp = 2000/np;
        for (int i=0; i<=nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            double pVLH = (T > driesner07_NaCl_VLH_T()) ? 0 : driesner07_H2O_NaCl_VLH_p_T(T);
            for (int j=0; j<=np; j++)
            {
                double p = (pmin + dp*j)*1e5;
                if (p<pVLH) continue;
                double Tm = driesner07_NaCl_LH_T_p(p);
                if (T>Tm) continue;
                double x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
           }
        }
    }
    id++;
    // VH
    {
        for (int i=0; i<=nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            if (T > driesner07_NaCl_VLH_T())
                continue;
            double VLH_p = driesner07_H2O_NaCl_VLH_p_T(T);
            for (int j=0; j<np; j++)
            {
                double p = (pmin + dp*j)*1e5;
                if (p > VLH_p) continue;
                double x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            }
            // double pVLH = driesner07_H2O_NaCl_VLH_p_T(T);
            // for (int j=0; j<=np; j++)
            // {
            //     double p = (pmin + dp*j)*1e5;
            //     if (p > pVLH) continue;
            //     double x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            //     fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            // }
        }
    }
    id++;
    // VL - xv
    {
        for (int i=0; i<=nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            // if (T > driesner07_NaCl_VLH_T())
            //     continue;
            double pc = driesner07_H2O_NaCl_pc_T(T);
            double pVLH = (T > driesner07_NaCl_VLH_T()) ? 0 : driesner07_H2O_NaCl_VLH_p_T(T);
            for (int j=0; j<=np; j++)
            {
                double p = (pmin + dp*j)*1e5;
                if (p < pVLH || p > pc) continue;
                if (p < if97_pc()*1e6 && p > if97_sat_p_T(T)*1e6)
                    continue;
                double x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
                if (x < 0) continue;
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            }
        }
    }
    id++;
    // VL - xl
    {
        for (int i=0; i<=nT; i++)
        {
            double T = min(Tmax+273.15, Tmin + dT*i + 273.15);
            // if (T > driesner07_NaCl_VLH_T())
            //     continue;
            double pc = driesner07_H2O_NaCl_pc_T(T);
            double pVLH = (T > driesner07_NaCl_VLH_T()) ? 0 : driesner07_H2O_NaCl_VLH_p_T(T);
            for (int j=0; j<=np; j++)
            {
                double p = (pmin + dp*j)*1e5;
                if (p < pVLH || p > pc) continue;
                double x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
                if (x < 0) continue;
                fprintf(fp, "%g, %g, %g, %d\n", p*1e-5, to_C(T), to_mass_frac(x), id);
            }
        }
    }
    fclose(fp);
}
#endif

UTEST(driesner07, NaCl1)
{
    const double rhoc = iapws95_get_rhoc();
    const double Tc = iapws95_get_Tc();

    //---------------------------------------------------------------
    // Driesner&Heinrich 2007: NaCl
    //---------------------------------------------------------------
    {
        printf("-> test driesner07_NaCl_LH_T_p(()\n");
        EXPECT_NEAR3(800.7, driesner07_NaCl_LH_T_p(0.0005e5)-273.15);
        EXPECT_NEAR3(900, driesner07_NaCl_LH_T_p(4016.02e5)-273.15);

        printf("-> test driesner07_NaCl_VH_p_T()\n");
        EXPECT_NEAR3(1.770379e-35, to_bar(driesner07_NaCl_VH_p_T(to_K(5))));
        EXPECT_NEAR3(4.918115e-04, to_bar(driesner07_NaCl_VH_p_T(to_K(800))));

        printf("-> test driesner07_NaCl_VL_p_T()\n");
        EXPECT_NEAR3(5.419375e-04, to_bar(driesner07_NaCl_VL_p_T(to_K(805))));
        EXPECT_NEAR3(2.762695e-03, to_bar(driesner07_NaCl_VL_p_T(to_K(900))));
    }
}

UTEST(driesner07, H2O_Nacl)
{
    //---------------------------------------------------------------
    // Driesner&Heinrich 2007: H2O-NaCl
    //---------------------------------------------------------------
    {
        printf("-> test driesner07_H2O_NaCl_pc_T()\n");
        EXPECT_NEAR3(2.206058e+01, driesner07_H2O_NaCl_pc_T(374+273.15)*1e-6);
        EXPECT_NEAR3(5.810101e+01, driesner07_H2O_NaCl_pc_T(500+273.15)*1e-6);
        EXPECT_NEAR3(9.166756e+01, driesner07_H2O_NaCl_pc_T(600+273.15)*1e-6);

        printf("-> test driesner07_H2O_NaCl_xc_T()\n");
        EXPECT_NEAR3(1.925758e-06, driesner07_H2O_NaCl_xc_T(374+273.15));
        EXPECT_NEAR3(4.570045e-02, driesner07_H2O_NaCl_xc_T(500+273.15));
        EXPECT_NEAR3(1.169455e-01, driesner07_H2O_NaCl_xc_T(800+273.15));

        printf("-> test driesner07_H2O_NaCl_LH_xl_pT()\n");
        EXPECT_NEAR3(1.006877e-01, driesner07_H2O_NaCl_LH_xl_pT(to_Pa(5e2), to_K(10)));
        EXPECT_NEAR3(1.578954e-01, driesner07_H2O_NaCl_LH_xl_pT(to_Pa(5e2), to_K(300)));
        EXPECT_NEAR3(4.685092e-01, driesner07_H2O_NaCl_LH_xl_pT(to_Pa(5e2), to_K(600)));

        printf("-> test driesner07_H2O_NaCl_VH_xv_pT()\n");
        EXPECT_NEAR3(1.931804e-07, driesner07_H2O_NaCl_VH_xv_pT(to_Pa(50), to_K(300)));
        EXPECT_NEAR3(9.746070e-09, driesner07_H2O_NaCl_VH_xv_pT(to_Pa(50), to_K(400)));
        EXPECT_NEAR3(2.505982e-07, driesner07_H2O_NaCl_VH_xv_pT(to_Pa(50), to_K(600)));
        EXPECT_NEAR3(5.028182e-06, driesner07_H2O_NaCl_VH_xv_pT(to_Pa(100), to_K(700)));

        printf("-> test driesner07_H2O_NaCl_VLH_p_T()\n");
        EXPECT_NEAR3(1.761260e+1, driesner07_H2O_NaCl_VLH_p_T(400+273.15)*1e-6);
        EXPECT_NEAR3(3.898861e+1, driesner07_H2O_NaCl_VLH_p_T(600+273.15)*1e-6);

        printf("-> test driesner07_H2O_NaCl_VLH_xl_T()\n");
        EXPECT_NEAR3(2.154955e-01, driesner07_H2O_NaCl_VLH_xl_T(to_K(400)));
        EXPECT_NEAR3(4.703249e-01, driesner07_H2O_NaCl_VLH_xl_T(to_K(600)));

        printf("-> test driesner07_H2O_NaCl_VLH_xv_T()\n");
        EXPECT_NEAR3(3.323791e-05, driesner07_H2O_NaCl_VLH_xv_T(to_K(400)));
        EXPECT_NEAR3(2.572002e-04, driesner07_H2O_NaCl_VLH_xv_T(to_K(600)));

        printf("-> test H2O_NaCl_VL vapor/liquid composition\n");
        double p, T;
        T = to_K(150);
        p = driesner07_H2O_NaCl_VLH_p_T(T); //p = 3.46654e5;
        EXPECT_NEAR3(1.151778e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(7.743987e-11, driesner07_H2O_NaCl_VL_xv_pT(p, T));
        p = to_Pa(4.731357e+00);
        EXPECT_NEAR3(2.560838e-03, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(1.481016e-11, driesner07_H2O_NaCl_VL_xv_pT(p, T));

        T = to_K(300);
        p = to_Pa(6.047119e+01);
        EXPECT_NEAR3(1.573207e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(5.708744e-07, driesner07_H2O_NaCl_VL_xv_pT(p, T));
        p = to_Pa(7.619851e+01);
        EXPECT_NEAR3(6.569883e-02, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(1.046986e-06, driesner07_H2O_NaCl_VL_xv_pT(p, T));

        T = to_K(400);
        p = driesner07_H2O_NaCl_VLH_p_T(T); // p = 1.761260e7;
        EXPECT_NEAR3(2.154955e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(3.323791e-05, driesner07_H2O_NaCl_VL_xv_pT(p, T));
        p = to_Pa(2.600932e+02);
        EXPECT_NEAR3(4.631438e-02, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(5.001350e-04, driesner07_H2O_NaCl_VL_xv_pT(p, T));

        T = to_K(600);
        p = driesner07_H2O_NaCl_VLH_p_T(T); //p = 3.898861e7;
        EXPECT_NEAR3(4.703249e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(2.572002e-04, driesner07_H2O_NaCl_VL_xv_pT(p, T));

        T = to_K(1000);
        p = to_Pa(1.524812e+01);
        EXPECT_NEAR3(9.910193e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(1.191219e-03, driesner07_H2O_NaCl_VL_xv_pT(p, T));
        p = to_Pa(1.057796e+03);
        EXPECT_NEAR3(5.365587e-01, driesner07_H2O_NaCl_VL_xl_pT(p, T));
        EXPECT_NEAR3(7.658204e-03, driesner07_H2O_NaCl_VL_xv_pT(p, T));
    }

//    return 0;
}

UTEST(driesner07, NaCl2)
{
    //---------------------------------------------------------------
    // Driesner 2007: NaCl
    //---------------------------------------------------------------
    {
        double p, T, x;
        printf("-> test driesner07_NaCl_H_rho_pT()\n");
        T = 300;
        p = 0.03e9;
        //printf("T=%g K, p=%g MPa\n", T, p*1e-6);
        EXPECT_NEAR0(1 / 0.4613e-3, driesner07_NaCl_H_rho_pT(p, T), 1e-3);
        T = 400;
        p = 0.2e9;
        //printf("T=%g K, p=%g MPa\n", T, p*1e-6);
        EXPECT_NEAR0(1 / 0.4636e-3, driesner07_NaCl_H_rho_pT(p, T), 1e-3);
        T = 800;
        p = 0.12e9;
        //printf("T=%g K, p=%g MPa\n", T, p*1e-6);
        EXPECT_NEAR0(1/0.4930e-3, driesner07_NaCl_H_rho_pT(p, T), 1e-3);

        //printf("# f(T)\n");
        //T = 300;
        //p = 0.27e9;
        //printf("T=%g K, p=%g MPa, exp=%g, calc=%g\n", T, p*1e-6, 1/0.4568e-3, driesner07_NaCl_H_rho_pT(p, T));
        //T = 400;
        //p = 0.2e9;
        //printf("T=%g K, p=%g MPa, exp=%g, calc=%g\n", T, p*1e-6, 1/0.4636e-3, driesner07_NaCl_H_rho_pT(p, T));
        //T = 500;
        //p = 0.27e9;
        //printf("T=%g K, p=%g MPa, exp=%g, calc=%g\n", T, p*1e-6, 1/0.4681e-3, driesner07_NaCl_H_rho_pT(p, T));


        printf("-> test driesner07_NaCl_L_rho_pT()\n");

        //8.130630e+02    5.000000e+02    1.000000e+00    1.572616e+03    1.110889e+06
        T = to_K(8.130630e+02);
        p = to_Pa(500);
        EXPECT_NEAR0(1.572616e+03, driesner07_NaCl_L_rho_pT(p, T), 2e-3);

        // 9.119670e+02    4.500000e+03    1.000000e+00    1.649377e+03    1.967562e+06
        T = to_K(9.119670e+02);
        p = to_Pa(4.500000e+03);
        EXPECT_NEAR0(1.649377e+03, driesner07_NaCl_L_rho_pT(p, T), 1e-3);

        printf("-> test driesner07_NaCl_H_cp_pT()\n");
        T = 300;
        p = to_Pa(1.);
        EXPECT_NEAR0(50.54*1./58.443e-3, driesner07_NaCl_H_cp_pT(p, T), 2e-2); //j/mol/K
        T = 500;
        p = to_Pa(1.);
        EXPECT_NEAR0(53.90*1. / 58.443e-3, driesner07_NaCl_H_cp_pT(p, T), 1e-2); //j/mol/K
        T = 800;
        p = to_Pa(1.);
        EXPECT_NEAR0(59.34*1. / 58.443e-3, driesner07_NaCl_H_cp_pT(p, T), 1e-2); //j/mol/K
        T = 1000;
        p = to_Pa(1.);
        EXPECT_NEAR0(64.84*1. / 58.443e-3, driesner07_NaCl_H_cp_pT(p, T), 1e-3); //j/mol/K
                                                                           //T = to_K(0.);
        //p = to_Pa(1.);
        //EXPECT_NEAR0(855, driesner07_NaCl_H_cp_pT(p, T), 1e-3);
        //T = to_K(400.);
        //p = to_Pa(1.);
        //EXPECT_NEAR0(975, driesner07_NaCl_H_cp_pT(p, T), 1e-3);
        //T = to_K(800.);
        //p = to_Pa(1.);
        //EXPECT_NEAR0(1090, driesner07_NaCl_H_cp_pT(p, T), 1e-3);

        printf("-> test driesner07_NaCl_H_h_pT()\n");
        T = 298.15;
        p = to_Pa(1.);
        double h0 = driesner07_NaCl_H_h_pT(p, T)*1e-3;
        //printf("h0=%g, h_tr=%g\n", h0, driesner07_NaCl_VLH_h_halite());
        T = 300;
        p = to_Pa(1.);
        EXPECT_NEAR0(h0 + 0.09*1. / 58.443e-3, driesner07_NaCl_H_h_pT(p, T)*1e-3, 1e-2);
        T = 500;
        p = to_Pa(1.);
        EXPECT_NEAR0(h0 + 10.56*1. / 58.443e-3, driesner07_NaCl_H_h_pT(p, T)*1e-3, 1e-2);
        T = 800;
        p = to_Pa(1.);
        EXPECT_NEAR0(h0 + 27.48*1. / 58.443e-3, driesner07_NaCl_H_h_pT(p, T)*1e-3, 1e-2);
        T = 1000;
        p = to_Pa(1.);
        EXPECT_NEAR0(h0 + 39.87*1. / 58.443e-3, driesner07_NaCl_H_h_pT(p, T)*1e-3, 1e-2);

    }
}

UTEST(driesner07, H2O_NaCl2)
{
    //---------------------------------------------------------------
    // Driesner 2007: H2O-NaCl
    //---------------------------------------------------------------
    {
        double p, T, x;
        printf("-> test driesner07_H2O_NaCl_Tv_pTx()\n");
        p = 1000.;
        T = 600;
        x = 0;
        EXPECT_NEAR0(T, driesner07_H2O_NaCl_Tv_pTx(p, T, x), 1e-3);
        x = to_mole_fraction(0.1);
        EXPECT_NEAR0(536, driesner07_H2O_NaCl_Tv_pTx(p, T, x), 1e-2);
        T = 400;
        x = to_mole_fraction(0.1);
        EXPECT_NEAR0(366, driesner07_H2O_NaCl_Tv_pTx(p, T, x), 1e-2);

        printf("-> test driesner07_H2O_NaCl_rho_pTx()\n");
#if 1
        // x=0
        T = to_K(15);
        p = to_Pa(1);
        x = 0;
        EXPECT_NEAR0(iapws95_rho_pT(p, T), driesner07_H2O_NaCl_rho_pTx(p, T, x), 1e-3);
        // x=1
        T = to_K(9.119670e+02);
        p = to_Pa(4.500000e+03);
        x = 1.0;
        // EXPECT_NEAR0(1.649377e+03, driesner07_prost_H2O_NaCl_rho_pTx(p, T, x), 1e-3);

        // EXPECT_EQ(driesner07_prost_H2O_NaCl_phase_type(p, T, x), driesner07_H2O_NaCl_phase_type(p, T, x));
        // printf("phase_type = %d\n", driesner07_prost_H2O_NaCl_phase_type(p, T, x));
        // {
        //     double xl = driesner07_prost_H2O_NaCl_LH_xl_pT(p, T);
        //     EXPECT_NEAR0(xl, driesner07_H2O_NaCl_LH_xl_pT(p, T), 1e-3);
        //     double mv = driesner07_prost_H2O_NaCl_vm_singlephase_pTx(p, T, xl);
        //     EXPECT_NEAR0(mv, driesner07_H2O_NaCl_vm_singlephase_pTx(p, T, xl), 1e-3);
        //     double rhol = driesner07_prost_H2O_NaCl_rho_singlephase_pTx(p, T, xl);
        //     EXPECT_NEAR0(rhol, driesner07_H2O_NaCl_rho_singlephase_pTx(p, T, xl), 1e-3);
        // }

        EXPECT_NEAR0(1.649377e+03, driesner07_H2O_NaCl_rho_pTx(p, T, x), 1e-3);

        // 0<x<1
        p = to_Pa(2.810849e+02);
        T = to_K(400);
        x = 6.832051e-03;
        EXPECT_NEAR0(4.324106e+02, driesner07_H2O_NaCl_rho_pTx(p, T, x), 1e-3);

        p = to_Pa(9.166756e+02);
        T = to_K(600);
        x = 7.777610e-02;
        EXPECT_NEAR0(6.248636e+02, driesner07_H2O_NaCl_rho_pTx(p, T, x), 1e-3);

        p = to_Pa(500);
        T = to_K(10);
        x = 1.006877e-01;
        EXPECT_NEAR0(1.223591e+03, driesner07_H2O_NaCl_rhol_pTx2(p, T, x), 1e-3);
        //printf("Tv=%g\n", driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T),x));
#endif
        // low-p,T, high x
        //printf("[exception: low-T]\n");
        p = to_Pa(3.466541);
        T = to_K(150);
        EXPECT_NEAR0(1.831081, driesner07_H2O_NaCl_VL_rhov_pT(p, T), 1e-3);
//        EXPECT_NEAR0(1.831081, driesner07_H2O_NaCl_rhov_pTx(p, T, x), 1e-3);

        // {
        //     double x_prost = driesner07_prost_H2O_NaCl_VL_xl_pT(p, T);
        //     double x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
        //     EXPECT_NEAR0(x_prost, x, 1e-3);
        //     double molar_v = driesner07_prost_H2O_NaCl_vm_singlephase_pTx(p, T, x_prost);
        //     EXPECT_NEAR0(molar_v, driesner07_H2O_NaCl_vm_singlephase_pTx(p, T, x), 1e-3);
        // }

        //TODO error>2% due to some difference in calculation of saturated properties between IAPS84 and PROST
        EXPECT_NEAR0(1.149964e+03, driesner07_H2O_NaCl_VL_rhol_pT(p, T), 3e-2);
        //printf("Tv=%g\n", driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T),x));
        // Prop* pliq = newProp('p', 'T', 0);
        // Prop* pvap = newProp('p', 'T', 0);
        // sat_p(p, pliq, pvap);
        // printf("Tsat=%g\n", to_C(pliq->T));
        // printf("x_l_sat=%g\n", driesner07_H2O_NaCl_VLH_xl_T(T));

        // p = to_Pa(5);
        // T = to_K(150);
        // x = 1.151778e-01; //liquid
        // EXPECT_NEAR0(1.149964e+03, driesner07_H2O_NaCl_rho_pTx(p, T, x), 1e-3);
        // sat_p(p, pliq, pvap);
        // T = pliq->T;
        // printf("Tsat=%g\n", to_C(pliq->T));
        // printf("Vl,sat=%g, Vg,sat=%g\n", 1./pliq->d*18.015e-3*1e6, 1./pvap->d*18.015e-3*1e6);
        // Prop* prop0 = newProp('p', 'T', 0);
        // T = to_K(driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T),x));
        // water_tp(T, p, 0, 1e-8, prop0);
        // //dumpProp(stdout, prop0);
        // printf("Tv=%g\n", to_C(T));
        // printf("Vv=%g\n", 1./prop0->d*18.015e-3*1e6);
        // // dumpProp(stdout, &pliq);
        // // dumpProp(stdout, &pvap);
        // freeProp(pliq);
        // freeProp(pvap);
        // freeProp(prop0);

#if 1
        // p<350bar, T>600C, high x
        //printf("[exception: high-T]\n");
        p = to_Pa(3.698947e+01);
        T = to_K(790);
        EXPECT_NEAR0(7.593799e+00, driesner07_H2O_NaCl_VL_rhov_pT(p, T), 1e-3);
        EXPECT_NEAR0(1.545472e+03, driesner07_H2O_NaCl_VL_rhol_pT(p, T), 1e-3);
        //printf("Tv=%g\n", driesner07_H2O_NaCl_Tv_pTx(to_bar(p), to_C(T),x));
#endif


        printf("-> test driesner07_H2O_NaCl_h_pTx()\n");
        p = to_Pa(2.810849e+02);
        T = to_K(400);
        x = 6.832051e-03;
        EXPECT_NEAR0(2.092815e+06, driesner07_H2O_NaCl_h_pTx(p, T, x), 1e-3);
        p = to_Pa(9.166756e+02);
        T = to_K(600);
        x = 7.777610e-02;
        EXPECT_NEAR0(2.352868e+06, driesner07_H2O_NaCl_h_pTx(p, T, x), 1e-3);

        printf("-> test driesner07_H2O_NaCl_cp_singlephase_pTx()\n");
        // 179 bar
        p = to_Pa(179);
        x = molality_to_mole_fraction(1.8601);
        T = to_K(50);
        EXPECT_NEAR0(3803, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        T = to_K(300);
        EXPECT_NEAR0(4588, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        x = molality_to_mole_fraction(2.9978);
        T = to_K(50);
        EXPECT_NEAR0(3640, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        T = to_K(300);
        EXPECT_NEAR0(4273, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        // 1000 bar
        x = molality_to_mole_fraction(1.);
        p = to_Pa(1000);
        T = to_K(50);
        EXPECT_NEAR0(3811, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        T = to_K(200);
        EXPECT_NEAR0(3928, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        x = molality_to_mole_fraction(5.);
        p = to_Pa(1000);
        T = to_K(50);
        EXPECT_NEAR0(3254, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
        T = to_K(200);
        EXPECT_NEAR0(3336, driesner07_H2O_NaCl_cp_singlephase_pTx(p, T, x), 1e-3);
    }
}

#if 0
UTEST(driesner07, test_300C)
{
    double T = to_K(300);
    double p, x;
    // printf("*pb(300C)=%g bar\n", to_bar(vapor_pressure(T)));
    // printf("*pc(300C)=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("########## T=%gC ##########\n", to_C(T));
    printf("[H2O]\n");
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_t(T, pliq, pvap);
    double p_boil = pliq->p;
    freeProp(pliq);
    freeProp(pvap);
    printf("*p_vapor=%g bar\n", to_bar(p_boil));
    printf("\n[NaCl]\n");
    printf("*p_vapor=%g bar\n", to_bar(driesner07_NaCl_sat_p_T(T)));
    printf("\n[H2O-NaCl]\n");
    printf("critical:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("*x=%g\n", driesner07_H2O_NaCl_xc_T(T));
    printf("VLH:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_VLH_p_T(T)));

    printf("LH:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 100, 200, 300, 400};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VH:\n");
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {1, 10, 20, 30, 40, 50, 60, to_bar(driesner07_H2O_NaCl_VLH_p_T(T))};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VL:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 65, 70, 75, 80, 82, 83, 84, 85, to_bar(p_boil)-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 65, 70, 75, 80, 82, 83, 84, 85, 85.5, to_bar(p_boil)-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    // printf("*(p, xl, xv, logK)\n");
    // {
    //     double pp[] = {61, 65, 70, 75, 80, 85};
    //     for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
    //     {
    //         p = to_Pa(pp[i]);
    //         printf("\t%g\t%g\t%g\t%g\n", to_bar(p),
    //             driesner07_H2O_NaCl_VL_xl_pT(p, T),
    //             driesner07_H2O_NaCl_VL_xv_pT(p, T),
    //             log10(driesner07_H2O_NaCl_VL_K(p,T)));
    //     }
    // }

    printf("K:\n*Kn=(p,pn,logKn,logKm,logK)\n");
    {
        double k;
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 65, 70, 75, 80, 82, 83, 84, 85, 85.5, to_bar(driesner07_H2O_NaCl_pc_T(T))};
        double p_NaCl = to_bar(driesner07_NaCl_sat_p_T(T));
        double pc = to_bar(driesner07_H2O_NaCl_pc_T(T));
        double x_pNaCl = driesner07_H2O_NaCl_LH_xl_pT(to_Pa(p_NaCl), T);
        printf("-> pNaCl        = %.3g bar\n", p_NaCl);
        printf("-> log x(pNaCl) = %.3g wt%% (%g)\n", log10(1e2*to_mass_frac(x_pNaCl)), to_mass_frac(x_pNaCl));
        printf("-> pc           = %.3g bar\n", pc);
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            double pn = (pp[i] - p_NaCl)/(pc - p_NaCl);
            pn = min(pn, 1.0);
            printf("\t%g\t%.3g\t%.3g\t%.3g\t%.3g\n",
                    pp[i], pn,
                    driesner07_H2O_NaCl_VL_log10_Kn(pn, to_C(T)),
                    driesner07_H2O_NaCl_VL_log10_Km(pp[i], T),
                    log10(driesner07_H2O_NaCl_VL_K(p,T)));
        }
    }

    printf("\n");
    //printf("###########################\n", to_C(T));

}

UTEST(driesner07, test_375C)
{
    double T = to_K(375);
    double p, x;
    // printf("*pb(300C)=%g bar\n", to_bar(vapor_pressure(T)));
    // printf("*pc(300C)=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("########## T=%gC ##########\n", to_C(T));
    printf("[H2O]\n");
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_t(T, pliq, pvap);
    double p_boil = pliq->p;
    freeProp(pliq);
    freeProp(pvap);
    printf("*p_vapor=%g bar\n", to_bar(p_boil));
    printf("\n[NaCl]\n");
    printf("*p_vapor=%g bar\n", to_bar(driesner07_NaCl_sat_p_T(T)));
    printf("\n[H2O-NaCl]\n");
    printf("critical:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("*x=%g\n", driesner07_H2O_NaCl_xc_T(T));
    printf("VLH:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_VLH_p_T(T)));

    printf("LH:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 200, 300, 400};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VH:\n");
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {0.7, 1, 1.936, 5, 10, 15, 20, 30, 40, 50, 60, 80, 90, 100, 110, 120, 130, to_bar(driesner07_H2O_NaCl_VLH_p_T(T))};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VL:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 150, 160, 170, 180, 190, 200, 210, 220, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 150, 160, 170, 180, 190, 200, 210, 220, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("\n");
    //printf("###########################\n", to_C(T));

}


UTEST(driesner07, test_500C)
{
    double T = to_K(500);
    double p, x;
    // printf("*pb(300C)=%g bar\n", to_bar(vapor_pressure(T)));
    // printf("*pc(300C)=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("########## T=%gC ##########\n", to_C(T));
    printf("[H2O]\n");
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_t(T, pliq, pvap);
    double p_boil = pliq->p;
    freeProp(pliq);
    freeProp(pvap);
    printf("*p_vapor=%g bar\n", to_bar(p_boil));
    printf("\n[NaCl]\n");
    printf("*p_vapor=%g bar\n", to_bar(driesner07_NaCl_sat_p_T(T)));
    printf("\n[H2O-NaCl]\n");
    printf("critical:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("*x=%g\n", driesner07_H2O_NaCl_xc_T(T));
    printf("VLH:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_VLH_p_T(T)));

    printf("LH:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 800};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VH:\n");
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {0.7, 1, 1.936, 5, 10, 15, 20, 30, 40, 50, 60, 80, 90, 100,
                        110, 120, 130, 200, 250, 300,
                        to_bar(driesner07_H2O_NaCl_VLH_p_T(T))};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("VL:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 350, 400, 450, 500, 550, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 350, 400, 450, 500, 550, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }

    printf("\n");
    //printf("###########################\n", to_C(T));

}


UTEST(driesner07, test_800C)
{
    double T = to_K(800);
    double p, x;
    // printf("*pb(300C)=%g bar\n", to_bar(vapor_pressure(T)));
    // printf("*pc(300C)=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("########## T=%gC ##########\n", to_C(T));
    printf("[H2O]\n");
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_t(T, pliq, pvap);
    double p_boil = pliq->p;
    freeProp(pliq);
    freeProp(pvap);
    printf("*p_vapor=%g bar\n", to_bar(p_boil));
    printf("\n[NaCl]\n");
    printf("*p_vapor=%g bar\n", to_bar(driesner07_NaCl_sat_p_T(T)));
    printf("\n[H2O-NaCl]\n");
    printf("critical:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("*x=%g\n", driesner07_H2O_NaCl_xc_T(T));
    printf("VLH:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_VLH_p_T(T)));

    printf("LH:\n");
    printf("*x_liquid=(log wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 2500};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
        printf("*x_liquid=(wt.%%, bar)\n");
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
        }
    }

    printf("VH:\n");
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {0.1, 0.5, 1, 2,
                        to_bar(driesner07_H2O_NaCl_VLH_p_T(T))};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
        printf("*x_vapor=(wt.%%, bar)\n");
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
            printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
        }
    }

    printf("VL:\n");
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 500, 1000, 1400, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
    }
    printf("*x_liquid=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 500, 1000, 1400, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
        }
    }
    printf("*x_vapor=(wt.%%, bar)\n");
    {
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 5, 10, 20, 40, 60, 80, 100,
                         500, 1000, 1400, to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
        printf("*x_vapor=(wt.%%, bar)\n");
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
        }
    }

    printf("\n");
    //printf("###########################\n", to_C(T));

}


UTEST(driesner07, test_1000C)
{
    double T = to_K(1000);
    double p, x;
    // printf("*pb(300C)=%g bar\n", to_bar(vapor_pressure(T)));
    // printf("*pc(300C)=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("########## T=%gC ##########\n", to_C(T));
    printf("[H2O]\n");
    Prop* pliq = newProp('p', 'T', 0);
    Prop* pvap = newProp('p', 'T', 0);
    sat_t(T, pliq, pvap);
    double p_boil = pliq->p;
    freeProp(pliq);
    freeProp(pvap);
    printf("*p_vapor=%g bar\n", to_bar(p_boil));
    printf("\n[NaCl]\n");
    printf("*p_vapor=%g bar\n", to_bar(driesner07_NaCl_sat_p_T(T)));
    printf("\n[H2O-NaCl]\n");
    printf("critical:\n");
    printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_pc_T(T)));
    printf("*x=%g\n", driesner07_H2O_NaCl_xc_T(T));
    printf("VLH:\n");
    //printf("*p=%g bar\n", to_bar(driesner07_H2O_NaCl_VLH_p_T(T)));

    printf("LH:\n");
    // printf("*x_liquid=(log wt.%%, bar)\n");
    // {
    //     double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)), 2500};
    //     for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
    //     {
    //         p = to_Pa(pp[i]);
    //         x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
    //         printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
    //     }
    //     printf("*x_liquid=(wt.%%, bar)\n");
    //     for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
    //     {
    //         p = to_Pa(pp[i]);
    //         x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
    //         printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
    //     }
    // }

    printf("VH:\n");
    // printf("*x_vapor=(wt.%%, bar)\n");
    // {
    //     double pp[] = {0.1, 0.5, 1, 2,
    //                     to_bar(driesner07_H2O_NaCl_VLH_p_T(T))};
    //     for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
    //     {
    //         p = to_Pa(pp[i]);
    //         x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
    //         printf("\t%g %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
    //     }
    //     printf("*x_vapor=(wt.%%, bar)\n");
    //     for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
    //     {
    //         p = to_Pa(pp[i]);
    //         x = driesner07_H2O_NaCl_VH_xv_pT(p, T);
    //         printf("\t%g %g\n", to_mass_frac(x), to_bar(p));
    //     }
    // }

    printf("VL:\n");
    printf("*x_liquid=(log10 wt.%%, bar)\n");

    if(0){
        double pp[] = {to_bar(driesner07_H2O_NaCl_VLH_p_T(T)),
                        1.524812e+01, 7.657447e+01, 100,
                        500, 1000, 1400, 1800, 2000, 2100,
                        to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%e %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
        printf("*x_liquid=(wt.%%, bar)\n");
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xl_pT(p, T);
            printf("\t%e %g\n", to_mass_frac(x), to_bar(p));
        }
        printf("*x_LH=(wt.%%, bar)\n");
        // for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        // {
            // p = driesner07_NaCl_sat_p_T(T);
            // x = driesner07_H2O_NaCl_LH_xl_pT(p, T);
            // printf("\tNaCl vapor p=%e\n", to_bar(p));
            // printf("\tNaCl melt. T=%e\n", to_C(driesner07_NaCl_LH_T_p((p)));
            // printf("\tx_LH,sat =%e\n", to_mass_frac(x));
            // x = driesner07_H2O_NaCl_LH_xl_pT(to_Pa(2e3), to_K(800));
            // printf("\tx_LH,sat =%e\n", x);

            // T = to_K(790);
            // p = to_Pa(3.698947e+01);
            // printf("\t%g %g %e %e\n", to_C(T), to_bar(p),
            //                         driesner07_H2O_NaCl_VL_xv_pT(p, T),
            //                         driesner07_H2O_NaCl_VL_xl_pT(p, T)
            //                         );
            // T = to_K(800);
            // p = to_Pa(3.358897e+01);
            // printf("\t%g %g %e %e\n", to_C(T), to_bar(p),
            //                         driesner07_H2O_NaCl_VL_xv_pT(p, T),
            //                         driesner07_H2O_NaCl_VL_xl_pT(p, T)
            //                         );

        printf("*debug\n");
            T = to_K(810);
            p = to_Pa(3.161988e+01);
            printf("\t%g %g %e %e\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T)
                                    );
            p = to_Pa(6.408160e+01);
            printf("\t%g %g %e %e\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T)
                                    );
            T = to_K(820);
            p = to_Pa(3.117224e+01);
            printf("\t%g %g %e %e\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T)
                                    );
            printf("\tT\tp\txv\t\txl\t\tlogK'\n");
            T = to_K(1000);
            p = to_Pa(15.24812);
            printf("\t%g\t%g\t%e\t%e\t%g\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T),
                                    driesner07_H2O_NaCl_VL_log10_Km(to_bar(p), T)
                                    );
            p = to_Pa(7.657447e+01);
            printf("\t%g\t%g\t%e\t%e\t%g\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T),
                                    driesner07_H2O_NaCl_VL_log10_Km(to_bar(p), T)
                                    );
            p = to_Pa(1.379008e+02);
            printf("\t%g\t%g\t%e\t%e\t%g\n", to_C(T), to_bar(p),
                                    driesner07_H2O_NaCl_VL_xv_pT(p, T),
                                    driesner07_H2O_NaCl_VL_xl_pT(p, T),
                                    driesner07_H2O_NaCl_VL_log10_Km(to_bar(p), T)
                                    );
            double p_NaCl = to_bar(driesner07_NaCl_VH_p_T(T));
            x = driesner07_H2O_NaCl_LH_xl_pT(to_Pa(p_NaCl), T);
            printf("\tx=%e\n", x);

            T = to_K(1000);
        // }
        //return;
    }
    printf("*x_vapor=(log10 wt.%%, bar)\n");
    {
        double pp[] = {0.1, 5, 10, 1.524812e+01, 20, 40, 60, 7.657447e+01, 80, 100,
                         200, 500, 1000, 1400, 1800, 2000, 2100,
                        to_bar(driesner07_H2O_NaCl_pc_T(T))-1e-12};
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%e %g\n", log10(1e2*to_mass_frac(x)), to_bar(p));
        }
        printf("*x_vapor=(wt.%%, bar)\n");
        for (int i=0; i<sizeof(pp)/sizeof(pp[0]); i++)
        {
            p = to_Pa(pp[i]);
            x = driesner07_H2O_NaCl_VL_xv_pT(p, T);
            printf("\t%e %g\n", to_mass_frac(x), to_bar(p));
        }
    }

    printf("\n");
    //printf("###########################\n", to_C(T));

}
#endif

UTEST(driesner07, bulk)
{
    double p, T, x;
    T = to_K(500); p = 50.e6;
    //
    {
        x = 0;
        double massfrac_l = driesner07_H2O_NaCl_mass_fraction_liquid(p, T, x);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        EXPECT_NEAR3(0.0, massfrac_l);
        EXPECT_NEAR3(0.0, volfrac_l);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR3(driesner07_H2O_NaCl_rho_singlephase_pTx(p, T, x), rho);
        EXPECT_NEAR3(driesner07_H2O_NaCl_h_singlephase_pTx(p, T, x), h);
    }

    {
        x = to_mole_fraction(0.5);
        double massfrac_l = driesner07_H2O_NaCl_mass_fraction_liquid(p, T, x);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        EXPECT_NEAR3(1.0, massfrac_l);
        EXPECT_NEAR3(1.0, volfrac_l);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR3(driesner07_H2O_NaCl_rho_singlephase_pTx(p, T, x), rho);
        EXPECT_NEAR3(driesner07_H2O_NaCl_h_singlephase_pTx(p, T, x), h);
    }

    //
    {
        x = to_mole_fraction(0.2);
        double rhol = driesner07_H2O_NaCl_VL_rhol_pT(p, T);
        double rhov = driesner07_H2O_NaCl_VL_rhov_pT(p, T);
        double xl = driesner07_H2O_NaCl_VL_xl_pT(p, T);
        double xv = driesner07_H2O_NaCl_VL_xv_pT(p, T);
        // printf("x=%g, xl=%g, xv=%g\n", x, xl, xv);
        // printf("rhol=%g, rhov=%g\n", rhol, rhov);
        double massfrac = driesner07_H2O_NaCl_mass_fraction_liquid(p, T, x);
        double volfrac = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        // printf("Ml=%g, satl=%g\n", massfrac, volfrac);
        EXPECT_NEAR3(0.48948, massfrac);
        EXPECT_NEAR3(0.23515, volfrac);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        double hl = driesner07_H2O_NaCl_h_singlephase_pTx(p, T, xl);
        double hv = driesner07_H2O_NaCl_h_singlephase_pTx(p, T, xv);
        // printf("hl=%g, hv=%g\n", hl*1e-3, hv*1e-3);

        EXPECT_NEAR3(410.45, rho);
        EXPECT_NEAR3(2213.36e3, h);
    }

    p = 20e6;
    {
        x = to_mole_fraction(0.5);
        double massfrac_l = driesner07_H2O_NaCl_mass_fraction_liquid(p, T, x);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        EXPECT_NEAR3(0.0, massfrac_l);
        EXPECT_NEAR3(0.0, volfrac_l);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        double xv = driesner07_H2O_NaCl_VH_xv_pT(p, T);
        EXPECT_NEAR3(driesner07_H2O_NaCl_rho_singlephase_pTx(p, T, xv), rho);
        EXPECT_NEAR3(driesner07_H2O_NaCl_h_singlephase_pTx(p, T, xv), h);
    }
}

UTEST(driesner07, inverse)
{
    double p, T, x, T0;
    T = to_K(500); p = 50.e6;
    T0 = T - 10.0;
    {
        x = to_mole_fraction(0.2);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        EXPECT_NEAR(volfrac_l, driesner07_H2O_NaCl_volume_fraction_liquid_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

    {
        x = to_mole_fraction(0.);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

    {
        x = to_mole_fraction(0.5);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR3(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

    p = 20e6;
    {
        x = to_mole_fraction(0.5);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR3(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

    p = 20e6;
    {
        x = to_mole_fraction(0.5);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR3(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

    T = to_K(300); p = 0.1e6;
    T0 = T - 10.0;
    {
        x = to_mole_fraction(0.0);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR3(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }
    {
        x = to_mole_fraction(0.03);
        double volfrac_l = driesner07_H2O_NaCl_volume_fraction_liquid(p, T, x);
        double rho = driesner07_H2O_NaCl_rho_pTx(p, T, x);
        EXPECT_NEAR3(T, driesner07_H2O_NaCl_T_rhopx(rho, p, x, T0));
        double h = driesner07_H2O_NaCl_h_pTx(p, T, x);
        EXPECT_NEAR(T, driesner07_H2O_NaCl_T_phx(p, h, x, T0));
    }

}

UTEST(driesner07, inverse_vl_p)
{
    double p, T, x, p0;
    T = to_K(300);
    p = 7.0e6;
    {
        p0 = p * 0.99;
        double xv = driesner07_H2O_NaCl_VL_xv_pT(p, T);
        double xl = driesner07_H2O_NaCl_VL_xl_pT(p, T);

        // for (int i=0; i<10; i++)
        // {
        //     double pp = 0.65e6 + 0.01e6 * i;
        //     printf("pMPa=%g, xv=%g\n", pp*1e-6, driesner07_H2O_NaCl_VL_xv_pT(pp, T));
        // }

        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txv(T, xv, p0));
        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txl(T, xl, p0));
    }
    x = 1e-5;
    {
        p0 = p * 0.99;
        double xv = driesner07_H2O_NaCl_VL_xv_pT(p, T);
        double xl = driesner07_H2O_NaCl_VL_xl_pT(p, T);
        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txv(T, xv, p0));
        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txl(T, xl, p0));
    }
    x = 0.1;
    {
        p = driesner07_H2O_NaCl_VL_p_Txl(T, x, p0);
        p0 = p * 0.99;
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txv(T, x, p0));
        EXPECT_NEAR3(x, driesner07_H2O_NaCl_VL_xl_pT(p, T));
    }
    x = 0.5;
    {
        p0 = p * 0.99;
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txv(T, x, p0));
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txl(T, x, p0));
    }

    T = to_K(500);
    p = 40.e6;
    {
        p0 = p * 0.99;
        double xv = driesner07_H2O_NaCl_VL_xv_pT(p, T);
        double xl = driesner07_H2O_NaCl_VL_xl_pT(p, T);
        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txv(T, xv, p0));
        EXPECT_NEAR3(p, driesner07_H2O_NaCl_VL_p_Txl(T, xl, p0));
    }
    x = 1e-5;
    {
        p0 = p * 0.99;
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txv(T, x, p0));
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txl(T, x, p0));
    }
    x = 0.01;
    {
        p = driesner07_H2O_NaCl_VL_p_Txv(T, x, p0);
        p0 = p * 0.99;
        EXPECT_NEAR3(x, driesner07_H2O_NaCl_VL_xv_pT(p, T));
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txl(T, x, p0));
    }
    x = 0.3;
    {
        p = driesner07_H2O_NaCl_VL_p_Txl(T, x, p0);
        p0 = p * 0.99;
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txv(T, x, p0));
        EXPECT_NEAR3(x, driesner07_H2O_NaCl_VL_xl_pT(p, T));
    }
    x = 0.8;
    {
        p0 = p * 0.99;
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txv(T, x, p0));
        EXPECT_NEAR3(-1, driesner07_H2O_NaCl_VL_p_Txl(T, x, p0));
    }
}

static double h2o_sat_p(double T)
{
    return iaps84_sat_p_T(T)*1e6;
}

UTEST(driesner07, sat_p_T300)
{
    double p, T, x, p0;
    double pp[4] = {0};
    T = to_K(300);
    x = 0;
    {
        p = h2o_sat_p(T);
        p0 = p*0.99;
        EXPECT_NEAR3(1, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(8.5838e6, pp[0]);
    }
    x = 1e-8;
    {
        p = h2o_sat_p(T);
        p0 = p*0.99;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(8.5815e+06, pp[0]);
        EXPECT_NEAR3(8.5838e6,   pp[1]);
    }
    x = 1e-6;
    {
        p = h2o_sat_p(T);
        p0 = 6.05e6; //p*0.99;
        EXPECT_NEAR3(4, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(6.0472e+06, pp[0]);
        EXPECT_NEAR3(7.3929e+06, pp[1]);
        EXPECT_NEAR3(8.0631e+06, pp[2]);
        EXPECT_NEAR3(8.5838e6,   pp[3]);
    }
    x = to_mole_fraction(3.2e-2);
    {
        p = h2o_sat_p(T);
        p0 = p*0.99;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(6.0472e6, pp[0]);
        EXPECT_NEAR3(8.4497e6, pp[1]);
    }
}

UTEST(driesner07, sat_p_T375)
{
    double p, T, x, p0;
    double pp[4] = {0};

    T = to_K(375.);
    p = h2o_sat_p(to_K(373.));
    p0 = p*0.99;
    x = 0;
    {
        EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    }
    x = 5e-5;
    {
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(14.143e6, pp[0]);
        EXPECT_NEAR3(18.449e6, pp[1]);
    }
    x = 1.5e-4;
    {
        // double pc = driesner07_H2O_NaCl_pc_T(T);
        // for (int i=0; i<10; i++)
        // {
        //     double px = pc-1e2*i;
        //     double xv = driesner07_H2O_NaCl_VL_xv_pT(px, T);
        //     double xl = driesner07_H2O_NaCl_VL_xl_pT(px, T);
        //     double p_NaCl = to_bar(driesner07_NaCl_VH_p_T(T));
        //     double pc = to_bar(driesner07_H2O_NaCl_pc_T(T));
        //     double pn = (to_bar(px) - p_NaCl)/(pc - p_NaCl);
        //     double log10Kn = driesner07_H2O_NaCl_VL_log10_Kn(pn, to_C(T));
        //     double log10Km = driesner07_H2O_NaCl_VL_log10_Km(to_bar(px), T);
        //     double K = driesner07_H2O_NaCl_VL_K(px, T);
        //     printf("%.12g %.5e %.5e %.5e %.5e %.5e\n", px*1e-6, xv, xl, log10Kn, log10Km, K);
        // }
        p0 = 22.2915e6;
        EXPECT_NEAR3(4, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(14.143e6, pp[0]);
        EXPECT_NEAR3(21.0269e6, pp[1]);
        EXPECT_NEAR3(22.2707e6, pp[2]);
        EXPECT_NEAR3(22.2953e6, pp[3]);
    }
    x = 0.1;
    {
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(14.143e6, pp[0]);
        EXPECT_NEAR3(17.6591e6, pp[1]);
    }
    x = 0.6;
    {
        EXPECT_NEAR3(1, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(14.143e6, pp[0]);
    }

}

UTEST(driesner07, sat_p_T500)
{
    double p, T, x, p0;
    double pp[4] = {0};

    T = to_K(500.);
    p0 = 50e6;
    x = 0;
    {
        EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    }
    x = 0.01;
    {
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(32.2003e6, pp[0]);
        EXPECT_NEAR3(54.2272e6, pp[1]);
    }
    x = 0.17;
    {
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(32.2003e6, pp[0]);
        EXPECT_NEAR3(46.7496e6, pp[1]);
    }
    x = 0.4;
    {
        EXPECT_NEAR3(1, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(32.2003e6, pp[0]);
    }

}


UTEST(driesner07, sat_p_T800)
{
    double p, T, x, p0;
    double pp[4] = {0};

    T = to_K(800.);
    p0 = 100e6;
    x = to_mole_fraction(1e-5);
    {
        EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    }
    x = to_mole_fraction(3e-4);
    {
        // for (int i=0; i<10; i++)
        // {
        //     double px = 10e6-1e6*i;
        //     double xv = driesner07_H2O_NaCl_VL_xv_pT(px, T);
        //     double xl = driesner07_H2O_NaCl_VL_xl_pT(px, T);
        //     printf("%.12g %.5e %.5e\n", px*1e-6, xv, xl);
        // }
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3( 5.9573e5, pp[0]);
        EXPECT_NEAR3(24.0527e6, pp[1]);
    }
    x = to_mole_fraction(1e-2);
    {
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3( 2.4641e5, pp[0]);
        EXPECT_NEAR3(79.8087e6, pp[1]);
    }

    x = to_mole_fraction(0.6);
    {
        p0 = 1e6;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(  2.4641e5, pp[0]);
        EXPECT_NEAR3(116.8677e6, pp[1]);
    }

}

UTEST(driesner07, sat_p_T900)
{
    double p, T, x, p0;
    double pp[4] = {0};

    T = to_K(900.);
    p0 = 100e6;
    // x = 0;
    // {
    //     EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    // }
    // x = to_mole_fraction(0.001);
    // {
    //     EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    // }
    x = to_mole_fraction(0.2);
    {
        // for (int i=0; i<10; i++)
        // {
        //     double px = 0.001e6-0.0001e6*i;
        //     double xv = driesner07_H2O_NaCl_VL_xv_pT(px, T);
        //     double xl = driesner07_H2O_NaCl_VL_xl_pT(px, T);
        //     printf("%.12g %.5e %.5e\n", px*1e-6, xv, xl);
        // }
        p0 = 1e4;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(4.46019e3, pp[0]);
        EXPECT_NEAR3(182.798e6, pp[1]);
    }

    // x = to_mole_fraction(0.6);
    // {
    //     p0 = 1e6;
    //     EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    //     EXPECT_NEAR3(  4.9587e3, pp[0]);
    //     EXPECT_NEAR3(196.1540e6, pp[1]);
    // }

}


UTEST(driesner07, sat_p_T1000)
{
    double p, T, x, p0;
    double pp[4] = {0};

    T = to_K(1000.);
    p0 = 100e6;
    x = 0;
    {
        EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    }
    x = to_mole_fraction(0.001);
    {
        EXPECT_NEAR3(0, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
    }
    x = to_mole_fraction(0.01);
    {
        // for (int i=0; i<10; i++)
        // {
        //     double px = 2e6-0.1e6*i;
        //     double xv = driesner07_H2O_NaCl_VL_xv_pT(px, T);
        //     double xl = driesner07_H2O_NaCl_VL_xl_pT(px, T);
        //     printf("%.12g %.5e %.5e\n", px*1e-6, xv, xl);
        // }
        p0 = 1e6;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(5.3167e5, pp[0]);
        EXPECT_NEAR3(71.791e6, pp[1]);
    }

    x = to_mole_fraction(0.6);
    {
        p0 = 1e6;
        EXPECT_NEAR3(2, driesner07_H2O_NaCl_sat_p_Tx(T, x, p0, pp));
        EXPECT_NEAR3(  4.9587e3, pp[0]);
        EXPECT_NEAR3(196.1540e6, pp[1]);
    }

}

#if 0
UTEST(driesner07, temp_write_sat)
{
    double dT = 10;
    int n = 900/dT;
    double satp[4];
    double x = to_mole_fraction(40.e-2);
    printf("Ts ps (x=%g)\n", to_mass_frac(x));
    for (int i=0; i<=n; i++)
    {
        double T = to_K(dT*i+100.);
        int np = driesner07_H2O_NaCl_sat_p_Tx(T, x, 1e6, satp);
        for (int j=0; j<np; j++)
        {
            printf("%g %g\n", to_C(T), to_bar(satp[j]));
        }
    }
}
#endif
