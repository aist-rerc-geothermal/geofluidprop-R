
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-SAT-92.h"
#include "model/iapws/IAPWS-Melt-11.h"
#include "model/iapws/IAPWS-Viscosity-08.h"
#include "model/iapws/IAPWS-ThermalConductivity-11.h"

#include "util/utility.h"

#include "utest.h"


UTEST(iapws95, double)
{
    // printf("-> test IAPWS-95\n");
#if 0
    {
        double p, h, rho, T;
#if 0
        printf("rho, T, p, h\n");
        rho = 1.0;
        for (int i = 0; i < 20; i++)
        {
            T = 273.16 + 5 * i;
            h = iapws95_h_rhoT(rho, T);
            p = iapws95_p_rhoT(rho, T);
            printf("%g, %g, %g, %g\n", rho, T, p, h*1e-3);
        }
#endif
        p = 1541.9493333333335;
        h = 63103.154430379742;
        rho = 4.3155138094697332;
        T = 286.59208251018418;
        iapws95_rhoT_ph_long(p, h, &rho, &T);
        p = 1541.95, h = 1.54126e+06;
        rho = 0.01; T = 282.51;
        //iapws95_rhoT_ph(p, h, &rho, &T);
        //p = 10.e6;
        //printf("h, rho, T\n");
        //for (int i = 0; i < 20; i++)
        //{
        //    h = 2000.e3 + 50.e3*i;
        //    iapws95_rhoT_ph(p, h, &rho, &T);
        //    printf("%g, %g, %g\n", h, rho, T);
        //}
        exit(0);
    }
#endif
    const double rhoc = iapws95_get_rhoc();
    const double Tc = iapws95_get_Tc();
    {
        double T=500, rho=838.025;
        double delta = rho/rhoc;
        double tau = Tc/T;
        EXPECT_NEAR(0.204797733e1, iapws95_phi0(delta, tau));
        EXPECT_NEAR(0.384236747, iapws95_dphi0_ddelta(delta, tau));
        EXPECT_NEAR(-0.147637878, iapws95_d2phi0_ddelta2(delta, tau));
        EXPECT_NEAR(0.904611106e1, iapws95_dphi0_dtau(delta, tau));
        EXPECT_NEAR(-0.193249185e1, iapws95_d2phi0_dtau2(delta, tau));
        EXPECT_NEAR(0.0, iapws95_d2phi0_ddelta_dtau(delta, tau));

        EXPECT_NEAR(-0.342693206e1, iapws95_phir(delta, tau));
        EXPECT_NEAR(-0.364366650, iapws95_dphir_ddelta(delta, tau));
        EXPECT_NEAR(0.856063701, iapws95_d2phir_ddelta2(delta, tau));
        EXPECT_NEAR(-0.581403435e1, iapws95_dphir_dtau(delta, tau));
        EXPECT_NEAR(-0.223440737e1, iapws95_d2phir_dtau2(delta, tau));
        EXPECT_NEAR(-0.112176915e1, iapws95_d2phir_ddelta_dtau(delta, tau));
    }
    {
        double T=647, rho=358;
        double delta = rho/rhoc;
        double tau = Tc/T;
        EXPECT_NEAR(-0.156319605e1, iapws95_phi0(delta, tau));
        EXPECT_NEAR(0.899441341, iapws95_dphi0_ddelta(delta, tau));
        EXPECT_NEAR(-0.808994726, iapws95_d2phi0_ddelta2(delta, tau));
        EXPECT_NEAR(0.980343918e1, iapws95_dphi0_dtau(delta, tau));
        EXPECT_NEAR(-0.343316334e1, iapws95_d2phi0_dtau2(delta, tau));
        EXPECT_NEAR(0.0, iapws95_d2phi0_ddelta_dtau(delta, tau));

        EXPECT_NEAR(-0.121202657e1, iapws95_phir(delta, tau));
        EXPECT_NEAR(-0.714012024, iapws95_dphir_ddelta(delta, tau));
        EXPECT_NEAR(0.475730696, iapws95_d2phir_ddelta2(delta, tau));
        EXPECT_NEAR(-0.321722501e1, iapws95_dphir_dtau(delta, tau));
        EXPECT_NEAR(-0.996029507e1, iapws95_d2phir_dtau2(delta, tau));
        EXPECT_NEAR(-0.133214720e1, iapws95_d2phir_ddelta_dtau(delta, tau));
    }

    {
        double T=300, rho=0.996556e3;
        EXPECT_NEAR(0.992418352e-1,     iapws95_p_rhoT(rho, T)*1e-6);
        EXPECT_NEAR(0.413018112e1,     iapws95_cv_rhoT(rho, T)*1e-3);
        EXPECT_NEAR(0.393062643,         iapws95_s_rhoT(rho, T)*1e-3);

        rho=0.1005308e4;
        EXPECT_NEAR(0.200022515e2,     iapws95_p_rhoT(rho, T)*1e-6);
        EXPECT_NEAR(0.406798347e1,     iapws95_cv_rhoT(rho, T)*1e-3);
        EXPECT_NEAR(0.387405401,         iapws95_s_rhoT(rho, T)*1e-3);

        T = 500;
        rho = 0.435;
        EXPECT_NEAR(0.999679423e-1,     iapws95_p_rhoT(rho, T)*1e-6);
        EXPECT_NEAR(0.150817541e1,     iapws95_cv_rhoT(rho, T)*1e-3);
        EXPECT_NEAR(0.794488271e1,     iapws95_s_rhoT(rho, T)*1e-3);

        T = 647;
        rho = 0.3580000e3;
        EXPECT_NEAR(0.220384756e2,     iapws95_p_rhoT(rho, T)*1e-6);
        EXPECT_NEAR(0.618315728e1,     iapws95_cv_rhoT(rho, T)*1e-3);
        EXPECT_NEAR(0.432092307e1,     iapws95_s_rhoT(rho, T)*1e-3);

        T = 900;
        rho = 0.870769e3;
        EXPECT_NEAR(0.700000006e3,     iapws95_p_rhoT(rho, T)*1e-6);
        EXPECT_NEAR(0.266422350e1,     iapws95_cv_rhoT(rho, T)*1e-3);
        EXPECT_NEAR(0.417223802e1,     iapws95_s_rhoT(rho, T)*1e-3);

        // T=300;
        // rho=996.534;
        // EXPECT_NEAR(0.05024,             iapws95_p_rhoT(rho, T)*1e-6);
        // EXPECT_NEAR(0.41303e1,         iapws95_cv_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(0.3931,             iapws95_s_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(112.557,             iapws95_u_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(112.608,             iapws95_h_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(4.1808,             iapws95_cp_rhoT(rho, T)*1e-3);

        // T = 400;
        // rho = 937.617;
        // EXPECT_NEAR(0.50058,            iapws95_p_rhoT(rho, T)*1e-6);
        // EXPECT_NEAR(3.6321,             iapws95_cv_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(1.6010,             iapws95_s_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(532.594,             iapws95_u_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(533.127,             iapws95_h_rhoT(rho, T)*1e-3);
        // EXPECT_NEAR(4.2548,             iapws95_cp_rhoT(rho, T)*1e-3);

    }
}

UTEST(iapws95, rho_pT)
{
    {
        double T = 300;
        double p = 0.05e6;
        EXPECT_NEAR(9.965338910044e+02,             iapws95_rho_pT(p, T));
        // double sat_p, sat_rhol, sat_rhov;
        // saturation(T, &sat_p, &sat_rhol, &sat_rhov);
        // printf("p_sat=%g, rhol_sat=%g, rhol_sat=%g\n", sat_p*1e-6, sat_rhol, sat_rhov);

        T = 400;
        p = 0.5e6;
        EXPECT_NEAR(9.376167009209e+02,             iapws95_rho_pT(p, T));
    }
}

UTEST(iapws95, T_rhop)
{
    {
        double T = 300;
        double p = 0.05e6;
        double rho = iapws95_rho_pT(p, T); //996.53338910044;
        EXPECT_NEAR(T,             iapws95_T_rhop(rho, p));

        T = 400;
        p = 0.5e6;
        rho = iapws95_rho_pT(p, T); //937.617
        EXPECT_NEAR(T,             iapws95_T_rhop(rho, p));
    }
}

UTEST(iapws95, rhoT_ph)
{
    {
        //printf("-> test rhoT_ph()\n");
        double Te = 300;
        double rhoe = 996.534;
        double p = 0.05024e6;
        double h = 112.608e3;
        double T = 0, rho = 0;
        iapws95_rhoT_ph(p, h, &rho, &T);
        EXPECT_NEAR0(Te, T, 1e-6);
        EXPECT_NEAR0(rhoe, rho, 1e-7);

        Te = 400;
        rhoe = 937.617;
        p = 0.50058e6;
        h = 533.127e3;
        iapws95_rhoT_ph(p, h, &rho, &T);
        EXPECT_NEAR0(Te, T, 1e-6);
        EXPECT_NEAR0(rhoe, rho, 1e-7);
    }
}

UTEST(iapws95, sat)
{
    {
        double T = 275;
        double sat_p = 1e6, sat_rhol = 1e3, sat_rhov=1e-3;
        iapws95_sat_prho_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.698451167e-3,     sat_p*1e-6);
        EXPECT_NEAR(0.999887406e3,     sat_rhol);
        EXPECT_NEAR(0.550664919e-2,     sat_rhov);

        T = 450;
        iapws95_sat_prho_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.932203564,         sat_p*1e-6);
        EXPECT_NEAR(0.890341250e3,     sat_rhol);
        EXPECT_NEAR(0.481200360e1,     sat_rhov);

        T = 625;
        iapws95_sat_prho_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.169082693e2,     sat_p*1e-6);
        EXPECT_NEAR(0.567090385e3,     sat_rhol);
        EXPECT_NEAR(0.118290280e3,     sat_rhov);
    }
}

#ifdef USE_LONGDOUBLE
UTEST(iapws95, sat_long)
{
    {
        long double T = 275;
        long double sat_p = 1e6, sat_rhol = 1e3, sat_rhov=1e-3;
        iapws95_sat_prho_T_long(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.698451167e-3,     sat_p*1e-6);
        EXPECT_NEAR(0.999887406e3,     sat_rhol);
        EXPECT_NEAR(0.550664919e-2,     sat_rhov);

        T = 450;
        sat_p = 1e6; sat_rhol = 890; sat_rhov = 4.8;
        iapws95_sat_prho_T_long(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.932203564,         sat_p*1e-6);
        EXPECT_NEAR(0.890341250e3,     sat_rhol);
        EXPECT_NEAR(0.481200360e1,     sat_rhov);

        T = 625;
        iapws95_sat_prho_T_long(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.169082693e2,     sat_p*1e-6);
        EXPECT_NEAR(0.567090385e3,     sat_rhol);
        EXPECT_NEAR(0.118290280e3,     sat_rhov);
    }
}
#endif

#ifdef USE_QUAD
UTEST(iapws95, sat_quad)
{
    {
        __float128 T = 275;
        __float128 sat_p = 0.698451167e3, sat_rhol = 0.999887406e3, sat_rhov=0.550664919e-2;
        iapws95_sat_prho_T_quad(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.698451167e-3,     sat_p*1e-6);
        EXPECT_NEAR(0.999887406e3,     sat_rhol);
        EXPECT_NEAR(0.550664919e-2,     sat_rhov);

        T = 450;
        sat_p = 1e6; sat_rhol = 890; sat_rhov = 4.8;
        iapws95_sat_prho_T_quad(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.932203564,         sat_p*1e-6);
        EXPECT_NEAR(0.890341250e3,     sat_rhol);
        EXPECT_NEAR(0.481200360e1,     sat_rhov);

        T = 625;
        iapws95_sat_prho_T_quad(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR(0.169082693e2,     sat_p*1e-6);
        EXPECT_NEAR(0.567090385e3,     sat_rhol);
        EXPECT_NEAR(0.118290280e3,     sat_rhov);
    }
}
#endif

// UTEST(iapws95, sat2)
// {
//     {
//         double T, sat_rhol, sat_rhov;
//         double p = 0.698451167e-3*1e6;
//         T=0; sat_rhol = 0; sat_rhov = 0;
//         iapws95_sat_Trho_p(p, &T, &sat_rhol, &sat_rhov);
//         EXPECT_NEAR(275.,                 T);
//         EXPECT_NEAR(0.999887406e3,     sat_rhol);
//         EXPECT_NEAR(0.550664919e-2,     sat_rhov);
//         //printf("%g, %g, %g\n", sat_p, sat_rhol, sat_rhov);
//         //abort();

//         p = 0.932203564*1e6;
//         //T=0; sat_rhol = 0; sat_rhov = 0;
//         T=450; sat_rhol = 1e3; sat_rhov = 10;
//         iapws95_sat_Trho_p(p, &T, &sat_rhol, &sat_rhov);
//         EXPECT_NEAR(450.,                 T);
//         EXPECT_NEAR(0.890341250e3,     sat_rhol);
//         EXPECT_NEAR(0.481200360e1,     sat_rhov);

//         p = 0.169082693e2*1e6;
//         //T=0; sat_rhol = 0; sat_rhov = 0;
//         T=625; sat_rhol = 600; sat_rhov = 120;
//         iapws95_sat_Trho_p(p, &T, &sat_rhol, &sat_rhov);
//         EXPECT_NEAR(625.,                 T);
//         EXPECT_NEAR(0.567090385e3,     sat_rhol);
//         EXPECT_NEAR(0.118290280e3,     sat_rhov);
//     }

// // #ifdef USE_QUAD
// //     {
// //         __float128 T, sat_rhol, sat_rhov;
// //         __float128 p = 0.698451167e-3*1e6;
// //         T=300;
// //         sat_rhol = 567;
// //         sat_rhov = 118;
// //         iapws95_sat_Trho_p_quad(p, &T, &sat_rhol, &sat_rhov);
// //         EXPECT_NEAR(275,                 (double)T);
// //         EXPECT_NEAR(0.999887406e3,     (double)sat_rhol);
// //         EXPECT_NEAR(0.550664919e-2,     (double)sat_rhov);
// //     }
// // #endif

// // //#ifdef USE_LONGDOUBLE
// //     {
// //         long double T, sat_rhol, sat_rhov;
// //         double p = 0.698451167e-3*1e6;
// //         T=300;
// //         sat_rhol = 567;
// //         sat_rhov = 118;
// //         iapws95_sat_Trho_p_long(p, &T, &sat_rhol, &sat_rhov);
// //         EXPECT_NEAR(275,                 T);
// //         EXPECT_NEAR(0.999887406e3,     sat_rhol);
// //         EXPECT_NEAR(0.550664919e-2,     sat_rhov);
// //     }
// // //#endif

// }

UTEST(iapws_sat92, 1)
{
    // printf("-> test IAPWS-sat-92\n");
    //---------------------------------------------------------------
    // IAPWS-SAT-92
    //---------------------------------------------------------------
    {
        double T=273.16;
        EXPECT_NEAR3(611.657, iapws92_sat_p_T(T));
        EXPECT_NEAR3(44.436693, iapws92_sat_dp_dT_T(T));
        EXPECT_NEAR3(999.789, iapws92_sat_rhol_T(T));
        EXPECT_NEAR3(0.00485426, iapws92_sat_rhov_T(T));
        EXPECT_NEAR3(-11.529101, iapws92_sat_alpha(T));
        EXPECT_NEAR3(0.611786, iapws92_sat_hl_T(T));
        EXPECT_NEAR3(2500.5e3, iapws92_sat_hv_T(T));
        EXPECT_NEAR0(-0.04, iapws92_sat_phi(T), 2e-1);
        EXPECT_NEAR3(0., iapws92_sat_sl_T(T));
        EXPECT_NEAR3(9.154e3, iapws92_sat_sg_T(T));

        T=373.1243;
        EXPECT_NEAR3(0.101325e6, iapws92_sat_p_T(T));
        EXPECT_NEAR3(3.616e3, iapws92_sat_dp_dT_T(T));
        EXPECT_NEAR3(958.365, iapws92_sat_rhol_T(T));
        EXPECT_NEAR3(0.597586, iapws92_sat_rhov_T(T));
        EXPECT_NEAR3(417.65e3, iapws92_sat_alpha(T));
        EXPECT_NEAR3(419.05e3, iapws92_sat_hl_T(T));
        EXPECT_NEAR3(2675.7e3, iapws92_sat_hv_T(T));
        EXPECT_NEAR3(1.303e3, iapws92_sat_phi(T));
        EXPECT_NEAR3(1.307e3, iapws92_sat_sl_T(T));
        EXPECT_NEAR3(7.355e3, iapws92_sat_sg_T(T));

        T=647.096;
        EXPECT_NEAR3(22.064e6, iapws92_sat_p_T(T));
        EXPECT_NEAR3(268e3, iapws92_sat_dp_dT_T(T));
        EXPECT_NEAR3(322., iapws92_sat_rhol_T(T));
        EXPECT_NEAR3(322., iapws92_sat_rhov_T(T));
        EXPECT_NEAR3(1548e3, iapws92_sat_alpha(T));
        EXPECT_NEAR3(2086.6e3, iapws92_sat_hl_T(T));
        EXPECT_NEAR3(2086.6e3, iapws92_sat_hv_T(T));
        EXPECT_NEAR3(3.578e3, iapws92_sat_phi(T));
        EXPECT_NEAR3(4.410e3, iapws92_sat_sl_T(T));
        EXPECT_NEAR3(4.410e3, iapws92_sat_sg_T(T));
    }
}

UTEST(iapws_melt11, 1)
{
    // printf("-> test IAPWS-Melt-2011\n");
    //---------------------------------------------------------------
    // IAPWS-Melt-2011
    //---------------------------------------------------------------
    {
        EXPECT_NEAR3(138.268e6, iapws11_melt_presssure_ice_Ih(260.));
        EXPECT_NEAR3(268.685e6, iapws11_melt_presssure_ice_III(254.));
        EXPECT_NEAR3(479.640e6, iapws11_melt_presssure_ice_V(265.));
        EXPECT_NEAR3(1356.76e6, iapws11_melt_presssure_ice_VI(320.));
        EXPECT_NEAR3(6308.71e6, iapws11_melt_presssure_ice_VII(550.));
        EXPECT_NEAR3(8.94735, iapws11_sublimation_pressure(230.));

        EXPECT_NEAR3(260., iapws11_melt_temperature(138.268e6));
        EXPECT_NEAR3(254., iapws11_melt_temperature(268.685e6));
        EXPECT_NEAR3(265., iapws11_melt_temperature(479.640e6));
        EXPECT_NEAR3(320., iapws11_melt_temperature(1356.76e6));
        EXPECT_NEAR3(550., iapws11_melt_temperature(6308.71e6));
        EXPECT_NEAR3(230., iapws11_melt_temperature(8.94735));


        // printf("Ih: %.9e ~ %.9e\n", melt_presssure_ice_Ih(273.16), melt_presssure_ice_Ih(251.165));
        // printf("III: %.9e ~ %.9e\n", melt_presssure_ice_III(251.165), melt_presssure_ice_III(256.164));
        // printf("V: %.9e ~ %.9e\n", melt_presssure_ice_V(256.164), melt_presssure_ice_V(273.31));
        // printf("VI: %.9e ~ %.9e\n", melt_presssure_ice_VI(273.31), melt_presssure_ice_VI(355));
        // printf("VII: %.9e ~ %.9e\n", melt_presssure_ice_VII(355), melt_presssure_ice_VII(715));
        // printf("s: %.9e ~ %.9e\n", sublimation_pressure(50), sublimation_pressure(273.16));
    }
}

UTEST(iapws_vis08, 1)
{
    // printf("-> test IAPWS-Viscosity-2008\n");
    //---------------------------------------------------------------
    // IAPWS-Viscosity-2008
    //---------------------------------------------------------------
    {
        EXPECT_NEAR3( 889.735100, iapws08_viscosity_rhoT( 998,     298.15)*1e6);
        EXPECT_NEAR3(1437.649467, iapws08_viscosity_rhoT(1200,     298.15)*1e6);
        EXPECT_NEAR3( 307.883622, iapws08_viscosity_rhoT(1000,     373.15)*1e6);
        EXPECT_NEAR3(  14.538324, iapws08_viscosity_rhoT(   1,     433.15)*1e6);
        EXPECT_NEAR3( 217.685358, iapws08_viscosity_rhoT(1000,     433.15)*1e6);
        EXPECT_NEAR3(  32.619287, iapws08_viscosity_rhoT(   1,    873.15)*1e6);
        EXPECT_NEAR3(  35.802262, iapws08_viscosity_rhoT( 100,     873.15)*1e6);
        EXPECT_NEAR3(  77.430195, iapws08_viscosity_rhoT( 600,     873.15)*1e6);
        EXPECT_NEAR3(  44.217245, iapws08_viscosity_rhoT(   1,     1173.15)*1e6);
        EXPECT_NEAR3(  47.640433, iapws08_viscosity_rhoT( 100,     1173.15)*1e6);
        EXPECT_NEAR3(  64.154608, iapws08_viscosity_rhoT( 400,     1173.15)*1e6);
    }
    {
        EXPECT_NEAR3( 25.520677, iapws08_viscosity_rhoT( 122,     647.35)*1e6);
        EXPECT_NEAR3( 31.337589, iapws08_viscosity_rhoT( 222,     647.35)*1e6);
        EXPECT_NEAR3( 36.228143, iapws08_viscosity_rhoT( 272,     647.35)*1e6);
        EXPECT_NEAR3( 42.961579, iapws08_viscosity_rhoT( 322,     647.35)*1e6);
        EXPECT_NEAR3( 45.688204, iapws08_viscosity_rhoT( 372,     647.35)*1e6);
        EXPECT_NEAR3( 49.436256, iapws08_viscosity_rhoT( 422,     647.35)*1e6);
    }
}

UTEST(iapws_thermal11, 1)
{
    // printf("-> test IAPWS-Thermal conductivity-2011\n");
    //---------------------------------------------------------------
    // IAPWS-Thermal conductivity-2011
    //---------------------------------------------------------------
    {// lambda2=0
        EXPECT_NEAR3(18.4341883, iapws11_thermal_conductivity_rhoT(   0, 298.15)*1e3);
        EXPECT_NEAR3(607.712868, iapws11_thermal_conductivity_rhoT( 998, 298.15)*1e3);
        EXPECT_NEAR3(799.038144, iapws11_thermal_conductivity_rhoT(1200, 298.15)*1e3);
        EXPECT_NEAR3(79.1034659, iapws11_thermal_conductivity_rhoT(   0, 873.15)*1e3);
    }

    {
        const double rhor = 322;
        const double Tr = 647.096;
        EXPECT_NEAR3(51.5764797, iapws11_lambda0bar(647.35/Tr));
        EXPECT_NEAR3(51.5764797, iapws11_lambda0bar(647.35/Tr));

        EXPECT_NEAR3(1.0068497, iapws11_lambda1bar(647.35/Tr,   1./rhor));
        EXPECT_NEAR3(2.1445173, iapws11_lambda1bar(647.35/Tr, 122./rhor));
        EXPECT_NEAR3(4.9681953, iapws11_lambda1bar(647.35/Tr, 322./rhor));

        EXPECT_NEAR3(0.00013000, iapws11_lambda2bar(647.35/Tr,   1./rhor));
        EXPECT_NEAR3(20.3162320, iapws11_lambda2bar(647.35/Tr, 122./rhor));
        EXPECT_NEAR3(1187.51354, iapws11_lambda2bar(647.35/Tr, 322./rhor));

        EXPECT_NEAR3(51.9298924, iapws11_thermal_conductivity_rhoT(   1, 647.35)*1e3);
        EXPECT_NEAR3(130.922885, iapws11_thermal_conductivity_rhoT( 122, 647.35)*1e3);
        EXPECT_NEAR3(367.787459, iapws11_thermal_conductivity_rhoT( 222, 647.35)*1e3);
        EXPECT_NEAR3(757.959776, iapws11_thermal_conductivity_rhoT( 272, 647.35)*1e3);
        EXPECT_NEAR3(1443.75556, iapws11_thermal_conductivity_rhoT( 322, 647.35)*1e3);
        EXPECT_NEAR3(650.319402, iapws11_thermal_conductivity_rhoT( 372, 647.35)*1e3);
        EXPECT_NEAR3(448.883487, iapws11_thermal_conductivity_rhoT( 422, 647.35)*1e3);
        EXPECT_NEAR3(600.961346, iapws11_thermal_conductivity_rhoT( 750, 647.35)*1e3);
    }
}
