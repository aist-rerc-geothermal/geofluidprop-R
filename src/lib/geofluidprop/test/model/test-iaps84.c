
#include <stdio.h>
#include <math.h>

#include "utest.h"

#include "util/utility.h"

#include "model/iapws/IAPS-84.h"

UTEST(iaps84, rho_pT)
{
    double p = 0.1; //MPa
    double T = 15.0;
    EXPECT_NEAR3(999.137, iaps84_rho_pT(p, to_K(T)));

    T = 300;
    p = 0.05;
    EXPECT_NEAR3(9.965397e2, iaps84_rho_pT(p, T));

    T = 400;
    p = 0.5;
    EXPECT_NEAR3(9.37667e2, iaps84_rho_pT(p, T));
}


UTEST(iaps84, T_rhop)
{
    {
        double T = 300;
        double p = 0.05;
        double rho = iaps84_rho_pT(p, T); //996.53338910044;
        EXPECT_NEAR3(T,             iaps84_T_rhop(rho, p));

        T = 400;
        p = 0.5;
        rho = iaps84_rho_pT(p, T); //937.617
        EXPECT_NEAR3(T,             iaps84_T_rhop(rho, p));
    }
}

UTEST(iaps84, rho_ph)
{
    {
        double rhoe = 996.534;
        double p = 0.05024;
        double h = 112.608e3;
        double rho = iaps84_rho_ph(p, h);
        EXPECT_NEAR3(rhoe, rho);

        rhoe = 937.617;
        p = 0.50058;
        h = 533.127e3;
        rho = iaps84_rho_ph(p, h);
        EXPECT_NEAR3(rhoe, rho);
    }
}

UTEST(iaps84, sat_prholv)
{
    {
        double T = 275;
        double sat_p = 1, sat_rhol = 1e3, sat_rhov=1e-3;
        iaps84_sat_prholv_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR3(0.698451167e-3,     sat_p);
        EXPECT_NEAR3(0.999887406e3,     sat_rhol);
        EXPECT_NEAR3(0.550664919e-2,     sat_rhov);

        T = 450;
        iaps84_sat_prholv_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR3(0.932203564,         sat_p);
        EXPECT_NEAR3(0.890341250e3,     sat_rhol);
        EXPECT_NEAR3(0.481200360e1,     sat_rhov);

        T = 625;
        iaps84_sat_prholv_T(T, &sat_p, &sat_rhol, &sat_rhov);
        EXPECT_NEAR3(0.169082693e2,     sat_p);
        EXPECT_NEAR3(0.567090385e3,     sat_rhol);
        EXPECT_NEAR3(0.118155691e3,     sat_rhov);
    }
}
