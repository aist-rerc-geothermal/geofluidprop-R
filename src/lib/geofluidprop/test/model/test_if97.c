
#include <stdio.h>
#include <math.h>

#include "util/utility.h"

#include "model/iapws/IAPWS-IF97.h"

#include "utest.h"

UTEST(if97, boundary23)
{
    // printf("-> test Auxiliary Equation for the Boundary between Regions 2 and 3\n");
    double T = 0.623150000e3;
    double p = 0.165291643e2;
    EXPECT_NEAR(p, if97_boundary23_p_T(T));
    EXPECT_NEAR(T, if97_boundary23_T_p(p));
}


UTEST(if97, sat)
{
    // printf("-> test saturation curve\n");
    double p, T;
    T = 300;
    p = 0.353658941e-2;
    EXPECT_NEAR(p, if97_sat_p_T(T));
    T = 500;
    p = 0.263889776e1;
    EXPECT_NEAR(p, if97_sat_p_T(T));
    T = 600;
    p = 0.123443146e2;
    EXPECT_NEAR(p, if97_sat_p_T(T));

    p = 0.1;
    T = 0.372755919e3;
    EXPECT_NEAR(T, if97_sat_T_p(p));
    p = 1.0;
    T = 0.453035632e3;
    EXPECT_NEAR(T, if97_sat_T_p(p));
    p = 10.0;
    T = 0.584149488e3;
    EXPECT_NEAR(T, if97_sat_T_p(p));
}

UTEST(if97, region1)
{
    // printf("-> test Region 1\n");
    double T = 300;
    double p = 3.0;
    double v = 0.100215168e-2;
    double h = 0.115331273e3;
    double cp = 0.417301218e1;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 300;
    p = 80.0;
    v = 0.971180894e-3;
    h = 0.184142828e3;
    cp = 0.401008987e1;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 500;
    p = 3.0;
    v = 0.120241800e-2;
    h = 0.975542239e3;
    cp = 0.465580682e1;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));
}

UTEST(if97, region2)
{
    // printf("-> test Region 2\n");
    double T = 300;
    double p = 0.0035;
    double v = 0.394913866e2;
    double h = 0.254991145e4;
    double cp = 0.191300162e1;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 700;
    p = 0.0035;
    v = 0.923015898e2;
    h = 0.333568375e4;
    cp = 0.208141274e1;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 700;
    p = 30.0;
    v = 0.542946619e-2;
    h = 0.263149474e4;
    cp = 0.103505092e2;
    EXPECT_NEAR(v, 1./if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));
}


UTEST(if97, region3)
{
    // printf("-> test Region 3\n");
    double T = 650;
    double rho = 500.0;
    double p = 0.255837018e2;
    double h = 0.186343019e4;
    double cp = 0.138935717e2;
    double cv = 0.319131787e1;
    EXPECT_NEAR(p, if97_p_vT_region3(1./rho, T));
    EXPECT_NEAR(rho, if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR0(cp, if97_cp_pT(p, T), 1e-3);
    EXPECT_NEAR(cv, if97_cv_rhoT_region3(rho, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 650;
    rho = 200.0;
    p = 0.222930643e2;
    h = 0.237512401e4;
    cp = 0.446579342e2;
    EXPECT_NEAR(p, if97_p_vT_region3(1./rho, T));
    EXPECT_NEAR0(rho, if97_rho_pT(p, T), 1.e-5);
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR0(cp, if97_cp_pT(p, T), 1.e-5);
    EXPECT_NEAR(T, if97_T_ph(p, h));

    T = 750;
    rho = 500.0;
    p = 0.783095639e2;
    h = 0.225868845e4;
    cp = 0.634165359e1;
    EXPECT_NEAR(p, if97_p_vT_region3(1./rho, T));
    EXPECT_NEAR(rho, if97_rho_pT(p, T));
    EXPECT_NEAR(h, if97_h_pT(p, T));
    EXPECT_NEAR(cp, if97_cp_pT(p, T));
    EXPECT_NEAR(T, if97_T_ph(p, h));
}

