
#include <stdio.h>
#include <stdlib.h>

#include "utest.h"

#include "eos.h"
#include "eos-iapws95.h"
#include "eos-linear.h"
#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-SAT-92.h"
#include "model/iapws/IAPWS-IF97.h"
#include "util/utility.h"

UTEST(eos, linear)
{
    {
        // constant case
        EOS* eos = eos_create(EOS_TYPE_LINEAR);
        ((EOS_LINEAR_VALUES*)eos->impl)->rho0 = 1e3;
        ((EOS_LINEAR_VALUES*)eos->impl)->beta_p = 0;
        ((EOS_LINEAR_VALUES*)eos->impl)->beta_T = 0;

        EOS_ARGS args;
        args.p = 1e6;
        EXPECT_NEAR3(1e3, eos_rho_pT(eos, &args));
        args.p = 2e6;
        EXPECT_NEAR3(1e3, eos_rho_pT(eos, &args));

        eos_free(eos);
    }
}


UTEST(eos, iapws95)
{
    {
        printf("-> test IAPWS95\n");
        EOS* eos = eos_create(EOS_TYPE_WATER_IAPWS95);
        EOS_ARGS args;

        //test
        args.p = 119780.;
        args.h = 2.75666e+06;
        EXPECT_NEAR3(0.634030804, eos_rho_ph(eos, &args));


        // liquid phase
        args.rho = 0.996556e3;
        args.T = 300.0;
        args.p = 0.992418352e-1*1e6;
        args.h = 112.652982e3;
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        args.p = 2.13969121e7;
        args.h = 1.68072000e+06;
        EXPECT_NEAR3(628.337803, eos_T_ph(eos, &args));
        EXPECT_NEAR3(586.959842, eos_rho_ph(eos, &args));

        // vapor phase
        args.rho = 0.01;
        args.T = 300.0;
        args.p = 0.00138362569*1e6;
        args.h = 2550.84999e3;
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        // supercrtical
        args.rho = 184.236786;
        args.T = 700.0;
        args.p = 30.e6;
        args.h = 2631.43982e3;
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        // two phase
        args.rho = 1.0e1;
        args.T = 300.0;
        args.p = 0.00353680675*1e6;
        args.h = 118.739374*1e3;
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(-1., eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        // saturated liquid
        args.rho = 996.513027;
        args.T = 300.0;
        args.p = 0.00353680675*1e6;
        args.h = 112.56486*1e3;
        EXPECT_NEAR3(args.p, eos_sat_p_T(eos, &args));
        EXPECT_NEAR3(args.T, eos_sat_T_p(eos, &args));
        EXPECT_NEAR3(args.rho, eos_sat_rhol_T(eos, &args));
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(-1., eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        args.T = 633.15;
        EXPECT_NEAR3(527.591629, eos_sat_rhol_T(eos, &args));

        // saturated vapor
        args.rho = 0.0255896737;
        args.T = 300.0;
        args.p = 0.00353680675*1e6;
        args.h = 2549.8541*1e3;
        EXPECT_NEAR3(args.p, eos_sat_p_T(eos, &args));
        EXPECT_NEAR3(args.T, eos_sat_T_p(eos, &args));
        EXPECT_NEAR3(args.rho, eos_sat_rhov_T(eos, &args));
        EXPECT_NEAR3(args.p, eos_p_rhoT(eos, &args));
        EXPECT_NEAR3(-1., eos_rho_pT(eos, &args));
        EXPECT_NEAR3(args.rho, eos_rho_ph(eos, &args));
        EXPECT_NEAR3(args.h, eos_h_rhoT(eos, &args));
        EXPECT_NEAR3(args.T, eos_T_ph(eos, &args));

        args.T = 633.15;
        EXPECT_NEAR3(143.898411, eos_sat_rhov_T(eos, &args));

        eos_free(eos);
    }
}

UTEST(eos, if97)
{
    {
        printf("-> test IF97\n");
        EOS* eos = eos_create(EOS_TYPE_WATER_IF97);
        EOS_ARGS args;
        args.p = 1.0e6;
        args.T = 400.0;
        EXPECT_NEAR3(937.9228368,  eos_rho_pT(eos, &args));
        // args.p = 700.;
        // args.T = 273.15;
        // EXPECT_NEAR3(937.9228368,  eos_h_pT(eos, &args));
        eos_free(eos);
    }
}

UTEST(eos, iaps84)
{
    {
        printf("-> test IAPS84\n");
        EOS* eos = eos_create(EOS_TYPE_WATER_IAPS84);
        EOS_ARGS args;
        args.rho = 937.9228368;
        args.p = 1.0e6;
        args.T = 400.0;
        EXPECT_NEAR3(937.9228368,  eos_rho_pT(eos, &args));
        EXPECT_NEAR3(1.0e6,     eos_p_rhoT(eos, &args));
        eos_free(eos);
    }
}

UTEST(eos, driesner07)
{
    {
        printf("-> test DRIESNER07\n");
        EOS* eos = eos_create(EOS_TYPE_H2ONaCl_DRIESNER07);
        EOS_ARGS args;
        args.p = 2.810849e7;
        args.T = 673.15;
        args.x1 = 6.832051e-03;
        EXPECT_NEAR3(4.324106e+02,  eos_rho_pT(eos, &args));
        eos_free(eos);
    }
}

#ifdef USE_FREESTEAM
UTEST(eos, if97_freesteam)
{
    EOS* eos = eos_create(EOS_TYPE_WATER_IF97_FREESTEAM);
    {
        EOS_ARGS args;
        args.p = 1.0e6;
        args.T = 400.0;
        EXPECT_NEAR3(937.9228368,  eos_rho_pT(eos, &args));
        // args.p = 700.;
        // args.T = 273.15;
        // EXPECT_NEAR3(937.9228368,  eos_h_pT(eos, &args));
    }
    {
        EOS_ARGS args;
        args.rho = 0.996556e3;
        args.T = 300.0;
        double p = eos_p_rhoT(eos, &args);
        double h = eos_h_rhoT(eos, &args);
        EXPECT_NEAR3(args.rho, eos_rho_pT(eos, &args));
        args.p = p;
        EXPECT_NEAR3(eos_h_pT(eos, &args), h);
    }
    {
        EOS_ARGS args;
        args.rho = 998;
        args.T = 298.15;
        EXPECT_NEAR3(8.897351028326e+02,  eos_vis_rhoT(eos, &args)*1e6);
    }
    eos_free(eos);
}
#endif

#ifdef USE_PROST
UTEST(eos, iaps84_prost)
{
    {
        printf("-> test IAPS84_PROST\n");
        EOS* eos = eos_create(EOS_TYPE_WATER_IAPS84_PROST);
        EOS_ARGS args;
        args.rho = 937.9228368;
        args.p = 1.0e6;
        args.T = 400.0;
        EXPECT_NEAR3(937.9228368,  eos_rho_pT(eos, &args));
        EXPECT_NEAR3(1.0,     eos_p_rhoT(eos, &args)*1e-6);
        eos_free(eos);
    }
}

UTEST(eos, driesner07_prost)
{
    {
        printf("-> test DRIESNER07_PROST\n");
        EOS* eos = eos_create(EOS_TYPE_H2ONaCl_DRIESNER07_PROST);
        EOS_ARGS args;
        args.p = 2.810849e7;
        args.T = 673.15;
        args.x1 = 6.832051e-03;
        EXPECT_NEAR3(4.324106e+02,  eos_rho_pT(eos, &args));
        eos_free(eos);
    }
}
#endif
