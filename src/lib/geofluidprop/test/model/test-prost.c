
#include <stdio.h>
#include <math.h>

#include <steam4.h>

#include "util/utility.h"
#include "utest.h"

#include "model/iapws/IAPS-84.h"

UTEST(prost, call)
{
    Prop* prop0 = newProp('p', 'T', 0);
    double p, T;
    p = 0.1e6;
    T = to_K(15.0);
    water_tp(T, p, 0, 1.e-8, prop0);
    EXPECT_NEAR0(999.137, prop0->d, 1e-2);
    freeProp(prop0);
}

static double prost_rho_pT(double p, double T)
{
    Prop* prop0 = newProp('p', 'T', 0);
    water_tp(T, p, 0, 1.e-8, prop0);
    double rho = prop0->d;
    freeProp(prop0);
    return rho;
}

static double prost_h_rhoT(double rho, double T)
{
    Prop* prop = newProp('d', 'T', 0);
    water_td(T, rho, prop);
    double h = prop->h;
    freeProp(prop);
    return h;
}


UTEST(prost, compare_with_iaps84)
{
    double p=30e6, T=400+273.15;
    EXPECT_NEAR0(prost_rho_pT(p, T), iaps84_rho_pT(p*1e-6,T), 1e-5);
    EXPECT_NEAR0(prost_h_rhoT(prost_rho_pT(p, T), T), iaps84_h_rhoT(iaps84_rho_pT(p*1e-6,T),T), 1e-5);
}
