
#include <stdio.h>
#include <assert.h>

#include "eos.h"

void main()
{
    EOS* iapws95 = eos_create(EOS_TYPE_WATER_IAPWS95);
    EOS_ARGS args;
    args.p = 1e6; // [Pa]
    args.T = 273.15 + 10; // [K]
    double iapws95_rho = eos_rho_pT(iapws95, &args); // [kg/m^3]
    printf("# Water density at %g Pa and %g K\n", args.p, args.T);
    printf("IAPWS95   = %10.8e kg/m^3\n", iapws95_rho);
    eos_free(iapws95);

    EOS* if97 = eos_create(EOS_TYPE_WATER_IF97);
    double if97_rho = eos_rho_pT(if97, &args); // [kg/m^3]
    printf("IF97      = %10.8e kg/m^3\n", if97_rho);
    eos_free(if97);

#ifdef USE_PROST
    EOS* iaps84 = eos_create(EOS_TYPE_WATER_IAPS84_PROST);
    double iaps84_rho = eos_rho_pT(iaps84, &args); // [kg/m^3]
    printf("IAPS84    = %10.8e kg/m^3\n", iaps84_rho);
    eos_free(iaps84);
#endif

}

