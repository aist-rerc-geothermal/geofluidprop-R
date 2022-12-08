#include <stdio.h>
#include "eos.h"

void main()
{
    EOS* eos = eos_create(EOS_TYPE_WATER_IAPWS95);

    EOS_ARGS args;
    args.p = 1e6; // [Pa]
    args.T = 273.15 + 10; // [K]
    double rho = eos_rho_pT(eos, &args); // [kg/m^3]

    printf("Water density at %g Pa and %g K is %g kg/m^3\n", args.p, args.T, rho);

    eos_free(eos);
}

