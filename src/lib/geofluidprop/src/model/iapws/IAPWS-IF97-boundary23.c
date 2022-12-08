// Auxiliary Equation for the Boundary between Regions 2 and 3
// p: 16.52916425 -100
// T: 623.15 - 863.15

#include "IAPWS-IF97.h"

#include <math.h>
#include <stdio.h>

const double B23N[5] =
{ 0.34805185628969E+3, -0.11671859879975E+1, 0.10192970039326E-2,
  0.57254459862746E+3, 0.13918839778870E+2};

double if97_boundary23_p_T(double temK)
{
    return (B23N[2] * temK + B23N[1]) * temK + B23N[0];
}

double if97_boundary23_T_p(double pres)
{
    return B23N[3] + sqrt((pres - B23N[4]) / B23N[2]);
}
