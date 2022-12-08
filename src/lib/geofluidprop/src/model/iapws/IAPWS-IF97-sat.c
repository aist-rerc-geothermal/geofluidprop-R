
#include "IAPWS-IF97.h"

#include <math.h>
#include <stdio.h>

const double SaAn[11] =
{
  0.0, // dummy
  0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02,
  0.12020824702470E+05, -0.32325550322333E+07,  0.14915108613530E+02,
 -0.48232657361591E+04,  0.40511340542057E+06, -0.23855557567849E+00,
  0.65017534844798E+03
};

// 273.15 - 1073.15K
double if97_sat_p_T(double temK)
{
    double teta = temK + SaAn[9] / (temK - SaAn[10]);

    double a = (teta + SaAn[1]) * teta + SaAn[2];
    double b = (SaAn[3] * teta + SaAn[4]) * teta + SaAn[5];
    double c = (SaAn[6] * teta + SaAn[7]) * teta + SaAn[8];

    double w = 2.0 * c / (-b + sqrt(b * b - 4.0 * a * c));

    return w * w * w * w;
}

// 0.001 - 100 MPa
double if97_sat_T_p(double pres)
{
    double beta = sqrt(sqrt(pres));

    double e = (beta + SaAn[3]) * beta + SaAn[6];
    double f = (SaAn[1] * beta + SaAn[4]) * beta + SaAn[7];
    double g = (SaAn[2] * beta + SaAn[5]) * beta + SaAn[8];
    double d = 2.0 * g / (-f - sqrt(f * f - 4.0 * e * g));

    return 0.5 * (SaAn[10] + d - sqrt((SaAn[10] + d) * (SaAn[10] + d) - 4.0 * (SaAn[9] + SaAn[10] * d)));
}
