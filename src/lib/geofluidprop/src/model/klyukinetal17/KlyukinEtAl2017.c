
#include "KlyukinEtAl2017.h"

#include <math.h>

#include "model/iapws/IAPWS-Viscosity-08.h"

static double get_equivalent_temperature(double T, double x)
{
    const double a1 = -35.9858;
    const double a2 = 0.80017;
    const double b1 = 1e-6;
    const double b2 = -0.05239;
    const double b3 = 1.32936;
    double x_a2 = pow(x, a2);
    double T_b2 = pow(T, b2);
    double e1 = a1 * x_a2;
    double e2 = 1. - b1*T_b2 - b3*x_a2*T_b2;
    double Ta = e1 + e2 * T;

    return Ta;
}

// x: mass fraction of NaCl
double klyukinetal2017_viscosity(double rho, double T, double x)
{
    return iapws08_viscosity_rhoT(rho, get_equivalent_temperature(T, x));
}
