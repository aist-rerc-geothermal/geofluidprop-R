
#include "model/klyukinetal17/KlyukinEtAl2017.h"

void R_klyukinetal2017_H2O_NaCl_viscosity_rhoTx(double*rho, double*T, double*x, double*out)
{
  *out = klyukinetal2017_viscosity(*rho, *T, *x);
}

