
#include "rwrap_iapws95.h"

#include "model/iapws/IAPWS-95.h"


void R_iapws95_pc(double*out)
{
  *out = iapws95_get_pc();
}

void R_iapws95_Tc(double*out)
{
  *out = iapws95_get_Tc();
}

void R_iapws95_rho_pT(double*p, double*T, double*out)
{
  *out = iapws95_rho_pT(*p, *T);
}

