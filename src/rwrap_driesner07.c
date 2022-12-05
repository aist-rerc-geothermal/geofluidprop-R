
#include "rwrap_driesner07.h"

#include "model/driesner07/Driesner2007_H2ONaCl.h"


void R_driesner07_H2O_NaCl_pc_T(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_pc_T(*T);
}

void R_driesner07_H2O_NaCl_pc_T2(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_pc_T2(*T, true);
}

void R_driesner07_H2O_NaCl_Tc_x(double*x, double*out)
{
  *out = driesner07_H2O_NaCl_Tc_x(*x);
}

void R_driesner07_H2O_NaCl_xc_T(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_xc_T(*T);
}

void R_driesner07_H2O_NaCl_phase_pTx(double*p, double*T, double*x, int*out)
{
  *out = driesner07_H2O_NaCl_phase_type(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_VL_xl_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VL_xl_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_VL_xv_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VL_xv_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_VL_rhov_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VL_rhov_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_VL_rhol_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VL_rhol_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_LH_xl_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_LH_xl_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_VH_xv_pT(double*p, double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VH_xv_pT(*p, *T);
}

void R_driesner07_H2O_NaCl_VLH_p_T(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VLH_p_T(*T);
}

void R_driesner07_H2O_NaCl_VLH_xl_T(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VLH_xl_T(*T);
}

void R_driesner07_H2O_NaCl_VLH_xv_T(double*T, double*out)
{
  *out = driesner07_H2O_NaCl_VLH_xv_T(*T);
}

void R_driesner07_H2O_NaCl_mass_fraction_liquid(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_mass_fraction_liquid(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_mass_fraction_vapor(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_mass_fraction_vapor(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_volume_fraction_liquid(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_volume_fraction_liquid(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_volume_fraction_vapor(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_volume_fraction_vapor(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_rho_pTx(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_rho_pTx(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_singlephase_rho_pTx(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_rho_singlephase_pTx(*p, *T, *x);
}

void R_driesner07_H2O_NaCl_singlephase_h_pTx(double*p, double*T, double*x, double*out)
{
  *out = driesner07_H2O_NaCl_h_singlephase_pTx(*p, *T, *x);
}
