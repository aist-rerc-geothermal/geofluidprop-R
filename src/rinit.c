
#include <stdlib.h> // for NULL

#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>

#include "rwrap_iapws95.h"
#include "rwrap_driesner07.h"
#include "rwrap_klyukinetal2017.h"


static R_NativePrimitiveArgType argtype_real1[] = {
      REALSXP
};
static R_NativePrimitiveArgType argtype_real2[] = {
      REALSXP, REALSXP
};
static R_NativePrimitiveArgType argtype_real3[] = {
      REALSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType argtype_real4[] = {
      REALSXP, REALSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType argtype_real3i1[] = {
      REALSXP, REALSXP, REALSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
  {"R_iapws95_pc_T",        (DL_FUNC) &R_iapws95_pc,        1, argtype_real1},
  {"R_iapws95_Tc_x",        (DL_FUNC) &R_iapws95_Tc,        1, argtype_real1},
  {"R_iapws95_rho_pT",        (DL_FUNC) &R_iapws95_rho_pT,        3, argtype_real3},
  {"R_driesner07_H2O_NaCl_pc_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_pc_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_pc_T2",        (DL_FUNC) &R_driesner07_H2O_NaCl_pc_T2,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_xc_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_xc_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_Tc_x",        (DL_FUNC) &R_driesner07_H2O_NaCl_Tc_x,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_phase_pTx",   (DL_FUNC) &R_driesner07_H2O_NaCl_phase_pTx,   4, argtype_real3i1},
  {"R_driesner07_H2O_NaCl_VL_xl_pT",    (DL_FUNC) &R_driesner07_H2O_NaCl_VL_xl_pT,    3, argtype_real3},
  {"R_driesner07_H2O_NaCl_VL_xv_pT",    (DL_FUNC) &R_driesner07_H2O_NaCl_VL_xv_pT,    3, argtype_real3},
  {"R_driesner07_H2O_NaCl_VL_rhol_pT",    (DL_FUNC) &R_driesner07_H2O_NaCl_VL_rhol_pT,    3, argtype_real3},
  {"R_driesner07_H2O_NaCl_VL_rhov_pT",    (DL_FUNC) &R_driesner07_H2O_NaCl_VL_rhov_pT,    3, argtype_real3},
  {"R_driesner07_H2O_NaCl_mass_fraction_vapor",    (DL_FUNC) &R_driesner07_H2O_NaCl_mass_fraction_vapor,    4, argtype_real4},
  {"R_driesner07_H2O_NaCl_volume_fraction_vapor",    (DL_FUNC) &R_driesner07_H2O_NaCl_volume_fraction_vapor,    4, argtype_real4},
  {"R_driesner07_H2O_NaCl_VLH_p_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_VLH_p_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_VLH_xl_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_VLH_xl_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_VLH_xv_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_VLH_xv_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_LH_xl_pT",        (DL_FUNC) &R_driesner07_H2O_NaCl_LH_xl_pT,        3, argtype_real3},
  {"R_driesner07_H2O_NaCl_rho_pTx",          (DL_FUNC) &R_driesner07_H2O_NaCl_rho_pTx,          4, argtype_real4},
  {"R_driesner07_H2O_NaCl_rho_singlephase_pTx",          (DL_FUNC) &R_driesner07_H2O_NaCl_rho_singlephase_pTx,          4, argtype_real4},
  {"R_driesner07_H2O_NaCl_singlephase_h_pTx",            (DL_FUNC) &R_driesner07_H2O_NaCl_singlephase_h_pTx,            4, argtype_real4},
  {"R_klyukinetal2017_H2O_NaCl_viscosity_rhoTx",  (DL_FUNC) &R_klyukinetal2017_H2O_NaCl_viscosity_rhoTx,  4, argtype_real4},
  {NULL, NULL, 0, NULL}
};

void R_init_geofluidprop(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
//  R_forceSymbols(dll, TRUE);
}
