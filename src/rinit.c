
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
{"R_iapws95_get_rhoc", (DL_FUNC) &R_iapws95_get_rhoc, 1, argtype_real1},
{"R_iapws95_get_Tc", (DL_FUNC) &R_iapws95_get_Tc, 1, argtype_real1},
{"R_iapws95_get_pc", (DL_FUNC) &R_iapws95_get_pc, 1, argtype_real1},
{"R_iapws95_get_hc", (DL_FUNC) &R_iapws95_get_hc, 1, argtype_real1},
{"R_iapws95_get_Tt", (DL_FUNC) &R_iapws95_get_Tt, 1, argtype_real1},
{"R_iapws95_get_rholt", (DL_FUNC) &R_iapws95_get_rholt, 1, argtype_real1},
{"R_iapws95_get_rhovt", (DL_FUNC) &R_iapws95_get_rhovt, 1, argtype_real1},
{"R_iapws95_get_hlt", (DL_FUNC) &R_iapws95_get_hlt, 1, argtype_real1},
{"R_iapws95_get_hvt", (DL_FUNC) &R_iapws95_get_hvt, 1, argtype_real1},
{"R_iapws95_f_rhoT", (DL_FUNC) &R_iapws95_f_rhoT, 1, argtype_real1},
{"R_iapws95_g_rhoT", (DL_FUNC) &R_iapws95_g_rhoT, 1, argtype_real1},
{"R_iapws95_u_rhoT", (DL_FUNC) &R_iapws95_u_rhoT, 1, argtype_real1},
{"R_iapws95_h_rhoT", (DL_FUNC) &R_iapws95_h_rhoT, 1, argtype_real1},
{"R_iapws95_h_pT", (DL_FUNC) &R_iapws95_h_pT, 1, argtype_real1},
{"R_iapws95_s_rhoT", (DL_FUNC) &R_iapws95_s_rhoT, 1, argtype_real1},
{"R_iapws95_p_rhoT", (DL_FUNC) &R_iapws95_p_rhoT, 1, argtype_real1},
{"R_iapws95_cp_rhoT", (DL_FUNC) &R_iapws95_cp_rhoT, 1, argtype_real1},
{"R_iapws95_cv_rhoT", (DL_FUNC) &R_iapws95_cv_rhoT, 1, argtype_real1},
{"R_iapws95_rho_pT", (DL_FUNC) &R_iapws95_rho_pT, 1, argtype_real1},
{"R_iapws95_rho_ph", (DL_FUNC) &R_iapws95_rho_ph, 1, argtype_real1},
{"R_iapws95_T_rhop", (DL_FUNC) &R_iapws95_T_rhop, 1, argtype_real1},
{"R_iapws95_T_ph", (DL_FUNC) &R_iapws95_T_ph, 1, argtype_real1},
{"R_iapws95_rhoT_ph", (DL_FUNC) &R_iapws95_rhoT_ph, 1, argtype_real1},
{"R_iapws95_wv_ph", (DL_FUNC) &R_iapws95_wv_ph, 1, argtype_real1},
{"R_iapws95_sat_p_T", (DL_FUNC) &R_iapws95_sat_p_T, 1, argtype_real1},
{"R_iapws95_sat_T_p", (DL_FUNC) &R_iapws95_sat_T_p, 1, argtype_real1},
{"R_iapws95_sat_rhol_T", (DL_FUNC) &R_iapws95_sat_rhol_T, 1, argtype_real1},
{"R_iapws95_sat_rhol_p", (DL_FUNC) &R_iapws95_sat_rhol_p, 1, argtype_real1},
{"R_iapws95_sat_rhov_T", (DL_FUNC) &R_iapws95_sat_rhov_T, 1, argtype_real1},
{"R_iapws95_sat_rhov_p", (DL_FUNC) &R_iapws95_sat_rhov_p, 1, argtype_real1},
{"R_iapws95_sat_hl_p", (DL_FUNC) &R_iapws95_sat_hl_p, 1, argtype_real1},
{"R_iapws95_sat_hl_T", (DL_FUNC) &R_iapws95_sat_hl_T, 1, argtype_real1},
{"R_iapws95_sat_hv_p", (DL_FUNC) &R_iapws95_sat_hv_p, 1, argtype_real1},
{"R_iapws95_sat_hv_T", (DL_FUNC) &R_iapws95_sat_hv_T, 1, argtype_real1},
{"R_iapws95_sat_prho_T", (DL_FUNC) &R_iapws95_sat_prho_T, 1, argtype_real1},
{"R_iapws95_sat_Trho_p", (DL_FUNC) &R_iapws95_sat_Trho_p, 1, argtype_real1},
  
  {"R_iapws95_rho_pT",        (DL_FUNC) &R_iapws95_rho_pT,        3, argtype_real3},
  {"R_driesner07_H2O_NaCl_pc_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_pc_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_pc_T2",        (DL_FUNC) &R_driesner07_H2O_NaCl_pc_T2,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_xc_T",        (DL_FUNC) &R_driesner07_H2O_NaCl_xc_T,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_Tc_x",        (DL_FUNC) &R_driesner07_H2O_NaCl_Tc_x,        2, argtype_real2},
  {"R_driesner07_H2O_NaCl_phase_type_pTx",   (DL_FUNC) &R_driesner07_H2O_NaCl_phase_type_pTx,   4, argtype_real3i1},
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
