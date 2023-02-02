
#ifndef _RWRAP_DRIESNER07_H_
#define _RWRAP_DRIESNER07_H_

void R_driesner07_H2O_NaCl_pc_T(double*T, double*out);
void R_driesner07_H2O_NaCl_pc_T2(double*T, double*out);
void R_driesner07_H2O_NaCl_Tc_x(double*x, double*out);
void R_driesner07_H2O_NaCl_xc_T(double*T, double*out);
void R_driesner07_H2O_NaCl_phase_pTx(double*p, double*T, double*x, int*out);
void R_driesner07_H2O_NaCl_VL_xl_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_VL_xv_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_VL_rhov_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_VL_rhol_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_LH_xl_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_VH_xv_pT(double*p, double*T, double*out);
void R_driesner07_H2O_NaCl_VLH_p_T(double*T, double*out);
void R_driesner07_H2O_NaCl_VLH_xl_T(double*T, double*out);
void R_driesner07_H2O_NaCl_VLH_xv_T(double*T, double*out);
void R_driesner07_H2O_NaCl_mass_fraction_liquid(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_mass_fraction_vapor(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_volume_fraction_liquid(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_volume_fraction_vapor(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_rho_pTx(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_rho_singlephase_pTx(double*p, double*T, double*x, double*out);
void R_driesner07_H2O_NaCl_singlephase_h_pTx(double*p, double*T, double*x, double*out);

void R_driesner07_NaCl_H_rho_pT(double* p_Pa, double* T_K, double* out);
void R_driesner07_NaCl_H_h_pT(double* p_Pa, double* T_K, double* out);
void R_driesner07_NaCl_H_cp_pT(double* p_Pa, double* T_K, double* out);

#endif
