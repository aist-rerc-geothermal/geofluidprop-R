#ifndef _EOS_IAPWS95_H_
#define _EOS_IAPWS95_H_

#include "eos_args.h"

void eos_iapws95_register();

void* eos_iapws95_create();

void eos_iapws95_free(void* eos);

double eos_iapws95_rho_pT(void* eos, EOS_ARGS* args);

double eos_iapws95_rho_ph(void* eos, EOS_ARGS* args);

double eos_iapws95_drho_dh(void* eos, EOS_ARGS* args);

double eos_iapws95_drho_dT(void* eos, EOS_ARGS* args);

double eos_iapws95_drho_dp(void* eos, EOS_ARGS* args);

double eos_iapws95_h_rhoT(void* eos, EOS_ARGS* args);

double eos_iapws95_h_pT(void* eos, EOS_ARGS* args);

double eos_iapws95_dh_dT(void* eos, EOS_ARGS* args);

double eos_iapws95_dh_dp(void* eos, EOS_ARGS* args);

double eos_iapws95_p_rhoT(void* eos, EOS_ARGS* args);

double eos_iapws95_T_ph(void* eos, EOS_ARGS* args);

double eos_iapws95_wv_ph(void* eos, EOS_ARGS* args);

double eos_iapws95_sat_p_T(void* eos, EOS_ARGS* args);
#ifdef USE_QUAD
double eos_iapws95_sat_p_T_quad(void* eos, EOS_ARGS* args);
#endif

double eos_iapws95_sat_T_p(void* eos, EOS_ARGS* args);

double eos_iapws95_sat_rhol_T(void* eos, EOS_ARGS* args);

double eos_iapws95_sat_rhov_T(void* eos, EOS_ARGS* args);

void eos_iapws95_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#ifdef USE_LONGDOUBLE
void eos_iapws95_sat_rholv_T_long(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#endif
#ifdef USE_QUAD
void eos_iapws95_sat_rholv_T_quad(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#endif

void eos_iapws95_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#ifdef USE_LONGDOUBLE
void eos_iapws95_sat_rholv_p_long(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#endif
#ifdef USE_QUAD
void eos_iapws95_sat_rholv_p_quad(void* eos, EOS_ARGS* args, double* rhol, double* rhov);
#endif

double eos_iapws95_sat_p_hl(void* eos, EOS_ARGS* args);

#ifdef USE_QUAD
double eos_iapws95_sat_p_hl_quad(void* eos, EOS_ARGS* args);
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_p_hl_long(void* eos, EOS_ARGS* args);

double eos_iapws95_sat_T_hl_long(void* eos, EOS_ARGS* args);
#endif

double eos_iapws95_sat_p_hv(void* eos, EOS_ARGS* args);

#ifdef USE_QUAD
double eos_iapws95_sat_p_hv_quad(void* eos, EOS_ARGS* args);
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_p_hv_long(void* eos, EOS_ARGS* args);
#endif

double eos_iapws95_sat_T_hv(void* eos, EOS_ARGS* args);

#ifdef USE_QUAD
double eos_iapws95_sat_T_hv_quad(void* eos, EOS_ARGS* args);
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_T_hv_long(void* eos, EOS_ARGS* args);
#endif

void eos_iapws95_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double*hv);
void eos_iapws95_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double*hv);

double eos_iapws95_h_crit(void* eos);

double eos_iapws95_hl_tr(void* eos);
double eos_iapws95_hv_tr(void* eos);
double eos_iapws95_rhol_tr(void* eos);
double eos_iapws95_rhov_tr(void* eos);

double eos_iapws95_vis_rhoT(void* eos, EOS_ARGS* args);


#endif // _EOS_IAPWS95_API_H_
