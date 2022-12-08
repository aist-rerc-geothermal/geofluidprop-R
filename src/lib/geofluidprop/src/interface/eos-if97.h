
#ifndef _EOS_IF97_H_
#define _EOS_IF97_H_

#include "eos_args.h"

void eos_if97_register();

void* eos_if97_create();

void eos_if97_free(void* eos);

double eos_if97_rho_ph(void* eos, EOS_ARGS* args);

double eos_if97_rho_pT(void* eos, EOS_ARGS* args);

double eos_if97_drho_dT(void* eos, EOS_ARGS* args);

double eos_if97_drho_dp(void* eos, EOS_ARGS* args);

double eos_if97_u_pT(void* eos, EOS_ARGS* args);

double eos_if97_h_rhoT(void* eos, EOS_ARGS* args);

double eos_if97_h_pT(void* eos, EOS_ARGS* args);

double eos_if97_dh_dT(void* eos, EOS_ARGS* args);

double eos_if97_dh_dp(void* eos, EOS_ARGS* args);

double eos_if97_p_rhoT(void* eos, EOS_ARGS* args);

double eos_if97_T_pu(void* eos, EOS_ARGS* args);

double eos_if97_T_ph(void* eos, EOS_ARGS* args);

double eos_if97_wv_ph(void* eos, EOS_ARGS* args);

double eos_if97_sat_p_T(void* eos, EOS_ARGS* args);

double eos_if97_sat_T_p(void* eos, EOS_ARGS* args);

double eos_if97_sat_rhol_T(void* eos, EOS_ARGS* args);

double eos_if97_sat_rhov_T(void* eos, EOS_ARGS* args);

void eos_if97_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_if97_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_if97_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double* hv);

void eos_if97_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double* hv);

double eos_if97_vis_rhoT(void* eos, EOS_ARGS* args);

#endif // _EOS_IF97_H_
