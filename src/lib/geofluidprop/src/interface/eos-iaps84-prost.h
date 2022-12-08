
#ifndef _EOS_IAPS84_PROST_H_
#define _EOS_IAPS84_PROST_H_

#include "eos_args.h"

void eos_iaps84_prost_register();

void* eos_iaps84_prost_create();

void eos_iaps84_prost_free(void* eos);

double eos_iaps84_prost_rho_pT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_rho_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_drho_dT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_drho_dp(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_h_rhoT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_h_pT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_dh_dT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_dh_dp(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_p_rhoT(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_T_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_p_T(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_T_p(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_rhol_T(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_rhov_T(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_hl_T(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_sat_hv_T(void* eos, EOS_ARGS* args);

void eos_iaps84_prost_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_iaps84_prost_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_iaps84_prost_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double* hv);

void eos_iaps84_prost_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double* hv);

double eos_iaps84_prost_wv_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_prost_vis_rhoT(void* eos, EOS_ARGS* args);

#endif // _EOS_IAPWS95_API_H_
