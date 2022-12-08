
#ifndef _EOS_IAPS84_H_
#define _EOS_IAPS84_H_

#include "eos_args.h"

void eos_iaps84_register();

void* eos_iaps84_create();

void eos_iaps84_free(void* eos);

double eos_iaps84_rho_pT(void* eos, EOS_ARGS* args);

double eos_iaps84_rho_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_drho_dT(void* eos, EOS_ARGS* args);

double eos_iaps84_drho_dp(void* eos, EOS_ARGS* args);

double eos_iaps84_h_rhoT(void* eos, EOS_ARGS* args);

double eos_iaps84_h_pT(void* eos, EOS_ARGS* args);

double eos_iaps84_dh_dT(void* eos, EOS_ARGS* args);

double eos_iaps84_dh_dp(void* eos, EOS_ARGS* args);

double eos_iaps84_p_rhoT(void* eos, EOS_ARGS* args);

double eos_iaps84_T_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_sat_rhol_T(void* eos, EOS_ARGS* args);

double eos_iaps84_sat_rhov_T(void* eos, EOS_ARGS* args);

double eos_iaps84_wv_ph(void* eos, EOS_ARGS* args);

double eos_iaps84_vis_rhoT(void* eos, EOS_ARGS* args);

#endif // _EOS_IAPS84_H_
