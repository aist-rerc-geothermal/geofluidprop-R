
#ifndef _EOS_LINEAR_H_
#define _EOS_LINEAR_H_

#include "eos_args.h"

typedef struct
{
    double rho0;
    double p0;
    double beta_p;
    double T0;
    double beta_T;
    double cp;
    double vis;
} EOS_LINEAR_VALUES;

void eos_linear_register();

void* eos_linear_create();

void eos_linear_free(void* eos);

double eos_linear_rho_pT(void* eos, EOS_ARGS* args);

double eos_linear_rho_ph(void* eos, EOS_ARGS* args);

double eos_linear_drho_dh(void* eos, EOS_ARGS* args);

double eos_linear_drho_dT(void* eos, EOS_ARGS* args);

double eos_linear_drho_dp(void* eos, EOS_ARGS* args);

double eos_linear_h_rhoT(void* eos, EOS_ARGS* args);

double eos_linear_h_pT(void* eos, EOS_ARGS* args);

double eos_linear_dh_dT(void* eos, EOS_ARGS* args);

double eos_linear_dh_dp(void* eos, EOS_ARGS* args);

double eos_linear_u_pT(void* eos, EOS_ARGS* args);

double eos_linear_p_rhoT(void* eos, EOS_ARGS* args);

double eos_linear_T_ph(void* eos, EOS_ARGS* args);

double eos_linear_sat_rhol_T(void* eos, EOS_ARGS* args);

double eos_linear_sat_rhov_T(void* eos, EOS_ARGS* args);

double eos_linear_wv_ph(void* eos, EOS_ARGS* args);

double eos_linear_vis_rhoT(void* eos, EOS_ARGS* args);

double eos_linear_vis_pT(void* eos, EOS_ARGS* args);

#endif // _EOS_LINEAR_H_
