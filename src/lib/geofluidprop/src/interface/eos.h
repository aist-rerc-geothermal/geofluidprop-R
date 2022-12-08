
#ifndef _EOS_H__
#define _EOS_H__

#include "eos_type.h"
#include "eos_args.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * create an EOS object with the given EoS type
 *
 * @param eos_type
 * @return EOS*
 */
EOS* eos_create(int eos_type);

/**
 * create an EOS object with the given EoS name
 *
 * @return EOS*
 */
EOS* eos_create_by_name(const char*);

/**
 * free an EoS object
 *
 */
void eos_free(EOS*);

////double eos_p(EOS* eos, EOS_ARGS* args);

double eos_rho_pT(EOS* eos, EOS_ARGS* args);

double eos_rho_ph(EOS* eos, EOS_ARGS* args);

double eos_drho_dT(EOS* eos, EOS_ARGS* args);

double eos_drho_dp(EOS* eos, EOS_ARGS* args);

double eos_drho_dh(EOS* eos, EOS_ARGS* args);

double eos_u_phrho(EOS* eos, EOS_ARGS* args);

double eos_u_pT(EOS* eos, EOS_ARGS* args);

double eos_h_rhoT(EOS* eos, EOS_ARGS* args);

double eos_h_pT(EOS* eos, EOS_ARGS* args);

double eos_dh_dT(EOS* eos, EOS_ARGS* args);

double eos_dh_dp(EOS* eos, EOS_ARGS* args);

double eos_p_rhoT(EOS* eos, EOS_ARGS* args);

double eos_T_pu(EOS* eos, EOS_ARGS* args);

double eos_T_ph(EOS* eos, EOS_ARGS* args);

double eos_wv_ph(EOS* eos, EOS_ARGS* args);

double eos_sat_p_T(EOS* eos, EOS_ARGS* args);

double eos_sat_T_p(EOS* eos, EOS_ARGS* args);

double eos_sat_rhol_T(EOS* eos, EOS_ARGS* args);

double eos_sat_rhov_T(EOS* eos, EOS_ARGS* args);

void eos_sat_rholv_T(EOS* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_sat_rholv_p(EOS* eos, EOS_ARGS* args, double* rhol, double* rhov);

void eos_sat_hlv_T(EOS* eos, EOS_ARGS* args, double* hl, double* hv);

void eos_sat_hlv_p(EOS* eos, EOS_ARGS* args, double* hl, double* hv);

double eos_vis_rhoT(EOS* eos, EOS_ARGS* args);

#ifdef __cplusplus
}
#endif

#endif

