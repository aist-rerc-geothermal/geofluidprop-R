#ifndef _EOS_IAPWS95SBTL_PH_H_
#define _EOS_IAPWS95SBTL_PH_H_

#include "eos_args.h"

void eos_iapws95sbtl_ph_register();

void* eos_iapws95sbtl_ph_create();

void eos_iapws95sbtl_ph_free(void* sbtl);

double eos_iapws95sbtl_ph_rho_ph(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_rho_pT(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_h_pT(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_T_ph(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_sat_p_T(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_sat_T_p(void* sbtl, EOS_ARGS* args);

void eos_iapws95sbtl_ph_sat_rholv_T(void* sbtl, EOS_ARGS* args, double* rhol, double* rhov);

void eos_iapws95sbtl_ph_sat_rholv_p(void* sbtl, EOS_ARGS* args, double* rhol, double* rhov);

void eos_iapws95sbtl_ph_sat_hlv_T(void* sbtl, EOS_ARGS* args, double* hl, double* hv);

void eos_iapws95sbtl_ph_sat_hlv_p(void* sbtl, EOS_ARGS* args, double* hl, double* hv);

double eos_iapws95sbtl_ph_wv_ph(void* sbtl, EOS_ARGS* args);

double eos_iapws95sbtl_ph_vis_rhoT(void* eos, EOS_ARGS* args);

#endif
