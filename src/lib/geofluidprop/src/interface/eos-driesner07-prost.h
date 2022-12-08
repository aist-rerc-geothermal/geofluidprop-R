
#ifndef _EOS_DRIESNER07_PROST_H_
#define _EOS_DRIESNER07_PROST_H_

#include "eos_args.h"

void eos_driesner07_prost_register();

void* eos_driesner07_prost_create();

void eos_driesner07_prost_free(void* eos);

double eos_driesner07_prost_rho_pTx(void* eos, EOS_ARGS* args);

double eos_driesner07_prost_rho_phx(void* eos, EOS_ARGS* args);

double eos_driesner07_prost_h_pTx(void* eos, EOS_ARGS* args);

double eos_driesner07_prost_T_phx(void* eos, EOS_ARGS* args);

double eos_driesner07_prost_wv_phx(void* eos, EOS_ARGS* args);

// double eos_driesner07_prost_sat_p_Tx(void* eos, EOS_ARGS* args);

double eos_driesner07_prost_vis_rhoTx(void* eos, EOS_ARGS* args);

#endif
