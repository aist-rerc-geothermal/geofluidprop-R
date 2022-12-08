
#ifndef _EOS_DRIESNER07_H_
#define _EOS_DRIESNER07_H_

#include "eos_args.h"

void eos_driesner07_register();

void* eos_driesner07_create();

void eos_driesner07_free(void* eos);

double eos_driesner07_rho_pTx(void* eos, EOS_ARGS* args);

double eos_driesner07_rho_phx(void* eos, EOS_ARGS* args);

double eos_driesner07_h_pTx(void* eos, EOS_ARGS* args);

double eos_driesner07_T_phx(void* eos, EOS_ARGS* args);

double eos_driesner07_wv_phx(void* eos, EOS_ARGS* args);

// double eos_driesner07_sat_p_Tx(void* eos, EOS_ARGS* args);

double eos_driesner07_vis_rhoTx(void* eos, EOS_ARGS* args);

#endif
