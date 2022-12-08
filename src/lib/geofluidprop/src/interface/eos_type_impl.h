
#ifndef _EOS_TYPE_IMPL_H__
#define _EOS_TYPE_IMPL_H__

#include <stdlib.h>

#include "eos_type.h"
#include "eos_args.h"

void eos_impl_init_register();

typedef void (*EOS_IMPL_REGISTER)();
typedef void* (*EOS_IMPL_CREATE)();
typedef void(*EOS_IMPL_FREE)(void*);
typedef double(*EOS_IMPL_CALC_RHO_PT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_RHO_PH)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_DRHO_DT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_DRHO_DH)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_DRHO_DP)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_U_PT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_H_RHOT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_H_PT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_DH_DT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_DH_DP)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_P_RHOT)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_T_PU)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_T_PH)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_WV_PH)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_SAT_P_T)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_SAT_T_P)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_SAT_RHOL_T)(void* eos, EOS_ARGS* args);
typedef double(*EOS_IMPL_CALC_SAT_RHOV_T)(void* eos, EOS_ARGS* args);
typedef void(*EOS_IMPL_CALC_SAT_RHOLV_T)(void* eos, EOS_ARGS* args, double*, double*);
typedef void(*EOS_IMPL_CALC_SAT_RHOLV_P)(void* eos, EOS_ARGS* args, double*, double*);
typedef void(*EOS_IMPL_CALC_SAT_HLV_T)(void* eos, EOS_ARGS* args, double*, double*);
typedef void(*EOS_IMPL_CALC_SAT_HLV_P)(void* eos, EOS_ARGS* args, double*, double*);
typedef double(*EOS_IMPL_CALC_VIS_RHOT)(void* eos, EOS_ARGS* args);

#define EOS_IMPL_ARRAY_MAX_SIZE 10

EOS_IMPL_REGISTER eos_impl_register[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CREATE eos_impl_create[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_FREE eos_impl_free[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_RHO_PT eos_impl_calc_rho_pT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_RHO_PH eos_impl_calc_rho_ph[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_DRHO_DT eos_impl_calc_drho_dh[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_DRHO_DT eos_impl_calc_drho_dT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_DRHO_DP eos_impl_calc_drho_dp[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_U_PT eos_impl_calc_u_pT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_H_RHOT eos_impl_calc_h_rhoT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_H_PT eos_impl_calc_h_pT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_DH_DT eos_impl_calc_dh_dT[EOS_IMPL_ARRAY_MAX_SIZE];
EOS_IMPL_CALC_DH_DP eos_impl_calc_dh_dp[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_P_RHOT eos_impl_calc_p_rhoT[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_T_PU eos_impl_calc_T_pu[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_T_PH eos_impl_calc_T_ph[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_WV_PH eos_impl_calc_wv_ph[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_P_T eos_impl_calc_sat_p_T[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_T_P eos_impl_calc_sat_T_p[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_RHOL_T eos_impl_calc_sat_rhol_T[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_RHOV_T eos_impl_calc_sat_rhov_T[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_RHOLV_T eos_impl_calc_sat_rholv_T[EOS_IMPL_ARRAY_MAX_SIZE];


EOS_IMPL_CALC_SAT_RHOLV_P eos_impl_calc_sat_rholv_p[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_RHOLV_T eos_impl_calc_sat_hlv_T[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_SAT_RHOLV_P eos_impl_calc_sat_hlv_p[EOS_IMPL_ARRAY_MAX_SIZE];

EOS_IMPL_CALC_VIS_RHOT eos_impl_calc_vis_rhoT[EOS_IMPL_ARRAY_MAX_SIZE];

#endif

