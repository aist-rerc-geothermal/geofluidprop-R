
#include "eos-iaps84.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#include "model/iapws/IAPS-84.h"

#include "eos_type_impl.h"

void eos_iaps84_register()
{
    const int index = EOS_TYPE_WATER_IAPS84;
    #define COMMAND(NAME)  eos_iaps84_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pT);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_ph);
    // eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    // eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    eos_impl_calc_h_rhoT[index] = COMMAND(h_rhoT);
    eos_impl_calc_h_pT[index] = COMMAND(h_pT);
    // eos_impl_calc_dh_dT[index] = COMMAND(dh_dT);
    // eos_impl_calc_dh_dp[index] = COMMAND(dh_dp);
    eos_impl_calc_p_rhoT[index] = COMMAND(p_rhoT);
    eos_impl_calc_T_ph[index] = COMMAND(T_ph);
    eos_impl_calc_wv_ph[index] = COMMAND(wv_ph);
    eos_impl_calc_sat_rhol_T[index] = COMMAND(sat_rhol_T);
    eos_impl_calc_sat_rhov_T[index] = COMMAND(sat_rhov_T);
    // eos_impl_calc_sat_rholv_T[index] = COMMAND(sat_rholv_T);
    // eos_impl_calc_sat_rholv_p[index] = COMMAND(sat_rholv_p);
    // eos_impl_calc_sat_hlv_T[index] = COMMAND(sat_hlv_T);
    // eos_impl_calc_sat_hlv_p[index] = COMMAND(sat_hlv_p);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoT);
}

void* eos_iaps84_create()
{
    return NULL;
}

void eos_iaps84_free(void* eos)
{
    free(eos);
}

double eos_iaps84_rho_pT(void* eos, EOS_ARGS* args)
{
    return iaps84_rho_pT(args->p*1e-6, args->T);
}

double eos_iaps84_rho_ph(void* eos, EOS_ARGS* args)
{
    return iaps84_rho_ph(args->p*1e-6, args->h);
}

double eos_iaps84_h_pT(void* eos, EOS_ARGS* args)
{
    double rho = eos_iaps84_rho_pT(eos, args);
    return iaps84_h_rhoT(rho, args->T);
}

// double eos_iaps84_drho_dT(void* eos, EOS_ARGS* args)
// {
//     assert(false);
// }

// double eos_iaps84_drho_dp(void* eos, EOS_ARGS* args)
// {
//     assert(false);
// }

double eos_iaps84_h_rhoT(void* eos, EOS_ARGS* args)
{
    return iaps84_h_rhoT(args->rho, args->T);
}

// double eos_iaps84_dh_dT(void* eos, EOS_ARGS* args)
// {
//     assert(false);
// }

// double eos_iaps84_dh_dp(void* eos, EOS_ARGS* args)
// {
//     assert(false);
// }

double eos_iaps84_p_rhoT(void* eos, EOS_ARGS* args)
{
    return iaps84_p_rhoT(args->rho, args->T)*1e6;
}

double eos_iaps84_T_ph(void* eos, EOS_ARGS* args)
{
    return iaps84_T_ph(args->p*1e-6, args->h);
}

double eos_iaps84_wv_ph(void* eos, EOS_ARGS* args)
{
    return iaps84_wv_ph(args->p*1e-6, args->h);
}

double eos_iaps84_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    return iaps84_sat_rhol_T(args->T);
}

double eos_iaps84_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    return iaps84_sat_rhov_T(args->T);
}

double eos_iaps84_vis_rhoT(void* eos, EOS_ARGS* args)
{
    return iaps84_viscosity_rhoT(args->rho, args->T);
}
