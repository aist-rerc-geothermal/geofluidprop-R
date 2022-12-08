
#include "eos-if97.h"

#include <stdlib.h>
#include <stdio.h>

#include "model/iapws/IAPWS-IF97.h"
#include "model/iapws/IAPWS-Viscosity-85.h"

#include "eos_type_impl.h"

void eos_if97_register()
{
    const int index = EOS_TYPE_WATER_IF97;
    #define COMMAND(NAME)  eos_if97_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pT);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_ph);
    eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    eos_impl_calc_h_rhoT[index] = COMMAND(h_rhoT);
    eos_impl_calc_h_pT[index] = COMMAND(h_pT);
    eos_impl_calc_dh_dT[index] = COMMAND(dh_dT);
    eos_impl_calc_dh_dp[index] = COMMAND(dh_dp);
    eos_impl_calc_p_rhoT[index] = COMMAND(p_rhoT);
    eos_impl_calc_T_ph[index] = COMMAND(T_ph);
    eos_impl_calc_wv_ph[index] = COMMAND(wv_ph);
    eos_impl_calc_sat_rhol_T[index] = COMMAND(sat_rhol_T);
    eos_impl_calc_sat_rhov_T[index] = COMMAND(sat_rhov_T);
    eos_impl_calc_sat_rholv_T[index] = COMMAND(sat_rholv_T);
    eos_impl_calc_sat_rholv_p[index] = COMMAND(sat_rholv_p);
    eos_impl_calc_sat_hlv_T[index] = COMMAND(sat_hlv_T);
    eos_impl_calc_sat_hlv_p[index] = COMMAND(sat_hlv_p);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoT);
}

void* eos_if97_create()
{
    if97_init();
    return NULL;
}

void eos_if97_free(void* eos)
{
    free(eos);
}

double eos_if97_rho_pT(void* eos, EOS_ARGS* args)
{
    return if97_rho_pT(args->p*1e-6, args->T);
}

double eos_if97_rho_ph(void* eos, EOS_ARGS* args)
{
    return if97_rho_ph(args->p*1e-6, args->h*1e-3);
}

double eos_if97_u_pT(void* eos, EOS_ARGS* args)
{
    return if97_u_pT(args->p*1e-6, args->T)*1e3;
}

double eos_if97_h_pT(void* eos, EOS_ARGS* args)
{
    return if97_h_pT(args->p*1e-6, args->T)*1e3;
}

double eos_if97_drho_dT(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_drho_dp(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_h_rhoT(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_dh_dT(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_dh_dp(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_p_rhoT(void* eos, EOS_ARGS* args)
{
    return -1.0;
}

double eos_if97_T_pu(void* eos, EOS_ARGS* args)
{
    return if97_T_pu(args->p*1e-6, args->u*1e-3);
}

double eos_if97_T_ph(void* eos, EOS_ARGS* args)
{
    return if97_T_ph(args->p*1e-6, args->h*1e-3);
}

double eos_if97_wv_ph(void* eos, EOS_ARGS* args)
{
    return if97_vapor_quality_ph(args->p*1e-6, args->h*1e-3);
}

double eos_if97_sat_p_T(void* eos, EOS_ARGS* args)
{
    return if97_sat_p_T(args->T)*1e6;
}

double eos_if97_sat_T_p(void* eos, EOS_ARGS* args)
{
    return if97_sat_T_p(args->p*1e-6);
}

double eos_if97_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_if97_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    return 0;
}

void eos_if97_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    if97_sat_rho_T(args->T, rhol, rhov);
}

void eos_if97_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    if97_sat_rho_p(args->p*1e-6, rhol, rhov);
}

void eos_if97_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    double pMPa = if97_sat_p_T(args->T);
    if97_sat_h_p(pMPa, hl, hv);
}

void eos_if97_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    if97_sat_h_p(args->p*1e-6, hl, hv);
}

double eos_if97_vis_rhoT(void* eos, EOS_ARGS* args)
{
    return iapws85_viscosity_rhoT(args->rho, args->T);
}
