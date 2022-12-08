
#include "eos-linear.h"

#include <stdlib.h>
#include <stdio.h>

#include "eos_type_impl.h"

void eos_linear_register()
{
    const int index = EOS_TYPE_LINEAR;
    #define COMMAND(NAME)  eos_linear_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pT);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_ph);
    eos_impl_calc_drho_dh[index] = COMMAND(drho_dh);
    eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    eos_impl_calc_h_rhoT[index] = COMMAND(h_rhoT);
    eos_impl_calc_h_pT[index] = COMMAND(h_pT);
    eos_impl_calc_dh_dT[index] = COMMAND(dh_dT);
    eos_impl_calc_dh_dp[index] = COMMAND(dh_dp);
    eos_impl_calc_p_rhoT[index] = COMMAND(p_rhoT);
    eos_impl_calc_T_ph[index] = COMMAND(T_ph);
    eos_impl_calc_wv_ph[index] = COMMAND(wv_ph);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoT);
}

void* eos_linear_create()
{
    EOS_LINEAR_VALUES* values = malloc(sizeof(EOS_LINEAR_VALUES));
    values->rho0 = 1e3;
    values->p0 = 1e6;
    values->beta_p = 4.5e-10;
    values->T0 = 20. + 273.15;
    values->beta_T = 207e-6;
    values->cp = 4180.0;
    values->vis = 1e-3;
    return values;
}

void eos_linear_free(void* data)
{
    free((EOS_LINEAR_VALUES*)data);
}

double eos_linear_rho_pT(void* eos, EOS_ARGS* args)
{
    EOS_LINEAR_VALUES* obj = (EOS_LINEAR_VALUES*)eos;
    return obj->rho0 * (1.0 + (args->p - obj->p0)*obj->beta_p - (args->T - obj->T0)*obj->beta_T);
}

double eos_linear_rho_ph(void* eos, EOS_ARGS* args)
{
    args->T = args->h /  ((EOS_LINEAR_VALUES*)eos)->cp;
    return eos_linear_rho_pT(eos, args);
}

double eos_linear_u_pT(void* eos, EOS_ARGS* args)
{
    EOS_LINEAR_VALUES* obj = (EOS_LINEAR_VALUES*)eos;
    double h = eos_linear_h_pT(eos, args);
    double rho = eos_linear_rho_pT(eos, args);
    return (h + args->p / rho);
}

double eos_linear_h_pT(void* eos, EOS_ARGS* args)
{
    return ((EOS_LINEAR_VALUES*)eos)->cp*args->T;
}

double eos_linear_drho_dh(void* eos, EOS_ARGS* args)
{
    double drho_dT = eos_linear_drho_dT(eos, args);
    double dh_dT = ((EOS_LINEAR_VALUES*)eos)->cp;
    return (drho_dT/dh_dT); //drho/dt /cp
}

double eos_linear_drho_dT(void* eos, EOS_ARGS* args)
{
    EOS_LINEAR_VALUES* obj = (EOS_LINEAR_VALUES*)eos;
    return (- obj->rho0 * obj->beta_T);
}

double eos_linear_drho_dp(void* eos, EOS_ARGS* args)
{
    EOS_LINEAR_VALUES* obj = (EOS_LINEAR_VALUES*)eos;
    return (obj->rho0 * obj->beta_p);
}

double eos_linear_h_rhoT(void* eos, EOS_ARGS* args)
{
    return ((EOS_LINEAR_VALUES*)eos)->cp*args->T;
}

double eos_linear_dh_dT(void* eos, EOS_ARGS* args)
{
    return ((EOS_LINEAR_VALUES*)eos)->cp;
}

double eos_linear_dh_dp(void* eos, EOS_ARGS* args)
{
    return 0;
}

double eos_linear_p_rhoT(void* eos, EOS_ARGS* args)
{
    EOS_LINEAR_VALUES* obj = (EOS_LINEAR_VALUES*)eos;
    //
    double p = obj->p0 + 1./obj->beta_p * (args->rho / obj->rho0 - 1.0 + (args->T - obj->T0)*obj->beta_T);
    return p;
}

double eos_linear_T_ph(void* eos, EOS_ARGS* args)
{
    return args->h / ((EOS_LINEAR_VALUES*)eos)->cp;
}

double eos_linear_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    return -1;
}

double eos_linear_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    return -1;
}

double eos_linear_wv_ph(void* eos, EOS_ARGS* args)
{
    return 0.0;
}

double eos_linear_vis_rhoT(void* eos, EOS_ARGS* args)
{
    return ((EOS_LINEAR_VALUES*)eos)->vis;
}

double eos_linear_vis_pT(void* eos, EOS_ARGS* args)
{
    return ((EOS_LINEAR_VALUES*)eos)->vis;
}
