
#include "eos-driesner07-prost.h"

#include <stdlib.h>

#include "model/driesner07/Driesner2007_H2ONaCl_prost.h"
#include "model/klyukinetal17/KlyukinEtAl2017.h"

#include "eos_type_impl.h"

void eos_driesner07_prost_register()
{
    const int index = EOS_TYPE_H2ONaCl_DRIESNER07_PROST;
    #define COMMAND(NAME)  eos_driesner07_prost_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pTx);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_phx);
    // eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    // eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    // eos_impl_calc_h_rhoT[index] = COMMAND(h_rhoT);
    eos_impl_calc_h_pT[index] = COMMAND(h_pTx);
    // eos_impl_calc_dh_dT[index] = COMMAND(dh_dT);
    // eos_impl_calc_dh_dp[index] = COMMAND(dh_dp);
    eos_impl_calc_T_ph[index] = COMMAND(T_phx);
    eos_impl_calc_wv_ph[index] = COMMAND(wv_phx);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoTx);
}

void* eos_driesner07_prost_create()
{
    return NULL;
}

void eos_driesner07_prost_free(void* eos)
{
    free(eos);
}

double eos_driesner07_prost_rho_pTx(void* eos, EOS_ARGS* args)
{
    return driesner07_prost_H2O_NaCl_rho_pTx(args->p, args->T, args->x1);
}

double eos_driesner07_prost_rho_phx(void* eos, EOS_ARGS* args)
{
    double T = eos_driesner07_prost_T_phx(eos, args);
    return driesner07_prost_H2O_NaCl_rho_pTx(args->p, T, args->x1);
}

double eos_driesner07_prost_h_pTx(void* eos, EOS_ARGS* args)
{
    return driesner07_prost_H2O_NaCl_h_pTx(args->p, args->T, args->x1);
}

double eos_driesner07_prost_T_phx(void* eos, EOS_ARGS* args)
{
    double T0 = 400.+273.15; //TODO
    return driesner07_prost_H2O_NaCl_T_phx(args->p, args->h, args->x1, T0);
}

double eos_driesner07_prost_wv_phx(void* eos, EOS_ARGS* args)
{
    double T = eos_driesner07_prost_T_phx(eos, args);
    double wl = driesner07_prost_H2O_NaCl_mass_fraction_liquid(args->p, T, args->x1);
    return (1.-wl);
}

// double eos_driesner07_prost_sat_p_Tx(void* eos, EOS_ARGS* args)
// {

// }

double eos_driesner07_prost_vis_rhoTx(void* eos, EOS_ARGS* args)
{
    return klyukinetal2017_viscosity(args->rho, args->T, args->x1);
}
