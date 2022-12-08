
#include "eos-if97-freesteam.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <freesteam/derivs.h>
#include <freesteam/steam_ph.h>
#include <freesteam/steam_pT.h>
#include <freesteam/steam_pv.h>
#include <freesteam/steam_uv.h>
#include <freesteam/steam_Tx.h>
#include <freesteam/region4.h>

//#include "model/iapws/IAPWS-Viscosity-85.h"

#include "eos_type_impl.h"

void eos_if97_freesteam_register()
{
    const int index = EOS_TYPE_WATER_IF97_FREESTEAM;
    #define COMMAND(NAME)  eos_if97_freesteam_ ## NAME
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
    eos_impl_calc_sat_p_T[index] = COMMAND(sat_p_T);
    eos_impl_calc_sat_T_p[index] = COMMAND(sat_T_p);
    eos_impl_calc_sat_rhol_T[index] = COMMAND(sat_rhol_T);
    eos_impl_calc_sat_rhov_T[index] = COMMAND(sat_rhov_T);
    eos_impl_calc_sat_rholv_T[index] = COMMAND(sat_rholv_T);
    eos_impl_calc_sat_rholv_p[index] = COMMAND(sat_rholv_p);
    eos_impl_calc_sat_hlv_T[index] = COMMAND(sat_hlv_T);
    eos_impl_calc_sat_hlv_p[index] = COMMAND(sat_hlv_p);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoT);
}

void* eos_if97_freesteam_create()
{
    return NULL;
}

void eos_if97_freesteam_free(void* eos)
{
    free(eos);
}

double eos_if97_freesteam_rho_pT(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_pT(args->p, args->T);
    double rho = freesteam_rho(S);
    return rho;
}

double eos_if97_freesteam_rho_ph(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_ph(args->p, args->h);
    double rho = freesteam_rho(S);
    return rho;
}

double eos_if97_freesteam_h_pT(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_pT(args->p, args->T);
    double h = freesteam_h(S);
    return h;
}

double eos_if97_freesteam_drho_dT(void* eos, EOS_ARGS* args)
{
    double h = eos_if97_freesteam_h_pT(eos, args);
    SteamState S = freesteam_set_ph(args->p, args->h);
    double dv_dT = freesteam_deriv(S,"Tpv");
    double v = freesteam_v(S);
    double drho_dT = -1./(v*v)*dv_dT;
    return drho_dT;
}

double eos_if97_freesteam_drho_dp(void* eos, EOS_ARGS* args)
{
    double h = eos_if97_freesteam_h_pT(eos, args);
    SteamState S = freesteam_set_ph(args->p, args->h);
    double dv_dp = freesteam_deriv(S,"pTv");
    double v = freesteam_v(S);
    double drho_dp = -1./(v*v)*dv_dp;
    return drho_dp;
}

double eos_if97_freesteam_h_rhoT(void* eos, EOS_ARGS* args)
{
    double p = eos_if97_freesteam_p_rhoT(eos, args);
    SteamState S = freesteam_set_pv(p, 1./args->rho);
    double h = freesteam_h(S);
    return h;
}

double eos_if97_freesteam_dh_dT(void* eos, EOS_ARGS* args)
{
    double h = eos_if97_freesteam_h_pT(eos, args);
    SteamState S = freesteam_set_ph(args->p, args->h);
    double dh_dT = freesteam_deriv(S,"Tph");
    return dh_dT;
}

double eos_if97_freesteam_dh_dp(void* eos, EOS_ARGS* args)
{
    double h = eos_if97_freesteam_h_pT(eos, args);
    SteamState S = freesteam_set_ph(args->p, args->h);
    double dh_dp = freesteam_deriv(S,"pTh");
    return dh_dp;
}

double eos_if97_freesteam_p_rhoT(void* eos, EOS_ARGS* args)
{
    // get initial guess
    double p0 = 0;
    if (args->rho >= eos_if97_freesteam_sat_rhol_T(eos, args))
        p0 = freesteam_region4_psat_T(args->T) * 1.2;
    else
        p0 = freesteam_region4_psat_T(args->T) * 0.8;

    // run Newton-Raphson
    double p = p0;
    double dp = 1e6*1e-3;
    EOS_ARGS args_tmp;
    args_tmp.rho = args->rho;
    for (int i=0; i<100; i++)
    {
        args_tmp.p = p;
        double T1 = eos_if97_freesteam_T_rhop(eos, &args_tmp);
        double r = T1 - args->T;
        // printf("%d: r=%g\n", i, r);
        if (fabs(r)<args->T*1e-8)
            return p;
        args_tmp.p = p + dp;
        double T2 = eos_if97_freesteam_T_rhop(eos, &args_tmp);
        double drdp = (T2-T1)/dp;
        p -= r/drdp;
    }

    assert(false);

    return -1.0;
}

double eos_if97_freesteam_T_ph(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_ph(args->p, args->h);
    double T = freesteam_T(S);
    return T;
}

double eos_if97_freesteam_T_rhop(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_pv(args->p, 1./args->rho);
    double T = freesteam_T(S);
    return T;
}

double eos_if97_freesteam_sat_p_T(void* eos, EOS_ARGS* args)
{
    return freesteam_region4_psat_T(args->T);
}

double eos_if97_freesteam_sat_T_p(void* eos, EOS_ARGS* args)
{
    return freesteam_region4_Tsat_p(args->p);
}

double eos_if97_freesteam_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_Tx(args->T, 0.);
    double rho = freesteam_rho(S);
    return rho;
}

double eos_if97_freesteam_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_Tx(args->T, 1.);
    double rho = freesteam_rho(S);
    return rho;
}

void eos_if97_freesteam_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    *rhol = eos_if97_freesteam_sat_rhol_T(eos, args);
    *rhov = eos_if97_freesteam_sat_rhov_T(eos, args);
}

void eos_if97_freesteam_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    args->T = eos_if97_freesteam_sat_T_p(eos, args);
    eos_if97_freesteam_sat_rholv_T(eos, args, rhol, rhov);
}

double eos_if97_freesteam_sat_hl_T(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_Tx(args->T, 0.);
    return freesteam_h(S);
}

double eos_if97_freesteam_sat_hv_T(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_Tx(args->T, 1.);
    return freesteam_h(S);
}

void eos_if97_freesteam_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    *hl = eos_if97_freesteam_sat_hl_T(eos, args);
    *hv = eos_if97_freesteam_sat_hv_T(eos, args);
}

void eos_if97_freesteam_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    args->T = eos_if97_freesteam_sat_T_p(eos, args);
    eos_if97_freesteam_sat_hlv_T(eos, args, hl, hv);
}

double eos_if97_freesteam_wv_ph(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_ph(args->p, args->h);
    double wv = freesteam_x(S);
    return wv;
}

double eos_if97_freesteam_vis_rhoT(void* eos, EOS_ARGS* args)
{
    // return iapws85_viscosity_rhoT(args->rho, args->T);
    double p = eos_if97_freesteam_p_rhoT(eos, args);
    SteamState S = freesteam_set_pv(p, 1./args->rho);
    double vis = freesteam_mu(S);
    return vis;
}

double eos_if97_freesteam_vis_pT(void* eos, EOS_ARGS* args)
{
    SteamState S = freesteam_set_pT(args->p, args->T);
    double vis = freesteam_mu(S);
    return vis;
}
