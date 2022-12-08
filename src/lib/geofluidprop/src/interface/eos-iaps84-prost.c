
#include "eos-iaps84-prost.h"

#include <stdlib.h>
#include <stdio.h>

#include <steam4.h>

#include "eos_type_impl.h"

void eos_iaps84_prost_register()
{
    const int index = EOS_TYPE_WATER_IAPS84_PROST;
    #define COMMAND(NAME)  eos_iaps84_prost_ ## NAME
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

void* eos_iaps84_prost_create()
{
    return NULL;
}

void eos_iaps84_prost_free(void* eos)
{
    free(eos);
}

double eos_iaps84_prost_rho_pT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'T', 0);
    water_tp(args->T, args->p, 1.e3, 1.e-8, prop);
    double rho = prop->d;
    freeProp(prop);
    return rho;
}

double eos_iaps84_prost_rho_ph(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'h', 0);
    water_ph(args->p, args->h, 300., 1000., 1.e-8, 1.e-8, prop);
    double rho = prop->d;
    freeProp(prop);
    return rho;
}

double eos_iaps84_prost_h_pT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'T', 0);
    double d=1000.0, dp=1.e-8;
    water_tp(args->T, args->p, d, dp, prop);
    double h = prop->h;
    freeProp(prop);
    return h;
}

double eos_iaps84_prost_drho_dT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'T', 0);
    double d=1000.0, dp=1.e-8;
    water_tp(args->T, args->p, d, dp, prop);
    double drho_dT = prop->dd->T_Cd;
    freeProp(prop);
    return drho_dT;
}

double eos_iaps84_prost_drho_dp(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'T', 0);
    double d=1000.0, dp=1.e-8;
    water_tp(args->T, args->p, d, dp, prop);
    double drho_dp = prop->dd->p_Ch;
    freeProp(prop);
    return drho_dp;
}

double eos_iaps84_prost_h_rhoT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('d', 'T', 0);
    water_td(args->T, args->rho, prop);
    double h = prop->h;
    freeProp(prop);
    return h;
}

double eos_iaps84_prost_dh_dT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('d', 'T', 0);
    water_td(args->T, args->rho, prop);
    double dhdT = prop->dh->T_Cd;
    freeProp(prop);
    return dhdT;
}

double eos_iaps84_prost_dh_dp(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('d', 'T', 0);
    water_td(args->T, args->rho, prop);
    double dhdp = prop->dh->p_Ch;
    freeProp(prop);
    return dhdp;
}

double eos_iaps84_prost_p_rhoT(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('T', 'd', 0);
    water_td(args->T, args->rho, prop);
    //dumpProp(stdout, prop);
    double p = prop->p;
    freeProp(prop);
    return p;
}

double eos_iaps84_prost_T_ph(void* eos, EOS_ARGS* args)
{
    Prop* prop = newProp('p', 'h', 0);
    water_ph(args->p, args->h, 300.0, 1000.0, 1e-8, 1e-8, prop);
    //dumpProp(stdout, prop);
    double T = prop->T;
    freeProp(prop);
    return T;
}


double eos_iaps84_prost_sat_p_T(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('p', 'p', 0);
    Prop* propg = newProp('p', 'p', 0);
    sat_t(args->T, propl, propg);
    double p = propl->p;
    freeProp(propl);
    freeProp(propg);
    return p;
}

double eos_iaps84_prost_sat_T_p(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('p', 'p', 0);
    Prop* propg = newProp('p', 'p', 0);
    sat_p(args->p, propl, propg);
    double T = propl->T;
    freeProp(propl);
    freeProp(propg);
    return T;
}

double eos_iaps84_prost_wv_ph(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('p', 'p', 0);
    Prop* propg = newProp('p', 'p', 0);
    sat_p(args->p, propl, propg);
    double hl = propl->h;
    double hv = propg->h;
    freeProp(propl);
    freeProp(propg);

    if (args->h <= hl)
        return 0.0;
    if (args->h >= hv)
        return 1.0;

    Prop* prop = newProp('p', 'h', 0);
    water_ph(args->p, args->h, 300.0, 1000.0, 1e-8, 1e-8, prop);
    //dumpProp(stdout, prop);
    double x = prop->x;
    freeProp(prop);
    return x;
}

double eos_iaps84_prost_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    double rho = propl->d;
    freeProp(propl);
    freeProp(propg);
    return rho;
}

double eos_iaps84_prost_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    double rho = propg->d;
    freeProp(propl);
    freeProp(propg);
    return rho;
}


void eos_iaps84_prost_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    *rhol = propl->d;
    *rhov = propg->d;
    freeProp(propl);
    freeProp(propg);
}

void eos_iaps84_prost_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_p(args->p, propl, propg);
    *rhol = propl->d;
    *rhov = propg->d;
    freeProp(propl);
    freeProp(propg);
}

double eos_iaps84_prost_sat_hl_T(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    double h = propl->h;
    freeProp(propl);
    freeProp(propg);
    return h;
}

double eos_iaps84_prost_sat_hv_T(void* eos, EOS_ARGS* args)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    double h = propg->h;
    freeProp(propl);
    freeProp(propg);
    return h;
}


void eos_iaps84_prost_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_t(args->T, propl, propg);
    *hl = propl->h;
    *hv = propg->h;
    freeProp(propl);
    freeProp(propg);
}

void eos_iaps84_prost_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double* hv)
{
    Prop* propl = newProp('T', 'T', 0);
    Prop* propg = newProp('T', 'T', 0);
    sat_p(args->p, propl, propg);
    *hl = propl->h;
    *hv = propg->h;
    freeProp(propl);
    freeProp(propg);
}

double eos_iaps84_prost_vis_rhoT(void* eos, EOS_ARGS* args)
{
    double p = eos_iaps84_prost_p_rhoT(eos, args);
    Prop* prop = newProp('d', 'p', 1);
    prop->p = p;
    prop->d = args->rho;
    prop->T = args->T;
    double vis = viscos(prop);
    freeProp(prop);
    return vis;
}
