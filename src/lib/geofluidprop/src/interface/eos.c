
#include "eos.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "eos_type_impl.h"

//#define str(s) #s
#ifdef NDEBUG
#define CHECK_IMPL(fname)
#else
#define CHECK_IMPL(fname) \
    if (!fname[eos->eos_type]) \
    { \
        printf("%s() is not implemented for %s\n", #fname, eos->eos_type_name); \
        assert(fname[eos->eos_type]); \
    }
#endif



void eos_initialization()
{
    static bool is_eos_registerd = false;
    if (is_eos_registerd) return;
    is_eos_registerd = true;

    eos_impl_init_register();
    for (int i=0; i<EOS_IMPL_ARRAY_MAX_SIZE; i++)
    {
        if (eos_impl_register[i] != NULL)
            eos_impl_register[i]();
    }
}

EOS* eos_create_by_name(const char* eos_type_name)
{
    eos_initialization();
    const int eos_type = eos_get_eos_type(eos_type_name);
    if (eos_type == EOS_TYPE_UNKNOWN)
    {
        printf("***error: The given EoS Type (%s) is not found in this build.\n", eos_type_name);
        return NULL;
    }
    if (eos_impl_register[eos_type]==NULL)
    {
        printf("***error: The given EoS Type (%s) is not supported in this build.\n", eos_type_name);
        return NULL;
    }
    EOS* eos = malloc(sizeof(EOS));
    eos->eos_type = eos_type;
    eos_get_eos_type_name(eos_type, eos->eos_type_name);
    eos->impl = eos_impl_create[eos_type]();
    return eos;
}

EOS* eos_create(int eos_type)
{
    eos_initialization();
    EOS* eos = malloc(sizeof(EOS));
    eos->eos_type = eos_type;
    eos_get_eos_type_name(eos_type, eos->eos_type_name);
    eos->impl = eos_impl_create[eos_type]();
    return eos;
}

void eos_free(EOS* eos)
{
    if (eos==NULL) return;
    CHECK_IMPL(eos_impl_free);
    eos_impl_free[eos->eos_type](eos->impl);
    free(eos);
}

double eos_rho_pT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_rho_pT);
    return eos_impl_calc_rho_pT[eos->eos_type](eos->impl, args);
}

double eos_rho_ph(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_rho_ph);
    return eos_impl_calc_rho_ph[eos->eos_type](eos->impl, args);
}

double eos_drho_dh(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_drho_dh);
    return eos_impl_calc_drho_dh[eos->eos_type](eos->impl, args);
}

double eos_drho_dT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_drho_dT);
    return eos_impl_calc_drho_dT[eos->eos_type](eos->impl, args);
}

double eos_drho_dp(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_drho_dp);
    return eos_impl_calc_drho_dp[eos->eos_type](eos->impl, args);
}

double eos_u_phrho(EOS* eos, EOS_ARGS* args)
{
    double u = args->h - args->p / args->rho;
    return u;
}

double eos_u_pT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_u_pT);
    return eos_impl_calc_u_pT[eos->eos_type](eos->impl, args);
}

double eos_h_rhoT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_h_rhoT);
    return eos_impl_calc_h_rhoT[eos->eos_type](eos->impl, args);
}

double eos_h_pT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_h_pT);
    return eos_impl_calc_h_pT[eos->eos_type](eos->impl, args);
}

double eos_dh_dT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_dh_dT);
    return eos_impl_calc_dh_dT[eos->eos_type](eos->impl, args);
}

double eos_dh_dp(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_dh_dp);
    return eos_impl_calc_dh_dp[eos->eos_type](eos->impl, args);
}

double eos_p_rhoT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_p_rhoT);
    return eos_impl_calc_p_rhoT[eos->eos_type](eos->impl, args);
}

double eos_T_pu(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_T_pu);
    return eos_impl_calc_T_pu[eos->eos_type](eos->impl, args);
}

double eos_T_ph(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_T_ph);
    return eos_impl_calc_T_ph[eos->eos_type](eos->impl, args);
}

double eos_wv_ph(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_wv_ph);
    return eos_impl_calc_wv_ph[eos->eos_type](eos->impl, args);
}

double eos_sat_p_T(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_sat_p_T);
    return eos_impl_calc_sat_p_T[eos->eos_type](eos->impl, args);
}

double eos_sat_T_p(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_sat_T_p);
    return eos_impl_calc_sat_T_p[eos->eos_type](eos->impl, args);
}

double eos_sat_rhol_T(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_sat_rhol_T);
    return eos_impl_calc_sat_rhol_T[eos->eos_type](eos->impl, args);
}

double eos_sat_rhov_T(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_sat_rhov_T);
    return eos_impl_calc_sat_rhov_T[eos->eos_type](eos->impl, args);
}

void eos_sat_rholv_T(EOS* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    CHECK_IMPL(eos_impl_calc_sat_rholv_T);
    eos_impl_calc_sat_rholv_T[eos->eos_type](eos->impl, args, rhol, rhov);
}

void eos_sat_rholv_p(EOS* eos, EOS_ARGS* args, double* rhol, double* rhov)
{
    CHECK_IMPL(eos_impl_calc_sat_rholv_p);
    eos_impl_calc_sat_rholv_p[eos->eos_type](eos->impl, args, rhol, rhov);
}

void eos_sat_hlv_p(EOS* eos, EOS_ARGS* args, double* hl, double* hv)
{
    CHECK_IMPL(eos_impl_calc_sat_hlv_p);
    eos_impl_calc_sat_hlv_p[eos->eos_type](eos->impl, args, hl, hv);
}

double eos_vis_rhoT(EOS* eos, EOS_ARGS* args)
{
    CHECK_IMPL(eos_impl_calc_vis_rhoT);
    return eos_impl_calc_vis_rhoT[eos->eos_type](eos->impl, args);
}
