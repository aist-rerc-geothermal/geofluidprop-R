#include "eos-iapws95sbtl-ph.h"

#include <assert.h>
#include <stdlib.h>

#include "sbtl/iapws_sbtl_ph/IAPWS95_SBTL_ph.h"
#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-Viscosity-85.h"
#include "model/iapws/IAPWS-Viscosity-08.h"
#include "model/iapws/IAPWS-IF97.h"

#include "util/utility.h"

#include "eos_type_impl.h"

void eos_iapws95sbtl_ph_register()
{
    const int index = EOS_TYPE_WATER_IAPWS95_SBTL;
    #define COMMAND(NAME)  eos_iapws95sbtl_ph_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pT);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_ph);
    // eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    // eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    //eos_impl_calc_u_pT[index] = COMMAND(u_pT);
    // eos_impl_calc_h_rhoT[index] = COMMAND(h_rhoT);
    eos_impl_calc_h_pT[index] = COMMAND(h_pT);
    // eos_impl_calc_dh_dT[index] = COMMAND(dh_dT);
    // eos_impl_calc_dh_dp[index] = COMMAND(dh_dp);
    // eos_impl_calc_p_rhoT[index] = COMMAND(p_rhoT);
    eos_impl_calc_T_ph[index] = COMMAND(T_ph);
    eos_impl_calc_wv_ph[index] = COMMAND(wv_ph);
    eos_impl_calc_sat_p_T[index] = COMMAND(sat_p_T);
    eos_impl_calc_sat_T_p[index] = COMMAND(sat_T_p);
    // eos_impl_calc_sat_rhol_T[index] = COMMAND(sat_rhol_T);
    // eos_impl_calc_sat_rhov_T[index] = COMMAND(sat_rhov_T);
    eos_impl_calc_sat_rholv_T[index] = COMMAND(sat_rholv_T);
    eos_impl_calc_sat_rholv_p[index] = COMMAND(sat_rholv_p);
    eos_impl_calc_sat_hlv_T[index] = COMMAND(sat_hlv_T);
    eos_impl_calc_sat_hlv_p[index] = COMMAND(sat_hlv_p);
    eos_impl_calc_vis_rhoT[index] = COMMAND(vis_rhoT);
}

void* eos_iapws95sbtl_ph_create()
{
    return iapws95sbtl_ph_create();
}

void eos_iapws95sbtl_ph_free(void* sbtl)
{
    iapws95sbtl_ph_free(sbtl);
}

double eos_iapws95sbtl_ph_rho_ph(void* sbtl, EOS_ARGS* args)
{
    double wv = 0, hl=0, hv=0;
    if (args->p < iapws95_get_pc())
    {
        hl = iapws95sbtl_ph_sat_hl_p(sbtl, args->p);
        hv = iapws95sbtl_ph_sat_hv_p(sbtl, args->p);
        wv = calc_vapor_quality(hl, hv, args->h);
    }

    if (wv==0.0 || wv==1.0)
        return iapws95sbtl_ph_rho_ph(sbtl, args->p, args->h);

    double rhol = iapws95sbtl_ph_rho_ph(sbtl, args->p, hl);
    double rhov = iapws95sbtl_ph_rho_ph(sbtl, args->p, hv);
    double v = wv/rhov + (1.-wv)/rhol;
    return 1./v;
}

double eos_iapws95sbtl_ph_rho_pT(void* sbtl, EOS_ARGS* args)
{
    args->h = eos_iapws95sbtl_ph_h_pT(sbtl, args);
    return eos_iapws95sbtl_ph_rho_ph(sbtl, args);
}

double eos_iapws95sbtl_ph_h_pT(void* sbtl, EOS_ARGS* args)
{
    return iapws95sbtl_ph_h_pT(sbtl, args->p, args->T);
}

double eos_iapws95sbtl_ph_T_ph(void* sbtl, EOS_ARGS* args)
{
    double wv = 0, hl=0, hv=0;
    if (args->p < iapws95_get_pc())
    {
        hl = iapws95sbtl_ph_sat_hl_p(sbtl, args->p);
        hv = iapws95sbtl_ph_sat_hv_p(sbtl, args->p);
        wv = calc_vapor_quality(hl, hv, args->h);
    }
    if (wv==0.0 || wv==1.0)
        return iapws95sbtl_ph_T_ph(sbtl, args->p, args->h);

    return iapws95sbtl_ph_sat_T_p(sbtl, args->p);
}

double eos_iapws95sbtl_ph_sat_p_T(void* sbtl, EOS_ARGS* args)
{
    return iapws95sbtl_ph_sat_T_p(sbtl, args->T);
}

double eos_iapws95sbtl_ph_sat_T_p(void* sbtl, EOS_ARGS* args)
{
    return iapws95sbtl_ph_sat_T_p(sbtl, args->p);
}

void eos_iapws95sbtl_ph_sat_rholv_p(void* sbtl, EOS_ARGS* args, double* rhol, double* rhov)
{
    iapws95sbtl_ph_sat_rho_p(sbtl, args->p, rhol, rhov);
}

void eos_iapws95sbtl_ph_sat_rholv_T(void* sbtl, EOS_ARGS* args, double* rhol, double* rhov)
{
    double ps = iapws95sbtl_ph_sat_p_T(sbtl, args->T);
    iapws95sbtl_ph_sat_rho_p(sbtl, ps, rhol, rhov);
}

void eos_iapws95sbtl_ph_sat_hlv_p(void* sbtl, EOS_ARGS* args, double* hl, double* hv)
{
    iapws95sbtl_ph_sat_h_p(sbtl, args->p, hl, hv);
}

void eos_iapws95sbtl_ph_sat_hlv_T(void* sbtl, EOS_ARGS* args, double* hl, double* hv)
{
    double ps = iapws95sbtl_ph_sat_p_T(sbtl, args->T);
    iapws95sbtl_ph_sat_h_p(sbtl, ps, hl, hv);
}

double eos_iapws95sbtl_ph_wv_ph(void* sbtl, EOS_ARGS* args)
{
    if (args->p >= iapws95_get_pc())
        return 0.0;

    double hl, hv;
    iapws95sbtl_ph_sat_h_p(sbtl, args->p, &hl, &hv);
    double wv = calc_vapor_quality(hl, hv, args->h);
    return wv;
}

double eos_iapws95sbtl_ph_vis_rhoT(void* eos, EOS_ARGS* args)
{
    return iapws85_viscosity_rhoT(args->rho, args->T);
//    return iapws08_viscosity_rhoT(args->rho, args->T);
}

