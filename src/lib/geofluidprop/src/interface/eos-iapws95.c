
#include "eos-iapws95.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#ifdef USE_QUAD
#include <quadmath.h>
#endif

#include "util/utility.h"

#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-Melt-11.h"
#include "model/iapws/IAPWS-SAT-92.h"
#include "model/iapws/IAPWS-Viscosity-85.h"
#include "model/iapws/IAPWS-Viscosity-08.h"
#include "model/iapws/IAPWS-IF97.h"

#include "eos_type_impl.h"

void eos_iapws95_register()
{
    const int index = EOS_TYPE_WATER_IAPWS95;
    #define COMMAND(NAME)  eos_iapws95_ ## NAME
    eos_impl_create[index] = COMMAND(create);
    eos_impl_free[index] = COMMAND(free);
    eos_impl_calc_rho_pT[index] = COMMAND(rho_pT);
    eos_impl_calc_rho_ph[index] = COMMAND(rho_ph);
    eos_impl_calc_drho_dT[index] = COMMAND(drho_dT);
    eos_impl_calc_drho_dp[index] = COMMAND(drho_dp);
    //eos_impl_calc_u_pT[index] = COMMAND(u_pT);
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

void* eos_iapws95_create()
{
    return NULL;
}

void eos_iapws95_free(void* eos)
{
    free(eos);
}

int eos_iapws95_check_validity_pT(void* eos, double p, double T, const char* fname)
{
    if (T < iapws95_get_Tt()) { //} || T > 1273) {
        printf("%s(): Outside limits of Ttr <= T <= 1275. Given T = %g\n", fname, T);
        abort();
        return 0;
    }
    if (p < iapws95_get_pt_approx()) { //} || p > 1000e6) {
        printf("%s(): Outside limits of pTr <= p <= 1000MPa. Given p = %g\n", fname, p);
        //abort();
        return 0;
    }
    // double Tmelt = melt_temperature(p);
    // if (T < TTr || T > 1273) {
    //     printf("%s(): Outside limits of Tmelt(%g) <= T <= 1275. Given T = %g\n", fname, Tmelt, T);
    //     return 0;
    // }
    // double pmelt = T>TTr ? melt_presssure(T) : sublimation_pressure(T);
    // if (p < pmelt || p > 1000e6) {
    //     printf("%s(): Outside limits of pmelt <= p <= 1000MPa. Given p = %g\n", fname, p);
    //     return 0;
    // }
    return 1;
}

int eos_iapws95_check_validity_ph(void* eos, double p, double h, const char* fname)
{
    if (p < 0) {
        printf("%s(): Outside limits of p. Given p = %g <= 0\n", fname, p);
        return 0;
    }
    if (p > iapws95_get_pt_approx() + 0.1)
    {
        double hmin_approx = if97_h_pT(p*1e-6, iapws95_get_Tt());
        if (h > hmin_approx*1.03) {
            return 1;
        }
    }
    EOS_ARGS args_tmp;
    args_tmp.p = p;
    args_tmp.T = iapws95_get_Tt();
    double hmin = eos_iapws95_h_pT(eos, &args_tmp);
    if (h < hmin) {
        printf("%s(): Outside limits of h. Given h = %g < hmin(p)=%g, p=%g\n", fname, h, hmin, p);
        return 0;
    }
    return 1;
}

double eos_iapws95_rho_pT(void* eos, EOS_ARGS* args)
{
    return iapws95_rho_pT(args->p, args->T);
}

double eos_iapws95_rho_ph(void* eos, EOS_ARGS* args)
{
    return iapws95_rho_ph(args->p, args->h);
}

double eos_iapws95_drho_dh(void* eos, EOS_ARGS* args)
{
    return iapws95_drho_dh_ph(args->p, args->h);
}


double eos_iapws95_drho_dT(void* eos, EOS_ARGS* args)
{
    return iapws95_drho_dT_pT(args->p, args->T, args->rho);
}

double eos_iapws95_drho_dp(void* eos, EOS_ARGS* args)
{
    return iapws95_drho_dp_pT(args->p, args->T, args->rho);
}

double eos_iapws95_h_rhoT(void* eos, EOS_ARGS* args)
{
    return iapws95_h_rhoT(args->rho, args->T);
}

double eos_iapws95_h_pT(void* eos, EOS_ARGS* args)
{
    return iapws95_h_pT(args->p, args->T);
}

double eos_iapws95_dh_dT(void* eos, EOS_ARGS* args)
{
    return iapws95_dh_dT_rhoT(args->rho, args->T);
}

double eos_iapws95_dh_dp(void* eos, EOS_ARGS* args)
{
    return iapws95_dh_dp_rhopT(args->rho, args->p, args->T);
}

double eos_iapws95_p_rhoT(void* eos, EOS_ARGS* args)
{
    return iapws95_p_rhoT(args->rho, args->T);
}

double eos_iapws95_T_ph(void* eos, EOS_ARGS* args)
{
    return iapws95_T_ph(args->p, args->h);
}

double eos_iapws95_wv_ph(void* eos, EOS_ARGS* args)
{
    return iapws95_wv_ph(args->p, args->h);
}

double eos_iapws95_sat_p_T(void* eos, EOS_ARGS* args)
{
    return iapws95_sat_p_T(args->T);
}

#ifdef USE_QUAD
double eos_iapws95_sat_p_T_quad(void* eos, EOS_ARGS* args)
{
    return iapws95_sat_p_T_quad(args->T);
}
#endif

double eos_iapws95_sat_T_p(void* eos, EOS_ARGS* args)
{
    return iapws95_sat_T_p(args->p);
}

void eos_iapws95_sat_Trho_p(void* eos, EOS_ARGS* args, double*T, double*rhol, double*rhov)
{
    iapws95_sat_Trho_p(args->p, T, rhol, rhov);
}

void eos_iapws95_sat_rholv_T(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    if (args->T == iapws95_get_Tc()) {
        *rhol = *rhov = iapws95_get_rhoc();
        return;
    }
    double sat_p, sat_rhol, sat_rhov;
    sat_p = if97_sat_p_T(args->T)*1e6;
    if97_sat_rho_p(sat_p*1e-6, &sat_rhol, &sat_rhov);
    if (args->T == iapws95_get_Tt())
        sat_p = iapws95_get_pt_approx();
    iapws95_sat_prho_T(args->T, &sat_p, &sat_rhol, &sat_rhov);
    *rhol = sat_rhol;
    *rhov = sat_rhov;
}

#ifdef USE_LONGDOUBLE
void eos_iapws95_sat_rholv_T_long(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    if (args->T == iapws95_get_Tc()) {
        *rhol = *rhov = iapws95_get_rhoc();
        return;
    }
    double sat_rhol, sat_rhov;
    if97_sat_rho_T(args->T, &sat_rhol, &sat_rhov);
    double sat_ps = if97_sat_p_T(args->T)*1e6;
    if (args->T == iapws95_get_Tt())
        sat_ps = iapws95_get_pt_approx();
    //printf("p=%g, rhol=%g, rhov=%g\n", sat_ps, sat_rhol, sat_rhov);
    long double pl = sat_ps;
    long double rholl, rhovl;
    rholl = sat_rhol; rhovl = sat_rhov;
    iapws95_sat_prho_T_long(args->T, &pl, &rholl, &rhovl);
    *rhol = rholl;
    *rhov = rhovl;
}
#endif

#ifdef USE_QUAD
void eos_iapws95_sat_rholv_T_quad(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    if (args->T == iapws95_get_Tc()) {
        *rhol = *rhov = iapws95_get_rhoc();
        return;
    }
    double sat_rhol, sat_rhov;
    if97_sat_rho_T(args->T, &sat_rhol, &sat_rhov);
    double sat_ps = if97_sat_p_T(args->T)*1e6;
    if (args->T == iapws95_get_Tt())
        sat_ps = iapws95_get_pt_approx();
    //printf("p=%g, rhol=%g, rhov=%g\n", sat_ps, sat_rhol, sat_rhov);
    __float128 pl = sat_ps;
    __float128 rholl, rhovl;
    rholl = sat_rhol; rhovl = sat_rhov;
    iapws95_sat_prho_T_quad(args->T, &pl, &rholl, &rhovl);
    *rhol = rholl;
    *rhov = rhovl;
}
#endif

void eos_iapws95_sat_rholv_p(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    double sat_T, sat_rhol, sat_rhov;
    sat_T = if97_sat_T_p(args->p*1e-6);
    if97_sat_rho_p(args->p*1e-6, &sat_rhol, &sat_rhov);
    iapws95_sat_Trho_p(args->p, &sat_T, &sat_rhol, &sat_rhov);
    *rhol = sat_rhol;
    *rhov = sat_rhov;
}

#ifdef USE_LONGDOUBLE
void eos_iapws95_sat_rholv_p_long(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    double sat_T, sat_rhol, sat_rhov;
    sat_T = if97_sat_T_p(args->p*1e-6);
    if97_sat_rho_p(args->p*1e-6, &sat_rhol, &sat_rhov);
    long double Tq = sat_T, rholq = sat_rhol, rhovq = sat_rhov;
    iapws95_sat_Trho_p_long(args->p, &Tq, &rholq, &rhovq);
    *rhol = rholq;
    *rhov = rhovq;
}
#endif

#ifdef USE_QUAD
void eos_iapws95_sat_rholv_p_quad(void* eos, EOS_ARGS* args, double* rhol, double*rhov)
{
    double sat_T, sat_rhol, sat_rhov;
    sat_T = if97_sat_T_p(args->p*1e-6);
    if97_sat_rho_p(args->p*1e-6, &sat_rhol, &sat_rhov);
    __float128 Tq = sat_T, rholq = sat_rhol, rhovq = sat_rhov;
    iapws95_sat_Trho_p_quad(args->p, &Tq, &rholq, &rhovq);
    *rhol = rholq;
    *rhov = rhovq;
}
#endif

double eos_iapws95_sat_rhol_T(void* eos, EOS_ARGS* args)
{
    return iapws95_sat_rhol_T(args->T);
}

double eos_iapws95_sat_rhov_T(void* eos, EOS_ARGS* args)
{
    return iapws95_sat_rhov_T(args->T);
}

double eos_iapws95_sat_p_hl(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    assert(args->h >= eos_iapws95_hl_tr(eos));
#if 0
    double pmin = iapws95_get_pt_approx();
    double pmax = iapws95_get_pc();
    double p1 = pmin, p2 = pmax;
    double r1 = iapws95_sat_hl_p(p1, eos) - args->h;
    double r2 = iapws95_sat_hl_p(p2, eos) - args->h;
    printf("r1=%.3e, r2=%.3e\n", r1, r2);
    for (int i=0; i<70; i++)
    {
        double pm = 0.5 * (p1 + p2);
        double hh = iapws95_sat_hl_p(pm, eos);
        double r = hh - args->h;
        printf("%d: p=%.10e, r=%.3e\n", i, pm, r/fabs(args->h));
        if (fabs(r)<1e-8*fabs(args->h)) {
            return pm;
        }

        if (r < 0)
            p1 = pm;
        else
            p2 = pm;
        if (fabs(p1-p2)<1.e-3)
            break;
    }
#else
    double p0 = 1.0e6;
    double p = p0;
    int inl = 0;
    for (inl=0; inl<15; inl++)
    {
        double hh = iapws95_sat_hl_p(p);
        double r = hh - args->h;
        //printf("%d: p=%g, h=%g, r=%g\n", inl, p, hh, fabs(r)/args->h);
        if (fabs(r)<1e-8*fabs(args->h)) {
            return p;
        }

        double d = p*1e-8;
        double drdp = 0;
        if (p+d <= iapws95_get_pc())
        {
            double hh1 = iapws95_sat_hl_p(p+d);
            drdp = (hh1 - hh)/d;
        } else {
            double hh0 = iapws95_sat_hl_p(p-d);
            drdp = (hh - hh0)/d;
        }
        p -= r/drdp;
        p = min(max(p, iapws95_get_pt_approx()), iapws95_get_pc());
    }
#endif
    printf("%s: did not converge at h=%.12g\n", __FUNCTION__, args->h);
    abort();
}

#ifdef USE_QUAD
double eos_iapws95_sat_p_hl_quad(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    double pmin = iapws95_get_pt_approx();
    double pmax = iapws95_get_pc();
    __float128 p1 = pmin, p2 = pmax;
    for (int i=0; i<70; i++)
    {
        __float128 pm = 0.5 * (p1 + p2);
        __float128 hh = iapws95_sat_hl_p_quad(pm);
        __float128 r = hh - args->h;
        //printf("%d: p=%g, r=%g\n", i, pm, r);
        if (fabsq(r)<1e-8*fabsq(args->h)) {
            return pm;
        }

        if (r < 0)
            p1 = pm;
        else
            p2 = pm;
    }
    printf("%s: did not converge at h=%g\n", __FUNCTION__, args->h);
    abort();
}
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_p_hl_long(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    double pmin = iapws95_get_pt_approx();
    double pmax = iapws95_get_pc();
    double p1 = pmin, p2 = pmax;
    for (int i=0; i<70; i++)
    {
        double pm = 0.5 * (p1 + p2);
        double hh = iapws95_sat_hl_p_long(pm);
        double r = hh - args->h;
        //printf("%d: p=%g, r=%g\n", i, pm, r);
        if (fabs(r)<1e-8*fabs(args->h)) {
            return pm;
        }

        if (r < 0)
            p1 = pm;
        else
            p2 = pm;
    }
#ifdef USE_QUAD
    return eos_iapws95_sat_p_hl_quad(eos, args);
#else
    printf("%s: did not converge at h=%.10e\n", __FUNCTION__, args->h);
    abort();
#endif
}
#endif

#ifdef USE_QUAD
double eos_iapws95_sat_T_hl_quad(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    __float128 Tmin = iapws95_get_Tt();
    __float128 Tmax = iapws95_get_Tc();
    __float128 T1 = Tmin, T2 = Tmax;
    for (int i=0; i<70; i++)
    {
        __float128 Tm = 0.5 * (T1 + T2);
        __float128 hh = iapws95_sat_hl_T_quad(Tm);
        __float128 r = hh - args->h;
        //printf("%d: T=%.10Lg, r=%.5Le\n", i, Tm, r/args->h);
        if (fabsq(r)<1e-8*fabs(args->h)) {
            return Tm;
        }

        if (r < 0)
            T1 = Tm;
        else
            T2 = Tm;
    }
    printf("%s: did not converge at h=%g\n", __FUNCTION__, args->h);
    abort();
}
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_T_hl_long(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    long double Tmin = iapws95_get_Tt();
    long double Tmax = iapws95_get_Tc();
#if 0
    long double T = (Tmin+Tmax)*0.5;
    for (int i=0; i<100; i++)
    {
        long double hh = iapws95_sat_hl_T_long(T, eos);
        long double r = hh - args->h;
        printf("%d: T=%.10Lg, r=%.5Le\n", i, T, r/args->h);
        if (fabsl(r)<1e-8*fabsl(args->h))
            return T;

        long double d = Tmax*1.e-8, drdT;
        if (T+d <= Tmax)
        {
            long double hh1 = iapws95_sat_hl_T_long(T+d, eos);
            drdT = (hh1-hh)/d;
        } else {
            long double hh1 = iapws95_sat_hl_T_long(T-d, eos);
            drdT = (hh-hh1)/d;
        }

        T -= r/drdT;
        T = maxl(minl(T, Tmax), Tmin);
    }
#else
    long double T1 = Tmin, T2 = Tmax;
    for (int i=0; i<70; i++)
    {
        long double Tm = 0.5 * (T1 + T2);
        long double hh = iapws95_sat_hl_T_long(Tm);
        long double r = hh - args->h;
        //printf("%d: T=%.10Lg, r=%.5Le\n", i, Tm, r/args->h);
        if (fabsl(r)<1e-8*fabsl(args->h)) {
            return Tm;
        }

        if (r < 0)
            T1 = Tm;
        else
            T2 = Tm;
    }
#endif
#ifdef USE_QUAD
    return eos_iapws95_sat_T_hl_quad(eos, args);
#else
    printf("%s: did not converge at h=%g\n", __FUNCTION__, args->h);
    abort();
#endif
}
#endif

double eos_iapws95_sat_p_hv(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    double pmin = iapws95_get_pt_approx();
    double pmax = iapws95_get_pc();
    double p1 = pmin, p2 = pmax;
    double r1 = iapws95_sat_hv_p(p1) - args->h;
    double r2 = iapws95_h_rhoT(iapws95_get_rhoc(), iapws95_get_Tc()) - args->h;
    for (int i=0; i<100; i++)
    {
        double pm = 0.5 * (p1 + p2);
        double hh = iapws95_sat_hv_p(pm);
        double r = hh - args->h;
        //printf("%d: p=%.10g, r=%g, p1=%.10g, p2=%.10g\n", i, pm, r, p1, p2);
        if (fabs(r)<1e-8*fabs(args->h)) {
            return pm;
        }

        if (r > 0)
            p1 = pm;
        else
            p2 = pm;
    }
    printf("%s: did not converge at h=%g\n", __FUNCTION__, args->h);
    abort();
}

#ifdef USE_QUAD
double eos_iapws95_sat_p_hv_quad(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    const __float128 pmin = iapws95_get_pt_approx();
    const __float128 pmax = iapws95_get_pc();
    const __float128 href = args->h;
#if 0
    const __float128 dp = pmax*1e-8;
    __float128 pm = pmax*0.99;
    for (int i=0; i<100; i++)
    {
        __float128 hh = iapws95_sat_hv_p_quad(pm, eos);
        __float128 r = hh - href;
        if (fabsq(r)<1e-8*href)
            return pm;

        __float128 hh1 = iapws95_sat_hv_p_quad(pm-dp, eos);
        __float128 hh2 = iapws95_sat_hv_p_quad(pm-2.*dp, eos);
        __float128 drdp = (hh-hh2)/(2*dp);
        printf("%d: p=%.12e, r=%.5e, hh=%.10e, jac=%.5e\n", i, (double)pm, (double)(r/href), (double)hh, (double)drdp);
        pm -= r/drdp;
        if (pm>pmax) pm = pmax*(1.q - 1.e-10q);
    }
#else
    __float128 p1 = pmin, p2 = pmax;
    __float128 r1 = eos_iapws95_hv_tr(eos) - href;
//    __float128 r2 = eos_iapws95_h_crit(eos) - href;
//    __float128 r2 = iapws95_h_rhoT_quad(iapws95_get_rhoc(), iapws95_get_Tc()) - href;
    __float128 r2 = iapws95_sat_hv_p_quad(p2) - href;
    printf("%s: pmin=%g, pmax%g, r_left=%g, r_right=%g\n", __FUNCTION__, (double)pmin, (double)pmax, (double)(r1/href), (double)(r2/href));
    if (r1*r2 > 0) {
        printf("%s: cannot apply bisecion because r_left=%g, r_right=%g\n", __FUNCTION__, (double)r1, (double)r2);
        abort();
    }
    if (r2<0) {
        __float128 tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    for (int i=0; i<100; i++)
    {
        __float128 pm = 0.5 * (p1 + p2);
        __float128 hh = iapws95_sat_hv_p_quad(pm);
        __float128 r = hh - href;
        printf("%d: p=%.12e, hh=%.5e, r=%.5e, p1=%.12e, p2=%.12e\n", i, (double)pm, (double)hh, (double)(r/href), (double)p1, (double)p2);
        if (fabsq(r)<1e-8*href) {
            return pm;
        }

        if (r > 0)
            p1 = pm;
        else
            p2 = pm;
            if (p2<pmax) break;

        if (fabsq(p1-p2)<1e-8) break;
    }
#endif
    printf("%s: did not converge at h=%.10g\n", __FUNCTION__, args->h);
    abort();
}
#endif

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_p_hv_long(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    double pmin = iapws95_get_pt_approx();
    double pmax = iapws95_get_pc();
    long double p1 = pmin, p2 = pmax;
//    double r1 = iapws95_sat_hv_p_long(p1, eos) - args->h;
//    double r2 = iapws95_h_rhoT(iapws95_get_rhoc(), iapws95_get_Tc()) - args->h;
    const long double h0 = args->h;
    const long double dh = h0*1e-8;
    for (int i=0; i<100; i++)
    {
        long double pm = 0.5 * (p1 + p2);
        long double hh = iapws95_sat_hv_p_long(pm);
        long double r = hh - h0;
        //printf("%d: p=%.10Lg, r=%Lg, p1=%.10Lg, p2=%.10Lg\n", i, pm, r, p1, p2);
        //if (fabs(r) < 1e-8*fabs(args->h)) {
        if (h0-dh<hh && hh<h0+dh) {
            return pm;
        }
        if (fabsl(p1 - p2) < p1*1e-12)
            break;

        if (r > 0)
            p1 = pm;
        else
            p2 = pm;
    }
#ifdef USE_QUAD
    // long double Ts = eos_iapws95_sat_T_hv_long(eos, args);
    // long double ps, rll, rlv;
    // ps = if97_sat_p_T(Ts);
    // double rld, rvd;
    // if97_sat_rho_T(Ts, &rld, &rvd);
    // rll = rld, rlv=rvd;
    // iapws95_sat_prho_T_long(Ts, &ps, &rll, &rlv);
    // return ps;
    __float128 Ts = eos_iapws95_sat_T_hv_quad(eos, args);
    double rld, rvd;
    if97_sat_rho_T(Ts, &rld, &rvd);
    __float128 ps, rll, rlv;
    ps = if97_sat_p_T(Ts);
    rll = rld, rlv=rvd;
    iapws95_sat_prho_T_quad(Ts, &ps, &rll, &rlv);
    return ps;
//    return eos_iapws95_sat_p_hv_quad(eos, args);
#else
    printf("%s: did not converge at h=%.10e\n", __FUNCTION__, args->h);
    abort();
#endif
}
#endif

double eos_iapws95_sat_T_hv(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    double Tmin = iapws95_get_Tt();
    double Tmax = iapws95_get_Tc();
    double T1 = Tmin, T2 = Tmax;
    for (int i = 0; i < 100; i++)
    {
        double Tm = 0.5 * (T1 + T2);
        double hh = iapws95_sat_hv_T(Tm);
        double r = hh - args->h;
        //printf("%d: p=%.10g, r=%g, p1=%.10g, p2=%.10g\n", i, pm, r, p1, p2);
        if (fabs(r) < 1e-8*fabs(args->h)) {
            return Tm;
        }

        if (r > 0)
            T1 = Tm;
        else
            T2 = Tm;
    }
    printf("%s: did not converge at h=%g\n", __FUNCTION__, args->h);
    abort();
}

#ifdef USE_LONGDOUBLE
double eos_iapws95_sat_T_hv_long(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    const long double h0 = args->h;
    long double Tmin = iapws95_get_Tt();
    long double Tmax = iapws95_get_Tc();
    long double T1 = Tmin, T2 = Tmax;
    long double r1 = iapws95_sat_hv_T_long(T1) - h0;
    long double r2 = iapws95_sat_hv_T_long(T2) - h0;
    printf("r1=%Lg, r2=%Lg\n", r1, r2);
    if (r1*r2>0) {
        printf("%s: r1*r2>0\n", __FUNCTION__);
        abort();
    }
    // for (int i = 0; i < 10; i++)
    // {
    //     long double T = Tmax - 1e-7L*i;
    //     long double hh = iapws95_sat_hv_T_long(T, eos);
    //     long double hq = iapws95_sat_hv_T_quad(T, eos);
    //     printf("T=%.14Lg, h=%.14Lg, r=%.14Lg, hq=%.14Lg, rq=%.14Lg\n",  T, hh, hh-h0, hq, hq-h0);
    // }
    const long double dh = h0 * 1e-8;
    for (int i = 0; i < 100; i++)
    {
        long double Tm = 0.5 * (T1 + T2);
        long double hh = iapws95_sat_hv_T_long(Tm);
        long double r = hh - h0;
        printf("%d: T=%.12Lg, hh=%.12Lg, r=%Lg, T1=%.12Lg, T2=%.12Lg\n", i, Tm, hh, r, T1, T2);
        //if (fabs(r) < 1e-8*fabs(args->h)) {
        if (h0 - dh < hh && hh < h0 + dh) {
            return Tm;
        }
        if (fabsl(T1 - T2) < T1*1e-12)
            break;

        if (r > 0)
            T1 = Tm;
        else
            T2 = Tm;
    }
    // printf("%s: did not converge at h=%.10g\n", __FUNCTION__, args->h);
    // abort();

#ifdef USE_QUAD
    return eos_iapws95_sat_T_hv_quad(eos, args);
#else
    return -1;
#endif
}
#endif

#ifdef USE_QUAD
double eos_iapws95_sat_T_hv_quad(void* eos, EOS_ARGS* args)
{
    //printf("%s\n", __FUNCTION__);
    const __float128 h0 = args->h;
    const __float128 Tmin = iapws95_get_Tt();
    const __float128 Tmax = iapws95_get_Tc();
#if 1
    const __float128 dT = Tmax * 1e-10q;
    const __float128 rtol = 1e-10q;
    __float128 T = (args->h < 2900e6) ? (Tmax-1.) : 0.5 * (Tmin + Tmax);
    __float128 du = 0.;
    for (int i = 0; i < 300; i++)
    {
        __float128 hh = iapws95_sat_hv_T_quad(T);
        __float128 r = hh - h0;
        printf("%d: T=%3.12f, dT=%.10e, hh=%3.12f, r=%.3e\n", i, (double)T, (double)du, (double)hh, (double)r/fabs(h0));
        if (fabs(r) < rtol*fabsl(h0))
            return T;

        __float128 drdT = (hh - iapws95_sat_hv_T_quad(T-dT))/dT;
        __float128 scaling = (args->h < 2900e6) ? 0.1q : 1.0q;
//        if (i>100) scaling = 0.1q;
        du = -r/drdT * scaling;
        T += du;
        T = max(Tmin, min(T, Tmax));
    }
#else
    __float128 T1 = Tmin, T2 = Tmax;
    if (h0 < 2090e3)
    {
        __float128 r2 = iapws95_sat_hv_T_quad(T2, eos) - h0;
        printf("r2=%g\n", (double)r2);
        int found = 0;
        for (int i=0; i<100; i++)
        {
            T1 = Tmax - 0.0001q*i;
            __float128 r1 = iapws95_sat_hv_T_quad(T1, eos) - h0;
            printf("T1=%g, r1=%g\n", (double)T1, (double)r1);
            if (r1*r2<0) {
                printf("r1=%g, r2=%g\n", (double)r1, (double)r2);
                found = 1;
                break;
            }
        }
        if (!found) {
            printf("%s: r1*r2>0 at h0=%.12g\n", __FUNCTION__, args->h);
            abort();
        }
    }
    const __float128 dh = h0 * 1e-8;
    for (int i = 0; i < 100; i++)
    {
        __float128 Tm = 0.5 * (T1 + T2);
        __float128 hh = iapws95_sat_hv_T_quad(Tm);
        __float128 r = hh - h0;
        printf("%d: T=%.12g, hh=%.12g, r=%g, T1=%.12g, T2=%.12g\n", i, (double)Tm, (double)hh, (double)r, (double)T1, (double)T2);
        //if (fabs(r) < 1e-8*fabs(args->h)) {
        if (h0 - dh < hh && hh < h0 + dh)
            return Tm;

        if (fabsq(T1 - T2) < T1*1e-12)
            break;

        if (r > 0)
            T1 = Tm;
        else
            T2 = Tm;
    }
#endif
    printf("%s: did not converge at h=%.10g\n", __FUNCTION__, args->h);
    abort();
}
#endif

void eos_iapws95_sat_hlv_T(void* eos, EOS_ARGS* args, double* hl, double*hv)
{
    double sat_rhol, sat_rhov;
    eos_iapws95_sat_rholv_T(eos, args, &sat_rhol, &sat_rhov);
    *hl = iapws95_h_rhoT(sat_rhol, args->T);
    *hv = iapws95_h_rhoT(sat_rhov, args->T);
}

void eos_iapws95_sat_hlv_p(void* eos, EOS_ARGS* args, double* hl, double*hv)
{
    args->T = eos_iapws95_sat_T_p(eos, args);
    eos_iapws95_sat_hlv_T(eos, args, hl, hv);
}


double eos_iapws95_h_crit(void* eos)
{
    return iapws95_get_hc();
}

double eos_iapws95_rhol_tr(void* eos)
{
    return iapws95_get_rholt();
}

double eos_iapws95_rhov_tr(void* eos)
{
    return iapws95_get_rhovt();
}

double eos_iapws95_hl_tr(void* eos)
{
    return iapws95_get_hlt();
}

double eos_iapws95_hv_tr(void* eos)
{
    return iapws95_get_hvt();
}

double eos_iapws95_vis_rhoT(void* eos, EOS_ARGS* args)
{
//    return iapws85_viscosity_rhoT(args->rho, args->T);
    return iapws08_viscosity_rhoT(args->rho, args->T);
}
