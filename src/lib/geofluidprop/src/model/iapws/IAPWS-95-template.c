// #include "IAPWS-95.h"
#if defined(REAL_AS_DOUBLE)
#define REAL_AS_DOUBLE_H
#elif defined(REAL_AS_LONG_DOUBLE)
#define REAL_AS_LONG_DOUBLE_H
#elif defined(REAL_AS_QUAD)
#define REAL_AS_QUAD_H
#else
#define REAL_TYPE_UNDEFINED
#endif

#ifndef REAL_TYPE_UNDEFINED

#include "IAPWS-95-template.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef REAL_AS_QUAD
#include <quadmath.h>
#endif

#include "util/utility.h"
#include "util/gauss3.h"

#include "IAPWS-IF97.h"

//---------------------------------------------------------
// floating point type setting
//---------------------------------------------------------
#if defined(REAL_AS_LONG_DOUBLE)
#define Real long double
#define RealSuffix long
#elif defined(REAL_AS_QUAD)
#define Real __float128
#define RealSuffix quad
#else
#ifndef REAL_AS_DOUBLE
#define REAL_AS_DOUBLE
#endif
#define Real double
#define RealSuffix
#endif

// #define DO_EXPAND(VAL)  VAL ## 1
// #define EXPAND(VAL)     DO_EXPAND(VAL)

// suffix to function name
#ifdef REAL_AS_DOUBLE
    #define CAT_I(a,b) a
#else
    #define CAT_I(a,b) a##_##b
#endif
#define CAT(a,b) CAT_I(a, b)

// suffix to float point constat
#if defined(REAL_AS_DOUBLE)
#define c(x) x
#elif defined(REAL_AS_LONG_DOUBLE)
#define c(x) x##l
#elif defined(REAL_AS_QUAD)
#define c(x) x##q
#endif
//   #define c(x) _Generic((x), long double: x##l, default: x)

//---------------------------------------------------------
// math functions
//---------------------------------------------------------
#if defined(REAL_AS_DOUBLE)
#define ABS(x) fabs(x)
#define EXP(x) exp(x)
#define LOG(x) log(x)
#define POW(a,b) pow(a,b)
#define GAUSS3(a,b,c) gauss3(a,b,c)

#elif defined(REAL_AS_LONG_DOUBLE)
#define ABS(x) fabsl(x)
#define EXP(x) expl(x)
#define LOG(x) logl(x)
#define POW(a,b) powl(a,b)
#define GAUSS3(a,b,c) gauss3_long(a,b,c)

#elif defined(REAL_AS_QUAD)
#define ABS(x) fabsq(x)
#define EXP(x) expq(x)
#define LOG(x) logq(x)
#define POW(a,b) powq(a,b)
#define GAUSS3(a,b,c) gauss3_quad(a,b,c)
#endif
// #define ABS(x) _Generic((x), long double: fabsl, __float128: fabsq, default: fabs)(x)

//---------------------------------------------------------
// internal parameters
//---------------------------------------------------------
// tolerance
#define RTOL_DOUBLE 1e-8
#ifdef _MSC_VER
#define RTOL_LONG RTOL_DOUBLE
#else
#define RTOL_LONG 1e-10l
#endif
#define RTOL_QUAD 1e-12q

#if defined(REAL_AS_DOUBLE)
#define RTOL RTOL_DOUBLE
#elif defined(REAL_AS_LONG_DOUBLE)
#define RTOL RTOL_LONG
#elif defined(REAL_AS_QUAD)
#define RTOL RTOL_QUAD
#endif


//-------------------------------------------------------------------
// Thermodynamic constants
//-------------------------------------------------------------------

#if defined(REAL_AS_DOUBLE)
#include "IAPWS-95-table.h"
#elif defined(REAL_AS_LONG_DOUBLE)
#include "IAPWS-95-tablel.h"
#elif defined(REAL_AS_QUAD)
#include "IAPWS-95-tableq.h"
#endif


// triple point
#define Tt            c(273.16)   // [K]
#define pt_approx     c(611.6547711)  // [Pa]
#define ht_approx     c(0.611782) // [J/kg]
#define rholt_approx  c(999.793)    // [kg/m^3]
#define rhovt_approx  c(0.00485458) // [kg/m^3]
// critical point
#define Tc    c(647.096)  // [K]
#define pc    c(22.064e6) // [Pa]
#define rhoc  c(322.0)  // [kg/m^3]
#define hc_approx    c(2084.25625591e3) // J/kg
// specific gas constant of water
#define R     c(0.46151805e3) // [J/kg/K]

#define iapws95_get_rhoc CAT(iapws95_get_rhoc, RealSuffix)
#define iapws95_get_Tc CAT(iapws95_get_Tc, RealSuffix)
#define iapws95_get_pc CAT(iapws95_get_pc, RealSuffix)
#define iapws95_get_Tt CAT(iapws95_get_Tt, RealSuffix)
#define iapws95_get_pt_approx CAT(iapws95_get_pt_approx, RealSuffix)

Real iapws95_get_rhoc() {return rhoc;}
Real iapws95_get_Tc() {return Tc;}
Real iapws95_get_pc() {return pc;}
Real iapws95_get_Tt() {return Tt;}
Real iapws95_get_pt_approx() {return pt_approx;}

//-------------------------------------------------------------------
// Function DELTA and its derivatives
//-------------------------------------------------------------------
#define iapws95_theta CAT(iapws95_theta, RealSuffix)
Real iapws95_theta(Real delta, Real tau, int i)
{
    return (c(1.)-tau) + A[i]*POW((delta-c(1.))*(delta-c(1.)), c(0.5)/beta[i]);
}

#define iapws95_DELTA CAT(iapws95_DELTA, RealSuffix)
Real iapws95_DELTA(Real delta, Real tau, int i)
{
    //Real th = theta(delta, tau, i);
    const Real dm1 = delta-c(1.);
    const Real delta_m1_2 = dm1*dm1;
    const Real th = (c(1.)-tau) + A[i]*POW(delta_m1_2, c(0.5)/beta[i]);
    Real DELTA = th*th + B[i]*POW(delta_m1_2, a[i]);
    return DELTA;
}

#define iapws95_dDELTA_ddelta CAT(iapws95_dDELTA_ddelta, RealSuffix)
Real iapws95_dDELTA_ddelta(Real delta, Real theta, int i)
{
    const Real dm1 = delta-c(1.);
    if (dm1==0.0) return 0;
    const Real delta_m1_2 = dm1*dm1;
    const Real ibeta = c(1.)/beta[i];
    Real ret = dm1*(A[i]*theta*c(2.)*ibeta*POW(delta_m1_2,c(0.5)*ibeta-c(1.)) + c(2.)*B[i]*a[i]*POW(delta_m1_2,a[i]-c(1.)));
    assert(!isnan(ret));
    return ret;
}

#define iapws95_d2DELTA_ddelta2 CAT(iapws95_d2DELTA_ddelta2, RealSuffix)
Real iapws95_d2DELTA_ddelta2(Real delta, Real theta, int i)
{
    //assert(delta!=1);
    const Real delta_m1 = delta-c(1.);
    if (delta_m1==0.0) return 0;
    const Real delta_m1_2 = delta_m1*delta_m1;
    const Real ai = a[i];
    const Real Ai = A[i];
    const Real ibetai = c(1.)/beta[i];
    Real dd = c(1.)/delta_m1*iapws95_dDELTA_ddelta(delta,theta,i)
                + delta_m1_2*(
                    c(4.)*B[i]*ai*(ai-c(1.))*POW(delta_m1_2,ai-c(2.))
                    + c(2.)*Ai*Ai*ibetai*ibetai*POW(POW(delta_m1_2,c(0.5)*ibetai-c(1.)),c(2.))
                    + Ai*theta*c(4.)*ibetai*(c(0.5)*ibetai-c(1.))*POW(delta_m1_2,c(0.5)*ibetai-c(2.))
                );
    assert(!isnan(dd));
    return dd;
}

//-------------------------------------------------------------------
// Function DELTA^b[i] and its derivatives
//-------------------------------------------------------------------

#define iapws95_DELTAb CAT(iapws95_DELTAb, RealSuffix)
Real iapws95_DELTAb(Real DELTA, int i)
{
    return POW(DELTA, b[i]);
}

#define iapws95_dDELTAb_ddelta CAT(iapws95_dDELTAb_ddelta, RealSuffix)
Real iapws95_dDELTAb_ddelta(Real delta, Real theta, Real DELTA, int i)
{
    if (DELTA==0.0) return 0;
    Real ret = b[i]*POW(DELTA, b[i]-c(1.))*iapws95_dDELTA_ddelta(delta, theta, i);
    assert(!isnan(ret));
    return ret;
}

#define iapws95_d2DELTAb_ddelta2 CAT(iapws95_d2DELTAb_ddelta2, RealSuffix)
Real iapws95_d2DELTAb_ddelta2(Real delta, Real theta, Real DELTA, int i)
{
    if (DELTA==0.0) return 0;
    Real ret = b[i]*(POW(DELTA, b[i]-c(1.))*iapws95_d2DELTA_ddelta2(delta, theta, i)
            + (b[i]-c(1.))*POW(DELTA,b[i]-c(2.))*POW(iapws95_dDELTA_ddelta(delta, theta, i),c(2.)));
    assert(!isnan(ret));
    return ret;
}

#define iapws95_dDELTAb_dtau CAT(iapws95_dDELTAb_dtau, RealSuffix)
Real iapws95_dDELTAb_dtau(Real delta, Real theta, Real DELTA, int i)
{
    if (DELTA==0.0) return 0;
    return -c(2.)*theta*b[i]*POW(DELTA,b[i]-c(1.));
}

#define iapws95_d2DELTAb_dtau2 CAT(iapws95_d2DELTAb_dtau2, RealSuffix)
Real iapws95_d2DELTAb_dtau2(Real delta, Real theta, Real DELTA, int i)
{
    if (DELTA==0.0) return 0;
    return (c(2.)*b[i]*POW(DELTA,b[i]-c(1.)) + c(4.)*theta*theta*b[i]*(b[i]-c(1.))*POW(DELTA,b[i]-c(2.)));
}

#define iapws95_d2DELTAb_ddelta_dtau CAT(iapws95_d2DELTAb_ddelta_dtau, RealSuffix)
Real iapws95_d2DELTAb_ddelta_dtau(Real delta, Real theta, Real DELTA, int i)
{
    if (DELTA==0.0) return 0;
    const Real delta_m1 = delta-c(1.);
    const Real delta_m1_2 = delta_m1*delta_m1;
    const Real bi = b[i];
    const Real ibetai = c(1.)/beta[i];
    const Real pow_DELTA_bi = POW(DELTA,bi-c(1.));
    Real dd = -A[i] * bi * c(2.) * ibetai * pow_DELTA_bi * delta_m1 *
                  POW(delta_m1_2, c(0.5) * ibetai - c(1.)) -
              c(2.) * theta * bi * (bi - c(1.)) * pow_DELTA_bi / DELTA *
                  iapws95_dDELTA_ddelta(delta, theta, i);
    // Real dd = -A[i]*bi*c(2.)*ibetai*POW(DELTA,bi-c(1.))*(delta-c(1.))*POW(delta_m1_2,
    // 0.5*ibetai-c(1.))
    //             -c(2.)*theta*bi*(bi-c(1.))*POW(DELTA,bi-2)*dDELTA_ddelta(delta,
    //             theta, i);
    assert(!isnan(dd));
    return dd;
}

//-------------------------------------------------------------------
// Function psi and its derivatives
//-------------------------------------------------------------------

#define iapws95_psi CAT(iapws95_psi, RealSuffix)
Real iapws95_psi(Real delta, Real tau, int i)
{
    const Real dm1 = delta-c(1.);
    const Real tm1 = tau-c(1.);
    return EXP(-C[i]*dm1*dm1-D[i]*tm1*tm1);
}

#define iapws95_dpsi_ddelta CAT(iapws95_dpsi_ddelta, RealSuffix)
Real iapws95_dpsi_ddelta(Real delta, Real psi, int i)
{
    return -c(2.)*C[i]*(delta-c(1.))*psi;
}

#define iapws95_d2psi_ddelta2 CAT(iapws95_d2psi_ddelta2, RealSuffix)
Real iapws95_d2psi_ddelta2(Real delta, Real psi, int i)
{
    const Real dm1 = delta-c(1.);
    return (c(2.)*C[i]*dm1*dm1 -c(1.))*c(2.)*C[i]*psi;
}

#define iapws95_dpsi_dtau CAT(iapws95_dpsi_dtau, RealSuffix)
Real iapws95_dpsi_dtau(Real tau, Real psi, int i)
{
    return -c(2.)*D[i]*(tau-c(1.))*psi;
}

#define iapws95_d2psi_dtau2 CAT(iapws95_d2psi_dtau2, RealSuffix)
Real iapws95_d2psi_dtau2(Real tau, Real psi, int i)
{
    return (c(2.)*D[i]*(tau-c(1.))*(tau-c(1.)) -c(1.))*c(2.)*D[i]*psi;
}

#define iapws95_d2psi_ddelta_dtau CAT(iapws95_d2psi_ddelta_dtau, RealSuffix)
Real iapws95_d2psi_ddelta_dtau(Real delta, Real tau, Real psi, int i)
{
    return c(4.)*C[i]*D[i]*(delta-c(1.))*(tau-c(1.))*psi;
}

//-------------------------------------------------------------------
// ideal gas part
//-------------------------------------------------------------------
#define iapws95_phi0 CAT(iapws95_phi0, RealSuffix)
Real iapws95_phi0(Real delta, Real tau)
{
    Real phi0 = LOG(delta) + n0[0] + n0[1]*tau + n0[2]*LOG(tau);
    for (int i=3; i<8; i++)
        phi0 += n0[i]*LOG(c(1.)-EXP(-gamma0[i]*tau));
    return phi0;
}

#define iapws95_dphi0_ddelta CAT(iapws95_dphi0_ddelta, RealSuffix)
Real iapws95_dphi0_ddelta(Real delta, Real tau)
{
    return c(1.)/delta;
}

#define iapws95_d2phi0_ddelta2 CAT(iapws95_d2phi0_ddelta2, RealSuffix)
Real iapws95_d2phi0_ddelta2(Real delta, Real tau)
{
    return -c(1.)/(delta*delta);
}

#define iapws95_dphi0_dtau CAT(iapws95_dphi0_dtau, RealSuffix)
Real iapws95_dphi0_dtau(Real delta, Real tau)
{
    Real dphi0 = n0[1] + n0[2]/tau;
    for (int i=3; i<8; i++)
        dphi0 += n0[i]*gamma0[i]*(c(1.)/(c(1.)-EXP(-gamma0[i]*tau))-c(1.));
    return dphi0;
}

#define iapws95_d2phi0_dtau2 CAT(iapws95_d2phi0_dtau2, RealSuffix)
Real iapws95_d2phi0_dtau2(Real delta, Real tau)
{
    Real dphi0 = - n0[2]/(tau*tau);
    for (int i=3; i<8; i++)
        dphi0 -= n0[i]*gamma0[i]*gamma0[i]*EXP(-gamma0[i]*tau)*POW(c(1.)-EXP(-gamma0[i]*tau), -c(2.));
    return dphi0;
}

#define iapws95_d2phi0_ddelta_dtau CAT(iapws95_d2phi0_ddelta_dtau, RealSuffix)
Real iapws95_d2phi0_ddelta_dtau(Real delta, Real tau)
{
    return 0.0;
}

//-------------------------------------------------------------------
// residual part
//-------------------------------------------------------------------

#define iapws95_phir CAT(iapws95_phir, RealSuffix)
Real iapws95_phir(Real delta, Real tau)
{
    // the residual part
#if 1
    Real pow_delta_i[17]; // from -1 to 15
    pow_delta_i[0] = 1./delta;
    for (int i=1; i<17; i++)
        pow_delta_i[i] = delta * pow_delta_i[i-1];

    Real phir = 0.0;
    for (int i=0; i<7; i++)
        phir += n[i]*pow_delta_i[d[i]+1]*POW(tau,t[i]);
    for (int i=7; i<51; i++)
        phir += n[i]*pow_delta_i[d[i]+1]*POW(tau,t[i])*EXP(-pow_delta_i[c[i]+1]);
    for (int i=51; i<54; i++)
        phir += n[i]*pow_delta_i[d[i]+1]*POW(tau,t[i])*EXP(-alpha[i]*(delta-eps[i])*(delta-eps[i])-beta[i]*(tau-gamma_[i])*(tau-gamma_[i]));
    const Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
    for (int i=54; i<56; i++)
    {
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        phir += n[i]*POW(DELTA,b[i])*delta*psi;
    }
#else
    Real phir = 0.0;
    for (int i=0; i<7; i++)
        phir += n[i]*POW(delta,d[i])*POW(tau,t[i]);
    for (int i=7; i<51; i++)
        phir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-POW(delta, c[i]));
    for (int i=51; i<54; i++)
        phir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)));
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        phir += n[i]*POW(DELTA,b[i])*delta*psi;
    }
#endif
    return phir;
}

#define iapws95_dphir_ddelta CAT(iapws95_dphir_ddelta, RealSuffix)
Real iapws95_dphir_ddelta(Real delta, Real tau)
{
#if 1
    Real pow_delta_i[17]; // from -1 to 15
    pow_delta_i[0] = c(1.)/delta;
    for (int i=1; i<17; i++)
        pow_delta_i[i] = delta * pow_delta_i[i-1];
    Real pow_tau_i[51]; // 0 to 50
    pow_tau_i[0] = 1.;
    for (int i = 1; i < 51; i++)
        pow_tau_i[i] = tau * pow_tau_i[i - 1];

    Real dphir = 0.0;
    for (int i=0; i<7; i++)
    {
        const int di = d[i];
        dphir += n[i]*di*pow_delta_i[di-1+1]*POW(tau,t[i]);
    }
    for (int i=7; i<51; i++)
    {
        const int ci = c[i];
        const int di = d[i];
        dphir += n[i]*EXP(-pow_delta_i[ci+1])*pow_delta_i[di-1+1]*POW(tau,t[i])*(di-ci*pow_delta_i[ci+1]);
    }
    for (int i=51; i<54; i++)
    {
        const Real de = delta - eps[i];
        const Real tg = tau - gamma_[i];
        dphir += n[i]*pow_delta_i[d[i]+1]*pow_tau_i[ti[i]]*EXP(-alpha[i]*de*de-beta[i]*tg*tg)
                *(d[i]/delta -c(2.)*alpha[i]*de);
    }
    const Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
    for (int i=54; i<56; i++)
    {
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        if (DELTA!=0.0)
        dphir += n[i]*(
                    POW(DELTA,b[i])*(psi+delta*iapws95_dpsi_ddelta(delta, psi, i))
                    +iapws95_dDELTAb_ddelta(delta, theta, DELTA, i)*delta*psi
                    );
    }
#else
    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*d[i]*POW(delta,d[i]-c(1.))*POW(tau,t[i]);
    for (int i=7; i<51; i++)
        dphir += n[i]*EXP(-POW(delta, c[i]))*POW(delta,d[i]-c(1.))*POW(tau,t[i])*(d[i]-c[i]*POW(delta,c[i]));
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*(delta-eps[i])*(delta-eps[i])-beta[i]*(tau-gamma_[i])*(tau-gamma_[i]))
                *(d[i]/delta -c(2.)*alpha[i]*(delta-eps[i]));
    }
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        dphir += n[i]*(
                    POW(DELTA,b[i])*(psi+delta*iapws95_dpsi_ddelta(delta, psi, i))
                    +iapws95_dDELTAb_ddelta(delta, theta, DELTA, i)*delta*psi
                    );
    }
#endif
    return dphir;
}


#define iapws95_d2phir_ddelta2 CAT(iapws95_d2phir_ddelta2, RealSuffix)
Real iapws95_d2phir_ddelta2(Real delta, Real tau)
{
#if 1
    Real pow_delta_i[17]; // from -1 to 15
    pow_delta_i[0] = c(1.)/delta;
    for (int i=1; i<17; i++)
        pow_delta_i[i] = delta * pow_delta_i[i-1];

    Real dphir = 0.0;
    for (int i=0; i<7; i++)
    {
        const int di = d[i];
        dphir += n[i]*di*(di-c(1.))*pow_delta_i[di-2+1]*POW(tau,t[i]);
    }
    for (int i=7; i<51; i++)
    {
        const int ci = c[i];
        const int di = d[i];
        dphir += n[i]*EXP(-pow_delta_i[ci+1])
                    *(
                        pow_delta_i[di-2+1]*POW(tau,t[i])
                        *(
                            (di-ci*pow_delta_i[ci+1])
                            *(di-1-ci*pow_delta_i[ci+1])
                            -ci*ci*pow_delta_i[ci+1]
                        )
                    );
    }
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)))
                *(-c(2.)*alpha[i]*pow_delta_i[d[i]+1] + c(4.)*alpha[i]*alpha[i]*pow_delta_i[d[i]+1]*POW(delta-eps[i],c(2.))
                  -c(4.)*d[i]*alpha[i]*pow_delta_i[d[i]-1+1]*(delta-eps[i]) + d[i]*(d[i]-c(1.))*pow_delta_i[d[i]-2+1]
                );
    }
    const Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
    for (int i=54; i<56; i++)
    {
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        dphir += n[i] *
                 (POW(DELTA, b[i]) *
                      (c(2.) * iapws95_dpsi_ddelta(delta, psi, i) +
                       delta * iapws95_d2psi_ddelta2(delta, psi, i)) +
                  c(2.) * iapws95_dDELTAb_ddelta(delta, theta, DELTA, i) *
                      (psi + delta * iapws95_dpsi_ddelta(delta, psi, i)) +
                  iapws95_d2DELTAb_ddelta2(delta, theta, DELTA, i) *
                      delta * psi);
    }
    return dphir;

#else
    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*d[i]*(d[i]-c(1.))*POW(delta,d[i]-2)*POW(tau,t[i]);
    for (int i=7; i<51; i++)
        dphir += n[i]*EXP(-POW(delta,c[i]))
                    *(
                        POW(delta,d[i]-2)*POW(tau,t[i])
                        *(
                            (d[i]-c[i]*POW(delta,c[i]))
                            *(d[i]-1-c[i]*POW(delta,c[i]))
                            -c[i]*c[i]*POW(delta,c[i])
                        )
                    );
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)))
                *(-c(2.)*alpha[i]*POW(delta,d[i]) + c(4.)*alpha[i]*alpha[i]*POW(delta,d[i])*POW(delta-eps[i],c(2.))
                  -4*d[i]*alpha[i]*POW(delta,d[i]-c(1.))*(delta-eps[i]) + d[i]*(d[i]-c(1.))*POW(delta, d[i]-2)
                );
    }
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        dphir += n[i]*(
                    POW(DELTA,b[i])*(c(2.)*dpsi_ddelta(delta, psi, i) + delta*d2psi_ddelta2(delta, psi, i))
                    +2*dDELTAb_ddelta(delta, theta, DELTA, i)*(psi + delta*dpsi_ddelta(delta, psi, i))
                    +d2DELTAb_ddelta2(delta, theta, DELTA, i)*delta*psi
                    );
    }

    return dphir;
#endif
}


#define iapws95_dphir_dtau CAT(iapws95_dphir_dtau, RealSuffix)
Real iapws95_dphir_dtau(Real delta, Real tau)
{
#if 1
    Real pow_delta_i[17]; // from -1 to 15 -> (0 to 16)
    pow_delta_i[0] = 1./delta;
    for (int i=1; i<17; i++)
        pow_delta_i[i] = delta * pow_delta_i[i-1];

    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*t[i]*pow_delta_i[d[i]+1]*POW(tau,t[i]-c(1.));
    for (int i=7; i<51; i++)
        dphir += n[i]*t[i]*pow_delta_i[d[i]+1]*POW(tau,t[i]-c(1.))*EXP(-pow_delta_i[c[i]+1]);
    for (int i=51; i<54; i++)
    {
        const Real tg = tau - gamma_[i];
        const Real de = delta - eps[i];
        dphir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*de*de-beta[i]*tg*tg)
                *(t[i]/tau -c(2.)*beta[i]*tg);
    }
    const Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
    for (int i=54; i<56; i++)
    {
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        if (DELTA!=0.0)
        dphir += n[i]*delta*(
                    iapws95_dDELTAb_dtau(delta, theta, DELTA, i)*psi + POW(DELTA,b[i])*iapws95_dpsi_dtau(tau, psi, i)
                    );
    }
#else
    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*t[i]*POW(delta,d[i])*POW(tau,t[i]-c(1.));
    for (int i=7; i<51; i++)
        dphir += n[i]*t[i]*POW(delta,d[i])*POW(tau,t[i]-c(1.))*EXP(-POW(delta, c[i]));
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)))
                *(t[i]/tau -c(2.)*beta[i]*(tau-gamma_[i]));
    }
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        dphir += n[i]*delta*(
                    dDELTAb_dtau(delta, theta, DELTA, i)*psi + POW(DELTA,b[i])*dpsi_dtau(tau, psi, i)
                    );
    }
#endif
    return dphir;
}


#define iapws95_d2phir_dtau2 CAT(iapws95_d2phir_dtau2, RealSuffix)
Real iapws95_d2phir_dtau2(Real delta, Real tau)
{
    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*t[i]*(t[i]-c(1.))*POW(delta,d[i])*POW(tau,t[i]-c(2.));
    for (int i=7; i<51; i++)
        dphir += n[i]*t[i]*(t[i]-c(1.))*POW(delta,d[i])*POW(tau,t[i]-c(2.))*EXP(-POW(delta, c[i]));
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)))
                *(POW(t[i]/tau -c(2.)*beta[i]*(tau-gamma_[i]),c(2.)) -t[i]/(tau*tau)-c(2.)*beta[i]);
    }
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        dphir += n[i]*delta*(
                    iapws95_d2DELTAb_dtau2(delta, theta, DELTA, i)*psi + c(2.)*iapws95_dDELTAb_dtau(delta, theta, DELTA, i)*iapws95_dpsi_dtau(tau, psi, i) + POW(DELTA,b[i])*iapws95_d2psi_dtau2(tau, psi, i)
                    );
    }

    return dphir;
}

#define iapws95_d2phir_ddelta_dtau CAT(iapws95_d2phir_ddelta_dtau, RealSuffix)
Real iapws95_d2phir_ddelta_dtau(Real delta, Real tau)
{
    Real dphir = 0.0;
    for (int i=0; i<7; i++)
        dphir += n[i]*d[i]*t[i]*POW(delta,d[i]-c(1.))*POW(tau,t[i]-c(1.));
    for (int i=7; i<51; i++)
        dphir += n[i]*t[i]*POW(delta,d[i]-c(1.))*POW(tau,t[i]-c(1.))*(d[i]-c[i]*POW(delta,c[i]))*EXP(-POW(delta, c[i]));
    for (int i=51; i<54; i++)
    {
        dphir += n[i]*POW(delta,d[i])*POW(tau,t[i])*EXP(-alpha[i]*POW(delta-eps[i],c(2.))-beta[i]*POW(tau-gamma_[i],c(2.)))
                *(d[i]/delta -c(2.)*alpha[i]*(delta-eps[i]))*(t[i]/tau -c(2.)*beta[i]*(tau-gamma_[i]));
    }
    for (int i=54; i<56; i++)
    {
        Real delta_m_1_2 = (delta-c(1.))*(delta-c(1.));
        Real theta = (c(1.)-tau) + A[i]*POW(delta_m_1_2, c(0.5)/beta[i]);
        Real DELTA = theta*theta + B[i]*POW(delta_m_1_2, a[i]);
        Real psi = EXP(-C[i]*delta_m_1_2-D[i]*(tau-c(1.))*(tau-c(1.)));
        Real tmp =
            (DELTA==0 ? 0 : POW(DELTA, b[i])) *
                (iapws95_dpsi_dtau(tau, psi, i) +
                 delta * iapws95_d2psi_ddelta_dtau(delta, tau, psi, i)) +
            delta * iapws95_dDELTAb_ddelta(delta, theta, DELTA, i) *
                iapws95_dpsi_dtau(tau, psi, i) +
            iapws95_dDELTAb_dtau(delta, theta, DELTA, i) *
                (psi + delta * iapws95_dpsi_ddelta(delta, psi, i)) +
            iapws95_d2DELTAb_ddelta_dtau(delta, theta, DELTA, i) * delta *
                psi;

        dphir +=
            n[i] * tmp;
            // (POW(DELTA, b[i]) *
            //      (iapws95_dpsi_dtau(tau, psi, i) +
            //       delta * iapws95_d2psi_ddelta_dtau(delta, tau, psi, i)) +
            //  delta * iapws95_dDELTAb_ddelta(delta, theta, DELTA, i) *
            //      iapws95_dpsi_dtau(tau, psi, i) +
            //  iapws95_dDELTAb_dtau(delta, theta, DELTA, i) *
            //      (psi + delta * iapws95_dpsi_ddelta(delta, psi, i)) +
            //  iapws95_d2DELTAb_ddelta_dtau(delta, theta, DELTA, i) * delta *
            //      psi);
    }
    assert(!isnan(dphir));

    return dphir;
}


//-------------------------------------------------------------------
// thermodynamic properties
//-------------------------------------------------------------------

#define iapws95_f_rhoT CAT(iapws95_f_rhoT, RealSuffix)
#define iapws95_g_rhoT CAT(iapws95_g_rhoT, RealSuffix)
#define iapws95_p_rhoT CAT(iapws95_p_rhoT, RealSuffix)
#define iapws95_p_rhoT_s CAT(iapws95_p_rhoT_s, RealSuffix)
#define iapws95_u_rhoT CAT(iapws95_u_rhoT, RealSuffix)
#define iapws95_s_rhoT CAT(iapws95_s_rhoT, RealSuffix)
#define iapws95_h_rhoT CAT(iapws95_h_rhoT, RealSuffix)
#define iapws95_h_pT CAT(iapws95_h_pT, RealSuffix)
#define iapws95_cv_rhoT CAT(iapws95_cv_rhoT, RealSuffix)
#define iapws95_cp_rhoT CAT(iapws95_cp_rhoT, RealSuffix)
#define iapws95_beta_s_rhoT CAT(iapws95_beta_s_rhoT, RealSuffix)

#define iapws95_rho_pT CAT(iapws95_rho_pT, RealSuffix)
#define iapws95_rho_ph CAT(iapws95_rho_ph, RealSuffix)
#define iapws95_T_rhop CAT(iapws95_T_rhop, RealSuffix)
#define iapws95_T_ph CAT(iapws95_T_ph, RealSuffix)

#define iapws95_rhoT_ph CAT(iapws95_rhoT_ph, RealSuffix)

#define iapws95_sat_p_T CAT(iapws95_sat_p_T, RealSuffix)
#define iapws95_sat_T_p CAT(iapws95_sat_T_p, RealSuffix)
#define iapws95_sat_rhol_T CAT(iapws95_sat_rhol_T, RealSuffix)
#define iapws95_sat_rhol_p CAT(iapws95_sat_rhol_p, RealSuffix)
#define iapws95_sat_rhov_T CAT(iapws95_sat_rhov_T, RealSuffix)
#define iapws95_sat_rhov_p CAT(iapws95_sat_rhov_p, RealSuffix)
#define iapws95_sat_hl_T CAT(iapws95_sat_hl_T, RealSuffix)
#define iapws95_sat_hl_p CAT(iapws95_sat_hl_p, RealSuffix)
#define iapws95_sat_hv_T CAT(iapws95_sat_hv_T, RealSuffix)
#define iapws95_sat_hv_p CAT(iapws95_sat_hv_p, RealSuffix)

#define iapws95_sat_prho_T CAT(iapws95_sat_prho_T, RealSuffix)
#define iapws95_sat_Trho_p CAT(iapws95_sat_Trho_p, RealSuffix)

#define iapws95_wv_ph CAT(iapws95_wv_ph, RealSuffix)
#define iapws95_wv_rhoT CAT(iapws95_wv_rhoT, RealSuffix)

#define iapws95_get_hc CAT(iapws95_get_hc, RealSuffix)
#define iapws95_get_pt CAT(iapws95_get_pt, RealSuffix)
#define iapws95_get_rholt CAT(iapws95_get_rholt, RealSuffix)
#define iapws95_get_rhovt CAT(iapws95_get_rhovt, RealSuffix)
#define iapws95_get_hlt CAT(iapws95_get_hlt, RealSuffix)
#define iapws95_get_hvt CAT(iapws95_get_hvt, RealSuffix)

#define iapws95_df_dT_rhoT CAT(iapws95_df_dT_rhoT, RealSuffix)
#define iapws95_du_drho_rhoT CAT(iapws95_du_drho_rhoT, RealSuffix)
#define iapws95_du_dT_rhoT CAT(iapws95_du_dT_rhoT, RealSuffix)
#define iapws95_drho_dp_pT CAT(iapws95_drho_dp_pT, RealSuffix)
#define iapws95_drho_dh_ph CAT(iapws95_drho_dh_ph, RealSuffix)
#define iapws95_drho_dT_pT CAT(iapws95_drho_dT_pT, RealSuffix)
#define iapws95_dp_drho_rhoT CAT(iapws95_dp_drho_rhoT, RealSuffix)
#define iapws95_dp_dT_rhoT CAT(iapws95_dp_dT_rhoT, RealSuffix)
#define iapws95_dh_drho_rhoT CAT(iapws95_dh_drho_rhoT, RealSuffix)
#define iapws95_dh_dT_rhoT CAT(iapws95_dh_dT_rhoT, RealSuffix)
#define iapws95_dh_dp_rhopT CAT(iapws95_dh_dp_rhopT, RealSuffix)

#define iapws95_r12 CAT(iapws95_r12, RealSuffix)
#define iapws95_r3 CAT(iapws95_r3, RealSuffix)
#define iapws95_dr12_dp CAT(iapws95_dr12_dp, RealSuffix)
#define iapws95_dr12_dT CAT(iapws95_dr12_dT, RealSuffix)
#define iapws95_dr12_drho CAT(iapws95_dr12_drho, RealSuffix)
#define iapws95_dr3_dp CAT(iapws95_dr3_dp, RealSuffix)
#define iapws95_dr3_dT CAT(iapws95_dr3_dT, RealSuffix)
#define iapws95_dr3_drhol CAT(iapws95_dr3_drhol, RealSuffix)
#define iapws95_dr3_drhov CAT(iapws95_dr3_drhov, RealSuffix)

Real iapws95_f_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    return (iapws95_phi0(delta, tau) + iapws95_phir(delta, tau)) * R * T;
}

Real iapws95_df_dT_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real df = iapws95_phi0(delta, tau) + iapws95_phir(delta, tau);
    df += - tau*(iapws95_dphi0_dtau(delta,tau) + iapws95_dphir_dtau(delta,tau));
    df *= R;
    return df;
}

Real iapws95_g_rhoT(Real rho, Real T)
{
//    return (iapws95_f_rhoT(rho,T) + iapws95_p_rhoT(rho,T)/rho);
    const Real delta = rho/rhoc;
    const Real tau = Tc/T;
    Real fRT = iapws95_phi0(delta, tau) + iapws95_phir(delta, tau);
    Real prhoRT = c(1.)+delta*iapws95_dphir_ddelta(delta, tau);
    Real gRT = fRT + prhoRT;
    return gRT*R*T;
}

Real iapws95_p_rhoT_s(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real p = c(1.)+delta*iapws95_dphir_ddelta(delta, tau);
    p *= rho*R*T;
    return p;
}

Real iapws95_p_rhoT(Real rho, Real T)
{
    assert(Tt <= T);
    if (T < Tc)
    {
        double rhold, rhovd;
        if97_sat_rho_T(T, &rhold, &rhovd);
        if (rho < rhold*1.1 && rho > rhovd*0.9)
        {
            Real ps, rhol, rhov;
            if (!iapws95_sat_prho_T(T, &ps, &rhol, &rhov))
                return -1;
            if (rhov < rho && rho < rhol)
                return ps;
        }
    }
    return iapws95_p_rhoT_s(rho, T);
}

Real iapws95_dp_drho_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real dp = c(1.) + c(2.)*delta*iapws95_dphir_ddelta(delta,tau);
    dp += delta*delta*iapws95_d2phir_ddelta2(delta,tau);
    // Real dp = 1+delta*dphir_ddelta(delta, tau);
    // dp += rho*(c(1.)/rhoc*dphir_ddelta(delta,tau)+delta/rhoc*d2phir_ddelta2(delta,tau));
    return dp*R*T;
}

Real iapws95_dp_dT_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real dp = c(1.)+delta*iapws95_dphir_ddelta(delta, tau);
    dp += -delta*tau*iapws95_d2phir_ddelta_dtau(delta,tau);
    dp *= rho*R;
    return dp;
}

Real iapws95_u_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;

    return tau*(iapws95_dphi0_dtau(delta, tau) + iapws95_dphir_dtau(delta, tau))*R*T;
}

Real iapws95_du_drho_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real du = iapws95_d2phi0_ddelta_dtau(delta, tau) + iapws95_d2phir_ddelta_dtau(delta, tau);
    du *= tau*R*T/rhoc;
    return du;
}

Real iapws95_du_dT_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real du = iapws95_d2phi0_dtau2(delta, tau) + iapws95_d2phir_dtau2(delta, tau);
    du *= -R*tau*tau;
    return du;
}

Real iapws95_s_rhoT(Real rho, Real T)
{
    return (iapws95_u_rhoT(rho,T)-iapws95_f_rhoT(rho,T))/T;
    // return (tau*(dphi0_dtau(delta, tau) + dphir_dtau(delta, tau))
    //         -phi0(delta,tau)-phir(delta,tau))*R;
}

Real iapws95_h_rhoT(Real rho, Real T)
{
    assert(Tt <= T);
    if (T < Tc)
    {
        Real ps, rhol, rhov;
        if (!iapws95_sat_prho_T(T, &ps, &rhol, &rhov))
            return -1.;
        if (rhov < rho && rho < rhol)
        {
            Real x = calc_vapor_quality(1./rhol, 1./rhov, 1./rho);
            Real sat_hl = iapws95_h_rhoT(rhol, T);
            Real sat_hv = iapws95_h_rhoT(rhov, T);
            return (sat_hl + x * (sat_hv - sat_hl));
        }
    }

    Real delta = rho / rhoc;
    Real tau = Tc / T;
    Real h = c(1.)+tau*(iapws95_dphi0_dtau(delta, tau) + iapws95_dphir_dtau(delta, tau))
             +  delta * iapws95_dphir_ddelta(delta,tau);
    h *= R * T;
    //Real h = iapws95_u_rhoT(rho,T) + iapws95_p_rhoT(rho,T)/rho;
    return h;
}

Real iapws95_h_pT(Real p, Real T)
{
    Real rho = iapws95_rho_pT(p, T);
    return iapws95_h_rhoT(rho, T);
}

Real iapws95_dh_drho_rhoT(Real rho, Real T)
{
    Real dh = iapws95_du_drho_rhoT(rho,T) - iapws95_p_rhoT(rho,T)/(rho*rho)+iapws95_dp_drho_rhoT(rho,T)/rho;
    return dh;
}

Real iapws95_dh_dT_rhoT(Real rho, Real T)
{
    Real dh = iapws95_du_dT_rhoT(rho,T) + iapws95_dp_dT_rhoT(rho,T)/rho;
    return dh;
}

Real iapws95_dh_dp_rhopT(Real rho, Real p, Real T)
{
    Real dh = iapws95_dh_drho_rhoT(rho, T)*iapws95_drho_dp_pT(p, T, rho);
    return dh;
}

Real iapws95_cv_rhoT(Real rho, Real T)
{
    return iapws95_du_dT_rhoT(rho,T);
}

Real iapws95_cp_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;

    Real v = POW(c(1.)+delta*iapws95_dphir_ddelta(delta,tau)-delta*tau*iapws95_d2phir_ddelta_dtau(delta,tau),c(2.));
    v /= (c(1.)+c(2.)*delta*iapws95_dphir_ddelta(delta,tau)+delta*delta*iapws95_d2phir_ddelta2(delta,tau));
    v *= R;
    v += iapws95_cv_rhoT(rho,T);
    return v;
    // return (cv(rho,T)
    //     + POW(1+delta*dphir_ddelta(delta,tau)-delta*tau*d2phir_ddelta_dtau(delta,tau),c(2.))
    //       /(1+2*delta*dphir_ddelta(delta,tau)+delta*delta*d2phir_ddelta2(delta,tau))*R);
}

Real iapws95_beta_s_rhoT(Real rho, Real T)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real v = (c(1.)+delta*iapws95_dphir_ddelta(delta,tau)-delta*tau*iapws95_d2phir_ddelta_dtau(delta,tau))
                / (
                    POW(c(1.)+delta*iapws95_dphir_ddelta(delta,tau)+delta*tau*iapws95_d2phir_ddelta_dtau(delta,tau),c(2.))
                    - tau*tau*(iapws95_d2phi0_dtau2(delta,tau)+iapws95_d2phir_dtau2(delta,tau))*(c(1.)+c(2.)*delta*iapws95_dphir_ddelta(delta,tau)+delta*delta*iapws95_d2phir_ddelta2(delta,tau))
                );
    return v;
}

//-------------------------------------------------------------------
// back-analysis
//-------------------------------------------------------------------

Real iapws95_rho_pT(Real p, Real T)
{
    if (T == iapws95_get_Tc() && p == iapws95_get_pc())
        return rhoc;
    if (T < iapws95_get_Tc() && p < iapws95_get_pc())
    {
        Real sat_T = iapws95_sat_T_p(p);
        if (ABS(T - sat_T) < c(1.e-8)*sat_T)
        {
            printf("***warning in %s(): two phase at p=%e MPa, T=%g C. cannot determine rho(p,T)\n", __FUNCTION__, (double)(p*1e-6), (double)(T-273.15));
            return -1; // two phase, impossible to determine rho
        }
    }

    const Real rho0 = if97_rho_pT(p*c(1e-6), T);
    Real rho = rho0;
    for (int i=0; i<50; i++)
    {
        Real pi = iapws95_p_rhoT(rho, T);
        Real r = pi - p;
#ifdef DEBUG_RHO_PT
        printf("\t %d: rho=%g, pi=%.5e, r=%.5e\n", i, (double)rho, (double)pi, (double)r);
#endif
        if (ABS(r) < p*c(1e-8))
            return rho;
        Real J = iapws95_dp_drho_rhoT(rho, T);
#ifdef DEBUG_RHO_PT
//        printf("\t J=%g\n", dp_drho);
        //printf("\t %d: rho=%.3e, pi=%.3e, r=%.3e, J=%.3e\n", i, rho, pi, r, dp_drho);
#endif
        rho -= r/J;
        rho = rho > 0 ? rho : c(1e-3);
    }
    printf("***ERROR: iteration did not converge in %s() at p=%g, T=%g, rho0=%g\n", __FUNCTION__, (double)p, (double)T, (double)rho0);
    abort();
#undef DEBUG_RHO_PT
    return -1;
}

Real iapws95_T_rhop(Real rho_, Real p_)
{
    // ad-hoc sampling for initial guess
    Real T = 300;
    Real T0 = T;
    Real r_min = 1e99;
    const Real dT = c(50.);
    const Real n = (c(1000.)-c(300.))/dT;
    for (int i=0; i<n; i++)
    {
        T = c(300.) + dT*i;
        Real r_p = iapws95_p_rhoT(rho_, T) - p_;
        if (ABS(r_p) < r_min) {
            T0 = T;
            r_min = ABS(r_p);
            //printf("\t T=%.3e, r=%.3e\n", T, r_min);
        }
    }
    // iteration
    T = T0;
    const Real rtol = RTOL;
    for (int i=0; i<100; i++)
    {
        Real pi = iapws95_p_rhoT(rho_, T);
        Real r_p = pi - p_;
        if (ABS(r_p) < p_*rtol)
        {
            return T;
        }
        Real J = iapws95_dp_dT_rhoT(rho_, T);
        //printf("\t %d: T=%.3e, pi=%.3e, r=%.3e, J=%.3e, err=%.3e\n", i, (double)T, (double)pi, (double)r_p, (double)J, (double)ABS(r_p/p_));
        T += -r_p/J;
        T = T >= 0 ? T : Tt; //TODO
    }
    printf("***ERROR: iteration did not converge in T_rhop()\n");
    return -1;
}


void iapws95_rhoT_ph(Real p_, Real h_, Real* o_rho, Real* o_T)
{
    if (p_ < pc)
    {
        // check two phase
        Real Ts, rhol, rhov;
        iapws95_sat_Trho_p(p_, &Ts, &rhol, &rhov);
        Real hl = iapws95_h_rhoT(rhol, Ts);
        Real hv = iapws95_h_rhoT(rhov, Ts);
        if (hl < h_ && h_ < hv)
        {
            Real wv = (h_ - hl) / (hv - hl);
            Real v = 1./rhol + wv * (1./rhov - 1./rhol);
            *o_rho = 1./v;
            *o_T = Ts;
            return;
        }
    }

    // single phase: solve rho, T which satisfy
    //  h=f(rho,T)
    //  p=f(rho,T)
    const Real rho0 = if97_rho_ph(p_*c(1e-6), h_*c(1e-3));
    const Real T0 = if97_T_ph(p_*c(1e-6), h_*c(1e-3));
    Real rho = rho0, T = T0;
    const Real tol = RTOL;
    // Real rp_hist[100], rh_hist[100];
    int i=0;
    for (i=0; i<100; i++)
    {
        Real r_p = iapws95_p_rhoT(rho, T) - p_;
        Real r_h = iapws95_h_rhoT(rho, T) - h_;

        // rp_hist[i] = r_p/p_;
        // rh_hist[i] = r_h/h_;
        if (ABS(r_p) < p_*tol && ABS(r_h) < h_*tol)
        {
            *o_rho = rho;
            *o_T = T;
            return;
        }

        Real Jp_rho = iapws95_dp_drho_rhoT(rho, T);
        Real Jp_T = iapws95_dp_dT_rhoT(rho, T);
        Real Jh_rho = iapws95_dh_drho_rhoT(rho, T);
        Real Jh_T = iapws95_dh_dT_rhoT(rho, T);
        Real drho = c(1.)/(Jh_rho-Jp_rho*Jh_T/Jp_T)*(r_h-Jh_T/Jp_T*r_p);
        Real dT = c(1.)/Jp_T*(r_p-Jp_rho*drho);
        //printf("\t %d: rp=%.3e, rh=%.3e, rho=%g, T=%g\n", i, (double)(r_p/p_), (double)(r_h/h_), (double)rho, (double)T);

        rho -= drho;
        T -= dT;
        rho = MAX(c(1e-5), rho);
        T = MAX(T, Tt);
    }

    *o_rho = -1;
    *o_T = -1;
    // printf("***ERROR: %s did not converge for p=%g, h=%g\n", __FUNCTION__, (double)p_, (double)h_);
    // //printf("***ERROR: rhoT_ph() did not converge for p=%g, h=%g, rho0=%g, T0=%g\n", p_d, h_d, rho0, T0);
    // for (int j=i-11; j<i; j++)
    //     printf("%d: rp=%g, rh=%g\n", j, (double)rp_hist[j], (double)rh_hist[j]);
    // abort();
}

Real iapws95_rho_ph(Real p, Real h)
{
    Real rho, T;
    iapws95_rhoT_ph(p, h, &rho, &T);
    return rho;
}

Real iapws95_T_ph(Real p, Real h)
{
    Real rho, T;
    iapws95_rhoT_ph(p, h, &rho, &T);
    return T;
}

Real iapws95_drho_dh_ph(Real p, Real h)
{
    Real d = h*c(1.e-8);
    Real rho1, rho2, T;
    iapws95_rhoT_ph(p, h+d, &rho1, &T);
    iapws95_rhoT_ph(p, h-d, &rho2, &T);
    Real drho = (rho1 - rho2)/(c(2.)*d);
    return drho;
}

Real iapws95_drho_dT_pT(Real p, Real T, Real rho0)
{
    Real dT = T*c(1.e-8);
    Real drho = (iapws95_rho_pT(p, T+dT) - iapws95_rho_pT(p, T-dT))/(c(2.)*dT);
    return drho;
}

Real iapws95_drho_dp_pT(Real p, Real T, Real rho0)
{
    Real dp = p*c(1.e-8);
    Real drho = (iapws95_rho_pT(p+dp, T) - iapws95_rho_pT(p-dp, T))/(c(2.)*dp);
    return drho;
}

//-------------------------------------------------------------------
// Phase boundaries and two-phase region
//-------------------------------------------------------------------

Real iapws95_r12(Real T, Real ps, Real rho)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    return (ps/(R*T*rho)-delta*iapws95_dphir_ddelta(delta,tau)-c(1.));
//    return (ps/(R*T*rho)-c(1.)-delta*iapws95_dphir_ddelta(delta,tau));
}

Real iapws95_dr12_dp(Real T, Real ps, Real rho)
{
    return c(1.)/(R*T*rho);
}

Real iapws95_dr12_dT(Real T, Real ps, Real rho)
{
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real dr = -ps/(R*T*T*rho);
    dr += delta*tau/T*iapws95_d2phir_ddelta_dtau(delta,tau);
    return dr;
}

Real iapws95_dr12_drho(Real T, Real ps, Real rho)
{
#if 0
    Real d = rho*c(1e-8);
    return (iapws95_r12(T,ps,rho+0.5*d)-iapws95_r12(T,ps,rho-0.5*d))/d;
#else
    Real delta = rho/rhoc;
    Real tau = Tc/T;
    Real dr = -( ps/(R*T*rho*rho) + c(1.)/rhoc*iapws95_dphir_ddelta(delta,tau));
    dr += - rho/(rhoc*rhoc) * iapws95_d2phir_ddelta2(delta,tau);
    return dr;
#endif
}

Real iapws95_r3(Real T, Real ps, Real rhol, Real rhov)
{
    Real deltal = rhol/rhoc;
    Real deltav = rhov/rhoc;
    Real tau = Tc/T;
    Real r = ps/(R*T)*(c(1.)/rhov-c(1.)/rhol) - LOG(rhol/rhov);
    r += -iapws95_phir(deltal,tau) + iapws95_phir(deltav,tau);
    return r;
}

Real iapws95_dr3_dp(Real T, Real ps, Real rhol, Real rhov)
{
//    return c(1.)/(R*T)*(c(1.)/rhov-c(1.)/rhol);
    // (1/rhov-1/rhol)/RT
    // = (rhol-rhov)/(rhov*rhol*R*T)
    return (rhol-rhov)/(rhol*rhov*R*T);
}

Real iapws95_dr3_dT(Real T, Real ps, Real rhol, Real rhov)
{
    Real deltal = rhol/rhoc;
    Real deltav = rhov/rhoc;
    Real tau = Tc/T;
    Real dr = -ps/(R*T*T)*(c(1.)/rhov-c(1.)/rhol);
    dr += tau/T*(iapws95_dphir_dtau(deltal,tau) - iapws95_dphir_dtau(deltav,tau));
    return dr;
}

Real iapws95_dr3_drhol(Real T, Real ps, Real rhol, Real rhov)
{
#if 0
    Real d = rhol*c(1e-8);
    return (iapws95_r3(T,ps,rhol+0.5*d,rhov)-iapws95_r3(T,ps,rhol-0.5*d,rhov))/d;
#else
    Real deltal = rhol/rhoc;
    Real tau = Tc/T;
    Real dr = ps/(R*T*rhol*rhol) - c(1.)/rhol;
    dr += -iapws95_dphir_ddelta(deltal,tau)/rhoc;
    return dr;
#endif
}

Real iapws95_dr3_drhov(Real T, Real ps, Real rhol, Real rhov)
{
#if 0
    Real d = rhov*c(1e-8);
    return (iapws95_r3(T,ps,rhol,rhov+0.5*d)-iapws95_r3(T,ps,rhol,rhov-0.5*d))/d;
#else
    Real deltav = rhov/rhoc;
    Real tau = Tc/T;
    Real dr = -ps/(R*T*rhov*rhov) + c(1.)/rhov;
    dr += iapws95_dphir_ddelta(deltav,tau)/rhoc;
    return dr;
#endif
}

int iapws95_sat_prho_T(Real T, Real *sat_p, Real* sat_rhol, Real* sat_rhov)
{
    if (T==Tc) {
        *sat_p = pc;
        *sat_rhol = *sat_rhov = rhoc;
        return 1;
    }
    Real p0 = if97_sat_p_T(T)*c(1e6); // (*sat_p>0 && *sat_p<1e9) ? *sat_p : 1e5;
    double rhold0, rhovd0;
    if97_sat_rho_T(T, &rhold0, &rhovd0);
    Real rhol0 = rhold0; // (rhoc <= *sat_rhol && *sat_rhol < 2.e3) ? *sat_rhol : 1.e3;
    Real rhov0 = rhovd0; // (1.e-5<*sat_rhov && *sat_rhov <= rhoc) ? *sat_rhov : 1.e-03;
    Real ps = p0;
    Real rhol = rhol0;
    Real rhov = rhov0;
    Real J[9], r[3];
    Real du[3]; // rhol, rhov, ps
    int converged = 0;
    const Real rtol = RTOL;
    // NR loop
    for (int i=0; i<20; i++)
    {
        r[0] = -iapws95_r12(T, ps, rhol);
        r[1] = -iapws95_r12(T, ps, rhov);
        r[2] = -iapws95_r3(T, ps, rhol, rhov);
        // Real pl_tmp = iapws95_p_rhoT(rhol, T);
        // Real pv_tmp = iapws95_p_rhoT(rhov, T);
        Real rpl = iapws95_p_rhoT_s(rhol, T) - ps;
        Real rpv = iapws95_p_rhoT_s(rhov, T) - ps;
        Real g1 = iapws95_g_rhoT(rhol, T);
        Real rg = g1 - iapws95_g_rhoT(rhov, T);
//#define DEBUG_SATP
#ifdef DEBUG_SATP
        printf("%d: u=(%.10e,%.10e,%.10e), g=%.3e, r=(%.3e,%.3e,%.3e), du=(%.3e,%.3e,%.3e)\n",
            i, (double)rhol, (double)rhov, (double)ps, (double)g1, (double)rpl/(double)ps, (double)(rpv/ps), (double)(rg/ABS(g1)),
            (double)du[0], (double)du[1], (double)du[2]);
        printf("\t r=(%.3e,%.3e,%.3e)\n", (double)r[0], (double)r[1], (double)r[2]);
//        printf("\t r=(%.3e,%.3e,%.2e), p=(%g, %g)\n", (double)rpl, (double)rpv, (double)rg, (double)pl_tmp, (double)pv_tmp);
//        printf("%d: u=(%.10e,%.10e,%.10e), r=(%g,%g,%g)\n", i, (double)rhol, (double)rhov, (double)ps, (double)r[0], (double)r[1], (double)r[2]);
#endif
//        if (ABS(r[0])<rtol && ABS(r[1])<rtol && ABS(r[2])<rtol)
        if (ABS(rpl)<rtol*ps && ABS(rpv)<rtol*ps && ABS(rg)<rtol*ABS(g1))
        {
            converged = 1;
            break;
        }
        // if (i==0 && ABS(r[0])<1e-6 && ABS(r[1])<1e-6 && ABS(r[2])<1e-6)
        //     pmin = MAX(pt, p0-c(100.0));

        for (int j=0;j<9;j++)
            J[j] = 0;
        J[0] = iapws95_dr12_drho(T, ps, rhol);
        J[2] = iapws95_dr12_dp(T, ps, rhol);
        J[3*1+1] = iapws95_dr12_drho(T, ps, rhov);
        J[3*1+2] = iapws95_dr12_dp(T, ps, rhov);
        J[3*2+0] = iapws95_dr3_drhol(T, ps, rhol, rhov);
        J[3*2+1] = iapws95_dr3_drhov(T, ps, rhol, rhov);
        J[3*2+2] = iapws95_dr3_dp(T, ps, rhol, rhov);

#ifdef DEBUG_SATP
        printf("\t J:\n");
        for (int j=0; j<3; j++)
            printf("\t%.5e\t%.5e\t%.5e\n", (double)J[j*3], (double)J[j*3+1], (double)J[j*3+2]);
#endif

        GAUSS3(J,r,du);
        Real scaling = 1.0 ; //(i<30) ? c(1.0) : c(0.1);
        rhol += scaling * du[0];
        rhov += scaling * du[1];
        ps += scaling * du[2];

#if 0
        if (rhol == rhov) {
            if (ps == pc && T==Tc)
                    break;
            rhol = rhoc + c(1e-4);
            rhov = rhoc - c(1e-4);
        }
#endif
        // rhol = max(rhol, rhoc);
        rhol = MAX(rhol, rhoc+c(1.e-10));
        rhov = MIN(MAX(rhov, c(1e-5)), rhoc-c(1.e-10));
        ps = MIN(pc, MAX(ps, pt_approx - c(1.e-8)));
    }
    if (!converged)
    {
        printf("***ERROR: NR loop did not converge in %s(): T=%.12g, p0=%.10e, rhol0=%e, rhov0=%e\n", __FUNCTION__, (double)T, (double)p0, (double)rhol0, (double)rhov0);
        *sat_p = *sat_rhol = *sat_rhov = -1;
        return 0;
        //abort();
    }

    *sat_p = ps;
    *sat_rhol = MAX(rhol, rhov);
    *sat_rhov = MIN(rhol, rhov);
    return 1;
}

int iapws95_sat_Trho_p(Real ps, Real *sat_T, Real* sat_rhol, Real* sat_rhov)
{
    if (ps > iapws95_get_pc()) {
        printf("%s(): Outside limits of p <= pCr. Given p = %g\n", __FUNCTION__, (double)ps);
        *sat_T = *sat_rhol = *sat_rhov = -1;
        return 0;
    }
    if (ps==pc) {
        *sat_T = Tc;
        *sat_rhol = *sat_rhov = rhoc;
        return 1;
    }

    Real T0 = if97_sat_T_p(ps*1e-6);
    double rhold, rhovd;
    if97_sat_rho_p(ps*1e-6, &rhold, &rhovd);
    Real rhol0 = rhold, rhov0 = rhovd;

#if 0
    if (T0==iapws95_get_Tc())
    {
        T0 -= c(1.e-5);
        rhol0 = iapws95_get_rhoc() * (c(1.)+c(1.e-5));
        rhov0 = iapws95_get_rhoc() * (c(1.)-c(1.e-5));
    }
    if (rhol0==rhov0)
    {
        rhol0 = iapws95_get_rhoc() * (c(1.)+c(1.e-5));
        rhov0 = iapws95_get_rhoc() * (c(1.)-c(1.e-5));
    }
#endif
    Real T = T0;
    Real rhol = rhol0;
    Real rhov = rhov0;
    Real J[9], r[3];
    Real du[3]; // rhol, rhov, T
    bool converged = false;
    const Real rtol = RTOL;
    Real Tmin = Tt;
    // NR loop
    for (int i=0; i<300; i++)
    {
        r[0] = -iapws95_r12(T, ps, rhol);
        r[1] = -iapws95_r12(T, ps, rhov);
        r[2] = -iapws95_r3(T, ps, rhol, rhov);
        Real rpl = iapws95_p_rhoT_s(rhol, T) - ps;
        Real rpv = iapws95_p_rhoT_s(rhov, T) - ps;
        Real g1 = iapws95_g_rhoT(rhol, T);
        Real rg = g1 - iapws95_g_rhoT(rhov, T);
//#define DEBUG_SAT
#ifdef DEBUG_SAT
        printf("%d: u=(%.8e,%.8e,%.8e), r=(%.3e,%.3e,%.3e)\n", i, (double)rhol, (double)rhov, (double)T, (double)rpl, (double)rpv, (double)rg);
        //printf("%d: u=(%.8e,%.8e,%.8e), r=(%.3e,%.3e,%.3e)\n", i, (double)rhol, (double)rhov, (double)T, (double)r[0], (double)r[1], (double)r[2]);
        // printf("%d: T=%g, r=(%s, %s, %s)\n", i, (double)T, buf1, buf2, buf3);
#endif
        const Real r12_ref = 1.0;// ps/(R*T*1.e3);
//        if (ABS(r[0])<rtol*r12_ref && ABS(r[1])<rtol*r12_ref && ABS(r[2])<rtol)
        if (ABS(rpl)<rtol*ABS(ps) && ABS(rpv)<rtol*ABS(ps) && ABS(rg)<rtol*ABS(g1))
        {
            converged = true;
            if (T==Tc)
                rhol = rhov = rhoc;
            break;
        }
        if (i==0 && ABS(r[0])<c(1e-6) && ABS(r[1])<c(1e-6) && ABS(r[2])<c(1e-6))
            Tmin = MAX(Tt, T0-c(10.0));

        for (int j=0;j<9;j++)
            J[j] = 0;
        J[0] = iapws95_dr12_drho(T, ps, rhol);
        J[2] = iapws95_dr12_dT(T, ps, rhol);
        J[3*1+1] = iapws95_dr12_drho(T, ps, rhov);
        J[3*1+2] = iapws95_dr12_dT(T, ps, rhov);
        J[3*2+0] = iapws95_dr3_drhol(T, ps, rhol, rhov);
        J[3*2+1] = iapws95_dr3_drhov(T, ps, rhol, rhov);
        J[3*2+2] = iapws95_dr3_dT(T, ps, rhol, rhov);

// #ifdef DEBUG_SAT
//         printf("\t J:\n");
//         for (int j=0; j<3; j++)
//         {
//             printf("\t%.5e\t%.5e\t%.5e\n", (double)J[j*3], (double)J[j*3+1], (double)J[j*3+2]);
//         }
// #endif

        GAUSS3(J,r,du);
// #ifdef DEBUG_SAT
//         printf("\t du: %.5e\t%.5e\t%.5e\n", (double)du[0], (double)du[1], (double)du[2]);
// #endif

        rhol += du[0];
        rhov += du[1];
        T += du[2];

        if (rhol<rhov) {
            Real tmp = rhol;
            rhol = rhov;
            rhov = tmp;
        }
        rhol = MAX(rhol, rhoc);
        rhov = MAX(rhov, c(1e-5));
        rhov = MIN(rhov, rhoc);
        T = MIN(Tc, MAX(T, Tmin));
        //if (T<0) T = Tt;
        if (isnan(rhol)) rhol = c(900.);
        if (isnan(T) || isnan(rhol) || isnan(rhov))
            break;
    }
    if (!converged)
    {
        // printf("***ERROR: NR loop did not converge in %s(): p=%.10Lg, T0=%.10Lg, rhol0=%.10Lg, rhov0=%.10Lg\n", __FUNCTION__, ps, T0, *sat_rhol, *sat_rhov);
        // printf("res: %g, %g, %g\n", (double)r[0], (double)r[1], (double)r[2]);
        *sat_T = *sat_rhol = *sat_rhov = 0;
        return 0;
        //abort();
    }

    *sat_T = T;
    *sat_rhol = MAX(rhol, rhov);
    *sat_rhov = MIN(rhol, rhov);
    return 1;
}

Real iapws95_sat_p_T(Real Ts)
{
    Real ps, rhol, rhov;
    if (!iapws95_sat_prho_T(Ts, &ps, &rhol, &rhov))
        return -1;
    return ps;
}

Real iapws95_sat_rhol_T(Real Ts)
{
    Real ps, rhol, rhov;
    if (!iapws95_sat_prho_T(Ts, &ps, &rhol, &rhov))
        return -1;
    return rhol;
}

Real iapws95_sat_rhov_T(Real Ts)
{
    Real ps, rhol, rhov;
    if (!iapws95_sat_prho_T(Ts, &ps, &rhol, &rhov))
        return -1;
    return rhov;
}

Real iapws95_sat_hl_T(Real Ts)
{
    Real ps, rhol, rhov;
    if (!iapws95_sat_prho_T(Ts, &ps, &rhol, &rhov))
        return -1;
    return iapws95_h_rhoT(rhol, Ts);
}

Real iapws95_sat_hv_T(Real Ts)
{
    Real ps, rhol, rhov;
    if (!iapws95_sat_prho_T(Ts, &ps, &rhol, &rhov))
        return -1;
    return iapws95_h_rhoT(rhov, Ts);
}

Real iapws95_sat_T_p(Real ps)
{
    Real Ts, rhol, rhov;
    if (!iapws95_sat_Trho_p(ps, &Ts, &rhol, &rhov))
        return -1;
    return Ts;
}

Real iapws95_sat_rhol_p(Real ps)
{
    Real Ts, rhol, rhov;
    if (!iapws95_sat_Trho_p(ps, &Ts, &rhol, &rhov))
        return -1;
    return rhol;
}

Real iapws95_sat_rhov_p(Real ps)
{
    Real Ts, rhol, rhov;
    if (!iapws95_sat_Trho_p(ps, &Ts, &rhol, &rhov))
        return -1;
    return rhov;
}

Real iapws95_sat_hl_p(Real ps)
{
    Real Ts, rhol, rhov;
    if (!iapws95_sat_Trho_p(ps, &Ts, &rhol, &rhov))
        return -1;
    return iapws95_h_rhoT(rhol, Ts);
}

Real iapws95_sat_hv_p(Real ps)
{
    Real Ts, rhol, rhov;
    if (!iapws95_sat_Trho_p(ps, &Ts, &rhol, &rhov))
        return -1;
    return iapws95_h_rhoT(rhov, Ts);
}

Real iapws95_wv_ph(Real p, Real h)
{
    if (p >= iapws95_get_pc())
        return 0.0;

    Real T = iapws95_T_ph(p, h);
    if (T >= iapws95_get_Tc())
        return 0.0;

    Real sat_T, sat_rhol, sat_rhov;
    if (!iapws95_sat_Trho_p(p, &sat_T, &sat_rhol, &sat_rhov))
        return -1;

    if (T < sat_T)
        return 0.0;
    if (T > sat_T)
        return 1.0;

    Real sat_hl = iapws95_h_rhoT(sat_rhol, sat_T);
    if (h <= sat_hl)
        return 0.0;

    Real sat_hv = iapws95_h_rhoT(sat_rhov, sat_T);
    if (sat_hv <= h)
        return 1.0;

    // two phase case
    return calc_vapor_quality(sat_hl, sat_hv, h);
}

Real iapws95_wv_rhoT(Real rho, Real T)
{
    if (T >= iapws95_get_Tc())
        return 0.0;

    {
        double rhol, rhov;
        if97_sat_rho_T(T, &rhol, &rhov);
        if (rho > 1.1*rhol) return 0.0;
        if (rho < 0.9*rhov) return 1.0;
    }

    Real sat_p, sat_rhol, sat_rhov;
    if (!iapws95_sat_prho_T(T, &sat_p, &sat_rhol, &sat_rhov))
        return -1;

    if (rho >= sat_rhol) return 0.0;
    if (sat_rhov <= rho) return 1.0;

    // two phase case
    return calc_vapor_quality(1./sat_rhol, 1./sat_rhov, 1./rho);
}

Real iapws95_get_hc()
{
    return iapws95_h_rhoT(rhoc, Tc);
}

Real iapws95_get_pt()
{
    return iapws95_sat_p_T(Tt);
}

Real iapws95_get_rholt()
{
    return iapws95_sat_rhol_T(Tt);
}

Real iapws95_get_rhovt()
{
    return iapws95_sat_rhov_T(Tt);
}

Real iapws95_get_hlt()
{
    return iapws95_sat_hl_T(Tt);
}

Real iapws95_get_hvt()
{
    return iapws95_sat_hv_T(Tt);
}

#endif
