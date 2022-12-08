//-------------------------------------------------------------------
// This code is written based on the following document:
// * Revised Release on the IAPWS Formulation 1995 for the
//   Thermodynamic Properties of Ordinary Water Substance for
//   General and Scientific Use (2016)
//-------------------------------------------------------------------

// Units: u [J/kg], rho [kg/m3], p [Pa], T [K], cp [J/kg/K]

//---------------------------------------------------------
// floating point type setting
//---------------------------------------------------------
#if defined(REAL_AS_LONG_DOUBLE_H)
  #define Real long double
  #define RealSuffix long
#elif defined(REAL_AS_QUAD_H)
  #define Real __float128
  #define RealSuffix quad
#else
  #ifndef REAL_AS_DOUBLE_H
    #define REAL_AS_DOUBLE_H
  #endif
  #define Real double
  #define RealSuffix
#endif

#ifdef REAL_AS_DOUBLE_H
    #define CAT_I(a,b) a
#else
    #define CAT_I(a,b) a##_##b
#endif
#define CAT(a,b) CAT_I(a, b)
#define FNAME(str) CAT(str, RealSuffix)


Real FNAME(iapws95_get_rhoc)();
Real FNAME(iapws95_get_Tc)();
Real FNAME(iapws95_get_pc)();
Real FNAME(iapws95_get_hc)();
//Real FNAME(iapws95_get_hc)();
Real FNAME(iapws95_get_Tt)();
Real FNAME(iapws95_get_pt_approx)();
Real FNAME(iapws95_get_rholt)();
Real FNAME(iapws95_get_rhovt)();
Real FNAME(iapws95_get_hlt)();
Real FNAME(iapws95_get_hvt)();

Real FNAME(iapws95_f_rhoT)(Real rho, Real T);
Real FNAME(iapws95_g_rhoT)(Real rho, Real T);
Real FNAME(iapws95_u_rhoT)(Real rho, Real T);
Real FNAME(iapws95_h_rhoT)(Real rho, Real T);
Real FNAME(iapws95_h_pT)(Real p, Real T);
Real FNAME(iapws95_s_rhoT)(Real rho, Real T);
Real FNAME(iapws95_p_rhoT)(Real rho, Real T);
Real FNAME(iapws95_cp_rhoT)(Real rho, Real T);
Real FNAME(iapws95_cv_rhoT)(Real rho, Real T);
Real FNAME(iapws95_rho_pT)(Real p, Real T);
Real FNAME(iapws95_rho_ph)(Real p, Real h);
Real FNAME(iapws95_T_rhop)(Real rho, Real p);
Real FNAME(iapws95_T_ph)(Real p, Real h);

void FNAME(iapws95_rhoT_ph)(Real p, Real h, Real* o_rho, Real* o_T);

Real FNAME(iapws95_wv_ph)(Real p, Real h);

Real FNAME(iapws95_sat_p_T)(Real T);
Real FNAME(iapws95_sat_T_p)(Real p);
Real FNAME(iapws95_sat_rhol_T)(Real T);
Real FNAME(iapws95_sat_rhol_p)(Real p);
Real FNAME(iapws95_sat_rhov_T)(Real T);
Real FNAME(iapws95_sat_rhov_p)(Real p);
Real FNAME(iapws95_sat_hl_p)(Real p);
Real FNAME(iapws95_sat_hl_T)(Real T);
Real FNAME(iapws95_sat_hv_p)(Real p);
Real FNAME(iapws95_sat_hv_T)(Real T);
int FNAME(iapws95_sat_prho_T)(Real T, Real *sat_p, Real* sat_rhol, Real* sat_rhov);
//void iapws95_sat_prho_T2)(Real T, Real *sat_p, Real* sat_rhol, Real* sat_rhov);
int FNAME(iapws95_sat_Trho_p)(Real p, Real *sat_T, Real* sat_rhol, Real* sat_rhov);

Real FNAME(iapws95_dh_drho_rhoT)(Real rho, Real T);
Real FNAME(iapws95_dh_dT_rhoT)(Real rho, Real T);
Real FNAME(iapws95_dh_dp_rhopT)(Real rho, Real p, Real T);
Real FNAME(iapws95_dp_drho_rhoT)(Real rho, Real T);
Real FNAME(iapws95_dp_dT_rhoT)(Real rho, Real T);
Real FNAME(iapws95_drho_dh_ph)(Real p, Real h);
Real FNAME(iapws95_drho_dT_pT)(Real p, Real T, Real rho0);
Real FNAME(iapws95_drho_dp_pT)(Real p, Real T, Real rho0);

// internal functions
Real FNAME(iapws95_phi0)(Real delta, Real tau);
Real FNAME(iapws95_dphi0_ddelta)(Real delta, Real tau);
Real FNAME(iapws95_d2phi0_ddelta2)(Real delta, Real tau);
Real FNAME(iapws95_dphi0_dtau)(Real delta, Real tau);
Real FNAME(iapws95_d2phi0_dtau2)(Real delta, Real tau);
Real FNAME(iapws95_d2phi0_ddelta_dtau)(Real delta, Real tau);
Real FNAME(iapws95_phir)(Real delta, Real tau);
Real FNAME(iapws95_dphir_ddelta)(Real delta, Real tau);
Real FNAME(iapws95_d2phir_ddelta2)(Real delta, Real tau);
Real FNAME(iapws95_dphir_dtau)(Real delta, Real tau);
Real FNAME(iapws95_d2phir_dtau2)(Real delta, Real tau);
Real FNAME(iapws95_d2phir_ddelta_dtau)(Real delta, Real tau);

#undef Real
#undef RealSuffix
#undef CAT_I
#undef CAT
#undef FNAME
