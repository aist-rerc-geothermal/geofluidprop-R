%module geofluidprop

%{
#define SWIG_FILE_WITH_INIT
#include "model/iapws/IAPWS-95.h"
#include "model/iapws/IAPWS-SAT-92.h"
#include "model/iapws/IAPWS-Melt-11.h"
#include "model/iapws/IAPWS-Viscosity-85.h"
#include "model/iapws/IAPWS-Viscosity-08.h"
#include "model/iapws/IAPWS-ThermalConductivity-11.h"
#include "model/driesner07/Driesner2007_NaCl.h"
#include "model/driesner07/Driesner2007_H2ONaCl.h"
#include "model/klyukinetal17/KlyukinEtAl2017.h"
#ifdef USE_PROST
#include "model/driesner07/Driesner2007_H2ONaCl_prost.h"
#endif
%}

double iapws92_sat_p_T(double T);
double iapws92_sat_T_p(double p);
double iapws92_sat_dp_dT_T(double T);
double iapws92_sat_rhol_T(double T);
double iapws92_sat_drhol_dT_T(double T);
double iapws92_sat_rhov_T(double T);
double iapws92_sat_alpha(double T);
double iapws92_sat_phi(double T);
double iapws92_sat_hl_T(double T);
double iapws92_sat_hv_T(double T);
double iapws92_sat_sl_T(double T);
double iapws92_sat_sg_T(double T);

double iapws95_get_rhoc();
double iapws95_get_Tc();
double iapws95_get_pc();
double iapws95_get_hc();
double iapws95_get_Tt();
double iapws95_get_pt_approx();
double iapws95_get_rholt();
double iapws95_get_rhovt();
double iapws95_get_hlt();
double iapws95_get_hvt();
double iapws95_f_rhoT(double rho, double T);
double iapws95_g_rhoT(double rho, double T);
double iapws95_u_rhoT(double rho, double T);
double iapws95_h_rhoT(double rho, double T);
double iapws95_h_pT(double p, double T);
double iapws95_s_rhoT(double rho, double T);
double iapws95_p_rhoT(double rho, double T);
double iapws95_cp_rhoT(double rho, double T);
double iapws95_cv_rhoT(double rho, double T);
double iapws95_rho_pT(double p, double T);
double iapws95_rho_ph(double p, double h);
double iapws95_T_rhop(double rho, double p);
double iapws95_T_ph(double p, double h);
double iapws95_wv_ph(double p, double h);
double iapws95_sat_p_T(double T);
double iapws95_sat_T_p(double p);
double iapws95_sat_rhol_T(double T);
double iapws95_sat_rhol_p(double p);
double iapws95_sat_rhov_T(double T);
double iapws95_sat_rhov_p(double p);
double iapws95_sat_hl_p(double p);
double iapws95_sat_hl_T(double T);
double iapws95_sat_hv_p(double p);
double iapws95_sat_hv_T(double T);

double iapws85_viscosity_rhoT(double rho, double T);

double iapws08_viscosity_rhoT(double rho, double T);
double iapws08_viscosity_rhoT_simplified(double rho, double T);

double iapws11_thermal_conductivity_rhoT(double rho, double T);

double driesner07_NaCl_VLH_p();
double driesner07_NaCl_VLH_T();
double driesner07_NaCl_VLH_TC();
double driesner07_NaCl_VLH_h_halite();
double driesner07_NaCl_LH_T_p(double p);
double driesner07_NaCl_VH_p_T(double T);
double driesner07_NaCl_VL_p_T(double T);
double driesner07_NaCl_sat_p_T(double T);
double driesner07_NaCl_H_rho_pT(double p_Pa, double T_K);
double driesner07_NaCl_H_h_pT(double p_Pa, double T_K);
double driesner07_NaCl_H_cp_pT(double p_Pa, double T_K);
double driesner07_NaCl_L_rho_pT(double p_Pa, double T_K);
double driesner07_NaCl_L_compressibility_pT(double T_K);
int driesner07_H2O_NaCl_phase_type(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_pc_T(double T);
double driesner07_H2O_NaCl_pc_T2(double T, bool use_satp);
double driesner07_H2O_NaCl_xc_T(double T);
double driesner07_H2O_NaCl_Tc_x(double x);
double driesner07_H2O_NaCl_VL_xl_pT(double p, double T);
double driesner07_H2O_NaCl_VL_xv_pT(double p, double T);
double driesner07_H2O_NaCl_VL_p_Txv(double T_K, double x, double p0);
double driesner07_H2O_NaCl_VL_p_Txl(double T_K, double x, double p0);
double driesner07_H2O_NaCl_VL_rhov_pT(double p_Pa, double T_K);
double driesner07_H2O_NaCl_VL_rhol_pT(double p_Pa, double T_K);
double driesner07_H2O_NaCl_LH_xl_pT(double p, double T);
double driesner07_H2O_NaCl_VH_xv_pT(double p, double T);
double driesner07_H2O_NaCl_VLH_p_T(double T);
double driesner07_H2O_NaCl_VLH_xl_T(double T);
double driesner07_H2O_NaCl_VLH_xv_T(double T);
double driesner07_H2O_NaCl_NaCl_solubility_xv_pT(double p, double T);
double driesner07_H2O_NaCl_NaCl_solubility_xl_pT(double p, double T);
double driesner07_H2O_NaCl_NaCl_solubility_x_pT(double p, double T);
double driesner07_H2O_NaCl_rho_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_vm_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_h_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_cp_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_mass_fraction_liquid(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_mass_fraction_vapor(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_volume_fraction_liquid(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_volume_fraction_vapor(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_volume_fraction_halite(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_T_rhopx(double rho, double p_Pa, double x, double TK0);
double driesner07_H2O_NaCl_T_phx(double p_Pa, double h, double x, double TK0);
double driesner07_H2O_NaCl_volume_fraction_liquid_rhopx(double rho, double p_Pa, double x, double TK0);
double driesner07_H2O_NaCl_rho_pTx(double p_Pa, double T_K, double x);
double driesner07_H2O_NaCl_h_pTx(double p_Pa, double T_K, double x);

double klyukinetal2017_viscosity(double rho, double T, double xNaCl);
