#ifndef _IAPS84_H_
#define _IAPS84_H_

double iaps84_Tt();
double iaps84_pt();
double iaps84_rholt();
double iaps84_rhovt();
double iaps84_hlt();
double iaps84_hvt();
double iaps84_Tc();
double iaps84_pc();
double iaps84_rhoc();
double iaps84_hc();
double iaps84_f_rhoT(double rhoS, double temp);
double iaps84_g_rhoT(double rhoS, double temp);
double iaps84_u_rhoT(double rhoS, double temp);
double iaps84_h_rhoT(double rhoS, double temp);
double iaps84_s_rhoT(double rhoS, double temp);
double iaps84_p_rhoT(double rhoS, double temp);
double iaps84_cp_rhoT(double rhoS, double temp);
double iaps84_cv_rhoT(double rhoS, double temp);
double iaps84_rho_pT(double pres, double temp);
double iaps84_rho_ph(double pres, double hJkg);
double iaps84_drho_dp_ph(double pres, double hJkg);
double iaps84_drho_dh_ph(double pres, double hJkg);
double iaps84_T_rhop(double rhoS, double pres);
double iaps84_T_ph(double pres, double hJkg);
double iaps84_dT_dp_ph(double pres, double hJkg);
double iaps84_dT_dh_ph(double pres, double hJkg);
double iaps84_viscosity_rhoT(double rhoS, double temp);
double iaps84_thermal_conductivity_rhoT(double rhoS, double temp);
double iaps84_wv_ph(double pres, double hJkg);
double iaps84_satv_ph(double pres, double hJkg);
double iaps84_sat_p_T(double temp);
double iaps84_sat_T_p(double pres);
double iaps84_sat_rhol_T(double temp);
double iaps84_sat_rhol_p(double pres);
double iaps84_sat_rhov_T(double temp);
double iaps84_sat_rhov_p(double pres);
double iaps84_sat_hl_T(double temp);
double iaps84_sat_hl_p(double pres);
double iaps84_sat_hv_T(double temp);
double iaps84_sat_hv_p(double pres);
int iaps84_sat_prholv_T(double temp, double* pres, double* rhoL, double* rhoV);
int iaps84_sat_Trholv_p(double pres, double* temp, double* rhoL, double* rhoV);

#endif
