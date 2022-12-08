
#ifndef _IAPWS_IF97_H_
#define _IAPWS_IF97_H_

double if97_pc();
double if97_Tc();
double if97_rhoc();

double if97_rho_pT(double p, double T);
double if97_rho_ph(double p, double h);
double if97_u_pT(double p, double T);
double if97_h_pT(double p, double T);
double if97_h_pu(double p, double u);
double if97_cp_pT(double p, double T);
double if97_T_pu(double p, double u);
double if97_T_ph(double p, double h);

double if97_sat_p_T(double T);
double if97_sat_T_p(double p);
void if97_sat_rho_p(double pres, double* rhol, double* rhov);
void if97_sat_rho_T(double T, double* rhol, double* rhov);
void if97_sat_h_p(double pres, double* hl, double* hv);
double if97_vapor_quality_ph(double p, double h);

int if97_find_region_ph(double p, double H);

//---------------------------------------------------------
// internal functions
//---------------------------------------------------------

void if97_init();

double if97_cp_pT_region1(double p, double T);
double if97_cp_pT_region2(double p, double T);
double if97_cp_rhoT_region3(double rho, double T);
double if97_cv_rhoT_region3(double rho, double T);

void if97_init_region1();
void if97_init_region2();
void if97_init_region3();

double if97_rho_pT_region1(double p, double T);
double if97_rho_pT_region2(double p, double T);
double if97_rho_pT_region3(double p, double T);

double if97_h_pT_region1(double p, double T);
double if97_h_pT_region2(double p, double T);
double if97_h_pT_region3(double p, double T);

double if97_h_rhoT_region3(double rho, double T);

double if97_T_ph_region1(double p, double h);
double if97_T_ph_region2(double p, double h);
double if97_T_ph_region30(double p, double h);
double if97_T_ph_region31(double p, double h);
double if97_T_ph_region32(double p, double h);

double if97_T_ph_region1_jsme(double p, double h);
double if97_T_ph_region2_jsme(double p, double h);

double if97_p_vT_region3(double v, double T);

double if97_sat_vL_T_region3(double T);
double if97_sat_vG_T_region3(double T);

double if97_boundary23_p_T(double T);
double if97_boundary23_T_p(double p);

#endif
