
#ifndef _RWRAP_IAPWS95_H_
#define _RWRAP_IAPWS95_H_


void R_iapws95_get_rhoc(double* out);
void R_iapws95_get_Tc(double* out);
void R_iapws95_get_pc(double* out);
void R_iapws95_get_hc(double* out);
void R_iapws95_get_Tt(double* out);
void R_iapws95_get_rholt(double* out);
void R_iapws95_get_rhovt(double* out);
void R_iapws95_get_hlt(double* out);
void R_iapws95_get_hvt(double* out);

void R_iapws95_f_rhoT(double* rho, double* T, double* out);
void R_iapws95_g_rhoT(double* rho, double* T, double* out);
void R_iapws95_u_rhoT(double* rho, double* T, double* out);
void R_iapws95_h_rhoT(double* rho, double* T, double* out);
void R_iapws95_h_pT(double* p, double* T, double* out);
void R_iapws95_s_rhoT(double* rho, double* T, double* out);
void R_iapws95_p_rhoT(double* rho, double* T, double* out);
void R_iapws95_cp_rhoT(double* rho, double* T, double* out);
void R_iapws95_cv_rhoT(double* rho, double* T, double* out);
void R_iapws95_rho_pT(double* p, double* T, double* out);
void R_iapws95_rho_ph(double* p, double* h, double* out);
void R_iapws95_T_rhop(double* rho, double* p, double* out);
void R_iapws95_T_ph(double* p, double* h, double* out);

void R_iapws95_rhoT_ph(double* p, double* h, double* o_rho, double* o_T);

void R_iapws95_wv_ph(double* p, double* h, double* out);

void R_iapws95_sat_p_T(double* T, double* out);
void R_iapws95_sat_T_p(double* p, double* out);
void R_iapws95_sat_rhol_T(double* T, double* out);
void R_iapws95_sat_rhol_p(double* p, double* out);
void R_iapws95_sat_rhov_T(double* T, double* out);
void R_iapws95_sat_rhov_p(double* p, double* out);
void R_iapws95_sat_hl_p(double* p, double* out);
void R_iapws95_sat_hl_T(double* T, double* out);
void R_iapws95_sat_hv_p(double* p, double* out);
void R_iapws95_sat_hv_T(double* T, double* out);
void R_iapws95_sat_prho_T(double* T, double*sat_p, double* sat_rhol, double* sat_rhov, int* error);
void R_iapws95_sat_Trho_p(double* p, double*sat_T, double* sat_rhol, double* sat_rhov, int* error);


#endif
