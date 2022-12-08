
#include "rwrap_iapws95.h"

#include "model/iapws/IAPWS-95.h"


void R_iapws95_get_rhoc(double* out)
{
  *out = iapws95_get_rhoc();
}

void R_iapws95_get_Tc(double*out)
{
  *out = iapws95_get_Tc();
}

void R_iapws95_get_pc(double*out)
{
  *out = iapws95_get_pc();
}

void R_iapws95_get_hc(double* out)
{
  *out = iapws95_get_hc();
}
void R_iapws95_get_Tt(double* out)
{
  *out = iapws95_get_Tt();
}
void R_iapws95_get_rholt(double* out)
{
  *out = iapws95_get_rholt();
}
void R_iapws95_get_rhovt(double* out)
{
  *out = iapws95_get_rhovt();

}
void R_iapws95_get_hlt(double* out)
{
  *out = iapws95_get_hlt();

}
void R_iapws95_get_hvt(double* out)
{
  *out = iapws95_get_hvt();
}

void R_iapws95_get_rho_pT(double*p, double*T, double*out)
{
  *out = iapws95_get_rho_pT(*p, *T);
}

void R_iapws95_f_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_f_rhoT(*rho, *T);
}

void R_iapws95_g_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_g_rhoT(*rho, *T);
}

void R_iapws95_u_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_u_rhoT(*rho, *T);
}

void R_iapws95_h_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_h_rhoT(*rho, *T);
}

void R_iapws95_h_pT(double* p, double* T, double* out)
{
  *out = iapws95_h_pT(*p, *T);
}

void R_iapws95_s_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_s_rhoT(*rho, *T);
}

void R_iapws95_p_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_p_rhoT(*rho, *T);
}

void R_iapws95_cp_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_cp_rhoT(*rho, *T);
}

void R_iapws95_cv_rhoT(double* rho, double* T, double* out)
{
  *out = iapws95_cv_rhoT(*rho, *T);
}

void R_iapws95_rho_pT(double* p, double* T, double* out)
{
  *out = iapws95_rho_pT(*p, *T);
}

void R_iapws95_rho_ph(double* p, double* h, double* out)
{
  *out = iapws95_rho_ph(*p, *h);
}

void R_iapws95_T_rhop(double* rho, double* p, double* out)
{
  *out = iapws95_T_rhop(*rho, *p);
}

void R_iapws95_T_ph(double* p, double* h, double* out)
{
  *out = iapws95_T_ph(*p, *h);
}


void R_iapws95_rhoT_ph(double* p, double* h, double* o_rho, double* o_T)
{
  iapws95_rhoT_ph(*p, *h, o_rho, o_T);
}

void R_iapws95_wv_ph(double* p, double* h, double* out)
{
  *out = iapws95_wv_ph(*p, *h);
}


void R_iapws95_sat_p_T(double* T, double* out)
{
  *out = iapws95_sat_p_T(*T);
}

void R_iapws95_sat_T_p(double* p, double* out)
{
  *out = iapws95_sat_T_p(*p);
}

void R_iapws95_sat_rhol_T(double* T, double* out)
{
  *out = iapws95_sat_rhol_T(*T);
}

void R_iapws95_sat_rhol_p(double* p, double* out)
{
  *out = iapws95_sat_rhol_p(*p);
}

void R_iapws95_sat_rhov_T(double* T, double* out)
{
  *out = iapws95_sat_rhov_T(*T);
}

void R_iapws95_sat_rhov_p(double* p, double* out)
{
  *out = iapws95_sat_rhov_p(*p);
}

void R_iapws95_sat_hl_p(double* p, double* out)
{
  *out = iapws95_sat_hl_p(*p);
}

void R_iapws95_sat_hl_T(double* T, double* out)
{
  *out = iapws95_sat_hl_T(*T);
}

void R_iapws95_sat_hv_p(double* p, double* out)
{
  *out = iapws95_sat_hv_p(*p);
}

void R_iapws95_sat_hv_T(double* T, double* out)
{
  *out = iapws95_sat_hv_T(*T);
}

void R_iapws95_sat_prho_T(double* T, double*sat_p, double* sat_rhol, double* sat_rhov, int* error)
{
  *error = iapws95_sat_prho_T(*T, sat_p, sat_rhol, sat_rhov);
}

void R_iapws95_sat_Trho_p(double* p, double*sat_T, double* sat_rhol, double* sat_rhov, int* error)
{
  *error = iapws95_sat_prho_T(*p, sat_T, sat_rhol, sat_rhov);
}
