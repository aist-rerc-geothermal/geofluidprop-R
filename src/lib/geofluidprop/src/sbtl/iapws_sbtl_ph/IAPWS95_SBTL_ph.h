#ifndef _MODEL_IAPWS95SBTL_PH_H_
#define _MODEL_IAPWS95SBTL_PH_H_

#include "sbtl/sbtl_utility.h"

typedef struct
{
    SbtlTable2d* tbl_rho_ph;
    SbtlTable2d* tbl_T_ph;
    SbtlTable2d* tbl_vis_ph;
    SbtlTable1d* tbl_sat_T_p;
    SbtlTable1d* tbl_sat_hl_p;
    SbtlTable1d* tbl_sat_hv_p;
} IAPWS95SBTL_PH;

void* iapws95sbtl_ph_create();

void iapws95sbtl_ph_free(IAPWS95SBTL_PH* sbtl);

double iapws95sbtl_ph_rho_pT(IAPWS95SBTL_PH* sbtl, double p, double h);

double iapws95sbtl_ph_rho_ph(IAPWS95SBTL_PH* sbtl, double p, double h);

double iapws95sbtl_ph_h_pT(IAPWS95SBTL_PH* sbtl, double p, double T);

double iapws95sbtl_ph_T_ph(IAPWS95SBTL_PH* sbtl, double p, double h);

double iapws95sbtl_ph_vis_ph(IAPWS95SBTL_PH* sbtl, double p, double h);

double iapws95sbtl_ph_sat_p_T(IAPWS95SBTL_PH* sbtl, double p);

double iapws95sbtl_ph_sat_T_p(IAPWS95SBTL_PH* sbtl, double T);

double iapws95sbtl_ph_sat_rhol_T(IAPWS95SBTL_PH* sbtl, double T);

double iapws95sbtl_ph_sat_rhov_T(IAPWS95SBTL_PH* sbtl, double T);

void iapws95sbtl_ph_sat_rho_p(IAPWS95SBTL_PH* sbtl, double p, double* rhol, double* rhov);

double iapws95sbtl_ph_sat_hl_p(IAPWS95SBTL_PH* sbtl, double p);

double iapws95sbtl_ph_sat_hv_p(IAPWS95SBTL_PH* sbtl, double p);

void iapws95sbtl_ph_sat_h_p(IAPWS95SBTL_PH* sbtl, double p, double* hl, double* hv);

#endif // _MODEL_IAPWS95SBTL_PH_H_

