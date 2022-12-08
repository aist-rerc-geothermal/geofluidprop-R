
#include "IAPWS-IF97.h"

#include <math.h>
#include <stdio.h>

//---------------------------------------------------------
// Thermodynamic constants used in IF97
//---------------------------------------------------------
#include "IAPWS-IF97-const.h"

//---------------------------------------------------------
// Internal (global) variables  <- to be removed
//---------------------------------------------------------
static double Pref;
static int is_if97_initialized = 0;

//---------------------------------------------------------
// Overall
//---------------------------------------------------------
void if97_init()
{
    if (is_if97_initialized) return;
    Pref = if97_sat_p_T(TKref);
    if97_init_region1();
    if97_init_region2();
    if97_init_region3();
    is_if97_initialized = 1;
}

double if97_pc() { return Pcrt; }
double if97_Tc() { return TKcrt; }
double if97_rhoc() { return Rhocrt; }

double if97_rho_pT(double pres, double tmpK)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres || tmpK < TKmin || TKmax < tmpK)
    {
        //printf("ERROR: out of range p=%g, T=%g\n", pres, tmpK);
        return -1;
    }

    double density = 0.0;
    if (pres <= Pref)
    {
        if (tmpK <= if97_sat_T_p(pres))
            density = if97_rho_pT_region1(pres, tmpK);
        else
            density = if97_rho_pT_region2(pres, tmpK);
    }
    else
    {
        if (tmpK <= TKref)
            density = if97_rho_pT_region1(pres, tmpK);
        else if (tmpK < if97_boundary23_T_p(pres))
            density = if97_rho_pT_region3(pres, tmpK);
        else
            density = if97_rho_pT_region2(pres, tmpK);
    }

    return density;
}

double if97_rho_ph(double p, double h)
{
    double wv = if97_vapor_quality_ph(p, h);
    if (wv < 0) return -1;
    if (wv==0.0 || wv==1.0)
    {
        double T = if97_T_ph(p, h);
        return if97_rho_pT(p, T);
    }

    double rhov, rhol;
    if97_sat_rho_p(p, &rhol, &rhov);

    double v = wv/rhov + (1.-wv)/rhol;
    return 1./v;
}

double if97_u_pT(double pMPa, double TK)
{
    // u = h + p/rho
    double h = if97_h_pT(pMPa, TK); // [kJ/kg]
    double rho = if97_rho_pT(pMPa, TK);
    double u = h + pMPa / rho * 1e-3;
    return u;
}

double if97_h_pT(double pres, double tmpK)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres || tmpK < TKmin || TKmax < tmpK)
    {
        //printf("ERROR: out of range p=%g, T=%g\n", pres, tmpK);
        return -1;
    }

    double enthalpy = 0.0;
    if (pres <= Pref)
    {
        if (tmpK <= if97_sat_T_p(pres))
            enthalpy = if97_h_pT_region1(pres, tmpK);
        else
            enthalpy = if97_h_pT_region2(pres, tmpK);
    }
    else
    {
        if (tmpK <= TKref)
            enthalpy = if97_h_pT_region1(pres, tmpK);
        else if (tmpK < if97_boundary23_T_p(pres))
            enthalpy = if97_h_pT_region3(pres, tmpK);
        else
            enthalpy = if97_h_pT_region2(pres, tmpK);
    }

    return enthalpy;
}


void if97_sat_h_p(double pres, double *hl, double *hv)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pcrt < pres)
    {
        //printf("ERROR: out of range p=%g, T=%g\n", pres, tmpK);
        *hl = *hv = -1;
        return;
    }

    double Ts = if97_sat_T_p(pres);
    if (Ts <= TKref)
    {
        *hl = if97_h_pT_region1(pres, Ts);
        *hv = if97_h_pT_region2(pres, Ts);
        return;
    }

    double rhol, rhov;
    if97_sat_rho_p(pres, &rhol, &rhov);
    *hl = if97_h_rhoT_region3(rhol, Ts);
    *hv = if97_h_rhoT_region3(rhov, Ts);
}

double if97_cp_pT(double pres, double tmpK)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres || tmpK < TKmin || TKmax < tmpK)
    {
        //printf("ERROR: out of range p=%g, T=%g\n", pres, tmpK);
        return -1;
    }

    double cp = 0.0;
    if (pres <= Pref)
    {
        if (tmpK <= if97_sat_T_p(pres))
            cp = if97_cp_pT_region1(pres, tmpK);
        else
            cp = if97_cp_pT_region2(pres, tmpK);
    }
    else
    {
        if (tmpK <= TKref)
            cp = if97_cp_pT_region1(pres, tmpK);
        else if (tmpK < if97_boundary23_T_p(pres))
        {
             double rho = if97_rho_pT_region3(pres, tmpK);
             cp = if97_cp_rhoT_region3(rho, tmpK);
        }
        else
             cp = if97_cp_pT_region2(pres, tmpK);
    }

    return cp;
}

double if97_h_pu(double pres, double u)
{
    double h = u;
    for (int i=0; i<50; i++)
    {
        double rho = if97_rho_ph(pres, h);
        if (rho < 0) return -1;
        double drhodh = (rho - if97_rho_ph(pres, h-h*1e-3))/(h*1e-3);
        double r = (h + pres/rho) - u;
        if (fabs(r) < 1e-3) {
            return h;
        }
        double drdh = 1. - pres/(rho*rho) * drhodh;
        h += -r/drdh;
    }
    printf("ERROR: dit not converged in if97_h_pu()\n");
    return -1;
}

double if97_T_pu(double pres, double u)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres)
    {
        //printf("ERROR: out of range p=%g, h=%g\n", pres, enth);
        return -1;
    }

    double h = if97_h_pu(pres, u);
    if (h<0) return -1;
    return if97_T_ph(pres, h);
}


double if97_T_ph(double pres, double enth)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres)
    {
        //printf("ERROR: out of range p=%g, h=%g\n", pres, enth);
        return -1;
    }

    int region = if97_find_region_ph(pres, enth);
    if (region == 1)
        return if97_T_ph_region1_jsme(pres, enth);
    if (region == 2)
        return if97_T_ph_region2_jsme(pres, enth);
    if (region == 30)
        return if97_T_ph_region30(pres, enth);
    if (region == 31)
        return if97_T_ph_region31(pres, enth);
    if (region == 32)
        return if97_T_ph_region32(pres, enth);

    return if97_sat_T_p(pres);
}


double if97_vapor_quality_ph(double pres, double enth)
{
    if (!is_if97_initialized) if97_init();
    const int regionID = if97_find_region_ph(pres, enth);
    if (regionID == 1 || regionID == 30 || regionID == 31)
        return 0.0;
    else if (regionID == 2 || regionID == 32)
        return 1.0;

    if (regionID != 12 && regionID != 33)
    {
        printf("ERROR: Out of IAPWS-IF97 applicable range! p=%g, h=%g, regionID=%d (in %s())\n", pres, enth, regionID, __FUNCTION__);
        return -1;
    }

    double TmpK = if97_sat_T_p(pres);
    double hL, hG;
    if (regionID == 12)
    {
        hL = if97_h_pT_region1(pres, TmpK);
        hG = if97_h_pT_region2(pres, TmpK);
    }
    else if (regionID == 33)
    {
        double dL = 1.0 / if97_sat_vL_T_region3(TmpK);
        double dG = 1.0 / if97_sat_vG_T_region3(TmpK);
        hL = if97_h_rhoT_region3(dL, TmpK);
        hG = if97_h_rhoT_region3(dG, TmpK);
    }
    double Dryn = (enth - hL) / (hG - hL);
    return Dryn;
}

void if97_sat_rho_p(double pres, double* rhol, double* rhov)
{
    if (!is_if97_initialized) if97_init();
    if (pres > Pcrt) {
        *rhol = *rhov = 0.0;
        return;
    }

    double TmpK = if97_sat_T_p(pres);
    if (pres <= Pref)
    {
        *rhol = if97_rho_pT_region1(pres, TmpK);
        *rhov = if97_rho_pT_region2(pres, TmpK);
    } else if (pres <= Pcrt) {
        *rhol = 1.0 / if97_sat_vL_T_region3(TmpK);
        *rhov = 1.0 / if97_sat_vG_T_region3(TmpK);
    }
}

void if97_sat_rho_T(double TmpK, double* rhol, double* rhov)
{
    if (!is_if97_initialized) if97_init();
    if (TmpK > TKcrt) {
        *rhol = *rhov = 0.0;
        return;
    }

    if (TmpK <= TKref)
    {
        double pres = if97_sat_p_T(TmpK);
        *rhol = if97_rho_pT_region1(pres, TmpK);
        *rhov = if97_rho_pT_region2(pres, TmpK);
    } else if (TmpK <= TKcrt) {
        *rhol = 1.0 / if97_sat_vL_T_region3(TmpK);
        *rhov = 1.0 / if97_sat_vG_T_region3(TmpK);
    }
}


int if97_find_region_ph(double pres, double enth)
{
    if (!is_if97_initialized) if97_init();
    if (pres < Pmin || Pmax < pres) return -1;
    if (enth > if97_h_pT_region2(pres, TKmax))
        return -1; // h > hmax
    if (enth < if97_h_pT_region1(pres, TKmin))
        return -1; // h < hin

    // Test blow saturation pressure at 350 degC
    if (pres <= Pref)
    {
        double satT = if97_sat_T_p(pres);
        if (enth >= if97_h_pT_region2(pres, satT))
            return 2;

        if (enth > if97_h_pT_region1(pres, satT))
            return 12;

        return 1;
    }
    // Test below critical pressure
    else if (pres <= Pcrt)
    {
        double bnd23_T = if97_boundary23_T_p(pres);
        if (enth >= if97_h_pT_region2(pres, bnd23_T))
            return 2;

        double sat_T = if97_sat_T_p(pres);
        double sat_rhoG = 1.0 / if97_sat_vG_T_region3(sat_T);
        if (enth >= if97_h_rhoT_region3(sat_rhoG, sat_T))
            return 32; // gas

        double sat_rhoL = 1.0 / if97_sat_vL_T_region3(sat_T);
        if (enth > if97_h_rhoT_region3(sat_rhoL, sat_T))
            return 33; // two phase

        if (enth > if97_h_pT_region1(pres, TKref))
            return 31; // liquid

        return 1;
    }
    // Test above critical pressure
    else
    {
        double bnd23_T = if97_boundary23_T_p(pres);
        if (enth >= if97_h_pT_region2(pres, bnd23_T))
            return 2;

        if (enth > if97_h_pT_region1(pres, TKref))
            return 30;

        return 1;
    }
}
