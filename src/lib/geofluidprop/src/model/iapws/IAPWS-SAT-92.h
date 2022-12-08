//-------------------------------------------------------------------
// This code is written based on the following document:
// * IAPWS SR1-86(1992) : Revised Supplementary Release on Saturation
//   Properties of Ordinary Water Substance September 1992
//
// Range of validity of the equations are
//   273.16K <= T <= 647.096K
//-------------------------------------------------------------------

#ifndef _IAPWS_SAT_92_H_
#define _IAPWS_SAT_92_H_

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

#endif // _IAPWS_SAT_92_H_

