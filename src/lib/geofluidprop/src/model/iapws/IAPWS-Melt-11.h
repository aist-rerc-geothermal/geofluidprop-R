//-------------------------------------------------------------------
// This code is written based on the following document:
// * Revised Release on the Pressure along the Melting and
//   Sublimation Curves of Ordinary Water Substance (2011)
//
// Range of validity of the equations are
//   273.16K <= T <= 647.096K
//-------------------------------------------------------------------

#ifndef _IAPWS_MELT_2011_H_
#define _IAPWS_MELT_2011_H_

//double iapws11_melt_presssure(double T);

// 251.165 ~ 273.16
double iapws11_melt_presssure_ice_Ih(double T);

// 251.165 ~ 256.164
double iapws11_melt_presssure_ice_III(double T);

// 256.164~ 273.31
double iapws11_melt_presssure_ice_V(double T);

// 273.31~355
double iapws11_melt_presssure_ice_VI(double T);

// 355~715
double iapws11_melt_presssure_ice_VII(double T);

// 50 ~ 273.16
double iapws11_sublimation_pressure(double T);

double iapws11_melt_temperature(double p);

#endif // _IAPWS_MELT_2011_H_

