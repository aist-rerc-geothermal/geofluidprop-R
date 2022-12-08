//-------------------------------------------------------------------
// This code is written based on the following document:
// * Klyukin YI, Lowell RP, Bodnar RJ (2017) A revised empirical
//   model to calculate the dynamic viscosity of H2O-NaCl fluids
//   at elevated temperatures and pressures (<=1000C, <=500 MPa,
//   0-100 wt% NaCl). Fluid Phase Equiibria 433, 193-205
//
//-------------------------------------------------------------------

#ifndef _KLYUKIN_ETAL_2017_H_
#define _KLYUKIN_ETAL_2017_H_

// x: mass fraction of NaCl
double klyukinetal2017_viscosity(double rho, double T, double xNaCl);

#endif // _KLYUKIN_ETAL_2017_H_

