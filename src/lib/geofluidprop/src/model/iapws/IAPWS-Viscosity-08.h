//-------------------------------------------------------------------
// This code is written based on the following document:
// * Release on the IAPWS Formulation 2008 for the Viscosity of
//   Ordinary Water Substance (September 2008)
//
// Range of validity of the equation are
//   0 < p < pt               and   273.16 K <= T <= 1173.15 K
//   pt <= p <= 300 MPa       and   Tm(p) <= T <= 1173.15 K
//   300 MPa < p <= 350 MPa   and   Tm(p) <= T <= 873.15 K
//   350 MPa < p <= 500 MPa   and   Tm(p) <= T <= 433.15 K
//   500 MPa < p <= 1000 MPa  and   Tm(p) <= T <= 373.15 K
//   (Tm is the pressure-dependent melting temperature and
//    pt is the triple-point pressure)
//-------------------------------------------------------------------

#ifndef _IAPWS_VISCOSITY_2008_H_
#define _IAPWS_VISCOSITY_2008_H_

double iapws08_viscosity_rhoT(double rho, double T);
double iapws08_viscosity_rhoT_simplified(double rho, double T);

#endif // _IAPWS_VISCOSITY_2008_H_

