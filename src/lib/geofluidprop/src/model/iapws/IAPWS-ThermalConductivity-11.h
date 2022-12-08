//-------------------------------------------------------------------
// This code is written based on the following document:
// * Release on the IAPWS Formulation 2011 for the Thermal
//   Conductivity of Ordinary Water Substance (September 2011)
//
// Range of validity of the equation are
//   0 < p < pt               and   273.16 K <= T <= 1173.15 Ks
//   pt <= p <= 100 MPa       and   Tm(p) <= T <= 1173.15 K
//   100 MPa < p <= 250 MPa   and   Tm(p) <= T <= 874 K
//   250 MPa < p <= 687 MPa   and   Tm(p) <= T <= 573 K
//   687 MPa < p <= 785 MPa   and   Tm(p) <= T <= 403 K
//   785 MPa < p <= 1000 MPa  and   Tm(p) <= T <= 348 K
//   (Tm is the pressure-dependent melting temperature and
//    pt is the triple-point pressure 611.657 Pa)
//-------------------------------------------------------------------

#ifndef _IAPWS_THERMAL_CONDUCTIVITY_2011_H_
#define _IAPWS_THERMAL_CONDUCTIVITY_2011_H_

double iapws11_thermal_conductivity_rhoT(double rho, double T);

double iapws11_lambda0bar(double Tbar);
double iapws11_lambda1bar(double Tbar, double rhobar);
double iapws11_lambda2bar(double Tbar, double rhobar);

#endif // _IAPWS_THERMAL_CONDUCTIVITY_2011_H_
