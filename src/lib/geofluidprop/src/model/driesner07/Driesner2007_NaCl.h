#ifndef _DRIESNER2007_NACL_H_
#define _DRIESNER2007_NACL_H_


//-------------------------------------------------------------------
// Pure NaCl properties based on
// * Driesner & Heinrich (2007) The system H2O–NaCl. Part I:
//   Correlation formulae for phase relations in temperature–pressure
//   –composition space from 0 to 1000 °C, 0 to 5000 bar, and 0 to 1 XNaCl.
//   Geochimica et Cosmochimica Acta 71(20), 4880-4901
// * Driesner (2007) The system H2O–NaCl. Part II: Correlations for
//   molar volume, enthalpy, and isobaric heat capacity from 0 to
//   1000 °C, 1 to 5000 bar, and 0 to 1 XNaCl.
//   Geochimica et Cosmochimica Acta 71(20), 4902-4919
//-------------------------------------------------------------------

// triple point
double driesner07_NaCl_VLH_p();
double driesner07_NaCl_VLH_T();
double driesner07_NaCl_VLH_TC();
double driesner07_NaCl_VLH_h_halite();

// melting temperature
double driesner07_NaCl_LH_T_p(double p);
// sublimation pressure
double driesner07_NaCl_VH_p_T(double T);
// boiling pressure
double driesner07_NaCl_VL_p_T(double T);
// phase equilibrium pressure
double driesner07_NaCl_sat_p_T(double T);

// halite properties
double driesner07_NaCl_H_rho_pT(double p_Pa, double T_K);
double driesner07_NaCl_H_h_pT(double p_Pa, double T_K);
double driesner07_NaCl_H_cp_pT(double p_Pa, double T_K);

// liquid properties
double driesner07_NaCl_L_rho_pT(double p_Pa, double T_K);
double driesner07_NaCl_L_compressibility_pT(double T_K);

// vapor properties (not supported in the papers)

#endif
