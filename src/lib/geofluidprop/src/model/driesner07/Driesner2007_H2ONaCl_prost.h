
#ifndef _DRIESNER2007_H2ONACL_PROST_H_
#define _DRIESNER2007_H2ONACL_PROST_H_

#include <stdbool.h>

#include "Driesner2007_H2ONaCl.h"

// H2O-NaCl
// * Driesner & Heinrich (2007) The system H2O–NaCl. Part I:
//   Correlation formulae for phase relations in temperature–pressure
//   –composition space from 0 to 1000 °C, 0 to 5000 bar, and 0 to 1 XNaCl.
//   Geochimica et Cosmochimica Acta 71(20), 4880-4901
// * Driesner (2007) The system H2O–NaCl. Part II: Correlations for
//   molar volume, enthalpy, and isobaric heat capacity from 0 to
//   1000 °C, 1 to 5000 bar, and 0 to 1 XNaCl.
//   Geochimica et Cosmochimica Acta 71(20), 4902-4919

// typedef enum
// {
//     H2O_NaCl_PhaseType_V = 0,
//     H2O_NaCl_PhaseType_VH = 1,
//     H2O_NaCl_PhaseType_L = 10,
//     H2O_NaCl_PhaseType_LH = 11,
//     H2O_NaCl_PhaseType_VL = 20,
//     H2O_NaCl_PhaseType_VLH = 21,
//     H2O_NaCl_PhaseType_SC = 30,
//     H2O_NaCl_PhaseType_SCH = 31,
//     H2O_NaCl_PhaseType_INVALID = -1
// } H2O_NaCl_PhaseType;

H2O_NaCl_PhaseType driesner07_prost_H2O_NaCl_phase_type(double p_Pa, double T_K, double x);

char* driesner07_prost_H2O_NaCl_phase_type_to_str(H2O_NaCl_PhaseType phase, char*);

//-------------------------------------------------------------------
// critical curve: the crest of the VL surface
//-------------------------------------------------------------------
double driesner07_prost_H2O_NaCl_pc_T(double T);
double driesner07_prost_H2O_NaCl_pc_T2(double T, bool use_satp);
double driesner07_prost_H2O_NaCl_xc_T(double T);
double driesner07_prost_H2O_NaCl_Tc_x(double x);

//-------------------------------------------------------------------
// vapor-liquid
//-------------------------------------------------------------------
// mole fraction of liquid phase in VL
double driesner07_prost_H2O_NaCl_VL_xl_pT(double p, double T);
// mole fraction of liquid phase in VL
double driesner07_prost_H2O_NaCl_VL_xv_pT(double p, double T);
// vapor pressure at the given temperatture and salinty
// Remark: the returned value can be different depending on the given initial guess
double driesner07_prost_H2O_NaCl_VL_p_Txv(double T_K, double x, double p0);
// find all vapor pressure at the given temperatture and salinty
// The function returns the number of vapor pressure
int driesner07_prost_H2O_NaCl_VL_p_all_Txv(double T_K, double x, double pmin, double pmax, double* p);
// liquid pressure at the given temperature and salinty
double driesner07_prost_H2O_NaCl_VL_p_Txl(double T_K, double x, double p0);

// density
double driesner07_prost_H2O_NaCl_VL_rhov_pT(double p_Pa, double T_K);
double driesner07_prost_H2O_NaCl_VL_rhol_pT(double p_Pa, double T_K);

//-------------------------------------------------------------------
// liquid-halite
//-------------------------------------------------------------------
// mole fraction of liquid phase in LH
double driesner07_prost_H2O_NaCl_LH_xl_pT(double p, double T);

//-------------------------------------------------------------------
// vapor-halite
//-------------------------------------------------------------------
// mole fraction of liquid phase in VH
double driesner07_prost_H2O_NaCl_VH_xv_pT(double p, double T);

//-------------------------------------------------------------------
// vapor-liquid-halite
//-------------------------------------------------------------------
// V-L-H pressure
double driesner07_prost_H2O_NaCl_VLH_p_T(double T);
// mole fraction of liquid at V-L-H
double driesner07_prost_H2O_NaCl_VLH_xl_T(double T);
// mole fraction of vapor at V-L-H
double driesner07_prost_H2O_NaCl_VLH_xv_T(double T);

//-------------------------------------------------------------------
// NaCl solubility in water
//-------------------------------------------------------------------
double driesner07_prost_H2O_NaCl_NaCl_solubility_xv_pT(double p, double T);
double driesner07_prost_H2O_NaCl_NaCl_solubility_xl_pT(double p, double T);
double driesner07_prost_H2O_NaCl_NaCl_solubility_x_pT(double p, double T);

//-------------------------------------------------------------------
// single phase property
//-------------------------------------------------------------------
double driesner07_prost_H2O_NaCl_rho_singlephase_pTx(double p_Pa, double T_K, double x);
// molar volume
double driesner07_prost_H2O_NaCl_vm_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_h_singlephase_pTx(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_cp_singlephase_pTx(double p_Pa, double T_K, double x);

//-------------------------------------------------------------------
// bulk properties
//-------------------------------------------------------------------
double driesner07_prost_H2O_NaCl_mass_fraction_liquid(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_mass_fraction_vapor(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_volume_fraction_liquid(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_volume_fraction_vapor(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_volume_fraction_halite(double p_Pa, double T_K, double x);

double driesner07_prost_H2O_NaCl_T_rhopx(double rho, double p_Pa, double x, double TK0);
double driesner07_prost_H2O_NaCl_T_phx(double p_Pa, double h, double x, double TK0);
double driesner07_prost_H2O_NaCl_volume_fraction_liquid_rhopx(double rho, double p_Pa, double x, double TK0);

double driesner07_prost_H2O_NaCl_sat_p_upper_Tx(double T_K, double x, double p0);
double driesner07_prost_H2O_NaCl_sat_p_lower_Tx(double T_K, double x, double p0);
int driesner07_prost_H2O_NaCl_sat_p_Tx(double T_K, double x, double p0, double* p);

// bulk density
double driesner07_prost_H2O_NaCl_rho_pTx(double p_Pa, double T_K, double x);
// bulk specific enthalpy
double driesner07_prost_H2O_NaCl_h_pTx(double p_Pa, double T_K, double x);

//-------------------------------------------------------------------
// internal functions
//-------------------------------------------------------------------
// vapor-liquid distribution coefficient
double driesner07_prost_H2O_NaCl_VL_K(double p, double T);
double driesner07_prost_H2O_NaCl_VL_log10_Kn(double pn, double T_C);
double driesner07_prost_H2O_NaCl_VL_log10_Km(double p_bar, double T_K);
// scaled temperature used for the molar volume
double driesner07_prost_H2O_NaCl_Tv_pTx(double p /*bar*/, double T /*C*/, double x);
double driesner07_prost_H2O_NaCl_vm_pTx0(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_vm_pTx1(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_vm_pTx2(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_NaCl_rhol_pTx2(double p_Pa, double T_K, double x);
double driesner07_prost_H2O_vm_pT(double p, double T);


#endif // _DRIESNER2007_H2ONACL_H_

