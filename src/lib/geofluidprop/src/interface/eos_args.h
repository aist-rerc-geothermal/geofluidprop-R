
#ifndef _EOS_ARGS_H__
#define _EOS_ARGS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Struct for storing input parameters required for EoS calculations
 *
 */
typedef struct
{
    double p; //!< pressure [Pa]
    double T; //!< temperature [K]
    double h; //!< specific enthalpy [J/kg]
    double u; //!< specific internal energy [J/kg]
    double rho; //!< density [kg/m^3]
    double x1; //!< concentration 1
    double x2; //!< concentration 2
} EOS_ARGS;

#define EOS_TYPE_NAME_MAX_SIZE 64

typedef struct
{
    int eos_type;
    char eos_type_name[EOS_TYPE_NAME_MAX_SIZE];
    void* impl;
} EOS;

#ifdef __cplusplus
}
#endif

#endif

