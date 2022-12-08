
#ifndef _EOS_TYPE_H__
#define _EOS_TYPE_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * List of supported EoS types
 *
 */
enum EOS_TYPE
{
    EOS_TYPE_LINEAR = 0 //!< Linear model
    ,EOS_TYPE_WATER_IAPWS95 //!< Water, IAPWS95
    ,EOS_TYPE_WATER_IF97 //!< Water, IF97
    ,EOS_TYPE_WATER_IAPS84 //!< Water, IAPS84
    ,EOS_TYPE_H2ONaCl_DRIESNER07 //!< H2O-NaCl, Driesner and Heinrich (2007), Driesner (2007)
    ,EOS_TYPE_WATER_IF97_FREESTEAM //!< Water, IF97 (freesteam)
    ,EOS_TYPE_WATER_IAPS84_PROST //!< Water, IAPS84 (PROST)
    ,EOS_TYPE_H2ONaCl_DRIESNER07_PROST //!< H2O-NaCl, Driesner and Heinrich (2007), Driesner (2007) using PROST
    ,EOS_TYPE_WATER_IAPWS95_SBTL //!< Water, IAPWS95 (SBTL)
    ,EOS_TYPE_UNKNOWN
};

/// get EoS type index
int eos_get_eos_type(const char* name);

/**
 * get EoS type name
 *
 * @param eos_type  EoS type index
 * @param name EoS name (output)
 * @return char*
 */
char* eos_get_eos_type_name(int eos_type, char* name);

#ifdef __cplusplus
}
#endif

#endif

