

#include "eos_type_impl.h"

#include "eos-linear.h"
#include "eos-iapws95.h"
#include "eos-if97.h"
#include "eos-iaps84.h"
#include "eos-driesner07.h"
#ifdef USE_FREESTEAM
#include "eos-if97-freesteam.h"
#endif
#ifdef USE_PROST
#include "eos-iaps84-prost.h"
#include "eos-driesner07-prost.h"
#endif
#include "eos-iapws95sbtl-ph.h"

void eos_impl_init_register()
{
    for (int i=0; i<EOS_IMPL_ARRAY_MAX_SIZE; i++)
        eos_impl_register[i] = NULL;

    eos_impl_register[EOS_TYPE_LINEAR] = eos_linear_register;
    eos_impl_register[EOS_TYPE_WATER_IAPWS95] = eos_iapws95_register;
    eos_impl_register[EOS_TYPE_WATER_IF97] = eos_if97_register;
    eos_impl_register[EOS_TYPE_WATER_IAPS84] = eos_iaps84_register;
    eos_impl_register[EOS_TYPE_H2ONaCl_DRIESNER07] = eos_driesner07_register;
#ifdef USE_FREESTEAM
    eos_impl_register[EOS_TYPE_WATER_IF97_FREESTEAM] = eos_if97_freesteam_register;
#endif
#ifdef USE_PROST
    eos_impl_register[EOS_TYPE_WATER_IAPS84_PROST] = eos_iaps84_prost_register;
    eos_impl_register[EOS_TYPE_H2ONaCl_DRIESNER07_PROST] = eos_driesner07_prost_register;
#endif
    eos_impl_register[EOS_TYPE_WATER_IAPWS95_SBTL] = eos_iapws95sbtl_ph_register;
}

