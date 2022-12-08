
#include "eos_type.h"

#include <string.h>

int eos_get_eos_type(const char* name)
{
    if (strcmp("LINEAR", name) == 0)
        return EOS_TYPE_LINEAR;
    if (strcmp("IAPWS95", name) == 0)
        return EOS_TYPE_WATER_IAPWS95;
    if (strcmp("IF97", name) == 0)
        return EOS_TYPE_WATER_IF97;
    if (strcmp("IAPS84", name) == 0)
        return EOS_TYPE_WATER_IAPS84;
    if (strcmp("DRIESNER07", name) == 0)
        return EOS_TYPE_H2ONaCl_DRIESNER07;
    if (strcmp("IF97_FREESTEAM", name) == 0)
        return EOS_TYPE_WATER_IF97_FREESTEAM;
    if (strcmp("IAPS84_PROST", name) == 0)
        return EOS_TYPE_WATER_IAPS84_PROST;
    if (strcmp("DRIESNER07_PROST", name) == 0)
        return EOS_TYPE_H2ONaCl_DRIESNER07_PROST;

    return EOS_TYPE_UNKNOWN;
}


char* eos_get_eos_type_name(int eos_type, char* name)
{
    if (eos_type == EOS_TYPE_LINEAR)
        strcpy(name, "LINEAR");
    else if (eos_type == EOS_TYPE_WATER_IAPWS95)
        strcpy(name, "IAPWS95");
    else if (eos_type == EOS_TYPE_WATER_IF97)
        strcpy(name, "IF97");
    else if (eos_type == EOS_TYPE_WATER_IAPS84)
        strcpy(name, "IAPS84");
    else if (eos_type == EOS_TYPE_H2ONaCl_DRIESNER07)
        strcpy(name, "DRIESNER07");
    else if (eos_type == EOS_TYPE_WATER_IF97_FREESTEAM)
        strcpy(name, "IF97_FREESTEAM");
    else if (eos_type == EOS_TYPE_WATER_IAPS84_PROST)
        strcpy(name, "IAPS84_PROST");
    else if (eos_type == EOS_TYPE_H2ONaCl_DRIESNER07_PROST)
        strcpy(name, "DRIESNER07_PROST");
    else
        strcpy(name, "UNKNOWN");

    return name;
}

