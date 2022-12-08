//-------------------------------------------------------------------
// This code is written based on the following document:
// * Revised Release on the IAPWS Formulation 1995 for the
//   Thermodynamic Properties of Ordinary Water Substance for
//   General and Scientific Use (2016)
//-------------------------------------------------------------------

#ifndef _IAPWS95_H_
#define _IAPWS95_H_

// #ifdef REAL_AS_DOUBLE
// #define REAL_AS_DOUBLE_PRESET
// #endif
// #ifdef REAL_AS_LONG_DOUBLE
// #define REAL_AS_LONG_PRESET
// #endif
// #ifdef REAL_AS_QUAD
// #define REAL_AS_QUAD_PRESET
// #endif
// #undef REAL_AS_DOUBLE
// #undef REAL_AS_LONG_DOUBLE
// #undef REAL_AS_QUAD

#define REAL_AS_DOUBLE_H
#include "IAPWS-95-template.h"
#undef REAL_AS_DOUBLE_H

#ifdef USE_LONGDOUBLE
#define REAL_AS_LONG_DOUBLE_H
#include "IAPWS-95-template.h"
#undef REAL_AS_LONG_DOUBLE_H
#endif

#ifdef USE_QUAD
#define REAL_AS_QUAD_H
#include "IAPWS-95-template.h"
#undef REAL_AS_QUAD_H
#endif

// #ifdef REAL_AS_DOUBLE_PRESET
// #define REAL_AS_DOUBLE
// #endif
// #ifdef REAL_AS_LONG_PRESET
// #define REAL_AS_LONG_DOUBLE
// #endif
// #ifdef REAL_AS_QUAD_PRESET
// #define REAL_AS_QUAD
// #endif

// #undef REAL_AS_DOUBLE_PRESET
// #undef REAL_AS_LONG_PRESET
// #undef REAL_AS_QUAD_PRESET

#endif

