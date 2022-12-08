
#ifndef _UTILITY_H__
#define _UTILITY_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef USE_QUAD
#include <quadmath.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double to_Pa(double bar) {return bar*1e5;}
inline double to_bar(double Pa) {return Pa*1e-5;}
inline double to_K(double C) {return C+273.15;}
inline double to_C(double K) {return K-273.15;}

double calc_vapor_quality(double sat_yl, double sat_yv, double y);

double calc_vapor_saturation(double vg, double vl, double wv);


#if !defined(MAX) && !defined(max)
#define MAX(a,b) (((a)>(b))?(a):(b))
#define max MAX
#endif
#if !defined(MAX) && defined(max)
#define MAX max
#endif

#if !defined(MIN) && !defined(min)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define min MIN
#endif
#if !defined(MIN) && defined(min)
#define MIN min
#endif

#endif //_UTILITY_H__

