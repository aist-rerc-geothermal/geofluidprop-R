
#include "utility.h"

#include <math.h>
#include <stdio.h>
#ifdef USE_QUAD
#include <quadmath.h>
#endif

extern double to_Pa(double bar);
extern double to_bar(double Pa);
extern double to_K(double C);
extern double to_C(double K);

// extern double max(double a, double b);
// extern double min(double a, double b);
//extern int sign(double a);

#ifdef USE_QUAD
extern __float128 maxq(__float128 a, __float128 b);
extern __float128 minq(__float128 a, __float128 b);
#endif

extern long double maxl(long double a, long double b);
extern long double minl(long double a, long double b);


// Calculate mass fraction of vapor in the mixture
// y can be specific enthalpy, specific entropy, specific volume
// or specific internal energy.
double calc_vapor_quality(double sat_yl, double sat_yv, double y)
{
    double x = 1.0;
    if (sat_yl < y && y < sat_yv)
        x = (y - sat_yl) / (sat_yv - sat_yl);
    else if (y <= sat_yl)
        x = 0.0;
    return x;
}

// Calculate volumetric fraction of vapor in the mixture
double calc_vapor_saturation(double vg, double vl, double wv)
{
    double satv = wv*vg/(wv*vg + (1.-wv)*vl);
    return satv;
}
