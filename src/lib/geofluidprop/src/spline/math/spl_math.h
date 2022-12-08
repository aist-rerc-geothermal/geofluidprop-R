
#ifndef _SPL_MATH_H__
#define _SPL_MATH_H__

#if !defined(MAX) && !defined(max)
#define MAX(a,b) (((a)>(b))?(a):(b))
#define max MAX
#endif

#if !defined(MIN) && !defined(min)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define min MIN
#endif

#if !defined(SIGN) && !defined(sign)
#define SIGN(a) ((0 < (a)) - ((a) < 0))
#define sign SIGN
#endif

#endif

