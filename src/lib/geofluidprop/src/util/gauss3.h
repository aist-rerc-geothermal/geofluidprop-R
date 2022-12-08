
#ifndef _GAUSS3_H__
#define _GAUSS3_H__

void gauss3(double* A, double*b, double*x);

void gauss3_long(long double* A, long double*b, long double*x);

#ifdef USE_QUAD
void gauss3_quad(__float128* A, __float128*b, __float128*x);
#endif

double det3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);
long double det3_long(long double a11, long double a12, long double a13, long double a21, long double a22, long double a23, long double a31, long double a32, long double a33);
#ifdef USE_QUAD
__float128 det3_quad(__float128 a11, __float128 a12, __float128 a13, __float128 a21, __float128 a22, __float128 a23, __float128 a31, __float128 a32, __float128 a33);
#endif

int cramer3(double* A, double*b, double*x);
int cramer3_long(long double* A, long double*b, long double*x);
#ifdef USE_QUAD
int cramer3_quad(__float128* A, __float128*b, __float128*x);
#endif

#endif
