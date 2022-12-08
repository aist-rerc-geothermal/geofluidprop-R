
#include "gauss3.h"

#include <math.h>
#include <stdio.h>

#undef Real
#undef ABS
#undef gauss3
#undef direct2
#undef det3
#undef cramer3

#if defined(REAL_AS_DOUBLE)
#define Real double
#define ABS(x) fabs(x)
#define gauss3 gauss3
#define direct2 direct2
#define det3 det3
#define cramer3 cramer3

#elif defined(REAL_AS_LONG_DOUBLE)
#define Real long double
#define ABS(x) fabsl(x)
#define gauss3 gauss3_long
#define direct2 direct2_long
#define det3 det3_long
#define cramer3 cramer3_long

#elif defined(REAL_AS_QUAD)
#include <quadmath.h>
#define Real __float128
#define ABS(x) fabsq(x)
#define gauss3 gauss3_quad
#define direct2 direct2_quad
#define det3 det3_quad
#define cramer3 cramer3_quad

#endif


void gauss3(Real* A, Real*b, Real*x)
{
#define n 3
    for (int i=0; i<n; i++)
    {
        // Search for maximum in this column
        Real maxEl = ABS(A[i*n+i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++)
        {
            if (ABS(A[k*n+i]) > maxEl)
            {
                maxEl = ABS(A[k*n+i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        if (i!=maxRow)
        {
            for (int k=i; k<n;k++)
            {
                Real tmp = A[maxRow*n+k];
                A[maxRow*n+k] = A[i*n+k];
                A[i*n+k] = tmp;
            }
            Real tmp = b[maxRow];
            b[maxRow] = b[i];
            b[i] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++)
        {
            Real c = - A[k*n+i] / A[i*n+i];
            for (int j=i+1; j<n; j++)
                A[k*n+j] += c * A[i*n+j];
            A[k*n+i] = 0;
            b[k] += c*b[i];
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--)
    {
        Real tmp = b[i];
        for (int k=n-1; k>i; k--)
            tmp -= A[i*n+k] * x[k];
        x[i] = tmp / A[i*n+i];
    }
#undef n
}

void direct2(Real* A, Real*b, Real*x)
{
    x[0] = 1./(A[1*2+0]-A[0*2+0]*A[1*2+1]/A[0*2+1])*(b[1]-A[1*2+1]/A[0*2+1]*b[0]);
    x[1] = 1./A[0*2+1]*(b[0] - A[0*2+0]*x[0]);
}

Real det3(Real a11, Real a12, Real a13, Real a21, Real a22, Real a23, Real a31, Real a32, Real a33)
{
    Real det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21* a32;
    det -= (a13 * a22 *a31 + a11 * a23 * a32 + a12 * a21 * a33);
    return det;
}

int cramer3(Real* A, Real*b, Real*x)
{
#define AA(i,j) A[i*3+j]
    Real D = det3(AA(0,0),AA(0,1),AA(0,2),AA(1,0),AA(1,1),AA(1,2),AA(2,0),AA(2,1),AA(2,2));
    if (D==0.) {
        printf("error in %s: det is zero\n",  __FUNCTION__);
        return 0;
    }
    Real Dx = det3(b[0], AA(0,1),AA(0,2),b[1],AA(1,1),AA(1,2),b[2],AA(2,1),AA(2,2));
    Real Dy = det3(AA(0,0),b[0],AA(0,2),AA(1,0),b[1],AA(1,2),AA(2,0),b[2],AA(2,2));
    Real Dz = det3(AA(0,0),AA(0,1),b[0],AA(1,0),AA(1,1),b[1],AA(2,0),AA(2,1),b[2]);
    x[0] = Dx/D;
    x[1] = Dy/D;
    x[2] = Dz/D;
    return 1;
#undef AA
}

