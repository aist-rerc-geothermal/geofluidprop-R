
#include "spl_darray.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spl_math.h"


double* darray_malloc(int n)
{
    return (double*)malloc(sizeof(double)*n);
}

void darray_set(double* x, int n, double v)
{
    for (int i=0; i<n; i++)
        x[i] = v;
}

void darray_copy(double* src, double* dest, int n)
{
    for (int i=0; i<n; i++)
        dest[i] = src[i];
}

double* darray_malloc_zero(int n)
{
    double* v = darray_malloc(n);
    if (v == NULL) {
        printf("ERROR : out-of-memory error in %s\n", __FUNCTION__);
        return NULL;
    }
    darray_set(v, n, 0.);
    return v;
}

double darray_norm2(double* x, int n)
{
    double v = 0;
    for (int i=0; i<n; i++)
        v += x[i]*x[i];
    v = sqrt(v);
    return v;
}

double darray_normi(double* x, int n)
{
    double vmax = 0;
    for (int i=0; i<n; i++)
        vmax = max(vmax, fabs(x[i]));
    return vmax;
}

void darray_print(double* x, int n)
{
    for (int i = 0; i < n; i++)
       printf("%g\n", x[i]);
    printf("\n");
}

void darray2d_print(double* A, int nrows, int ncols)
{
#if 0
    for (int i = 0; i < nrows; i++)
    {
       for (int j = 0; j < ncols; j++)
       {
           if (A[i*ncols + j]!=0.0)
                printf("%d, %d: %g\n", i, j, A[i*ncols + j]);
       }
    }
#else
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
            printf("%g ", A[i*ncols + j]);
        printf("\n");
    }
#endif
}


void darray_gauss(int n, double* A, double*b, double*x)
{
    for (int i=0; i<n; i++)
    {
        // Search for maximum in this column
        double maxEl = fabs(A[i*n+i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++)
        {
            if (fabs(A[k*n+i]) > maxEl)
            {
                maxEl = fabs(A[k*n+i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        if (i!=maxRow)
        {
            for (int k=i; k<n;k++)
            {
                double tmp = A[maxRow*n+k];
                A[maxRow*n+k] = A[i*n+k];
                A[i*n+k] = tmp;
            }
            double tmp = b[maxRow];
            b[maxRow] = b[i];
            b[i] = tmp;
        }

        assert(A[i * n + i] != 0.0);

        //printf("i:%d\n", i);
        //darray2d_print(A, n, n);

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++)
        {
            double c = - A[k*n+i] / A[i*n+i];
            for (int j=i+1; j<n; j++)
                A[k*n+j] += c * A[i*n+j];
            A[k*n+i] = 0;
            b[k] += c*b[i];
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--)
    {
        double tmp = b[i];
        for (int k=n-1; k>i; k--)
            tmp -= A[i*n+k] * x[k];
        x[i] = tmp / A[i*n+i];
    }
}
