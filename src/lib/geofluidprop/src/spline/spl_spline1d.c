
#include "spl_spline1d.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "math/spl_darray.h"
#include "math/spl_math.h"


double spl_spline1d_f_local(SplSpline1D* tbl, int node_index, double x)
{
    double dx = x - tbl->node_x[node_index];
    double* a = &tbl->a[node_index*3];
    double f = a[0] + a[1]*dx + a[2]*dx*dx;
    return f;
}

double spl_spline1d_df_local(SplSpline1D* tbl, int node_index, double x)
{
    double dx = x - tbl->node_x[node_index];
    double* a = &tbl->a[node_index*3];
    double dfdx = a[1] + 2.*a[2]*dx;
    return dfdx;
}

double spl_spline1d_d2f_local(SplSpline1D* tbl, int node_index, double x)
{
    double* a = &tbl->a[node_index * 3];
    double d2fdx2 = 2.*a[2];
    return d2fdx2;
}

int spl_spline1d_getNodeIndex(SplSpline1D* tbl, double x)
{
    int ibegin = 0;
    int iend = tbl->nnodes_x;
    int i = 0;
    while (1)
    {
        i = (iend + ibegin) / 2;
        if (x < tbl->node_x[i])
            iend = i;
        else if (x > tbl->node_x[i+1])
            ibegin = i;
        else
            break;
    }
    return i;
}

double spl_spline1d_interpolate(SplSpline1D* tbl, double x)
{
    int i = spl_spline1d_getNodeIndex(tbl, x);
    double f = spl_spline1d_f_local(tbl, i, x);
    return f;
}


void spl_spline1d_setSplineCoefficients(SplSpline1D* tbl)
{
    const int nseg = tbl->nnodes_x - 1;
    const int na = nseg * 3;
    darray_set(tbl->a, na, 1.);

    double* J = darray_malloc_zero(na*na);
    double* du = darray_malloc_zero(na);
    double* r = darray_malloc_zero(na);

    for (int inl=0; inl<3; inl++)
    {
        //-----------------------------------------------------------
        // calc r
        //-----------------------------------------------------------
        // match exact at the left and right edges of the segment
        for (int i = 0; i < nseg; i++)
        {
            r[i * 3] = spl_spline1d_f_local(tbl, i, tbl->node_x[i]) - tbl->node_f[i];
            r[i * 3 + 1] = spl_spline1d_f_local(tbl, i, tbl->node_x[i+1]) - tbl->node_f[i+1];
        }
        // continuity of 1st deriv at the shared nodes
        for (int i = 0; i < nseg-1; i++)
        {
            r[i * 3 + 2] = spl_spline1d_df_local(tbl, i, tbl->node_x[i+1])
                            - spl_spline1d_df_local(tbl, i+1, tbl->node_x[i+1]);
        }
        // f''(xmin) = 0
        //r[na-1] = spl_spline1d_d2f_local(tbl, 0, tbl->node_x[0]);
        r[na - 1] = spl_spline1d_d2f_local(tbl, nseg-1, tbl->node_x[tbl->nnodes_x-1]);

        //printf("r=\n");
        //vec_print(r, na);

        // check
        double norm_r = darray_norm2(r,na);
        printf("%d: |r|=%g\n", inl, norm_r);
        if (norm_r < 1.e-12)
            break;

        //-----------------------------------------------------------
        // eval J=dr/da
        //-----------------------------------------------------------
        // match exact at the left and right edges of the segment
        for (int i = 0; i < nseg; i++)
        {
            // fi(x0) = f(x0)
            J[(i * 3) * na + i * 3] = 1.;
            // fi(x1) = f(x1)
            double dx = tbl->node_x[i+1] - tbl->node_x[i];
            J[(i * 3 + 1) * na + i * 3] = 1;
            J[(i * 3 + 1) * na + i * 3 + 1] = dx;
            J[(i * 3 + 1) * na + i * 3 + 2] = dx*dx;
        }

        // continuity of 1st deriv at the shared nodes
        for (int i = 0; i < nseg-1; i++)
        {
            double dx = tbl->node_x[i+1] - tbl->node_x[i];
            // r = a12 + 2*a13*dx - a22
            J[(i * 3 + 2) * na + i * 3 + 1] = 1.;
            J[(i * 3 + 2) * na + i * 3 + 2] = 2*dx;
            J[(i * 3 + 2) * na + (i + 1) * 3 + 1] = -1.;
        }
        // f''(xmin) = 0
        J[(na-1)*na + na-1] = 2.;
        //J[(na - 1)*na + 2] = 2.;

        //printf("J=\n");
        //mat_print(J, na, na);

        //-----------------------------------------------------------
        // update solution
        //-----------------------------------------------------------
        darray_gauss(na, J, r, du);

        for (int i=0; i<na; i++)
            tbl->a[i] -= du[i];
    }

    free(J);
    free(du);
    free(r);
}

SplSpline1D* spl_spline1d_malloc(int n)
{
    SplSpline1D* tbl = (SplSpline1D*)malloc(sizeof(SplSpline1D));
    tbl->nnodes_x = n;
    tbl->node_x = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->node_f = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->a = (double*)malloc(sizeof(double)*(tbl->nnodes_x-1) * 3);
    return tbl;
}

SplSpline1D* spl_spline1d_create(double* x, double* y, int n)
{
    SplSpline1D* tbl = spl_spline1d_malloc(n);
    darray_copy(x, tbl->node_x, n);
    darray_copy(y, tbl->node_f, n);

    spl_spline1d_setSplineCoefficients(tbl);
    return tbl;
}

void spl_spline1d_free(SplSpline1D* tbl)
{
    free(tbl->node_x);
    free(tbl->node_f);
    free(tbl->a);
    free(tbl);
}
