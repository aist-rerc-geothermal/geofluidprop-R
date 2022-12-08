
#include "spl_splinek1d.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "math/spl_math.h"

SplSplinek1D* spl_splinek1d_malloc(int nsegments)
{
    SplSplinek1D* tbl = malloc(sizeof(SplSplinek1D));
    tbl->nnodes_x = nsegments+1;
    tbl->node_x = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->node_f = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->a = (double*)malloc(sizeof(double)*tbl->nnodes_x *3);
    tbl->nknots_x = tbl->nnodes_x+1;
    tbl->knot_x = (double*)malloc(sizeof(double)*tbl->nknots_x);
    tbl->inverse_tbl = NULL;

    return tbl;
}

void spl_splinek1d_free(SplSplinek1D* tbl)
{
    free(tbl->node_x);
    free(tbl->node_f);
    free(tbl->a);
    free(tbl->knot_x);
    tbl->inverse_tbl = 0;
    tbl->grid_x = 0;
    free(tbl);
}

void spl_splinek1d_write(SplSplinek1D* tbl, FILE* file)
{
    fwrite(&tbl->nnodes_x, sizeof(int), 1, file);
    fwrite(&tbl->nknots_x, sizeof(int), 1, file);
    fwrite(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    fwrite(tbl->node_f, sizeof(double), tbl->nnodes_x, file);
    fwrite(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    fwrite(tbl->a, sizeof(double), tbl->nnodes_x*3, file);
    spl_grid_write(tbl->grid_x, file);
    if (!tbl->inverse_tbl)
    {
        int dummy = 0;
        fwrite(&dummy, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek1d_write(tbl->inverse_tbl, file);
    }
}

SplSplinek1D* spl_splinek1d_read(FILE* file)
{
    int nnodes = 0;
    size_t ret;
    ret = fread(&nnodes, sizeof(int), 1, file);
    SplSplinek1D* tbl = spl_splinek1d_malloc(nnodes-1);
    tbl->nnodes_x = nnodes;
    ret = fread(&tbl->nknots_x, sizeof(int), 1, file);
    ret = fread(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    ret = fread(tbl->node_f, sizeof(double), tbl->nnodes_x, file);
    ret = fread(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    ret = fread(tbl->a, sizeof(double), tbl->nnodes_x*3, file);
    tbl->grid_x = spl_grid_read(file);
    int dummy = 0;
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl = spl_splinek1d_read(file);

    return tbl;
}


void spl_splinek1d_setInverse(SplSplinek1D* tbl, SplSplinek1D* inv_tbl)
{
    tbl->inverse_tbl = (SplSplinek1D*)inv_tbl;
}


double spl_splinek1d_spline_local(SplSplinek1D* tbl, int node_index, double x)
{
    double dx = x - tbl->node_x[node_index];
    double* a = &tbl->a[node_index*3];
    double f = a[0] + a[1]*dx + a[2]*dx*dx;
    return f;
}

double spl_splinek1d_dspline_local(SplSplinek1D* tbl, int node_index, double x)
{
    double dx = x - tbl->node_x[node_index];
    double* a = &tbl->a[node_index*3];
    double dfdx = a[1] + 2.*a[2]*dx;
    return dfdx;
}

double spl_splinek1d_d2spline_local(SplSplinek1D* tbl, int node_index, double x)
{
    double dx = x - tbl->node_x[node_index];
    double* a = &tbl->a[node_index * 3];
    double d2fdx2 = 2.*a[2];
    return d2fdx2;
}

int spl_splinek1d_getNodeIndex(SplSplinek1D* tbl, double x)
{
    return spl_grid_getNodeIndex(tbl->grid_x, x);
}

double spl_splinek1d_interpolate(SplSplinek1D* tbl, double x)
{
    int i = spl_splinek1d_getNodeIndex(tbl, x);
    double f = spl_splinek1d_spline_local(tbl, i, x);
    return f;
}

double spl_splinek1d_dx(SplSplinek1D* tbl, double x)
{
    int i = spl_splinek1d_getNodeIndex(tbl, x);
    double dfdx = spl_splinek1d_dspline_local(tbl, i, x);
    return dfdx;
}

double spl_splinek1d_tryInverse(SplSplinek1D* tbl, int i, double val)
{
    double A = tbl->a[i*3 + 2];
    double B = tbl->a[i*3 + 1];
    double C = tbl->a[i*3 + 0] - val;
    assert(A!=0.0 || B!=0.0);
    if (B==0.0)
        return (sqrt(-C/A)+ tbl->node_x[i]);
    if (fabs(A) < 1.e-10)
        return (-C/B + tbl->node_x[i]);

    if (B*B - 4 * A*C <= 0)
        return DBL_MAX;
    double xbar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];
    double xbar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];

    double df = spl_splinek1d_dspline_local(tbl, i, tbl->node_x[i]);
    double d2f = spl_splinek1d_d2spline_local(tbl, i, tbl->node_x[i]);
    if (sign(A)*df*d2f < 0)
        return xbar1;
    else
        return xbar2;

    //TODO non-monotonic function
}

double spl_splinek1d_inverse(SplSplinek1D* tbl, double val)
{
    int i = spl_splinek1d_getNodeIndex(tbl->inverse_tbl, val);

    for (int ii=0; ii<tbl->nknots_x; ii++)
    {
        double x = spl_splinek1d_tryInverse(tbl, i, val);
        if (x < tbl->knot_x[i])
            i--;
        else if (x > tbl->knot_x[i+1])
            i++;
        else {
            printf("found i after %d attempts\n", ii);
            return x;
        }
    }

    return nan("");

#if 0
    int i = 0;
    for (int ii=0; ii<tbl->nnodes_x; ii++)
    {
        double x = SBTL_tryInverse1D(tbl, i, val);
        //printf("i=%d: x=%g\n", i, x);
        if (tbl->knot_x[i] <= x && x <= tbl->knot_x[i + 1])
            return x;
        i++;
    }
    return nan(NULL);
#endif
}

