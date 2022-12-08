
#include "spl_splinek2d.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "math/spl_math.h"


SplSplinek2D* spl_splinek2d_malloc(int nsegments_x, int nsegments_y)
{
    SplSplinek2D* tbl = malloc(sizeof(SplSplinek2D));
    tbl->nnodes_x = nsegments_x + 1;
    tbl->nnodes_y = nsegments_y + 1;
    tbl->node_x = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->node_y = (double*)malloc(sizeof(double)*tbl->nnodes_y);

    tbl->a = (double*)malloc(sizeof(double)*tbl->nnodes_x*tbl->nnodes_y*9);

    tbl->nknots_x = tbl->nnodes_x + 1;
    tbl->nknots_y = tbl->nnodes_y + 1;
    tbl->knot_x = (double*)malloc(sizeof(double)*tbl->nknots_x);
    tbl->knot_y = (double*)malloc(sizeof(double)*tbl->nknots_y);
    tbl->inverse_tbl = NULL;
    tbl->inverse_tbl2 = NULL;
    tbl->algorithm = SPL_SPLINEK2D_EX1;

    return tbl;
}


void spl_splinek2d_free(SplSplinek2D* tbl)
{
    free(tbl->node_x);
    free(tbl->node_y);
    free(tbl->a);
    free(tbl->knot_x);
    free(tbl->knot_y);
    tbl->inverse_tbl = 0;
    tbl->inverse_tbl2 = 0;
    tbl->grid_x = 0;
    tbl->grid_y = 0;
    free(tbl);
}

void spl_splinek2d_write(SplSplinek2D* tbl, FILE* file)
{
    fwrite(&tbl->nnodes_x, sizeof(int), 1, file);
    fwrite(&tbl->nnodes_y, sizeof(int), 1, file);
    fwrite(&tbl->nknots_x, sizeof(int), 1, file);
    fwrite(&tbl->nknots_y, sizeof(int), 1, file);
    fwrite(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    fwrite(tbl->node_y, sizeof(double), tbl->nnodes_y, file);
    fwrite(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    fwrite(tbl->knot_y, sizeof(double), tbl->nknots_y, file);
    fwrite(tbl->a, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*9, file);
    spl_grid_write(tbl->grid_x, file);
    spl_grid_write(tbl->grid_y, file);
    if (!tbl->inverse_tbl)
    {
        fwrite(&tbl->inverse_tbl, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek2d_write(tbl->inverse_tbl, file);
    }
    if (!tbl->inverse_tbl2)
    {
        fwrite(&tbl->inverse_tbl2, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek2d_write(tbl->inverse_tbl2, file);
    }
    fwrite(&tbl->algorithm, sizeof(int), 1, file);

}

SplSplinek2D* spl_splinek2d_read(FILE* file)
{
    int nnodesx = 0, nnodesy = 0;
    size_t ret;
    ret = fread(&nnodesx, sizeof(int), 1, file);
    ret = fread(&nnodesy, sizeof(int), 1, file);
    SplSplinek2D* tbl = spl_splinek2d_malloc(nnodesx-1, nnodesy-1);
    tbl->nnodes_x = nnodesx;
    tbl->nnodes_y = nnodesy;
    ret = fread(&tbl->nknots_x, sizeof(int), 1, file);
    ret = fread(&tbl->nknots_y, sizeof(int), 1, file);
    ret = fread(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    ret = fread(tbl->node_y, sizeof(double), tbl->nnodes_y, file);
    ret = fread(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    ret = fread(tbl->knot_y, sizeof(double), tbl->nknots_y, file);
    ret = fread(tbl->a, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*9, file);
    tbl->grid_x = spl_grid_read(file);
    tbl->grid_y = spl_grid_read(file);
    int dummy = 0;
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl = spl_splinek2d_read(file);
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl2 = spl_splinek2d_read(file);
    ret = fread(&tbl->algorithm, sizeof(int), 1, file);

    return tbl;
}


void spl_splinek2d_setInverse(SplSplinek2D* tbl, SplSplinek2D* inv_tbl)
{
    tbl->inverse_tbl = (SplSplinek2D*)inv_tbl;
}

void spl_splinek2d_setInverseY(SplSplinek2D* tbl, SplSplinek2D* inv_tbl)
{
    tbl->inverse_tbl2 = (SplSplinek2D*)inv_tbl;
}

static int spl_splinek2d_getNodeIndexX(SplSplinek2D* tbl, int cell_index)
{
    return cell_index % tbl->nnodes_x;
}

static int spl_splinek2d_getNodeIndexY(SplSplinek2D* tbl, int cell_index)
{
    return cell_index / tbl->nnodes_x;
}

static int spl_splinek2d_getCellIndex(SplSplinek2D* tbl, double x, double y)
{
    int i = spl_grid_getNodeIndex(tbl->grid_x, x);
    int j = spl_grid_getNodeIndex(tbl->grid_y, y);
    if (i < 0 || j < 0) return -1;
    return (j*tbl->nnodes_x + i);
}

double spl_splinek2d_spline_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index*9];
    double f = a[0]
               + a[1]*dx + a[2]*dx*dx
               + a[3]*dy + a[4]*dy*dy
               + a[5]*dx*dy + a[6]*dx*dx*dy + a[7]*dx*dy*dy
               + a[8]*dx*dx*dy*dy;
    return f;
}

double spl_splinek2d_dsplinedx_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index*9];
    double dfdx = a[1] + 2.*a[2]*dx
               + a[5]*dy + 2.*a[6]*dx*dy + a[7]*dy*dy
               + 2.*a[8]*dx*dy*dy;
    return dfdx;
}

double spl_splinek2d_dsplinedy_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index*9];
    double dfdy = a[3] + 2.*a[4]*dy
               + a[5]*dx + a[6]*dx*dx + 2.*a[7]*dx*dy
               + 2.*a[8]*dx*dx*dy;
    return dfdy;
}

double spl_splinek2d_d2splinedxdy_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index * 9];
    double d2fdxdy = a[5] + 2.*a[6] * dx + 2.*a[7] * dy + 4.*a[8] * dx*dy;
    return d2fdxdy;
}

double spl_splinek2d_d2splinedx2_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index*9];
    double d2fdx2 = 2.*a[2]
               + 2.*a[6]*dy
               + 2.*a[8]*dy*dy;
    return d2fdx2;
}

double spl_splinek2d_d2splinedy2_local(SplSplinek2D* tbl, int cell_index, double x, double y)
{
    int node_index_x = spl_splinek2d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek2d_getNodeIndexY(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double* a = &tbl->a[cell_index*9];
    double d2fdy2 = 2.*a[4]
               + 2.*a[7]*dx
               + 2.*a[8]*dx*dx;
    return d2fdy2;
}

double spl_splinek2d_interpolate(SplSplinek2D* tbl, double x, double y)
{
    int ic = spl_splinek2d_getCellIndex(tbl, x, y);
    if (ic<0) return nan("");
    //printf("[spl_splinek2d_spline] ic=%d\n", ic);
    double f = spl_splinek2d_spline_local(tbl, ic, x, y);
    return f;
}

double spl_splinek2d_tryInverseX(SplSplinek2D* tbl, int ic, double y, double val)
{
    int i = spl_splinek2d_getNodeIndexX(tbl, ic);
    int j = spl_splinek2d_getNodeIndexY(tbl, ic);
    double dy = y - tbl->node_y[j];
    double A = tbl->a[ic*9 + 2] + dy*(tbl->a[ic*9 + 6] + dy*tbl->a[ic*9 + 8]);
    double B = tbl->a[ic*9 + 1] + dy*(tbl->a[ic*9 + 5] + dy*tbl->a[ic*9 + 7]);
    double C = tbl->a[ic*9 + 0] + dy*(tbl->a[ic*9 + 3] + dy*tbl->a[ic*9 + 4]) - val;
    if (fabs(A) < 1.e-10) {
        double x = -C/B + tbl->node_x[i];
        return x;
    }
    if (B*B - 4 * A*C <= 0)
        return nan("");
    double xbar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];
    double xbar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];

    double dfdx = spl_splinek2d_dsplinedx_local(tbl, ic, tbl->node_x[i], tbl->node_y[j]);
    double d2fdx2 = spl_splinek2d_d2splinedx2_local(tbl, ic, tbl->node_x[i], tbl->node_y[j]);
    if (sign(A)*dfdx*d2fdx2 < 0)
        return xbar1;
    else
        return xbar2;

    //TODO non-monotonic function
}

double spl_splinek2d_tryInverseY(SplSplinek2D* tbl, int ic, double x, double val)
{
    int i = spl_splinek2d_getNodeIndexX(tbl, ic);
    int j = spl_splinek2d_getNodeIndexY(tbl, ic);
    double dx = x - tbl->node_x[i];
    // printf("  i=%d, j=%d\n", i, j);
    // printf("  tbl->node_x[i]=%g\n", tbl->node_x[i]);
    // printf("  tbl->node_y[j]=%g\n", tbl->node_y[j]);
    double A = tbl->a[ic*9 + 4] + dx*(tbl->a[ic*9 + 7] + dx*tbl->a[ic*9 + 8]);
    double B = tbl->a[ic*9 + 3] + dx*(tbl->a[ic*9 + 5] + dx*tbl->a[ic*9 + 6]);
    double C = tbl->a[ic*9 + 0] + dx*(tbl->a[ic*9 + 1] + dx*tbl->a[ic*9 + 2]) - val;
    if (fabs(A)<1e-9) {
        double y = -C/B + tbl->node_y[j];
        return y;
    }
    if (B*B - 4 * A*C <= 0)
        return nan("");
    double ybar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_y[j];
    double ybar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_y[j];

    double dfdy = spl_splinek2d_dsplinedy_local(tbl, ic, tbl->node_x[i], tbl->node_y[j]);
    double d2fdy2 = spl_splinek2d_d2splinedy2_local(tbl, ic, tbl->node_x[i], tbl->node_y[j]);
    if (sign(A)*dfdy*d2fdy2 < 0)
        return ybar1;
    else
        return ybar2;

    //TODO non-monotonic function
}

int spl_splinek2d_get_adjacent_cellids_y(SplSplinek2D* tbl, int ic, double y, int* search_ic, int n)
{
    int is = 0;
    const int jy = spl_grid_getNodeIndex(tbl->grid_y, y);
    if (jy < 0) return -1;
    int j = spl_splinek2d_getNodeIndexY(tbl, ic);
    if (j == jy)
    {
        search_ic[is++] = ic;
        int i = spl_splinek2d_getNodeIndexX(tbl, ic);
        if (i > 0)
            search_ic[is++] = ic - 1;
        if (i < tbl->nnodes_x - 1)
            search_ic[is++] = ic + 1;
    }
    return is;
}

int spl_splinek2d_get_adjacent_cellids_x(SplSplinek2D* tbl, int ic, double x, int* search_ic, int n)
{
    int is = 0;
    int i = spl_splinek2d_getNodeIndexX(tbl, ic);
    const int ix = spl_grid_getNodeIndex(tbl->grid_x, x);
    if (ix < 0) return -1;
    if (i == ix)
    {
        search_ic[is++] = ic;
        int j = spl_splinek2d_getNodeIndexY(tbl, ic);
        if (j > 0)
            search_ic[is++] = ic - tbl->nnodes_x;
        if (j < tbl->nnodes_y - 1)
            search_ic[is++] = ic + tbl->nnodes_x;
    }
    return is;
}

double spl_splinek2d_inverseX(SplSplinek2D* tbl, double y, double val)
{
    //1st guess
    double x0 = 0.0;
    if (tbl->inverse_tbl)
    {
        x0 = spl_splinek2d_interpolate(tbl->inverse_tbl, y, val);
    } else {
        x0 = 0.5*(tbl->grid_x->intervals[0] + tbl->grid_x->intervals[1]);
    }
    return spl_splinek2d_inverseX2(tbl, y, val, x0);
}

double spl_splinek2d_inverseX2(SplSplinek2D* tbl, double y, double val, double x0)
{
    //printf("x0 = %g\n", x0);
    int ic0 = spl_splinek2d_getCellIndex(tbl, x0, y);
    const int nc = tbl->nnodes_x*tbl->nnodes_y;

    int search_ic[3];
    int nis = spl_splinek2d_get_adjacent_cellids_y(tbl, ic0, y, search_ic, 3);

    for (int ii=0; ii<nis; ii++)
    {
        int ic = search_ic[ii];
        //printf("%d: ic=%d\n", ii, ic);
        double x = spl_splinek2d_tryInverseX(tbl, ic, y, val);
        //assert(!isnan(x));
        if (!isnan(x))
        {
            int i = spl_splinek2d_getNodeIndexX(tbl, ic);
            if (tbl->knot_x[i] <= x && x <= tbl->knot_x[i+1])
            {
                // printf("\tfound ic after %d attempts\n", ii);
                return x;
            }
        }
    }

    const int jy = spl_grid_getNodeIndex(tbl->grid_y, y);
    if (jy < 0) return nan("");
    for (int i=0; i<tbl->nnodes_x; i++)
    {
        int ic = jy * tbl->nnodes_x + i;
        double x = spl_splinek2d_tryInverseX(tbl, ic, y, val);
        //assert(!isnan(x));
        if (!isnan(x))
        {
            if (tbl->knot_x[i] <= x && x <= tbl->knot_x[i+1])
            {
                // printf("\tfound ic after %d attempts\n", ic);
                return x;
            }
        }
    }

    return nan("");
}

double spl_splinek2d_inverseY2(SplSplinek2D* tbl, double x, double val, double y0)
{
    // printf("x=%g, val=%g, y0 = %g\n", x, val, y0);
    int ic0 = spl_splinek2d_getCellIndex(tbl, x, y0);
    const int nc = tbl->nnodes_x*tbl->nnodes_y;
    int search_ic[3];
    int nis = spl_splinek2d_get_adjacent_cellids_x(tbl, ic0, x, search_ic, 3);

    // printf("ic0 = %d\n", ic0);
    // printf("nc = %d\n", nc);
    // printf("nis = %d\n", nis);

    for (int ii=0; ii<nis; ii++)
    {
        int ic = search_ic[ii];
        // printf("%d: ic=%d\n", ii, ic);
        double y = spl_splinek2d_tryInverseY(tbl, ic, x, val);
        //assert(!isnan(y));
        if (!isnan(y))
        {
            int j = spl_splinek2d_getNodeIndexY(tbl, ic);
            if (tbl->knot_y[j] <= y && y <= tbl->knot_y[j+1])
            {
                // printf("\tfound ic after %d attempts\n", ii);
                return y;
            }
        }
    }

    const int ix = spl_grid_getNodeIndex(tbl->grid_x, x);
    if (ix < 0) return nan("");
    // printf("ix=%d\n", ix);
    for (int j=0; j<tbl->nnodes_y; j++)
    {
        int ic = j * tbl->nnodes_x + ix;
        // printf("j=%d: ic=%d\n", j, ic);
        double y = spl_splinek2d_tryInverseY(tbl, ic, x, val);
        if (!isnan(y))
        {
            if (tbl->knot_y[j] <= y && y <= tbl->knot_y[j+1])
            {
                // printf("\tfound ic after %d attempts\n", ic);
                return y;
            }
        }
    }

    return nan("");
}

int ij_to_cell(SplSplinek2D* tbl, int i, int j)
{
    return (j*tbl->nnodes_x + i);
}

double spl_splinek2d_inverseY(SplSplinek2D* tbl, double x, double val)
{
    // printf("x=%g, val=%g\n", x, val);
    // i = 305;
    // j = 166;
    // ic = ij_to_cell(tbl, i,j);
    // printf("ic=%d\n", ic);
    // printf("node_x[i-1]=%g, nnode_x[i]=%g, node_x[i+1]=%g\n", tbl->node_x[i-1], tbl->node_x[i], tbl->node_x[i+1]);
    // printf("node_y[j-1]=%g, node_y[j]=%g, node_y[j+1]=%g\n", tbl->node_y[j-1], tbl->node_y[j], tbl->node_y[j+1]);
    // printf("knot_x[i]=%g, knot_x[i+1]=%g\n", tbl->knot_x[i], tbl->knot_x[i+1]);
    // printf("knot_y[j]=%g, knot_y[j+1]=%g\n", tbl->knot_y[j], tbl->knot_y[j+1]);
    // printf("T(i,j)=%g\n", spl_splinek2d_spline_local(tbl, ic, tbl->knot_x[i], tbl->knot_y[j]));
    // printf("T(i,j+1)=%g\n", spl_splinek2d_spline_local(tbl, ic, tbl->knot_x[i], tbl->knot_y[j+1]));
    // printf("T(i+1,j)=%g\n", spl_splinek2d_spline_local(tbl, ic, tbl->knot_x[i+1], tbl->knot_y[j]));
    // printf("T(i+1,j+1)=%g\n", spl_splinek2d_spline_local(tbl, ic, tbl->knot_x[i+1], tbl->knot_y[j+1]));
    // printf("y=%g\n", spl_splinek2d_tryInverseY(tbl, ic, x, val));
    // return -1;

    //1st guess
    double y0 = 0.0;
    if (tbl->inverse_tbl2)
    {
        y0 = spl_splinek2d_interpolate(tbl->inverse_tbl2, x, val);
    } else {
        // int i = spl_grid_getNodeIndex(tbl->grid_x, x);
        y0 = 0.5*(tbl->grid_y->intervals[0] + tbl->grid_y->intervals[1]);
        // printf("tbl->grid_y->intervals[0,1] = %g, %g\n", tbl->grid_y->intervals[0], tbl->grid_y->intervals[1]);
    }
    return spl_splinek2d_inverseY2(tbl, x, val, y0);
}
