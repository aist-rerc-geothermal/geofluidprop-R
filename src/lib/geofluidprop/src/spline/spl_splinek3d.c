
#include "spl_splinek3d.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "math/spl_math.h"

static double nr_eps = 1.e-4;

SplSplinek3D* spl_splinek3d_malloc(int nsegments_x, int nsegments_y, int nsegments_z)
{
    SplSplinek3D* tbl = malloc(sizeof(SplSplinek3D));
    tbl->nnodes_x = nsegments_x + 1;
    tbl->nnodes_y = nsegments_y + 1;
    tbl->nnodes_z = nsegments_z + 1;
    tbl->node_x = (double*)malloc(sizeof(double)*tbl->nnodes_x);
    tbl->node_y = (double*)malloc(sizeof(double)*tbl->nnodes_y);
    tbl->node_z = (double*)malloc(sizeof(double)*tbl->nnodes_z);

    tbl->node_f = (double*)malloc(sizeof(double)*tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z);
    tbl->a = (double*)malloc(sizeof(double)*tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z*27);

    tbl->nknots_x = tbl->nnodes_x + 1;
    tbl->nknots_y = tbl->nnodes_y + 1;
    tbl->nknots_z = tbl->nnodes_z + 1;
    tbl->knot_x = (double*)malloc(sizeof(double)*tbl->nknots_x);
    tbl->knot_y = (double*)malloc(sizeof(double)*tbl->nknots_y);
    tbl->knot_z = (double*)malloc(sizeof(double)*tbl->nknots_z);
    tbl->inverse_tbl = NULL;
    tbl->inverse_tbl2 = NULL;
    tbl->inverse_tbl3 = NULL;

    return tbl;
}


void spl_splinek3d_free(SplSplinek3D* tbl)
{
    free(tbl->node_x);
    free(tbl->node_y);
    free(tbl->node_z);
    free(tbl->node_f);
    free(tbl->a);
    free(tbl->knot_x);
    free(tbl->knot_y);
    free(tbl->knot_z);
    tbl->inverse_tbl = 0;
    tbl->inverse_tbl2 = 0;
    tbl->inverse_tbl3 = 0;
    tbl->grid_x = 0;
    tbl->grid_y = 0;
    tbl->grid_z = 0;
    free(tbl);
}

void spl_splinek3d_write(SplSplinek3D* tbl, FILE* file)
{
    fwrite(&tbl->nnodes_x, sizeof(int), 1, file);
    fwrite(&tbl->nnodes_y, sizeof(int), 1, file);
    fwrite(&tbl->nnodes_z, sizeof(int), 1, file);
    fwrite(&tbl->nknots_x, sizeof(int), 1, file);
    fwrite(&tbl->nknots_y, sizeof(int), 1, file);
    fwrite(&tbl->nknots_z, sizeof(int), 1, file);
    fwrite(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    fwrite(tbl->node_y, sizeof(double), tbl->nnodes_y, file);
    fwrite(tbl->node_z, sizeof(double), tbl->nnodes_z, file);
    fwrite(tbl->node_f, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z, file);
    fwrite(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    fwrite(tbl->knot_y, sizeof(double), tbl->nknots_y, file);
    fwrite(tbl->knot_z, sizeof(double), tbl->nknots_z, file);
    fwrite(tbl->a, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z*27, file);
    spl_grid_write(tbl->grid_x, file);
    spl_grid_write(tbl->grid_y, file);
    spl_grid_write(tbl->grid_z, file);
    if (!tbl->inverse_tbl)
    {
        fwrite(&tbl->inverse_tbl, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek3d_write(tbl->inverse_tbl, file);
    }
    if (!tbl->inverse_tbl2)
    {
        fwrite(&tbl->inverse_tbl2, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek3d_write(tbl->inverse_tbl2, file);
    }
    if (!tbl->inverse_tbl3)
    {
        fwrite(&tbl->inverse_tbl3, sizeof(int), 1, file);
    }
    else
    {
        int dummy = 1;
        fwrite(&dummy, sizeof(int), 1, file);
        spl_splinek3d_write(tbl->inverse_tbl3, file);
    }

}

SplSplinek3D* spl_splinek3d_read(FILE* file)
{
    int nnodesx = 0, nnodesy = 0, nnodesz = 0;
    size_t ret;
    ret = fread(&nnodesx, sizeof(int), 1, file);
    ret = fread(&nnodesy, sizeof(int), 1, file);
    ret = fread(&nnodesz, sizeof(int), 1, file);
    SplSplinek3D* tbl = spl_splinek3d_malloc(nnodesx-1, nnodesy-1, nnodesz-1);
    tbl->nnodes_x = nnodesx;
    tbl->nnodes_y = nnodesy;
    tbl->nnodes_z = nnodesz;
    ret = fread(&tbl->nknots_x, sizeof(int), 1, file);
    ret = fread(&tbl->nknots_y, sizeof(int), 1, file);
    ret = fread(&tbl->nknots_z, sizeof(int), 1, file);
    ret = fread(tbl->node_x, sizeof(double), tbl->nnodes_x, file);
    ret = fread(tbl->node_y, sizeof(double), tbl->nnodes_y, file);
    ret = fread(tbl->node_z, sizeof(double), tbl->nnodes_z, file);
    ret = fread(tbl->node_f, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z, file);
    ret = fread(tbl->knot_x, sizeof(double), tbl->nknots_x, file);
    ret = fread(tbl->knot_y, sizeof(double), tbl->nknots_y, file);
    ret = fread(tbl->knot_z, sizeof(double), tbl->nknots_z, file);
    ret = fread(tbl->a, sizeof(double), tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z*27, file);
    tbl->grid_x = spl_grid_read(file);
    tbl->grid_y = spl_grid_read(file);
    tbl->grid_z = spl_grid_read(file);
    int dummy = 0;
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl = spl_splinek3d_read(file);
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl2 = spl_splinek3d_read(file);
    ret = fread(&dummy, sizeof(int), 1, file);
    if (dummy)
        tbl->inverse_tbl3 = spl_splinek3d_read(file);

    return tbl;
}

void spl_splinek3d_setInverse(SplSplinek3D* tbl, SplSplinek3D* inv_tbl)
{
    tbl->inverse_tbl = (SplSplinek3D*)inv_tbl;
}

void spl_splinek3d_setInverseY(SplSplinek3D* tbl, SplSplinek3D* inv_tbl)
{
    tbl->inverse_tbl2 = (SplSplinek3D*)inv_tbl;
}

void spl_splinek3d_setInverseZ(SplSplinek3D* tbl, SplSplinek3D* inv_tbl)
{
    tbl->inverse_tbl3 = (SplSplinek3D*)inv_tbl;
}

int spl_splinek3d_getNodeIndexX(SplSplinek3D* tbl, int cell_index)
{
    return (cell_index % (tbl->nnodes_x * tbl->nnodes_y)) % tbl->nnodes_x;
}

int spl_splinek3d_getNodeIndexY(SplSplinek3D* tbl, int cell_index)
{
    return (cell_index % (tbl->nnodes_x * tbl->nnodes_y)) / tbl->nnodes_x;
}

int spl_splinek3d_getNodeIndexZ(SplSplinek3D* tbl, int cell_index)
{
    return cell_index / (tbl->nnodes_x * tbl->nnodes_y);
}

static int spl_splinek3d_getCellIndex(SplSplinek3D* tbl, double x, double y, double z)
{
    int i = spl_grid_getNodeIndex(tbl->grid_x, x);
    int j = spl_grid_getNodeIndex(tbl->grid_y, y);
    int k = spl_grid_getNodeIndex(tbl->grid_z, z);
    return (k*tbl->nnodes_x*tbl->nnodes_y + j*tbl->nnodes_x + i);
}

double spl_splinek3d_spline_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double f = 0;
    for (int n=0; n<3; n++)
        for (int m=0; m<3; m++)
            for (int l=0; l<3; l++)
                f += a[l+m*3+n*3*3]*pow(dx,l)*pow(dy,m)*pow(dz,n);
    return f;
}

double spl_splinek3d_dsplinedx_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double df = 0;
    for (int n=0; n<3; n++)
        for (int m=0; m<3; m++)
            for (int l=1; l<3; l++)
                df += a[l+m*3+n*3*3]*l*pow(dx,l-1)*pow(dy,m)*pow(dz,n);
    return df;
}

double spl_splinek3d_dsplinedy_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double df = 0;
    for (int n=0; n<3; n++)
        for (int m=1; m<3; m++)
            for (int l=0; l<3; l++)
                df += a[l+m*3+n*3*3]*m*pow(dx,l)*pow(dy,m-1)*pow(dz,n);
    return df;
}

double spl_splinek3d_dsplinedz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double df = 0;
    for (int n=1; n<3; n++)
        for (int m=0; m<3; m++)
            for (int l=0; l<3; l++)
                df += a[l+m*3+n*3*3]*n*pow(dx,l)*pow(dy,m)*pow(dz,n-1);
    return df;
}

double spl_splinek3d_d2splinedxdy_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double d2f = 0;
    for (int n=0; n<3; n++)
        for (int m=1; m<3; m++)
            for (int l=1; l<3; l++)
                d2f += a[l+m*3+n*3*3]*l*m*pow(dx,l-1)*pow(dy,m-1)*pow(dz,n);
    return d2f;
}

double spl_splinek3d_d2splinedxdz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double d2f = 0;
    for (int n=1; n<3; n++)
        for (int m=0; m<3; m++)
            for (int l=1; l<3; l++)
                d2f += a[l+m*3+n*3*3]*l*n*pow(dx,l-1)*pow(dy,m)*pow(dz,n-1);
    return d2f;
}

double spl_splinek3d_d2splinedydz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double d2f = 0;
    for (int n=1; n<3; n++)
        for (int m=1; m<3; m++)
            for (int l=0; l<3; l++)
                d2f += a[l+m*3+n*3*3]*m*n*pow(dx,l)*pow(dy,m-1)*pow(dz,n-1);
    return d2f;
}

double spl_splinek3d_d3splinedxdydz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
{
    int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
    int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
    int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
    double dx = x - tbl->node_x[node_index_x];
    double dy = y - tbl->node_y[node_index_y];
    double dz = z - tbl->node_z[node_index_z];
    double* a = &tbl->a[cell_index*27];
    double d3f = 0;
    for (int n=1; n<3; n++)
        for (int m=1; m<3; m++)
            for (int l=1; l<3; l++)
                d3f += a[l+m*3+n*3*3]*l*m*n*pow(dx,l-1)*pow(dy,m-1)*pow(dz,n-1);
    return d3f;
}

 double spl_splinek3d_d2splinedx2_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
 {
     int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
     int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
     int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
     double dx = x - tbl->node_x[node_index_x];
     double dy = y - tbl->node_y[node_index_y];
     double dz = z - tbl->node_z[node_index_z];
     double* a = &tbl->a[cell_index * 27];
     double d2fdx2 = 2.*a[2] + 2.*a[5]*dy + 2.*a[8]*dy*dy
                + 2.*a[11]*dz + 2*a[14]*dy*dz + 2*a[17]*dy*dy*dz
                + 2.*a[20]*dz*dz + 2*a[23]*dy*dz*dz + 2*a[26]*dy*dy*dz*dz;
     return d2fdx2;
 }

 double spl_splinek3d_d2splinedy2_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
 {
     int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
     int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
     int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
     double dx = x - tbl->node_x[node_index_x];
     double dy = y - tbl->node_y[node_index_y];
     double dz = z - tbl->node_z[node_index_z];
     double* a = &tbl->a[cell_index * 27];
     double d2fdy2 = 2.*a[6] + 2.*a[7]*dx + 2.*a[8]*dx*dx
                + 2.*a[15]*dz + 2*a[16]*dx*dz + 2*a[17]*dx*dx*dz
                + 2.*a[24]*dz*dz + 2*a[25]*dx*dz*dz + 2*a[26]*dx*dx*dz*dz;
     return d2fdy2;
 }

 double spl_splinek3d_d2splinedz2_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z)
 {
     int node_index_x = spl_splinek3d_getNodeIndexX(tbl, cell_index);
     int node_index_y = spl_splinek3d_getNodeIndexY(tbl, cell_index);
     int node_index_z = spl_splinek3d_getNodeIndexZ(tbl, cell_index);
     double dx = x - tbl->node_x[node_index_x];
     double dy = y - tbl->node_y[node_index_y];
     double dz = z - tbl->node_z[node_index_z];
     double* a = &tbl->a[cell_index * 27];
     double d2fdz2 = 2.*a[18] + 2.*a[19]*dx + 2.*a[20]*dx*dx
                + 2.*a[21]*dy + 2*a[22]*dx*dy + 2*a[23]*dx*dx*dy
                + 2.*a[24]*dy*dy + 2*a[25]*dx*dy*dy + 2*a[26]*dx*dx*dy*dy;
     return d2fdz2;
 }

double spl_splinek3d_interpolate(SplSplinek3D* tbl, double x, double y, double z)
{
    int ic = spl_splinek3d_getCellIndex(tbl, x, y, z);
    //printf("[spl_splinek3d_spline] ic=%d\n", ic);
    double f = spl_splinek3d_spline_local(tbl, ic, x, y, z);
    return f;
}

double spl_splinek3d_tryInverseX(SplSplinek3D* tbl, int ic, double y, double z, double val)
{
    int i = spl_splinek3d_getNodeIndexX(tbl, ic);
    int j = spl_splinek3d_getNodeIndexY(tbl, ic);
    int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
    double dy = y - tbl->node_y[j];
    double dz = z - tbl->node_z[k];
    double A = tbl->a[ic*27 + 2] + dy*tbl->a[ic*27 + 5] + dy*dy*tbl->a[ic*27 + 8]
               + dz*tbl->a[ic*27 + 11] + dy*dz*tbl->a[ic*27 + 14] + dy*dy*dz*tbl->a[ic*27 + 17]
               + dz*dz*tbl->a[ic*27 + 20] + dy*dz*dz*tbl->a[ic*27 + 23] + dy*dy*dz*dz*tbl->a[ic*27 + 26];
    double B = tbl->a[ic*27 + 1] + dy*tbl->a[ic*27 + 4] + dy*dy*tbl->a[ic*27 + 7]
               + dz*tbl->a[ic*27 + 10] + dy*dz*tbl->a[ic*27 + 13] + dy*dy*dz*tbl->a[ic*27 + 16]
               + dz*dz*tbl->a[ic*27 + 19] + dy*dz*dz*tbl->a[ic*27 + 22] + dy*dy*dz*dz*tbl->a[ic*27 + 25];
    double C = tbl->a[ic*27 + 0] + dy*tbl->a[ic*27 + 3] + dy*dy*tbl->a[ic*27 + 6]
               + dz*tbl->a[ic*27 + 9] + dy*dz*tbl->a[ic*27 + 12] + dy*dy*dz*tbl->a[ic*27 + 15]
               + dz*dz*tbl->a[ic*27 + 18] + dy*dz*dz*tbl->a[ic*27 + 21] + dy*dy*dz*dz*tbl->a[ic*27 + 24]
                - val;
    if (fabs(A)<1.e-10) {
        double x = -C/B + tbl->node_x[i];
        return x;
    }
    if (B*B - 4 * A*C <= 0)
        return nan("");
    double xbar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];
    double xbar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_x[i];

    double dfdx = spl_splinek3d_dsplinedx_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    double d2fdx2 = spl_splinek3d_d2splinedx2_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    if (sign(A)*dfdx*d2fdx2 < 0)
        return xbar1;
    else
        return xbar2;

    //TODO non-monotonic function
}

double spl_splinek3d_tryInverseY(SplSplinek3D* tbl, int ic, double x, double z, double val)
{
    int i = spl_splinek3d_getNodeIndexX(tbl, ic);
    int j = spl_splinek3d_getNodeIndexY(tbl, ic);
    int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
    double dx = x - tbl->node_x[i];
    double dz = z - tbl->node_z[k];
    double A = tbl->a[ic*27 + 6] + dx*tbl->a[ic*27 + 7] + dx*dx*tbl->a[ic*27 + 8]
               + dz*tbl->a[ic*27 + 15] + dx*dz*tbl->a[ic*27 + 16] + dx*dx*dz*tbl->a[ic*27 + 17]
               + dz*dz*tbl->a[ic*27 + 24] + dx*dz*dz*tbl->a[ic*27 + 25] + dx*dx*dz*dz*tbl->a[ic*27 + 26];
    double B = tbl->a[ic*27 + 3] + dx*tbl->a[ic*27 + 4] + dx*dx*tbl->a[ic*27 + 5]
               + dz*tbl->a[ic*27 + 12] + dx*dz*tbl->a[ic*27 + 13] + dx*dx*dz*tbl->a[ic*27 + 14]
               + dz*dz*tbl->a[ic*27 + 21] + dx*dz*dz*tbl->a[ic*27 + 22] + dx*dx*dz*dz*tbl->a[ic*27 + 23];
    double C = tbl->a[ic*27 + 0] + dx*tbl->a[ic*27 + 1] + dx*dx*tbl->a[ic*27 + 2]
               + dz*tbl->a[ic*27 + 9] + dx*dz*tbl->a[ic*27 + 10] + dx*dx*dz*tbl->a[ic*27 + 11]
               + dz*dz*tbl->a[ic*27 + 18] + dx*dz*dz*tbl->a[ic*27 + 19] + dx*dx*dz*dz*tbl->a[ic*27 + 20]
               - val;
    if (fabs(A)<1.e-10) {
        double y = -C/B + tbl->node_y[j];
        return y;
    }
    if (B*B - 4 * A*C <= 0)
        return nan("");
    double ybar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_y[j];
    double ybar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_y[j];

    double dfdy = spl_splinek3d_dsplinedy_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    double d2fdy2 = spl_splinek3d_d2splinedy2_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    if (sign(A)*dfdy*d2fdy2 < 0)
        return ybar1;
    else
        return ybar2;

    //TODO non-monotonic function
}

double spl_splinek3d_tryInverseZ(SplSplinek3D* tbl, int ic, double x, double y, double val)
{
    int i = spl_splinek3d_getNodeIndexX(tbl, ic);
    int j = spl_splinek3d_getNodeIndexY(tbl, ic);
    int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
    double dx = x - tbl->node_x[i];
    double dy = y - tbl->node_y[j];
    double A = tbl->a[ic*27 + 18] + dx*tbl->a[ic*27 + 19] + dx*dx*tbl->a[ic*27 + 20]
               + dy*tbl->a[ic*27 + 21] + dx*dy*tbl->a[ic*27 + 22] + dx*dx*dy*tbl->a[ic*27 + 23]
               + dy*dy*tbl->a[ic*27 + 24] + dx*dy*dy*tbl->a[ic*27 + 25] + dx*dx*dy*dy*tbl->a[ic*27 + 26];
    double B = tbl->a[ic*27 + 9] + dx*tbl->a[ic*27 + 10] + dx*dx*tbl->a[ic*27 + 11]
               + dy*tbl->a[ic*27 + 12] + dx*dy*tbl->a[ic*27 + 13] + dx*dx*dy*tbl->a[ic*27 + 14]
               + dy*dy*tbl->a[ic*27 + 15] + dx*dy*dy*tbl->a[ic*27 + 16] + dx*dx*dy*dy*tbl->a[ic*27 + 17];
    double C = tbl->a[ic*27 + 0] + dx*tbl->a[ic*27 + 1] + dx*dx*tbl->a[ic*27 + 2]
               + dy*tbl->a[ic*27 + 3] + dx*dy*tbl->a[ic*27 + 4] + dx*dx*dy*tbl->a[ic*27 + 5]
               + dy*dy*tbl->a[ic*27 + 6] + dx*dy*dy*tbl->a[ic*27 + 7] + dx*dx*dy*dy*tbl->a[ic*27 + 8]
               - val;
    if (fabs(A)<1.e-10) {
        double z = -C/B + tbl->node_z[k];
        return z;
    }
    if (B*B - 4 * A*C <= 0)
        return nan("");
    double zbar1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_z[k];
    double zbar2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A) + tbl->node_z[k];

    double dfdz = spl_splinek3d_dsplinedz_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    double d2fdz2 = spl_splinek3d_d2splinedz2_local(tbl, ic, tbl->node_x[i], tbl->node_y[j], tbl->node_z[k]);
    if (sign(A)*dfdz*d2fdz2 < 0)
        return zbar1;
    else
        return zbar2;

    //TODO non-monotonic function
}

int spl_splinek3d_get_adjacent_cellids_x(SplSplinek3D* tbl, int ic, double x, int* search_ic, int n)
{
    int is = 0;
    int i = spl_splinek3d_getNodeIndexX(tbl, ic);
    const int ix = spl_grid_getNodeIndex(tbl->grid_x, x);
    if (i == ix)
    {
        search_ic[is++] = ic;
        int j = spl_splinek3d_getNodeIndexY(tbl, ic);
        if (j > 0)
            search_ic[is++] = ic - tbl->nnodes_x;
        if (j < tbl->nnodes_y - 1)
            search_ic[is++] = ic + tbl->nnodes_x;
        int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
        if (k > 0)
            search_ic[is++] = ic - tbl->nnodes_x*tbl->nnodes_y;
        if (k < tbl->nnodes_z - 1)
            search_ic[is++] = ic + tbl->nnodes_x*tbl->nnodes_y;
    }
    return is;
}

int spl_splinek3d_get_adjacent_cellids_y(SplSplinek3D* tbl, int ic, double y, int* search_ic, int n)
{
    int is = 0;
    const int jy = spl_grid_getNodeIndex(tbl->grid_y, y);
    int j = spl_splinek3d_getNodeIndexY(tbl, ic);
    if (j == jy)
    {
        search_ic[is++] = ic;
        int i = spl_splinek3d_getNodeIndexX(tbl, ic);
        if (i > 0)
            search_ic[is++] = ic - 1;
        if (i < tbl->nnodes_x - 1)
            search_ic[is++] = ic + 1;
        int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
        if (k > 0)
            search_ic[is++] = ic - tbl->nnodes_x*tbl->nnodes_y;
        if (k < tbl->nnodes_z - 1)
            search_ic[is++] = ic + tbl->nnodes_x*tbl->nnodes_y;
    }
    return is;
}

int spl_splinek3d_get_adjacent_cellids_z(SplSplinek3D* tbl, int ic, double z, int* search_ic, int n)
{
    int is = 0;
    const int kz = spl_grid_getNodeIndex(tbl->grid_z, z);
    int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
    if (k == kz)
    {
        search_ic[is++] = ic;
        int i = spl_splinek3d_getNodeIndexX(tbl, ic);
        if (i > 0)
            search_ic[is++] = ic - 1;
        if (i < tbl->nnodes_x - 1)
            search_ic[is++] = ic + 1;
        int j = spl_splinek3d_getNodeIndexY(tbl, ic);
        if (j > 0)
            search_ic[is++] = ic - tbl->nnodes_x;
        if (j < tbl->nnodes_y - 1)
            search_ic[is++] = ic + tbl->nnodes_x;
    }
    return is;
}



double spl_splinek3d_inverseX(SplSplinek3D* tbl, double y, double z, double val)
{
    //1st guess
    double x0 = spl_splinek3d_interpolate(tbl->inverse_tbl, y, z, val);
    //printf("x0 = %g\n", x0);
    int ic0 = spl_splinek3d_getCellIndex(tbl, x0, y, z);
    const int nc = tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z;

    int search_ic[6];
    int nis = spl_splinek3d_get_adjacent_cellids_y(tbl, ic0, y, search_ic, 6);

    for (int ii=0; ii<nis; ii++)
    {
        int ic = search_ic[ii];
        int ic_j = spl_splinek3d_getNodeIndexY(tbl, ic);
        int ic_k = spl_splinek3d_getNodeIndexZ(tbl, ic);
        if (y < tbl->knot_y[ic_j] || tbl->knot_y[ic_j+1] < y)
            continue;
        if (z < tbl->knot_z[ic_k] || tbl->knot_z[ic_k+1] < z)
            continue;

        //printf("%d: ic=%d\n", ii, ic);
        double x = spl_splinek3d_tryInverseX(tbl, ic, y, z, val);
        //assert(!isnan(x));
        if (!isnan(x))
        {
            int i = spl_splinek3d_getNodeIndexX(tbl, ic);
            if (tbl->knot_x[i] <= x && x <= tbl->knot_x[i+1])
            {
                printf("\tfound ic after %d attempts\n", ii);
                return x;
            }
        }
    }

    const int jy = spl_grid_getNodeIndex(tbl->grid_y, y);
    const int kz = spl_grid_getNodeIndex(tbl->grid_z, z);
    for (int i=0; i<tbl->nnodes_x; i++)
    {
        int ic = kz * tbl->nnodes_x * tbl->nnodes_y + jy * tbl->nnodes_x + i;
        double x = spl_splinek3d_tryInverseX(tbl, ic, y, z, val);
        //assert(!isnan(x));
        if (!isnan(x))
        {
            if (tbl->knot_x[i] <= x && x <= tbl->knot_x[i+1])
            {
                printf("\tfound ic after %d attempts\n", ic);
                return x;
            }
        }
    }

    return nan("");
}

double spl_splinek3d_inverseY(SplSplinek3D* tbl, double x, double z, double val)
{
    //1st guess
    double y0 = spl_splinek3d_interpolate(tbl->inverse_tbl2, x, z, val);
   // printf("y0 = %g\n", y0);
    int ic0 = spl_splinek3d_getCellIndex(tbl, x, y0, z);
    const int nc = tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z;
    int search_ic[6];
    int nis = spl_splinek3d_get_adjacent_cellids_x(tbl, ic0, x, search_ic, 6);

    for (int ii=0; ii<nis; ii++)
    {
        int ic = search_ic[ii];
        int ic_i = spl_splinek3d_getNodeIndexX(tbl, ic);
        int ic_k = spl_splinek3d_getNodeIndexZ(tbl, ic);
        if (x < tbl->knot_x[ic_i] || tbl->knot_x[ic_i+1] < x)
            continue;
        if (z < tbl->knot_z[ic_k] || tbl->knot_z[ic_k+1] < z)
            continue;

        //printf("%d: ic=%d\n", ii, ic);
        double y = spl_splinek3d_tryInverseY(tbl, ic, x, z, val);
        //assert(!isnan(y));
        if (!isnan(y))
        {
            int j = spl_splinek3d_getNodeIndexY(tbl, ic);
            if (tbl->knot_y[j] <= y && y <= tbl->knot_y[j+1])
            {
                printf("\tfound ic after %d attempts\n", ii);
                return y;
            }
        }
    }

    const int ix = spl_grid_getNodeIndex(tbl->grid_x, x);
    const int kz = spl_grid_getNodeIndex(tbl->grid_z, z);
    for (int j=0; j<tbl->nnodes_y; j++)
    {
        int ic = kz * tbl->nnodes_x * tbl->nnodes_y + j * tbl->nnodes_x + ix;
        double y = spl_splinek3d_tryInverseY(tbl, ic, x, z, val);
        if (!isnan(y))
        {
            if (tbl->knot_y[j] <= y && y <= tbl->knot_y[j+1])
            {
                printf("\tfound ic after %d attempts\n", ic);
                return y;
            }
        }
    }

    return nan("");
}

double spl_splinek3d_inverseZ(SplSplinek3D* tbl, double x, double y, double val)
{
    //1st guess
    double z0 = spl_splinek3d_interpolate(tbl->inverse_tbl3, x, y, val);
   // printf("y0 = %g\n", y0);
    int ic0 = spl_splinek3d_getCellIndex(tbl, x, y, z0);
    if (ic0 < 0) ic0 = 0;
    const int nc = tbl->nnodes_x*tbl->nnodes_y*tbl->nnodes_z;
    int search_ic[6];
    int nis = spl_splinek3d_get_adjacent_cellids_x(tbl, ic0, x, search_ic, 6);

    for (int ii=0; ii<nis; ii++)
    {
        int ic = search_ic[ii];
        int ic_i = spl_splinek3d_getNodeIndexX(tbl, ic);
        int ic_j = spl_splinek3d_getNodeIndexY(tbl, ic);
        if (x < tbl->knot_x[ic_i] || tbl->knot_x[ic_i+1] < x)
            continue;
        if (y < tbl->knot_y[ic_j] || tbl->knot_y[ic_j+1] < y)
            continue;

        //printf("%d: ic=%d\n", ii, ic);
        double z = spl_splinek3d_tryInverseZ(tbl, ic, x, y, val);
        //assert(!isnan(y));
        if (!isnan(z))
        {
            int k = spl_splinek3d_getNodeIndexZ(tbl, ic);
            if (tbl->knot_z[k] <= z && z <= tbl->knot_z[k+1])
            {
                printf("\tfound ic after %d attempts\n", ii);
                return z;
            }
        }
    }

    const int ix = spl_grid_getNodeIndex(tbl->grid_x, x);
    const int jy = spl_grid_getNodeIndex(tbl->grid_y, y);
    for (int k=0; k<tbl->nnodes_z; k++)
    {
        int ic = k * tbl->nnodes_x * tbl->nnodes_y + jy * tbl->nnodes_x + ix;
        double z = spl_splinek3d_tryInverseZ(tbl, ic, x, y, val);
        if (!isnan(z))
        {
            if (tbl->knot_y[k] <= z && z <= tbl->knot_y[k+1])
            {
                printf("\tfound ic after %d attempts\n", ic);
                return z;
            }
        }
    }

    return nan("");
}


