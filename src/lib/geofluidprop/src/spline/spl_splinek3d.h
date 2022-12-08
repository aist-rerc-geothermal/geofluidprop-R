
#ifndef _SPL_SPLINEK3D_H__
#define _SPL_SPLINEK3D_H__

#include "spl_grid.h"

typedef struct
{
    int nnodes_x;
    int nnodes_y;
    int nnodes_z;
    int nknots_x;
    int nknots_y;
    int nknots_z;
    double* node_x;
    double* node_y;
    double* node_z;
    double* node_f;
    double* knot_x;
    double* knot_y;
    double* knot_z;
    double* a;
    void* inverse_tbl;
    void* inverse_tbl2;
    void* inverse_tbl3;
    SplGrid* grid_x;
    SplGrid* grid_y;
    SplGrid* grid_z;
} SplSplinek3D;

SplSplinek3D* spl_splinek3d_malloc(int nsegments_x, int nsegments_y, int nsegments_z);

void spl_splinek3d_free(SplSplinek3D* tbl);

void spl_splinek3d_write(SplSplinek3D* tbl, FILE* file);

SplSplinek3D* spl_splinek3d_read(FILE* file);

double spl_splinek3d_interpolate(SplSplinek3D* tbl, double x, double y, double z);

double spl_splinek3d_spline_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_dsplinedx_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_dsplinedy_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_dsplinedz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_d2splinedxdy_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_d2splinedxdz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_d2splinedydz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

double spl_splinek3d_d3splinedxdydz_local(SplSplinek3D* tbl, int cell_index, double x, double y, double z);

void spl_splinek3d_setInverse(SplSplinek3D* tbl, SplSplinek3D* inv_tbl);

void spl_splinek3d_setInverseY(SplSplinek3D* tbl, SplSplinek3D* inv_tbl);

void spl_splinek3d_setInverseZ(SplSplinek3D* tbl, SplSplinek3D* inv_tbl);

double spl_splinek3d_inverseX(SplSplinek3D* tbl, double y, double z, double val);

double spl_splinek3d_inverseY(SplSplinek3D* tbl, double x, double z, double val);

double spl_splinek3d_inverseZ(SplSplinek3D* tbl, double x, double y, double val);

int spl_splinek3d_getNodeIndexX(SplSplinek3D* tbl, int cell_index);
int spl_splinek3d_getNodeIndexY(SplSplinek3D* tbl, int cell_index);
int spl_splinek3d_getNodeIndexZ(SplSplinek3D* tbl, int cell_index);

#endif
