
#ifndef _SPL_SPLINEK2D_H__
#define _SPL_SPLINEK2D_H__

#include "spl_grid.h"

enum SPL_SPLINEK2D_ALGORITHM
{
    SPL_SPLINEK2D_SPAETH = 0,
    SPL_SPLINEK2D_EX1 = 1
};

typedef struct
{
    int nnodes_x;
    int nnodes_y;
    int nknots_x;
    int nknots_y;
    double* node_x;
    double* node_y;
    double* knot_x;
    double* knot_y;
    double* a;
    void* inverse_tbl;
    void* inverse_tbl2;
    SplGrid* grid_x;
    SplGrid* grid_y;
    int algorithm;
} SplSplinek2D;

SplSplinek2D* spl_splinek2d_malloc(int nsegments_x, int nsegments_y);

void spl_splinek2d_free(SplSplinek2D* tbl);

void spl_splinek2d_write(SplSplinek2D* tbl, FILE* file);

SplSplinek2D* spl_splinek2d_read(FILE* file);

double spl_splinek2d_interpolate(SplSplinek2D* tbl, double x, double y);

double spl_splinek2d_spline_local(SplSplinek2D* tbl, int cell_index, double x, double y);

double spl_splinek2d_dsplinedx_local(SplSplinek2D* tbl, int cell_index, double x, double y);

double spl_splinek2d_dsplinedy_local(SplSplinek2D* tbl, int cell_index, double x, double y);

double spl_splinek2d_d2splinedxdy_local(SplSplinek2D* tbl, int cell_index, double x, double y);

void spl_splinek2d_setInverse(SplSplinek2D* tbl, SplSplinek2D* inv_tbl);

void spl_splinek2d_setInverseY(SplSplinek2D* tbl, SplSplinek2D* inv_tbl);

double spl_splinek2d_inverseX(SplSplinek2D* tbl, double y, double val);

double spl_splinek2d_inverseX2(SplSplinek2D* tbl, double y, double val, double x0);

double spl_splinek2d_inverseY(SplSplinek2D* tbl, double x, double val);

double spl_splinek2d_inverseY2(SplSplinek2D* tbl, double x, double val, double y0);

#endif
