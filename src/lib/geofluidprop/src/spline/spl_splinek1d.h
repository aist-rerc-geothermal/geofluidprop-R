
#ifndef _SPL_SPLINEK1D_H__
#define _SPL_SPLINEK1D_H__

#include "spl_grid.h"

typedef struct
{
    int nnodes_x;
    int nknots_x;
    double* node_x;
    double* node_f;
    double* knot_x;
    double* a;
    void* inverse_tbl;
    SplGrid* grid_x;
} SplSplinek1D;

SplSplinek1D* spl_splinek1d_malloc(int nsegments);

void spl_splinek1d_free(SplSplinek1D* tbl);

SplSplinek1D* spl_splinek1d_load(char const* filename);

double spl_splinek1d_interpolate(SplSplinek1D* tbl, double x);

double spl_splinek1d_dx(SplSplinek1D* tbl, double x);

void spl_splinek1d_setInverse(SplSplinek1D* tbl, SplSplinek1D* inv_tbl);

double spl_splinek1d_inverse(SplSplinek1D* tbl, double val);

void spl_splinek1d_write(SplSplinek1D* tbl, FILE* file);

SplSplinek1D* spl_splinek1d_read(FILE* file);

#endif
