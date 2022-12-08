
#ifndef _SPL_spl_spline1d_H__
#define _SPL_spl_spline1d_H__

typedef struct
{
    int nnodes_x;
    double* node_x;
    double* node_f;
    double* a;
} SplSpline1D;

SplSpline1D* spl_spline1d_create(double* x, double* y, int n);

double spl_spline1d_interpolate(SplSpline1D* obj, double x);

void spl_spline1d_free(SplSpline1D* obj);

#endif
