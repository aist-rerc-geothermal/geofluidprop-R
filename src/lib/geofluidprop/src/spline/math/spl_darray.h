
#ifndef _SPL_DARRAY_H__
#define _SPL_DARRAY_H__

double* darray_malloc(int n);
double* darray_malloc_zero(int n);
void darray_set(double* x, int n, double v);
void darray_copy(double* src, double* dest, int n);
double darray_norm2(double* x, int n);
double darray_normi(double* x, int n);
void darray_print(double* x, int n);
void darray2d_print(double* A, int nrows, int ncols);

void darray_gauss(int n, double* A, double*b, double*x);

#endif
