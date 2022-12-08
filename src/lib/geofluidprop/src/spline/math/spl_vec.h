
#ifndef _SPL_VEC_H__
#define _SPL_VEC_H__

#include "spl_spmat.h"

typedef struct
{
    int n;
    double* v;
} Vector;

Vector* vec_alloc_zero(int n);
void vec_set(Vector* x, int i, double v);
void vec_set_all(Vector* x, double v);
void vec_add(Vector* x, int i, double v);
double vec_get(Vector* x, int i);
double vec_norm2(Vector* x);
double vec_normi(Vector* x);
void vec_print(Vector* x);
void vec_free(Vector* x);

typedef struct
{
    int nrows;
    int ncols;
    SparseMatrix* A;
} Matrix;

Matrix* mat_alloc_zero(int nrows, int ncols);
double mat_get(Matrix* A, int i, int j);
void mat_set(Matrix* A, int i, int j, double v);
void mat_assemble(Matrix* A);
void mat_print(Matrix* A);
void mat_reset(Matrix* A);
void mat_free(Matrix* A);
void mat_matvec(Matrix* A, Vector* x, Vector* y);
void mat_matTvec(Matrix* A, Vector* x, Vector* y);
void mat_matmat(Matrix* A, Matrix* B, Matrix* C);
void mat_matTmat(Matrix* A, Matrix* B, Matrix* C);

#endif
