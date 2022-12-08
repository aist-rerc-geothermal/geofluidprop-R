
#include "spl_vec.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "spl_math.h"
#include "spl_darray.h"


Vector* vec_alloc_zero(int n)
{
    Vector* x = malloc(sizeof(Vector));
    x->n = n;
    x->v = darray_malloc_zero(n);
    return x;
}

void vec_set(Vector* x, int i, double v)
{
    assert(i < x->n);
    x->v[i] = v;
}

void vec_add(Vector* x, int i, double v)
{
    assert(i < x->n);
    x->v[i] += v;
}

void vec_set_all(Vector* x, double val)
{
    darray_set(x->v, x->n, val);
}

double vec_get(Vector* x, int i)
{
    assert(i < x->n);
    return x->v[i];
}

double vec_norm2(Vector* x)
{
    double norm = darray_norm2(x->v, x->n);
    return norm;
}

double vec_normi(Vector* x)
{
    double norm = darray_normi(x->v, x->n);
    return norm;
}

void vec_print(Vector* x)
{
    darray_print(x->v, x->n);
}

void vec_free(Vector* x)
{
    free(x->v);
    free(x);
}

Matrix* mat_alloc_zero(int nrows, int ncols)
{
    Matrix* mat = malloc(sizeof(Matrix));
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->A = spmat_alloc(nrows, ncols);
    return mat;
}

double mat_get(Matrix* A, int i, int j)
{
    return spmat_get(A->A, i, j);
}

void mat_set(Matrix* A, int i, int j, double v)
{
    spmat_set(A->A, i, j, v);
}

void mat_assemble(Matrix* A)
{
    spmat_assemble(A->A);
}

void mat_print(Matrix* A)
{
    spmat_print(A->A);
}

void mat_reset(Matrix* A)
{
    spmat_reset(A->A);
}

void mat_free(Matrix* A)
{
    spmat_free(A->A);
    free(A);
}

void mat_matvec(Matrix* A, Vector* x, Vector* y)
{
    vec_set_all(y, 0);
    spmat_matvec(A->A, x->v, y->v);
}


void mat_matTvec(Matrix* A, Vector* x, Vector* y)
{
    vec_set_all(y, 0);
    spmat_matTvec(A->A, x->v, y->v);
}

void mat_matmat(Matrix* A, Matrix* B, Matrix* C)
{
    spmat_matmat(A->A, B->A, C->A);
}

void mat_matTmat(Matrix* A, Matrix* B, Matrix* C)
{
    spmat_matTmat(A->A, B->A, C->A);
}
