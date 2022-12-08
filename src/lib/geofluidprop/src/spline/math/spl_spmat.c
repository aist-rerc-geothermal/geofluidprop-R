
#include "spl_spmat.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spl_math.h"
#include "spl_darray.h"

#define NROWCOLS_SIZE 10

int rowcol_entry_compare(const void* aa, const void* bb)
{
    rowcol_entry* a = (rowcol_entry*)aa;
    rowcol_entry* b = (rowcol_entry*)bb;
    if (a->col_index < b->col_index) return -1;
    if (a->col_index > b->col_index) return 1;
    return 0;
}

int rowcol_entry_find_pos(int n, rowcol_entry* array, int index)
{
    //if (n == 0) return -1;
    //if (n == 1) return (array[0].col_index == index) ? 0 : -1;

    //if (index < array[0].col_index || array[n - 1].col_index < index) return -1;


    rowcol_entry* result = bsearch(&index, array, n, sizeof(rowcol_entry), rowcol_entry_compare);
    if (result == NULL)
        return -1;
    return (int)(result - &array[0]);
}

SparseMatrix* spmat_alloc(int nrows, int ncols)
{
    SparseMatrix* mat = malloc(sizeof(SparseMatrix));
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->nrowcols = malloc(sizeof(int) * nrows);
    mat->nrowcols_reserved = malloc(sizeof(int) * nrows);
    mat->rowcol_entries = malloc(sizeof(rowcol_entry*) * nrows);
    int nrowcols_reserved = ncols < NROWCOLS_SIZE ? ncols : NROWCOLS_SIZE;
    for (int i = 0; i < nrows; i++)
    {
        mat->nrowcols[i] = 0;
        mat->nrowcols_reserved[i] = nrowcols_reserved;
        mat->rowcol_entries[i] = malloc(sizeof(rowcol_entry) * nrowcols_reserved);
    }
    mat->is_assembled = 0;
    return mat;
}

void spmat_free(SparseMatrix* A)
{
    for (int i = 0; i < A->nrows; i++)
        free(A->rowcol_entries[i]);
    free(A->rowcol_entries);
    free(A->nrowcols);
    free(A->nrowcols_reserved);
    free(A);
}

void spmat_set(SparseMatrix* A, int i, int j, double v)
{
    assert(i < A->nrows && j < A->ncols);
    A->is_assembled = 0;
    rowcol_entry* rowptr = A->rowcol_entries[i];
    // overwrite
    const int A_nrowcols = A->nrowcols[i];
    for (int jj = 0; jj < A_nrowcols; jj++)
    {
        if (rowptr[jj].col_index == j)
        {
            rowptr[jj].value = v;
            return;
        }
    }

    // add
    if (v == 0.0) return;
    int jj = A->nrowcols[i];
    if (jj == A->nrowcols_reserved[i])
    {
        int dsize = NROWCOLS_SIZE;
        A->rowcol_entries[i] = realloc(A->rowcol_entries[i], sizeof(rowcol_entry) * (A->nrowcols_reserved[i] + dsize));
        rowptr = A->rowcol_entries[i];
        assert(A->rowcol_entries[i] != NULL);
        A->nrowcols_reserved[i] += dsize;
    }
    rowptr[jj].col_index = j;
    rowptr[jj].value = v;
    A->nrowcols[i] += 1;
}


void spmat_add(SparseMatrix* A, int i, int j, double v)
{
    assert(i < A->nrows && j < A->ncols);
    A->is_assembled = 0;
    rowcol_entry* rowptr = A->rowcol_entries[i];
    // existing
    for (int jj = 0; jj < A->nrowcols[i]; jj++)
    {
        if (rowptr[jj].col_index == j)
        {
            rowptr[jj].value += v;
            return;
        }
    }

    // insert
    if (v == 0.0) return;
    int jj = A->nrowcols[i];
    if (jj == A->nrowcols_reserved[i])
    {
        int dsize = NROWCOLS_SIZE;
        A->rowcol_entries[i] = realloc(A->rowcol_entries[i], sizeof(rowcol_entry) * (A->nrowcols_reserved[i] + dsize));
        rowptr = A->rowcol_entries[i];
        assert(A->rowcol_entries[i] != NULL);
        A->nrowcols_reserved[i] += dsize;
    }
    rowptr[jj].col_index = j;
    rowptr[jj].value = v;
    A->nrowcols[i] += 1;
}

void spmat_assemble(SparseMatrix* A)
{
    if (A->is_assembled == 1)
        return;
    // check zero entries
    for (int i = 0; i < A->nrows; i++)
    {
        rowcol_entry* rowptr = A->rowcol_entries[i];
        for (int jj = 0; jj < A->nrowcols[i]; jj++)
        {
            if (rowptr[jj].value != 0.)
                continue;
            rowptr[jj].col_index = INT_MAX;
        }
    }
    // sort
    for (int i = 0; i < A->nrows; i++)
    {
        qsort(A->rowcol_entries[i], A->nrowcols[i], sizeof(rowcol_entry), rowcol_entry_compare);
        //shrink col size if zero exists
        rowcol_entry* rowptr = A->rowcol_entries[i];
        int jj = 0;
        for (; jj < A->nrowcols[i]; jj++)
        {
            if (rowptr[jj].col_index == INT_MAX)
                break;
        }
        A->nrowcols[i] = jj;
    }
    A->is_assembled = 1;
}

void spmat_reset(SparseMatrix* A)
{
    for (int i = 0; i < A->nrows; i++)
        A->nrowcols[i] = 0;
    A->is_assembled = 0;
}

void spmat_print(SparseMatrix* A)
{
    for (int i = 0; i < A->nrows; i++)
        for (int jj = 0; jj < A->nrowcols[i]; jj++)
            printf("%d %d: %g\n", i, A->rowcol_entries[i][jj].col_index, A->rowcol_entries[i][jj].value);
}

double spmat_get(SparseMatrix* A, int i, int j)
{
    rowcol_entry* rowptr = A->rowcol_entries[i];
    if (A->is_assembled)
    {
        int jj = rowcol_entry_find_pos(A->nrowcols[i], rowptr, j);
        return jj < 0 ? 0. : rowptr[jj].value;
    }

    // unsorted. check one by one
    for (int jj = 0; jj < A->nrowcols[i]; jj++)
        if (rowptr[jj].col_index == j)
            return rowptr[jj].value;

    return 0.0;
}

void spmat_matvec(SparseMatrix* A, double* x, double* y)
{
    assert(A->is_assembled);
    for (int i = 0; i < A->nrows; i++)
        for (int j = 0; j < A->ncols; j++)
            y[i] += spmat_get(A, i, j) * x[j];
}

void spmat_matTvec(SparseMatrix* A, double* x, double* y)
{
    assert(A->is_assembled);
    for (int j = 0; j < A->ncols; j++)
        for (int i = 0; i < A->nrows; i++)
            y[j] += spmat_get(A, i, j) * x[i];
}

void spmat_matmat(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C)
{
    assert(A->is_assembled && B->is_assembled);
    for (int i = 0; i < A->nrows; i++)
    {
        for (int l = 0; l < B->ncols; l++)
        {
            double v = 0;
            rowcol_entry* coldata = A->rowcol_entries[i];
            for (int jj = 0; jj < A->nrowcols[i]; jj++)
            {
                int j = coldata[jj].col_index;
                double a = coldata[jj].value;
                if (a == 0.0) continue;
                double b = spmat_get(B, j, l);
                v += a * b;
            }
            spmat_add(C, i, l, v);
        }
    }
}

void spmat_matTmat(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C)
{
    assert(A->is_assembled && B->is_assembled);

    // AT in row major
    SparseMatrix* AT = spmat_alloc(A->ncols, A->nrows);
    for (int i = 0; i < A->nrows; i++)
    {
        const int nrowcols = A->nrowcols[i];
        const rowcol_entry* rowcolent = A->rowcol_entries[i];
        for (int jj = 0; jj < nrowcols; jj++)
            spmat_set(AT, rowcolent[jj].col_index, i, rowcolent[jj].value);
    }
    spmat_assemble(AT);
    // printf("AT=\n"); spmat_print(AT);

    // BT in row major
    SparseMatrix* BT = spmat_alloc(B->ncols, B->nrows);
    for (int i = 0; i < B->nrows; i++)
    {
        const int nrowcols = B->nrowcols[i];
        const rowcol_entry* rowcolent = B->rowcol_entries[i];
        for (int jj = 0; jj < nrowcols; jj++)
            spmat_set(BT, rowcolent[jj].col_index, i, rowcolent[jj].value);
    }
    spmat_assemble(BT);
    // printf("BT=\n"); spmat_print(BT);

    // C = AT*(BT)^T
    int i = 0;
    const int AT_nrows = AT->nrows;
    const int BT_nrows = BT->nrows;
    #pragma omp parallel for private(i)
    for (i = 0; i < AT_nrows; i++)
    {
        const int AT_nrowcols = AT->nrowcols[i];
        if (AT_nrowcols == 0) continue;
        const rowcol_entry* AT_coldata = AT->rowcol_entries[i];
        for (int l = 0; l < BT_nrows; l++)
        {
            const int BT_nrowcols = BT->nrowcols[l];
            if (BT_nrowcols == 0) continue;
            const rowcol_entry* BT_coldata = BT->rowcol_entries[l];
            double v = 0;
            int kk = 0;
            int k = BT_coldata[kk].col_index;
            for (int jj = 0; jj < AT_nrowcols; jj++)
            {
                const int j = AT_coldata[jj].col_index;
                while (k < j)
                {
                    kk++;
                    if (kk >= BT_nrowcols)
                        break;
                    assert(kk < BT_nrowcols);
                    k = BT_coldata[kk].col_index;
                }
                if (kk >= BT_nrowcols) continue;
                if (j == k)
                    v += AT_coldata[jj].value * BT_coldata[kk].value;
            }
            spmat_set(C, i, l, v);
        }
    }
    spmat_assemble(C);

    spmat_free(AT);
    spmat_free(BT);
}


void spmat_get_crs(SparseMatrix* A, int* nnz, int** row_ptr_, int** col_index_, double** value_)
{
    assert(A->is_assembled);
    int* row_ptr = malloc(sizeof(int) * (A->nrows + 1));
    int n_nonzero = 0;
    for (int i = 0; i < A->nrows; i++)
        n_nonzero += A->nrowcols[i];
    int* col_index = malloc(sizeof(int) * n_nonzero);
    double* value = darray_malloc_zero(n_nonzero);

    int counter_row = 0;
    int counter_ptr = 0;
    for (int i = 0; i < A->nrows; i++)
    {
        row_ptr[counter_row++] = counter_ptr;         // starting point of the row
        for (int jj = 0; jj < A->nrowcols[i]; jj++)
        {
            col_index[counter_ptr] = A->rowcol_entries[i][jj].col_index;
            value[counter_ptr] = A->rowcol_entries[i][jj].value;
            counter_ptr++;
        }
    }
    row_ptr[A->nrows] = counter_ptr;

    *nnz = n_nonzero;
    *row_ptr_ = row_ptr;
    *col_index_ = col_index;
    *value_ = value;
}

void spmat_get_dense(SparseMatrix* A, double** val)
{
    assert(A->is_assembled);
    double* value = darray_malloc_zero(A->nrows * A->ncols);
    for (int i = 0; i < A->nrows; i++)
        for (int j = 0; j < A->ncols; j++)
            value[i*A->ncols + j] = spmat_get(A, i, j);

    *val = value;
}
