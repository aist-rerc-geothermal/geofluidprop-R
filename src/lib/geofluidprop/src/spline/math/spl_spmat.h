
#ifndef _SPL_SPMAT_H__
#define _SPL_SPMAT_H__

typedef struct
{
    int col_index;
    double value;
} rowcol_entry;

int rowcol_entry_compare(const void* a, const void* b);

typedef struct
{
    int nrows;
    int ncols;
    int* nrowcols;
    int* nrowcols_reserved;
    rowcol_entry** rowcol_entries;
    int is_assembled;
} SparseMatrix;

SparseMatrix* spmat_alloc(int nrows, int ncols);
void spmat_free(SparseMatrix*);
double spmat_get(SparseMatrix* A, int i, int j);
void spmat_set(SparseMatrix* A, int i, int j, double v);
void spmat_add(SparseMatrix* A, int i, int j, double v);
void spmat_assemble(SparseMatrix* A);
void spmat_reset(SparseMatrix* A);
void spmat_print(SparseMatrix* A);
void spmat_matvec(SparseMatrix* A, double* x, double* y);
void spmat_matTvec(SparseMatrix* A, double* x, double* y);
void spmat_matmat(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C);
void spmat_matTmat(SparseMatrix* A, SparseMatrix* B, SparseMatrix* C);

void spmat_get_crs(SparseMatrix* A, int* nnz, int** row_ptr, int** col_indices, double** val);
void spmat_get_dense(SparseMatrix* A, double** val);

#endif
