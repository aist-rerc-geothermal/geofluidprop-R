
#include "spl_lsolver.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "spl_math.h"
#include "spl_darray.h"




void mat_solve_linear_dense(Matrix* A, Vector* b, Vector* x)
{
    double* dense_values;
    spmat_get_dense(A->A, &dense_values);
    darray_gauss(A->nrows, dense_values, b->v, x->v);
    free(dense_values);
}


static int mat_linear_solver_type = LINEAR_SOLVER_TYPE_GAUSS;

void set_linear_solver_type(int solver_type)
{
    mat_linear_solver_type = solver_type;
}

void solve_linear(Matrix* A, Vector* b, Vector* x)
{
    if (mat_linear_solver_type == LINEAR_SOLVER_TYPE_GAUSS)
        mat_solve_linear_dense(A, b, x);
    else
    {
        printf("***error: the linear solver type (%d) is not supported\n", mat_linear_solver_type);
    }

}

