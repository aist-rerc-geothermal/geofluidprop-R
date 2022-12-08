
#ifndef _SPL_LSOLVER_H__
#define _SPL_LSOLVER_H__

#include "spl_vec.h"

enum LINEAR_SOLVER_TYPE
{
    LINEAR_SOLVER_TYPE_GAUSS = 0
};

void set_linear_solver_type(int solver_type);
void solve_linear(Matrix* A, Vector* b, Vector* x);

#endif
