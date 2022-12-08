
#include "spl_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

SplGrid* spl_grid_malloc()
{
    SplGrid* grid = (SplGrid*)malloc(sizeof(SplGrid));
    grid->intervals = NULL;
    grid->nnodes = NULL;
    grid->nintervals = 0;
    grid->first_node_ids = NULL;
    grid->dx = NULL;
    grid->nallnodes = 0;
    return grid;
}

void spl_grid_free(SplGrid* grid)
{
    free(grid->intervals);
    free(grid->nnodes);
    free(grid->first_node_ids);
    free(grid->dx);
    free(grid);
}

void spl_grid_write(SplGrid* grid, FILE* stream)
{
    fwrite(&grid->nintervals, sizeof(int), 1, stream);
    fwrite(grid->intervals, sizeof(double), grid->nintervals, stream);
    fwrite(grid->nnodes, sizeof(int), grid->nintervals-1, stream);
    fwrite(grid->first_node_ids, sizeof(int), grid->nintervals-1, stream);
    fwrite(grid->dx, sizeof(double), grid->nintervals-1, stream);
    fwrite(&grid->nallnodes, sizeof(int), 1, stream);
}

SplGrid* spl_grid_read(FILE* stream)
{
    SplGrid* grid = spl_grid_malloc();
    size_t ret;
    ret = fread(&grid->nintervals, sizeof(int), 1, stream);
    grid->intervals = malloc(sizeof(double)*grid->nintervals);
    grid->nnodes = malloc(sizeof(int)*(grid->nintervals-1));
    grid->first_node_ids = malloc(sizeof(int)*(grid->nintervals-1));
    grid->dx = malloc(sizeof(double)*(grid->nintervals-1));
    ret = fread(grid->intervals, sizeof(double), grid->nintervals, stream);
    ret = fread(grid->nnodes, sizeof(int), grid->nintervals-1, stream);
    ret = fread(grid->first_node_ids, sizeof(int), grid->nintervals-1, stream);
    ret = fread(grid->dx, sizeof(double), grid->nintervals-1, stream);
    ret = fread(&grid->nallnodes, sizeof(int), 1, stream);
    return grid;
}

SplGrid* spl_grid_createUniform(double vmin, double vmax, int nnodes)
{
    SplGrid* grid = spl_grid_malloc();
    grid->nintervals = 2;
    grid->intervals = malloc(sizeof(double)*grid->nintervals);
    grid->intervals[0] = vmin;
    grid->intervals[1] = vmax;
    grid->nnodes = malloc(sizeof(int)*1);
    grid->nnodes[0] = nnodes;
    grid->first_node_ids = malloc(sizeof(int)*1);
    grid->first_node_ids[0] = 0;
    grid->dx = malloc(sizeof(double)*1);
    grid->dx[0] = (vmax-vmin)/(nnodes-1);
    grid->nallnodes = nnodes;
    return grid;
}

SplGrid* spl_grid_create(int nintervals, double* intervals, int* nnodes)
{
    SplGrid* grid = spl_grid_malloc();
    grid->nintervals = nintervals;
    grid->intervals = malloc(sizeof(double)*grid->nintervals);
    for (int i=0; i<nintervals; i++)
        grid->intervals[i] = intervals[i];
    grid->nnodes = malloc(sizeof(int)*(nintervals-1));
    for (int i=0; i<nintervals-1; i++)
        grid->nnodes[i] = nnodes[i];
    grid->first_node_ids = malloc(sizeof(int)*(nintervals-1));
    grid->first_node_ids[0] = 0;
    for (int i=1; i<nintervals-1; i++)
        grid->first_node_ids[i] = grid->first_node_ids[i - 1] + grid->nnodes[i - 1] - 1;
    grid->dx = malloc(sizeof(double)*(nintervals-1));
    for (int i=0; i<nintervals-1; i++)
        grid->dx[i] = (intervals[i+1]-intervals[i])/(nnodes[i]-1);

    grid->nallnodes = grid->nnodes[0];
    for (int i=1; i<grid->nintervals-1; i++)
        grid->nallnodes += grid->nnodes[i] - 1;

    return grid;
}

int spl_grid_getNNodes(SplGrid* grid)
{
    return grid->nallnodes;
}

int spl_grid_getInterval(SplGrid* grid, double x)
{
    for (int i=1; i<grid->nintervals; i++)
        if (x <= grid->intervals[i])
            return i-1;
    //assert(0);
    return -1;
}

int spl_grid_get_intervali(SplGrid* grid, int node_index)
{
    int nsections = grid->nintervals - 1;
    for (int i=nsections-1; i>-1; i--)
        if (grid->first_node_ids[i] <= node_index)
            return i;
    assert(0);
    return -1;
}

int spl_grid_getNodeIndex(SplGrid* grid, double x)
{
    int ii = spl_grid_getInterval(grid, x);
    if (ii<0) {
        printf("spl_grid_getInterval() failed: x=%g\n", x);
        return -1;
    }
    double x0 = grid->intervals[ii] - 0.5*grid->dx[ii];
    if (x0 > x) {
        printf("x0=%g, x=%g\n", x0, x);
        return -1;
    }
    assert(x0 <= x);
    int i_node = (int)floor((x - x0) / grid->dx[ii]);
    i_node += grid->first_node_ids[ii];
    // printf("spl_grid_getNodeIndex()\n");
    // printf("  x=%g\n", x);
    // printf("  ii=%d\n", ii);
    // printf("  grid->intervals[ii]=%g\n", grid->intervals[ii]);
    // printf("  grid->dx[ii]=%g\n", grid->dx[ii]);
    return i_node;
}

double spl_grid_get_nodex(SplGrid* grid, int node_index)
{
    assert(node_index < grid->nallnodes);
    const int ii = spl_grid_get_intervali(grid, node_index);

    double x = 0.;
    int nsections = grid->nintervals - 1;
    if (ii < nsections-1 && node_index == grid->first_node_ids[ii+1])
        x = grid->intervals[ii+1];
    else
        x = grid->intervals[ii] + (node_index - grid->first_node_ids[ii]) * grid->dx[ii];

    return x;
}
