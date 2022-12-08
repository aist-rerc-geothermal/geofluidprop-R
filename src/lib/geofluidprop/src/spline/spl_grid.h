
#ifndef _SPL_GRID_H__
#define _SPL_GRID_H__

#include <stdio.h>

// Structure containing 1-d grid information
typedef struct
{
    int nintervals;         // the number of intervals
    double* intervals;      // positions of the intervals
    int* nnodes;            // the number of nodes in each interval
    int* first_node_ids;    // the first node index of each interval
    double* dx;             // cell size of each interval
    int nallnodes;          // the total number of nodes in all intervals
} SplGrid;


// create a grid object of a single interval
SplGrid* spl_grid_createUniform(double vmin, double vmax, int nnodes);

// create a grid object
SplGrid* spl_grid_create(int nintervals, double* intervals, int* nnodes);

// free a grid object
void spl_grid_free(SplGrid* grid);

// write a grid object into a file
void spl_grid_write(SplGrid* grid, FILE* stream);

// read a grid object from a file
SplGrid* spl_grid_read(FILE* stream);

// get the total number of nodes in the given grid
int spl_grid_getNNodes(SplGrid* grid);

// get interval index at the given position
int spl_grid_getInterval(SplGrid* grid, double x);

// get interval index containing the given node
int spl_grid_get_intervali(SplGrid* grid, int node_index);

// get the node index at the fiven position
int spl_grid_getNodeIndex(SplGrid* grid, double x);

// get the node position
double spl_grid_get_nodex(SplGrid* grid, int node_index);

#endif
