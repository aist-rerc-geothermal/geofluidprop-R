#ifndef _SBTL_UTILITY_H__
#define _SBTL_UTILITY_H__

#include "spline/spl_splinek1d.h"
#include "spline/spl_splinek2d.h"
#include "spline/spl_splinek3d.h"

// Enumerate type of transformation function
enum TF_TYPE {
    TF_TYPE_NORMAL = 0,
    TF_TYPE_LOG,
    TF_TYPE_LOG_MPa
};

// Sbtl object for 1-d spline table
typedef struct
{
    SplSplinek1D* tbl;
    enum TF_TYPE tf;

} SbtlTable1d;

// Sbtl object for 2-d spline table
typedef struct
{
    SplSplinek2D* tbl;
    enum TF_TYPE tf1;
    enum TF_TYPE tf2;

} SbtlTable2d;

// Sbtl object for 3-d spline table
typedef struct
{
    SplSplinek3D* tbl;
    enum TF_TYPE tf1;
    enum TF_TYPE tf2;
    enum TF_TYPE tf3;

} SbtlTable3d;

double tf_normal(double v);
double tf_log(double v);
double tf_log_MPa(double v);
//double tf_pv(double v, double p);
double tf_inv_normal(double v_tf);
double tf_inv_log(double v_tf);
double tf_inv_log_MPa(double v);
//double tf_inv_pv(double v_tf, double p);

typedef double(*SBTL_TRANSFORM)(double x);
typedef double(*SBTL_INV_TRANSFORM)(double x);

// Function array of transformation function
extern SBTL_TRANSFORM sbtl_transform[4];

// Function array of inverse transformation function
extern SBTL_INV_TRANSFORM sbtl_inv_transform[4];

SbtlTable1d* sbtl_splinek1d_read(const char* dir_path, const char* name);
SbtlTable2d* sbtl_splinek2d_read(const char* dir_path, const char* name);
SbtlTable3d* sbtl_splinek3d_read(const char* dir_path, const char* name);
void sbtl_table1d_free(SbtlTable1d* table1d);
void sbtl_table2d_free(SbtlTable2d* table2d);
void sbtl_table3d_free(SbtlTable3d* table2d);

void sbtl_read_error(char* key);

#endif