#include "sbtl_utility.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef _MAX_PATH
#define _MAX_PATH 1024
#endif

// This is dummy function for non-transformation.
double tf_normal(double v)
{
    return v;
}

// This function is for logarithmic transformation.
double tf_log(double v)
{
    return log(v);
}

double tf_log_MPa(double v)
{
    return log(v*1.e-2);
}

// This is inverse function for non-transformation.
double tf_inv_normal(double v_tf)
{
    return v_tf;
}

// This function is for logarithmic inverse transformation.
double tf_inv_log(double v_tf)
{
    return exp(v_tf);
}

double tf_inv_log_MPa(double v_tf)
{
    return exp(v_tf)*1e2;
}

// Function arrays for parameter transformation
SBTL_TRANSFORM sbtl_transform[4] = {
    tf_normal,
    tf_log,
    tf_log_MPa,
    NULL
};

// Function arrays for parameter inverse transformation
SBTL_INV_TRANSFORM sbtl_inv_transform[4] = {
    tf_inv_normal,
    tf_inv_log,
    tf_inv_log_MPa,
    NULL
};

// This function reads 1-d spline table.
SplSplinek1D* sbtl_splinek1d_table_read(const char* path)
{
    FILE *fp;
    fp = fopen(path, "rb");

    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    SplSplinek1D* tbl = spl_splinek1d_read(fp);

    fclose(fp);

    return tbl;
}

// This function reads 2-d spline table.
SplSplinek2D* sbtl_splinek2d_table_read(const char* path)
{
    FILE *fp;
    fp = fopen(path, "rb");

    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    SplSplinek2D* tbl = spl_splinek2d_read(fp);

    fclose(fp);

    return tbl;
}

// This function reads 3-d spline table.
SplSplinek3D* sbtl_splinek3d_table_read(const char* path)
{
    FILE *fp;
    fp = fopen(path, "rb");

    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    SplSplinek3D* tbl = spl_splinek3d_read(fp);

    fclose(fp);

    return tbl;
}

// This function reads auxiliary data for 1-d sbtl.
static void sbtl_splinek1d_info_read(const char* path, SbtlTable1d* tbl)
{
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    size_t ret = fread(&tbl->tf, sizeof(int), 1, fp);

    fclose(fp);
}

// This function reads auxiliary data for 2-d sbtl.
static void sbtl_splinek2d_info_read(const char* path, SbtlTable2d* tbl)
{
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    size_t ret = fread(&tbl->tf1, sizeof(int), 1, fp);
    ret = fread(&tbl->tf2, sizeof(int), 1, fp);

    fclose(fp);
}

// This function reads auxiliary data for 3-d sbtl.
static void sbtl_splinek3d_info_read(const char* path, SbtlTable3d* tbl)
{
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        printf("ERROR: cannot open %s\n", path);
        exit(1);
    }

    size_t ret = fread(&tbl->tf1, sizeof(int), 1, fp);
    ret = fread(&tbl->tf2, sizeof(int), 1, fp);
    ret = fread(&tbl->tf3, sizeof(int), 1, fp);

    fclose(fp);
}

// This function cerate 1-d sbtl table object and reads data from database file.
SbtlTable1d* sbtl_splinek1d_read(const char* dir_path, const char* name)
{
    SbtlTable1d* tbl = malloc(sizeof(SbtlTable1d));

    // check path
    char dir_path_tmp[_MAX_PATH];

    size_t len = strlen(dir_path);
    if (dir_path[len - 1] == '/' || dir_path[len - 1] == '\\') {
        sprintf(dir_path_tmp, "%s", dir_path);
    }
    else {
        sprintf(dir_path_tmp, "%s/", dir_path);
    }

    // read spline table
    char table_path[_MAX_PATH];
    sprintf(table_path, "%s%s.bin", dir_path_tmp, name);
    tbl->tbl = sbtl_splinek1d_table_read(table_path);

    // read other information
    char info_path[_MAX_PATH];
    sprintf(info_path, "%s%s_info.bin", dir_path_tmp, name);
    sbtl_splinek1d_info_read(info_path, tbl);

    return tbl;
}

// This function cerate 2-d sbtl table object and reads data from database file.
SbtlTable2d* sbtl_splinek2d_read(const char* dir_path, const char* name)
{
    SbtlTable2d* tbl = malloc(sizeof(SbtlTable2d));

    // check path
    char dir_path_tmp[_MAX_PATH];

    size_t len = strlen(dir_path);
    if (dir_path[len - 1] == '/' || dir_path[len - 1] == '\\') {
        sprintf(dir_path_tmp, "%s", dir_path);
    }
    else {
        sprintf(dir_path_tmp, "%s/", dir_path);
    }

    // read spline table
    char table_path[_MAX_PATH];
    sprintf(table_path, "%s%s.bin", dir_path_tmp, name);
    tbl->tbl = sbtl_splinek2d_table_read(table_path);

    // read other information
    char info_path[_MAX_PATH];
    sprintf(info_path, "%s%s_info.bin", dir_path_tmp, name);
    sbtl_splinek2d_info_read(info_path, tbl);

    return tbl;
}

// This function cerate 3-d sbtl table object and reads data from database file.
SbtlTable3d* sbtl_splinek3d_read(const char* dir_path, const char* name)
{
    SbtlTable3d* tbl = malloc(sizeof(SbtlTable3d));

    // check path
    char dir_path_tmp[_MAX_PATH];

    size_t len = strlen(dir_path);
    if (dir_path[len - 1] == '/' || dir_path[len - 1] == '\\') {
        sprintf(dir_path_tmp, "%s", dir_path);
    }
    else {
        sprintf(dir_path_tmp, "%s/", dir_path);
    }

    // read spline table
    char table_path[_MAX_PATH];
    sprintf(table_path, "%s%s.bin", dir_path_tmp, name);
    tbl->tbl = sbtl_splinek3d_table_read(table_path);

    // read other information
    char info_path[_MAX_PATH];
    sprintf(info_path, "%s%s_info.bin", dir_path_tmp, name);
    sbtl_splinek3d_info_read(info_path, tbl);

    return tbl;
}

// This function releases memory of 1-d sbtl table object.
void sbtl_table1d_free(SbtlTable1d* table1d)
{
    if (table1d)
        spl_splinek1d_free(table1d->tbl);
}

// This function releases memory of 2-d sbtl table object.
void sbtl_table2d_free(SbtlTable2d* table2d)
{
    if (table2d)
        spl_splinek2d_free(table2d->tbl);
}

// This function releases memory of 3-d sbtl table object.
void sbtl_table3d_free(SbtlTable3d* table3d)
{
    if (table3d)
        spl_splinek3d_free(table3d->tbl);
}

// This function is for error display in reading sbtl file list and spline table.
void sbtl_read_error(char* key)
{
    printf("READ_ERROR : @%s\n", key);
    exit(1);
}
