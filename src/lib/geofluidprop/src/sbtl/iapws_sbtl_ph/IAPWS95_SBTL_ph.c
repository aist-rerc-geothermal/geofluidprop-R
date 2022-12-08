#include "IAPWS95_SBTL_ph.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#ifndef _MAX_PATH
#define _MAX_PATH 1024
#endif

#ifndef _MAX_FNAME
#define _MAX_FNAME 1024
#endif

#ifndef BUFSIZ
#define BUFSIZ 1024
#endif

#define Tc    647.096  // [K]
#define pc    22.064e6  // [Pa]
#define rhoc  322.0  // [kg/m^3]

#define DEFAULT_FILENAME_TABLELIST_IAPWS95 "iapws95_sbtl_ph_table_list.txt"


// This function reads spline tables according to external file list.
static void iapws95sbtl_ph_read(IAPWS95SBTL_PH* sbtl)
{
//    printf("Reading iapws95-sbtl database...\n");

    // set root dir
    char *root_dir = NULL;
#ifndef _MSC_VER
    char root_dir_tmp[_MAX_PATH];
    root_dir = getenv("SBTL_IAPWS95_DIR");
    if (root_dir)
    {
        size_t len = strlen(root_dir);
        if (root_dir[len - 1] == '/' || root_dir[len - 1] == '\\')
            strcpy(root_dir_tmp, root_dir);
        else
            sprintf(root_dir_tmp, "%s/", root_dir);
        root_dir = root_dir_tmp;
//        printf("-> Root dir path for IAPWS95 SBTL tables: %s\n", root_dir);
    }
#endif

    // set path of the list file
    char list_path[_MAX_PATH];
    if (root_dir)
        sprintf(list_path, "%s%s", root_dir, DEFAULT_FILENAME_TABLELIST_IAPWS95);
    else
        strcpy(list_path, DEFAULT_FILENAME_TABLELIST_IAPWS95);

    // reading
    FILE *fp;
    fp = fopen(list_path, "r");

    if (fp == NULL) {
        printf("ERROR: cannot open '%s'\n", list_path);
#ifndef _MSC_VER
        printf("please check SBTL_IAPWS95_DIR\n");
#endif
        exit(1);
    }

    char key[BUFSIZ];
    char readline[BUFSIZ];
    char *ptr;
    char name[_MAX_FNAME];
    char dir_path[_MAX_PATH] = "./";
    char tmp_buff[_MAX_PATH];

    while (fgets(readline, BUFSIZ, fp) != NULL)
    {
        //printf("%s\n", readline);
        ptr = strtok(readline, " :");
        strcpy(key, ptr);

        if (strcmp(key, "*dir_path") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", dir_path) != 1) {
                sbtl_read_error(key);
            }
            const size_t len = strlen(dir_path);
            if (dir_path[len - 1] != '/' && dir_path[len - 1] != '\\')
                strcat(dir_path, "/");

            if (root_dir) {
                strcpy(tmp_buff, dir_path);
                sprintf(dir_path, "%s%s", root_dir, tmp_buff);
            }
        }
        else if (strcmp(key, "*sat_T_p") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_sat_T_p = sbtl_splinek1d_read(dir_path, name);
        }
        else if (strcmp(key, "*sat_hv_p") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_sat_hv_p = sbtl_splinek1d_read(dir_path, name);
        }
        else if (strcmp(key, "*sat_hl_p") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_sat_hl_p = sbtl_splinek1d_read(dir_path, name);
        }
        else if (strcmp(key, "*rho_ph") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_rho_ph = sbtl_splinek2d_read(dir_path, name);
        }
        else if (strcmp(key, "*T_ph") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_T_ph = sbtl_splinek2d_read(dir_path, name);
        }
        else if (strcmp(key, "*vis_ph") == 0) {
            ptr = strtok(NULL, " ");
            if (sscanf(ptr, "%s", name) != 1) {
                sbtl_read_error(key);
            }
            sbtl->tbl_vis_ph = sbtl_splinek2d_read(dir_path, name);
        }
        else {

        }
    }

    fclose(fp);

    if (sbtl->tbl_sat_T_p == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_sat_T_p is not set\n");
    if (sbtl->tbl_sat_hv_p == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_sat_hv_p is not set\n");
    if (sbtl->tbl_sat_hl_p == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_sat_hl_p is not set\n");
    if (sbtl->tbl_rho_ph == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_rho_ph is not set\n");
    if (sbtl->tbl_T_ph == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_T_ph is not set\n");
    if (sbtl->tbl_vis_ph == NULL)
        printf("***warning: IAPWS95_SBTL_ph: tbl_vis_ph is not set\n");
//    printf("Complete.\n");
}

// This function create IAPWS95SBTL_PH object and read spline tables.
void* iapws95sbtl_ph_create()
{
    IAPWS95SBTL_PH* sbtl = malloc(sizeof(IAPWS95SBTL_PH));
    sbtl->tbl_rho_ph = NULL;
    sbtl->tbl_T_ph = NULL;
    sbtl->tbl_vis_ph = NULL;
    sbtl->tbl_sat_T_p = NULL;
    sbtl->tbl_sat_hl_p = NULL;
    sbtl->tbl_sat_hv_p = NULL;

    iapws95sbtl_ph_read(sbtl);

    return sbtl;
}

void iapws95sbtl_ph_free(IAPWS95SBTL_PH* sbtl)
{
    sbtl_table2d_free(sbtl->tbl_rho_ph);
    sbtl_table2d_free(sbtl->tbl_T_ph);
    sbtl_table2d_free(sbtl->tbl_vis_ph);
    sbtl_table1d_free(sbtl->tbl_sat_T_p);
    sbtl_table1d_free(sbtl->tbl_sat_hv_p);
    sbtl_table1d_free(sbtl->tbl_sat_hl_p);

    free(sbtl);
}

double iapws95sbtl_ph_rho_pT(IAPWS95SBTL_PH* sbtl, double p, double T)
{
    double h = iapws95sbtl_ph_h_pT(sbtl, p, T);
    double rho = iapws95sbtl_ph_rho_ph(sbtl, p, h);
    return rho;
}

double iapws95sbtl_ph_rho_ph(IAPWS95SBTL_PH* sbtl, double p, double h)
{
    double p_tf = sbtl_transform[sbtl->tbl_rho_ph->tf1](p);
    double h_tf = sbtl_transform[sbtl->tbl_rho_ph->tf2](h);
    SplSplinek2D* tbl = sbtl->tbl_rho_ph->tbl;

    return spl_splinek2d_interpolate(tbl, p_tf, h_tf);
}

double iapws95sbtl_ph_h_pT(IAPWS95SBTL_PH* sbtl, double p, double T)
{
    double p_tf = sbtl_transform[sbtl->tbl_T_ph->tf1](p);
    double h_tf = spl_splinek2d_inverseY(sbtl->tbl_T_ph->tbl, p_tf, T);
    double h = sbtl_inv_transform[sbtl->tbl_T_ph->tf2](h_tf);
    return h;
}

double iapws95sbtl_ph_T_ph(IAPWS95SBTL_PH* sbtl, double p, double h)
{
    double p_tf = sbtl_transform[sbtl->tbl_T_ph->tf1](p);
    double h_tf = sbtl_transform[sbtl->tbl_T_ph->tf2](h);
    SplSplinek2D* tbl = sbtl->tbl_T_ph->tbl;

    return spl_splinek2d_interpolate(tbl, p_tf, h_tf);
}

double iapws95sbtl_ph_vis_ph(IAPWS95SBTL_PH* sbtl, double p, double h)
{
    double p_tf = sbtl_transform[sbtl->tbl_vis_ph->tf1](p);
    double h_tf = sbtl_transform[sbtl->tbl_vis_ph->tf2](h);
    SplSplinek2D* tbl = sbtl->tbl_vis_ph->tbl;

    return 1.e-6*spl_splinek2d_interpolate(tbl, p_tf, h_tf);
}

double iapws95sbtl_ph_sat_p_T(IAPWS95SBTL_PH* sbtl, double T)
{
    if (T > Tc) {
        printf("ERROR : Invalid temperature %g > Tc @ %s\n", T, __FUNCTION__);
        exit(0);
    }

    SplSplinek1D* tbl_Tp = sbtl->tbl_sat_T_p->tbl;
    double p_tf = spl_splinek1d_inverse(tbl_Tp, T);
    double p = sbtl_transform[sbtl->tbl_sat_T_p->tf](p_tf);
    return p;
}


double iapws95sbtl_ph_sat_T_p(IAPWS95SBTL_PH* sbtl, double p)
{
    if (p > pc) {
        printf("ERROR : Invalid pressure %g @ > pc %s\n", p, __FUNCTION__);
        // exit(0);
        return -1;
    }
    double p_tf = sbtl_transform[sbtl->tbl_sat_T_p->tf](p);
    SplSplinek1D* tbl = sbtl->tbl_sat_T_p->tbl;
    double T = spl_splinek1d_interpolate(tbl, p_tf);
    return T;
}

void iapws95sbtl_ph_sat_rho_p(IAPWS95SBTL_PH* sbtl, double p, double* rhol, double* rhov)
{
    if (p > pc) {
        printf("ERROR : Invalid p %g > pc @ %s\n", p, __FUNCTION__);
        *rhol = *rhov = -1;
        return;
    }
    double hl, hv;
    iapws95sbtl_ph_sat_h_p(sbtl, p, &hl, &hv);
    *rhol = iapws95sbtl_ph_rho_ph(sbtl, p, hl);
    *rhov = iapws95sbtl_ph_rho_ph(sbtl, p, hv);
}

double iapws95sbtl_ph_sat_hl_p(IAPWS95SBTL_PH* sbtl, double p)
{
    if (p > pc) {
        printf("ERROR : Invalid p %g > pc @ %s\n", p, __FUNCTION__);
        return -1;
    }
    double h0;
    {
        double p_tf = sbtl_transform[sbtl->tbl_sat_hl_p->tf](p);
        SplSplinek1D* tbl = sbtl->tbl_sat_hl_p->tbl;
        h0 = spl_splinek1d_interpolate(tbl, p_tf);
    }

    double T = iapws95sbtl_ph_sat_T_p(sbtl, p);
    double p_tf = sbtl_transform[sbtl->tbl_T_ph->tf1](p);
    double h_tf = spl_splinek2d_inverseY2(sbtl->tbl_T_ph->tbl, p_tf, T, h0);
    double h = sbtl_inv_transform[sbtl->tbl_T_ph->tf2](h_tf);
    return h;
}

double iapws95sbtl_ph_sat_hv_p(IAPWS95SBTL_PH* sbtl, double p)
{
    if (p > pc) {
        printf("ERROR : Invalid p %g > pc @ %s\n", p, __FUNCTION__);
        return -1;
    }
    double h0;
    {
        double p_tf = sbtl_transform[sbtl->tbl_sat_hv_p->tf](p);
        SplSplinek1D* tbl = sbtl->tbl_sat_hv_p->tbl;
        h0 = spl_splinek1d_interpolate(tbl, p_tf);
    }

    double T = iapws95sbtl_ph_sat_T_p(sbtl, p);
    double p_tf = sbtl_transform[sbtl->tbl_T_ph->tf1](p);
    double h_tf = spl_splinek2d_inverseY2(sbtl->tbl_T_ph->tbl, p_tf, T, h0);
    double h = sbtl_inv_transform[sbtl->tbl_T_ph->tf2](h_tf);
    return h;
}

void iapws95sbtl_ph_sat_h_p(IAPWS95SBTL_PH* sbtl, double p, double* hl, double* hv)
{
    *hl = iapws95sbtl_ph_sat_hl_p(sbtl, p);
    *hv = iapws95sbtl_ph_sat_hv_p(sbtl, p);
}
