The geofluidprop R package: a wrapper for the geofluidprop C library
================

<!-- badges: start -->
[![R-CMD-check](https://github.com/aist-rerc-geothermal/geofluidprop-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aist-rerc-geothermal/geofluidprop-R/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
  
This is an R wrapper for [`geofluidprop`](https://github.com/aist-rerc-geothermal/geofluidprop). geofluidprop is a C library for computing physical properties of geological fluids.
Please note that this R package is under development and it supports limited functionalities of 'geofluidprop'.

## Installation

```r
devtools::install_github('aist-rerc-geothermal/geofluidprop-R')
```


## Examples 

### Water density [kg/m^3] at 1 MPa and 20 degree C using IAPWS-95 formulation:

```r
iapws95_rho_pT(p=1e6, T=293.15)
## [1] 998.6184
```

### H<sub>2</sub>O-NaCl fluid properties at 10 MPa, 200 degree C, 0.1 wt% NaCl based on Driesner & Heinrich (2007) and Driesner (2007):

```r
prop = driesner07_H2O_NaCl_get_properties_TpX(TK=273.15+200, p=10e6, X=0.1e-2)
str(prop)
## List of 7
##  $ given_condition:List of 3
##   ..$ TdegC        : num 200
##   ..$ pMPa         : num 10
##   ..$ wtp_NaCl_bulk: num 0.1
##  $ critical_point :List of 2
##   ..$ TdegC: num 377
##   ..$ pMPa : num 22.7
##  $ phase          : chr "L"
##  $ density        : num 872
##  $ h              : num 855
##  $ H              : num 745609
##  $ viscosity      : num 0.000137
```


### H<sub>2</sub>O-NaCl fluid electrical conductivity at 10 MPa, 200 degree C, 0.1 wt% NaCl using the formulation in Watanabe et al. (2021) :

```r
H2ONaCl_ec_WatanabeEtAl2021_Tpm(TK=273.15+200, pMPa=10, m=H2ONaCl_massfrac_to_b(0.1e-2))
## [1] 0.9172785
```

