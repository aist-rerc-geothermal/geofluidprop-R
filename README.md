# R package for the `geofluidprop` library
R package for [`geofluidprop`](https://github.com/aist-rerc-geothermal/geofluidprop), that is a C/C++ library for computing physical properties of geological fluids (currently we support only water and H<sub>2</sub>O-NaCl fluids)

## Installation

```r
devtools::install_github('aist-rerc-geothermal/geofluidprop-R')
```

Please note that, at the moment, this package works only on Windows x64.
For other environments, one needs to build the codes by themselves.
See details in [`geofluidprop`](https://github.com/aist-rerc-geothermal/geofluidprop).


## Examples 

### Water density at 1 MPa and 20 degree C using IAPWS-95 formulation:

```r
iapws95_rho_pT(p=1e6, T=293.15)
## [1] 998.6184
```

### H<sub>2</sub>O-NaCl fluid density at 10 MPa, 200 degree C, 0.1 wt% NaCl:

```r
driesner07_H2O_NaCl_rho_pTx(p_Pa=10e6, T_K=273.15+200, x=H2ONaCl_massfrac_to_x(0.1e-2))
```

### H<sub>2</sub>O-NaCl fluid phase at 10 MPa, 400 degree C, 0.1 wt% NaCl:

```r
driesner07_H2O_NaCl_phase_name_Tpx(TK=273.15+400, p=10e6, x=H2ONaCl_massfrac_to_x(0.1e-2))
## [1] "VH"
```

### H<sub>2</sub>O-NaCl fluid electrical conductivity at 10 MPa, 200 degree C, 0.1 wt% NaCl using the formulation in Watanabe et al. (2021) :

```r
H2ONaCl_ec_WatanabeEtAl2021_Tpm(TK=273.15+200, pMPa=10, m=H2ONaCl_massfrac_to_b(0.1e-2))
## [1] 0.9172785
```


### To know more:
geofluidprop: https://github.com/aist-rerc-geothermal/geofluidprop
R-package for geofluidprop: https://github.com/aist-rerc-geothermal/geofluidprop-R
