R interface to geofluidprop
====

## How to use
Unfortunately we are not yey ready to publish a R-package for geofluidprop.  
Nevertheless, one can use some functions of geofluidprop from R with the following steps:

1. Build a shared library for R  
This can be done by passing `-DBUILD_R_LIB=ON` to the cmake arguments. 
Internally, we use SWIG to generate R interface.  
After running `make`, the following files will be created under `build/lib/r` directory.
    * `geofluidprop.R`
    * `geofluidprop.so` or `geofluidprop.dll`

2. Copy all the R related files to your work directory  
Copy the above files and all the R script files under `wrapper/R` to your work directory for R.

3. Load the wrapper script into R  
In R, move to the work directory and type 
    ```
    source("geofluidprop_init.R")
    ```  
    This should load the shared library and related wrapper functions.
    Then, one can type, e.g.,
    ```
    water_density <- iapws95_rho_pT(p=1e6, T=293.15)
    ```  

## Additional functions in the R scripts
`wrapper/R` directory includes the following R scripts, that implement additional functions based on  `geofluidprop` functionality: 
- [H2ONaCl_electrical_conductivity.R](H2ONaCl_electrical_conductivity.R): Implementing several formulations to calculate electrical conductivity of H2O-NaCl fluids.
