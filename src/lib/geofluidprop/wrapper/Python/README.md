Python interface to geofluidprop
====

## How to use
One can use some functions of geofluidprop from Python with the following steps:

1. Build a shared library for Python  
This can be done by passing `-DBUILD_PYTHON_LIB=ON` to the cmake arguments. 
Internally, we use SWIG to generate Python interface.  
After running `make`, the following files will be created under `build/lib/python` directory.
    * `geofluidprop.py`
    * `_geofluidprop_py.so` or `_geofluidprop_py.dll`

2. Copy the above two files to your work directory  

3. Call the geofluidprop functions in Python, e.g., 
    ```python
    from geofluidprop import *
    water_density = iapws95_rho_pT(p=1e6, T=293.15)
    ```  
