# pycaz - a collection of analysis functions
`pycaz`, a play on sounds on package, is a python package gathering the analysis tools 
developed during modelling with SCHISM, but now expanded to contains all sorts of analysis 
and modelling functionalities. The original package was called (vary unoriginally) pyschism.
However following the development of the toolbox by the SCHISM/NOAA called [pyschism](https://github.com/schism-dev/pyschism), the name is released for their use. The name `pycaz` comes from combination of `py` (short of python) and `caz` (from Bengali word "কাজ", which means work).

The package try to focuses on the application procedure, using a functional pattern - which means that there will be no side-effect on the data that are being worked on. Indeed some dataclasses will get in the way from time to time. The package will rely heavily on the pandas-scipy-numpy stack.

Currently the inherited packages has several modules which is listed below - 
* core : provides core data structures
* io: input/output functionality based on core data structures
* plot: plotting routines
* preprocess: routines to help preprocessing tasks based on core and io
* postprocess: routines to help postprocessing tasks based on core and io
* tide: routines to help with tidal analysis
* cyclone: routines to generate analytical cyclone fields
* model: wrapper routines to build schism model

Evidently, this will be continuously removed in the future releases. The intended structure will follow a collection of 
functions per file basis. For example, schism.py may include function to read mesh, create bctides, process discharge etc.

# Required packages
* python>=3.7
* numpy
* scipy 
* matplotlib
* cartopy
* cmocean
* pandas 
* xarray 
* netCDF4 
* jupyter notebook  

# Usage
The module is imported using the following command - 

```
from pycaz import schism # to import schism related functions into the workspace
```

# Contact
Feel free to use the scripts and if you find any bug report to https://github.com/jamal919/pycaz