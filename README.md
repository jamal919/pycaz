[![PyPI download total](https://img.shields.io/pypi/dt/pyschism.svg)](https://pypi.python.org/pypi/pyschism/)

# pyschism - Python SCHISM Toolbox
Python SCHISM Toolbox (pyschism) is a collections of modules and classes written in python to build, edit, analyze input files, forcings for SCHISM modelling system. This package is currently under development and a lot of changes are made everyday. 

The package several modules which is listed below - 
* core : provides core data structures
* io: input/output functionality based on core data structures
* plot: plotting routines
* preprocess: routines to help preprocessing tasks based on core and io
* postprocess: routines to help postprocessing tasks based on core and io
* tide: routines to help with tidal analysis
* cyclone: routines to generate analytical cyclone fields
* model: wrapper routines to build schism model

Currently we are developing scripts to tackle individual tasks. Please check out the scripts directory to see developed scripts which eventually got ported to the package.


# Usage
The module is imported using the following command - 

```
import schism
```

# Contact
Feel free to use the scripts and if you find any bug report to https://github.com/jamal919/pyschism