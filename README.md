# pycaz - a collection of analysis functions
pycaz - pronounced as 'package' - is a collections of functions and datastructure written in python for the analysis of geophysical data. The package only focuses on the application procedure, using a functional pattern - which means that there will be no side-effect on the data that are being worked on. The package will rely heavily on the pandas-scipy-numpy stack.

The current package structure is inherited from the original package which was built for handling SCHISM model outputs. However, as there is now a dedicted project for handling SCHISM model, and with the request from the SCHISM developers the original name of this package 'pyschism' is relased from pip and conda repository for their use. Consequently, the name of this repository is changed along with its intended function. 

Currently the inherited packages has several modules which is listed below - 
* core : provides core data structures
* io: input/output functionality based on core data structures
* plot: plotting routines
* preprocess: routines to help preprocessing tasks based on core and io
* postprocess: routines to help postprocessing tasks based on core and io
* tide: routines to help with tidal analysis
* cyclone: routines to generate analytical cyclone fields
* model: wrapper routines to build schism model

Evidently, this will be continuously removed in the future releases. The intended structure will follow a collection of functions per file basis. For example, schism.py may include function to read mesh, create bctides, process discharge etc.

# Usage
The module is imported using the following command - 

```
from pycaz import schism # to import schism related functions into the workspace
```

# Contact
Feel free to use the scripts and if you find any bug report to https://github.com/jamal919/pycaz