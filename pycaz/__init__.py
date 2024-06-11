#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`pycaz`, a play on sounds on package, is a python package gathering the analysis tools 
developed during modelling with SCHISM, but now expanded to contains all sorts of analysis 
and modelling functionalities. The original package was called (vary unoriginally) pyschism.
However following the development of the toolbox by the SCHISM/NOAA called pyschism, the
name is released for their use. The name `pycaz` comes from combination of py (short of python)
and caz (from Bengali word "কাজ", which means work).

.. danger::
    This toolbox is still in its major version 0. The interface might change during its 
    development, caution is adviced. Please direct your advice, contribution, 
    bug reports to the [github repository](https://github.com/jamal919/pycaz)

# Structure of the package
The package is divided into several modules - 

- `schism` : Modules to prepare schism input/output, handle model config.
- `tide` : Modules to handle tide data, preset for tidal analysis, supplementary function to utide.
- `webdata` : Modules to download data from the web portals - NOAA, Copernicus etc.
- `cyclone` : Modules to read cyclone tracks, generate wind and pressure field.
- `core` : Data classes, currently not used and expected to be depricated.
- `convert` : A conversion function between various units, has very limited functionalities.

Each module can be imported using the standard python import commands - `from pycaz import <module>`

# Installation
The `pycaz` package depends on standard scientific python libraries - `numpy`, `scipy`, 
`matplotlib`. For map plotting `cartopy` is used. Tidal anlaysis is mostly based on `utide`.

# How to use
Since pycaz is a collection of heterogeneous modules, the documentation for each module
is provided under the Submodule itself. The submodules are listed in the left pane under
the header **Submodules**.


# License
This toolbox is licensed under Apache License 2.0 - A permissive license whose main 
conditions require preservation of copyright and license notices. Contributors provide 
an express grant of patent rights. Licensed works, modifications, and larger works may 
be distributed under different terms and without source code.

"""
__version__ = '0.2'