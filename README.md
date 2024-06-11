# pycaz - a collection of analysis functions
`pycaz`, a play on sounds on package, is a python package gathering the analysis tools developed during modelling with 
SCHISM, but now expanded to contains all sorts of analysis  and modelling functionalities. The original package was 
called (vary unoriginally) pyschism. However following the development of the toolbox by the SCHISM/NOAA called 
[pyschism](https://github.com/schism-dev/pyschism), the name is released for their use. The name `pycaz` comes from combination of `py` (short of python) 
and `caz` (from Bengali word "কাজ", which means work).

Being a collection of different types of analysis, pre- and post-processing methods, some modules follows a more
object oriented pattern, and some follows a more functional/procedural pattern. The documentation of the modules are 
not yet well-developed, and will be incorporated in a best-effort basis.

# Installation
For best experience, the package should be installed in a conda environment, which can be obtained through any conda
distribution such a [anaconda](https://www.anaconda.com/download) or [miniconda](https://docs.anaconda.com/free/miniconda/).

```shell
conda create -n pycaz -c conda-forge python numpy scipy matplotlib xarray netcdf4 utide cmocean rioxarray tqdm ipykernel pyproj cartopy geopandas shapely jupyterlab jupyter notebook 
```

Then the toolbox can be installed using - 

```shell
conda activate pycaz
pip install .
```

# Contact
Feel free to use the scripts and if you find any bug report in the [repo issues](https://github.com/jamal919/pycaz/issues).