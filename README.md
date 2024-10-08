# pycaz - a package for doing work!
`pycaz`, a play on sounds on package, is a python package gathering the analysis tools developed during modelling with 
SCHISM, but now expanded to contains all sorts of analysis  and modelling functionalities. The original package was 
called (vary unoriginally) pyschism. However, following the development of the toolbox by the SCHISM/NOAA called 
[pyschism](https://github.com/schism-dev/pyschism), the name is released for their use and `pycaz` was born. The name `pycaz` comes from combination of 
`py` (obviously short form of python) and `caz` (from Bengali word "কাজ", which means work). It is also intended as a pun
since `pycaz` can also be read like 'package' (an overarching package to collect processing, analysis, tools). 

Being a collection of different types of analysis, pre- and post-processing methods, some modules follows a more
object oriented pattern, and some follows a more functional/procedural pattern. The documentation of the modules are 
not yet well-developed, and will be incorporated in a best-effort basis. 

# Installation
For best experience, the package should be installed in a conda environment, which can be obtained through any conda
distribution such a [anaconda](https://www.anaconda.com/download) or [miniconda](https://docs.anaconda.com/free/miniconda/).

```shell
conda create -n pycaz -c conda-forge python numpy scipy matplotlib pandas openpyxl xarray netcdf4 utide cmocean rioxarray tqdm ipykernel pyproj cartopy geopandas shapely beautifulsoup4 lxml jupyterlab jupyter notebook=6.4.12 ipywidgets jupyter_contrib_nbextensions hydromt hydromt_sfincs
```

Then the toolbox can be installed using - 

```shell
git clone https://github.com/jamal919/pycaz
cd pycaz
conda activate pycaz
pip install .
```

If you do not want to install the toolbox, or want to use it in a development mode then do the following the in the
beginning of your python script/notebook - 

```python
import sys
sys.path.append('/path/to/pycaz') # the outer directory downloaded from git

# now the modules/functions can be loaded as necessary
from pycaz.cyclone.jtwc import read_jtwc
```

# Contact
If you find any bug, please report in the [repo issues](https://github.com/jamal919/pycaz/issues).