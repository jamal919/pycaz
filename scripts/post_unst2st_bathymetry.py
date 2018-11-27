#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Bathymetry interpolation from SCHISM grid to structured grid.
Created on Wed May 23 14:12:23 2018

@author: Jamal Khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import time
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import LinearNDInterpolator

# input dataset
ml = Dataset('./schout_1_junk.nc')

x = ml.variables['SCHISM_hgrid_node_x'][:]
y = ml.variables['SCHISM_hgrid_node_y'][:]
depth = ml.variables['depth'][:]

# Create netcdf file
indata = Dataset('depth.nc', 'w', format='NETCDF4_CLASSIC')
lonrange = np.arange(86.19563, 93.11962, 0.006)
latrange = np.arange(20.02389, 24.24789, 0.006)

# Dimensions
dlon = indata.createDimension(dimname='lon', size=len(lonrange))
dlat = indata.createDimension(dimname='lat', size=len(latrange))

# Variables
vlats = indata.createVariable(varname='lat', datatype=np.float32, dimensions=('lat'))
vlats.units = 'degrees_north'
vlats.long_name = 'latitude'
vlats[:] = latrange

vlons = indata.createVariable(varname='lon', datatype=np.float32, dimensions=('lon'))
vlons.units = 'degrees_east'
vlons.long_name = 'longitude'
vlons[:] = lonrange

vdepth = indata.createVariable(varname='depth', datatype=np.float32, dimensions=('lat', 'lon'))
vdepth.units = 'm'
vdepth.long_name = 'Bathymetry from MSL (Positive downward)'

# Global attribute
indata.description = 'Bathymetry interpolated from SCHISM tide experiment'
indata.history = 'Created ' + time.ctime(time.time())
indata.source = 'SCHISM v5.6.1 - Mesh v21.6'

# output grid to interpolate
grid_x, grid_y = np.mgrid[86.19563:93.11962:1154j, 20.02389:24.24789:704j]

# interpolate bathymetry
grid_z0 = LinearNDInterpolator(list(zip(x, y)), depth, fill_value=np.nan)
depth_z0 = grid_z0(grid_x, grid_y)
vdepth[:, :] = np.transpose(depth_z0)

# Closing files
indata.close()