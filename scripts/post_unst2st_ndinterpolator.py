# -*- coding: utf-8 -*-
"""
Water level interpolation from SCHISM grid to structured grid.
Created on Wed May 23 14:12:23 2018

@author: Jamal Khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import time
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from scipy.interpolate import LinearNDInterpolator

# input dataset
ml = Dataset('./schout.nc')

intime = ml.variables['time'][:]
x = ml.variables['SCHISM_hgrid_node_x'][:]
y = ml.variables['SCHISM_hgrid_node_y'][:]
depth = ml.variables['depth'][:]

# Create netcdf file
indata = Dataset('interpolated.nc', 'w', format='NETCDF4_CLASSIC')
lonrange = np.arange(86.19563, 93.11962, 0.006)
latrange = np.arange(20.02389, 24.24789, 0.006)

# Dimensions
dlon = indata.createDimension(dimname='lon', size=len(lonrange))
dlat = indata.createDimension(dimname='lat', size=len(latrange))
dtime =indata.createDimension(dimname='time', size=None)

# Variables
vtimes = indata.createVariable(varname='time', datatype=np.float64, dimensions=('time'))
vtimes.units = 'hours since 1970-01-01 00:00:00'
vtimes.calendar = 'standard'

vlats = indata.createVariable(varname='lat', datatype=np.float32, dimensions=('lat'))
vlats.units = 'degrees_north'
vlats.long_name = 'latitude'
vlats[:] = latrange

vlons = indata.createVariable(varname='lon', datatype=np.float32, dimensions=('lon'))
vlons.units = 'degrees_east'
vlons.long_name = 'longitude'
vlons[:] = lonrange

velev = indata.createVariable(varname='elev', datatype=np.float32, dimensions=('time', 'lat', 'lon'), chunksizes=(1, len(latrange), len(lonrange)))
velev.units = 'm'
velev.long_name = 'Water level above MSL'

# Global attribute
indata.description = 'Water level interpolated from SCHISM tide experiment'
indata.history = 'Created ' + time.ctime(time.time())
indata.source = 'SCHISM v5.6.1 - Mesh v21.6'


# put proper time on the output netcdf
intime = num2date(times=intime, units=ml.variables['time'].units, calendar='standard')
vtimes[:] = date2num(dates=intime, units=vtimes.units, calendar=vtimes.calendar)

# output grid to interpolate
grid_x, grid_y = np.mgrid[86.19563:93.11962:1154j, 20.02389:24.24789:704j]

# interpolate depth to mask out the land value
grid_z0 = LinearNDInterpolator(list(zip(x, y)), depth, fill_value=np.nan)
depth_z0 = grid_z0(grid_x, grid_y)
depth_z0_mask = depth_z0 <= 0

# interpolate over time
for i in np.arange(len(intime)):
    elev = ml.variables['elev'][i, :]
    
    grid_z0 = LinearNDInterpolator(list(zip(x, y)), elev, fill_value=np.nan)
    grid_z0_val = grid_z0(grid_x, grid_y)
    grid_z0_val[depth_z0_mask] = np.nan
    
    velev[i, :, :] = np.transpose(grid_z0_val)
    
    if i%100==0:
        indata.sync()
        print(i)

# Closing files
indata.close()
