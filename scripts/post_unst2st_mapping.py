# -*- coding: utf-8 -*-
"""
Creating a structured netcdf from unstructured SCHISM netcdf output.
This script use nearby point mapping for interpolation using scipy.spatial.

This script is an alternative version of unstrucutred2structured_ndinterpolator.py.

In future the script will be added as a part of schism python toolbox module. 

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
from netCDF4 import Dataset, num2date, date2num
import time
import numpy as np
from scipy.spatial import distance
from pyproj import Proj

def closest_node(node, nodes):
    dist_mat = distance.cdist([node], nodes)
    closest_dist = dist_mat.min()
    closest_ind = dist_mat.argmin()
    return closest_dist, closest_ind
    
# Define projection
btm = Proj('+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +units=m +no_defs')

# Loading data
ml = Dataset('./schout_1_junk.nc')

intime = ml.variables['time'][:]
x = ml.variables['SCHISM_hgrid_node_x'][:]
y = ml.variables['SCHISM_hgrid_node_y'][:]
depth = ml.variables['depth'][:]
x_btm, y_btm = btm(x, y)
xy_btm = np.array([x_btm, y_btm]).transpose()

# output grid
grid_x, grid_y = np.mgrid[86.19563:93.11962:1154j, 20.02389:24.24789:704j]
grid_x_val = grid_x.flatten()
grid_y_val = grid_y.flatten()
grid_x_btm, grid_y_btm = btm(grid_x_val, grid_y_val)
grid_xy_btm = np.array([grid_x_btm, grid_y_btm]).transpose()

# find nearest position
print('Finding the nearest position. This may take a while...')
grid_position = np.array([closest_node(grid_xy_btm[i, :], xy_btm) for i in np.arange(0, len(grid_xy_btm))])
np.savetxt(fname='nearest_position.csv', X=grid_position, fmt='%i', delimiter=',')

grid_dist = np.array(grid_position[:, 0].reshape(grid_x.shape), dtype='int16')
grid_pos = np.array(grid_position[:, 1].reshape(grid_x.shape), dtype='int16')

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

velev = indata.createVariable(varname='elev', datatype=np.float32, dimensions=('time', 'lat', 'lon'))
velev.units = 'm'
velev.long_name = 'Water level above MSL'

# Global attribute
indata.description = 'Water level interpolated from SCHISM tide experiment'
indata.history = 'Created ' + time.ctime(time.time())
indata.source = 'SCHISM v5.6.1 - Mesh v21.6'

# put proper time on the output netcdf
intime = num2date(times=intime, units=ml.variables['time'].units, calendar='standard')
vtimes[:] = date2num(dates=intime, units=vtimes.units, calendar=vtimes.calendar)

# interpolate over time
for i in np.arange(len(intime)):
    elev = ml.variables['elev'][i, :]
    elev[grid_dist >= 7000] = np.NaN 
    elevint = elev[grid_pos]
    
    velev[i, :, :] = np.transpose(elevint)
    
    if i%100==0:
        indata.sync()
        print(i)

# Closing files
indata.close()
print('Dataset is written to file! Process completed.')