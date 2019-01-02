# -*- coding: utf-8 -*-
"""
Plot Tidal Atlas Amplitude.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
import os

# Datafile loading
dloc = './Atlas'
dconst = 'M2'
units = '(m)'
dfile = '{wave:s}-elevation-atlas.nc'.format(wave=dconst)
dpath = os.path.join(dloc, dfile)
ds = Dataset(dpath)
xs = ds['lon'][:]
ys = ds['lat'][:]
pvalue = ds['elevation_G'][:]

# Triangulations
triangles = ds['element'][:]
triangles = triangles - 1
tri = mtri.Triangulation(np.asarray(xs), np.asarray(ys), triangles)

# Plotting the triangulation points
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_title('{wave:s} Amplitude {units:s}'.format(wave=dconst, units=units))
cmap = cmocean.cm.phase
dplot = ax.tricontourf(tri, pvalue, cmap=cmap)

# Colorbar hack for equal axis plot
daxis = make_axes_locatable(ax)
cax = daxis.append_axes('right', size='5%', pad=0.1)
cbar = plt.colorbar(dplot, cax=cax, extend='both')

# Showing plot
plt.show()