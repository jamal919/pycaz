# -*- coding: utf-8 -*-
"""
Plot or extract values at selected nodes from schout file. Currently plotting 
is only implemented.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
import os

# Datafile loading
dloc = './'
dfile = 'schout.nc'
dpath = os.path.join(dloc, dfile)
ds = Dataset(dpath)
xs = ds['SCHISM_hgrid_node_x'][:]
ys = ds['SCHISM_hgrid_node_y'][:]
depth = ds['depth'][:]
elev = ds['elev']
time = ds['time'][:]
time = num2date(times=time, units=ds.variables['time'].units, calendar='standard')
triangles = ds['SCHISM_hgrid_face_nodes'][:]

# Triangulations
triangles = triangles[:, 0:3]
triangles = triangles - 1
tri = mtri.Triangulation(np.asarray(xs), np.asarray(ys), triangles)

# Plotting the triangulation points
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_title('Click on the point to extract')
cmap = cmocean.cm.deep
dplot = ax.tricontourf(tri, depth, cmap=cmap)

# Colorbar hack for equal axis plot
daxis = make_axes_locatable(ax)
cax = daxis.append_axes('right', size='5%', pad=0.1)
cbar = plt.colorbar(dplot, cax=cax, extend='both')

# 2 points tolerance line plotting for future selection
dpoints = ax.triplot(tri, 'ko-', linewidth=0.3, markersize=0.3)
line, = ax.plot(xs, ys, 'ko', markersize=0.3, picker=2)
ax.set_xlabel('Lon (degrees)')
ax.set_ylabel('Lat (degrees)')

# Data holder
selpoints = []

# Function definitions
def onpick(event):
    if event.artist != line:
        print('Not same type data')
        return True
    
    N = len(event.ind)
    if not N:
        print('No points found')
        return True

    figi = plt.figure()
    for num, ind in enumerate(event.ind):
        # Add point to the map
        ax.plot(xs[ind], ys[ind], '-ro')
        fig.canvas.draw()
        print(ind)
        selpoints.append([xs[ind], ys[ind]])
        
        # Plot the time series
        axi = figi.add_subplot(N, 1, num+1)
        axi.plot(time, elev[:, ind])
        
    figi.show()
    return True

# Connecting function
fig.canvas.mpl_connect('pick_event', onpick)

# Showing plot
plt.show()