# -*- coding: utf-8 -*-
"""
Plotting gr3 format file in a Map Background.

This scripts reads a gr3 format file and plot it using maptplotlib.

Author: Jamal Uddin Khan
Date : 17-05-2017
"""
import timeit
starttime = timeit.default_timer()

import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# Plotting for various settings
fig = plt.figure(figsize=(9,6), dpi=300)

graticulespace = [0.3, 0.2]
plotbound=[89.2, 90.7, 21.4, 22.2]

# Plot 01
indir = "./tide_surge_yann_561_implicit"
ax = fig.add_subplot(221)
ax.set_title("SCHISM_v561_Implicit")
# Reading gr3
with open(os.path.join(indir, 'maxelev.gr3')) as f:
    ds = f.readlines()
    grname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)

    dsval = np.genfromtxt(fname=ds[2:nnode+2])
    dselm = np.genfromtxt(fname=ds[nnode+2:nnode+nelem+2])

    # Reading mesh file
with open(os.path.join(indir, 'hgrid.gr3')) as f:
    ds = f.readlines()
    mshname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)
    
    nodes = np.genfromtxt(fname=ds[2:nnode+2])
    
    # Setting up the variables
    plotx = dsval[:, 1]
    ploty = dsval[:, 2]
    gr3val = dsval[:, 3]
    depthval = nodes[:, 3]
        
    plotval = (gr3val)*(depthval > 0) + (gr3val + depthval)*(depthval <= 0)

# Triangulation
triang = dselm[:, 2:5]
triang = triang-1

# Map plotting with builtin coastline at 'full' resolution
m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1],\
            urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg=4326)
m.drawcoastlines(linewidth=2)
#m.fillcontinents(zorder=0)
m.drawcountries()
#m.bluemarble()
m.arcgisimage(service='ESRI_Imagery_World_2D')
m.drawparallels(circles=np.arange(plotbound[2], plotbound[3], graticulespace[1]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
m.drawmeridians(meridians=np.arange(plotbound[0], plotbound[1], graticulespace[0]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
clevs = [0, 1, 2, 3, 4, 5, 6]
cs = m.contourf(x=plotx, y=ploty, data=plotval, levels=clevs,\
                tri=True, latlon=True, triangles=triang)

# Plot 2
ax = fig.add_subplot(222)
ax.set_title("SCHISM_v561_Explicit")
indir="./tide_surge_yann_561_explicit"
# Reading gr3
with open(os.path.join(indir, 'maxelev.gr3')) as f:
    ds = f.readlines()
    grname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)

    dsval = np.genfromtxt(fname=ds[2:nnode+2])
    dselm = np.genfromtxt(fname=ds[nnode+2:nnode+nelem+2])

# Reading mesh file
with open(os.path.join(indir, 'hgrid.gr3')) as f:
    ds = f.readlines()
    mshname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)
    
    nodes = np.genfromtxt(fname=ds[2:nnode+2])
    
    # Setting up the variables
    plotx = dsval[:, 1]
    ploty = dsval[:, 2]
    gr3val = dsval[:, 3]
    depthval = nodes[:, 3]
    
    plotval = (gr3val)*(depthval > 0) + (gr3val + depthval)*(depthval <= 0)

# Triangulation
triang = dselm[:, 2:5]
triang = triang-1

# Map plotting with builtin coastline at 'full' resolution
m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1],\
            urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg=4326)
m.drawcoastlines(linewidth=2)
#m.fillcontinents(zorder=0)
m.drawcountries()
#m.bluemarble()
m.arcgisimage(service='ESRI_Imagery_World_2D')
m.drawparallels(circles=np.arange(plotbound[2], plotbound[3], graticulespace[1]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
m.drawmeridians(meridians=np.arange(plotbound[0], plotbound[1], graticulespace[0]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
clevs = [0, 1, 2, 3, 4, 5, 6]
cs = m.contourf(x=plotx, y=ploty, data=plotval, levels=clevs,\
                tri=True, latlon=True, triangles=triang)

ax = fig.add_subplot(223)
ax.set_title("SCHISM_v561_Implicit_All_600")
indir="./tide_surge_yann_561_all_implicit_600sec"
# Reading gr3
with open(os.path.join(indir, 'maxelev.gr3')) as f:
    ds = f.readlines()
    grname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)

    dsval = np.genfromtxt(fname=ds[2:nnode+2])
    dselm = np.genfromtxt(fname=ds[nnode+2:nnode+nelem+2])

# Reading mesh file
with open(os.path.join(indir, 'hgrid.gr3')) as f:
    ds = f.readlines()
    mshname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)
    
    nodes = np.genfromtxt(fname=ds[2:nnode+2])
    
    # Setting up the variables
    plotx = dsval[:, 1]
    ploty = dsval[:, 2]
    gr3val = dsval[:, 3]
    depthval = nodes[:, 3]
    
    plotval = (gr3val)*(depthval > 0) + (gr3val + depthval)*(depthval <= 0)

# Triangulation
triang = dselm[:, 2:5]
triang = triang-1

# Map plotting with builtin coastline at 'full' resolution
m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1],\
            urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg=4326)
m.drawcoastlines(linewidth=2)
#m.fillcontinents(zorder=0)
m.drawcountries()
#m.bluemarble()
m.arcgisimage(service='ESRI_Imagery_World_2D')
m.drawparallels(circles=np.arange(plotbound[2], plotbound[3], graticulespace[1]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
m.drawmeridians(meridians=np.arange(plotbound[0], plotbound[1], graticulespace[0]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
clevs = [0, 1, 2, 3, 4, 5, 6]
cs = m.contourf(x=plotx, y=ploty, data=plotval, levels=clevs,\
                tri=True, latlon=True, triangles=triang)

ax = fig.add_subplot(224)
ax.set_title("SCHISM_v561_Hydro")
indir="./tide_surge_yann_561_hydro"
# Reading gr3
with open(os.path.join(indir, 'maxelev.gr3')) as f:
    ds = f.readlines()
    grname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)

    dsval = np.genfromtxt(fname=ds[2:nnode+2])
    dselm = np.genfromtxt(fname=ds[nnode+2:nnode+nelem+2])

# Reading mesh file
with open(os.path.join(indir, 'hgrid.gr3')) as f:
    ds = f.readlines()
    mshname = ds[0].split()
    nelem, nnode = ds[1].split()
    nelem = int(nelem)
    nnode = int(nnode)
    
    nodes = np.genfromtxt(fname=ds[2:nnode+2])
    
    # Setting up the variables
    plotx = dsval[:, 1]
    ploty = dsval[:, 2]
    gr3val = dsval[:, 3]
    depthval = nodes[:, 3]
    
    plotval = (gr3val)*(depthval > 0) + (gr3val + depthval)*(depthval <= 0)

# Triangulation
triang = dselm[:, 2:5]
triang = triang-1

# Map plotting with builtin coastline at 'full' resolution
m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1],\
            urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg=4326)
m.drawcoastlines(linewidth=2)
#m.fillcontinents(zorder=0)
m.drawcountries()
#m.bluemarble()
m.arcgisimage(service='ESRI_Imagery_World_2D')
m.drawparallels(circles=np.arange(plotbound[2], plotbound[3], graticulespace[1]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
m.drawmeridians(meridians=np.arange(plotbound[0], plotbound[1], graticulespace[0]),\
                labels=[True, False, False, True], dashes=[2, 2], fmt='%.01f')
clevs = [0, 1, 2, 3, 4, 5, 6]
cs = m.contourf(x=plotx, y=ploty, data=plotval, levels=clevs,\
                tri=True, latlon=True, triangles=triang)


# Colorbar axis
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(cs, cax=cbar_ax)

# Saving file
# plt.tight_layout()
plt.savefig('maxelev_traiangulation.png')

stoptime = timeit.default_timer()

print(stoptime - starttime)