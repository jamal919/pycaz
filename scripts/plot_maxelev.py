# -*- coding: utf-8 -*-
"""
Plotting gr3 format file in a Map Background.

This scripts reads a gr3 format file and plot it using maptplotlib.

Author: Jamal Uddin Khan
Date : 17-05-2017
"""
import timeit
starttime = timeit.default_timer()

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

# Reading .gr3 file
with open('maxelev.gr3') as f:
    ds = f.readlines()
    grname = ds[0]
    nelem, nnode = ds[1].split()[0:2]
    nelem = int(nelem)
    nnode = int(nnode)
    
    dsval = np.genfromtxt(fname=ds[2:nnode+2])
    dselm = np.genfromtxt(fname=ds[nnode+2:nnode+nelem+2])

# Reading mesh file
with open('hgrid.gr3') as f:
    ds = f.readlines()
    mshname = ds[0].split()
    nelem, nnode = ds[1].split()[0:2]
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
plotbound = [89.2, 90.7, 21.4, 22.2]
plotbound = [87, 93, 20, 24] # Whole delta
graticulespace = [0.5, 0.5]
plt.figure(figsize=(11,8), dpi=300)
m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[2], urcrnrlon=plotbound[1],\
            urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg=4326)
m.drawcoastlines(linewidth=2)
#m.fillcontinents(zorder=0)
m.drawcountries()
#m.bluemarble()
m.arcgisimage(service='ESRI_Imagery_World_2D')
m.drawparallels(circles=np.arange(plotbound[2], plotbound[3], graticulespace[1]),\
                labels=[True, False, False, True], dashes=[2, 2])
m.drawmeridians(meridians=np.arange(plotbound[0], plotbound[1], graticulespace[0]),\
                labels=[True, False, False, True], dashes=[2, 2])
clevs = [0, 1, 2, 3, 4, 5, 6]
cs = m.contourf(x=plotx, y=ploty, data=plotval, levels=clevs,\
                tri=True, latlon=True, triangles=triang)
cbar = m.colorbar(cs, location='right')
cbar.ax.set_title('(m)')
plt.gca().set_title('Maximum WL for Sidr')
plt.savefig('maxelev_traiangulation.png')

stoptime = timeit.default_timer()
print(stoptime - starttime)