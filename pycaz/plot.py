#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plotting options for schism output. This is intended for all the outputs except
gr3. The map plotting option from gr3 is implemented in Gr3 object.

Output types - 
1. BP format
2. netcdf format
3. station format

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


class Gr3(object):
    """
    SCHISM .gr3 object. This object loads the .gr3 like file, breaks it in
    it's subcomponent, and store it for the model object to use.
     
    The intended way for this model is to provide the original file from SMS.
    """
    def __init__(self, path = None):
        # path options
        self.path = path
        if os.access(self.path, os.F_OK):
            self.read()
        else:
            print('No file is given!')
        
    def read(self): 
        """ Read gr3 format """
        with open(self.path) as f:
            self.ds = f.readlines()
            self.grname = ds[0].split()
            self.nelem, self.nnode = ds[1].split()
            self.nelem = int(self.nelem)
            self.nnode = int(self.nnode)

            self.dnodes = np.genfromtxt(fname=self.ds[2:self.nnode+2])
            self.delems = np.genfromtxt(fname=self.ds[self.nnode+2:self.nnode+self.nelem+2])
            self.triang = self.delems[:, 2:5]
            self.triang = self.triang - 1
            
    def plot(self, ):
        """ Plot the gr3 in map """

        self.fig = plt.figure(figsize=(9,6), dpi=300)
        self.ax = self.fig.add_subplot(111)
        
        self.x = dsval[:, 1]
        self.y = dsval[:, 2]
        self.val = dsval[:, 3]

        self.plotbound = [np.min(self.x), np.min(self.y), np.max(self.x), np.max(self.y)]

        self.m = Basemap(llcrnrlon=plotbound[0], llcrnrlat=plotbound[1], urcrnrlon=plotbound[2], urcrnrlat=plotbound[3], projection='merc', resolution='f', epsg='4326')
        self.m.drawcoastlines(linewidth=2)
        self.m.drawcountries()
        self.m.bluemarble()
        self.m.controurf(x=self.x, y=self.y, data=self.val, tri=True, latlon=True, triangles=self.triang)
        return(self.fig)
