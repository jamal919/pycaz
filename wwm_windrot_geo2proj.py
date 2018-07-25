# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:55:28 2018

@author: khan
"""

from schism import io
grid = io.Gr3()
grid.readfromfile(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3')

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0

# Writing the grid file
grid.writetofile(path='./windrot_geo2proj.gr3', overwrite=True, writebounds=False, nodevalfmt='%4i')