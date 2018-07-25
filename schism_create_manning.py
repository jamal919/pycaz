# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 15:04:20 2018

@author: khan
"""

from SCHISM import schismio
grid = schismio.Gr3()
grid.readfromfile(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3')

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0.02

# Writing the grid file
grid.writetofile(path='./manning.gr3', overwrite=True, writebounds=False, nodevalfmt='%4i')