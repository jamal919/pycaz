# -*- coding: utf-8 -*-
"""
Sample script to create the windrot_geo2proj file for SCHISM using schism toolbox.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import schism

grid = schism.io.Gr3()
grid.readfromfile(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3')

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0

# Writing the grid file
grid.writetofile(path='./windrot_geo2proj.gr3', overwrite=True, writebounds=False, nodevalfmt='%4i')