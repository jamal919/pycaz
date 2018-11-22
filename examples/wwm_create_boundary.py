# -*- coding: utf-8 -*-
"""
Sample script to create the WWM boundary for SCHISM using schism toolbox.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import schism

grid = schism.io.Gr3()
grid.readfromfile(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3')

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0

# All boundaries are exterior boundary (1)
for boundary in grid.openbnd.boundaries:
    grid.dnodes[[i - 1 for i in boundary.nodes], 3] = 1
    
for boundary in grid.landbnd.boundaries:
    grid.dnodes[[i - 1 for i in boundary.nodes], 3] = 1

# Writing the grid file
grid.writetofile(path='./wwmbnd.gr3', overwrite=True, nodevalfmt='%4i')