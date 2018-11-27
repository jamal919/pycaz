#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sample script to create the WWM boundary for SCHISM using schism toolbox.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import sys
sys.path.append('/home/khan/MEGA/Codes/SCHISMMB')

import os
from schism import io as schismio

prjpath = '/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure'
grid = schismio.Gr3()
grid.readfromfile(path=os.path.join(prjpath, 'hgrid.gr3'))

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0

# All boundaries are exterior boundary (1)
for boundary in grid.openbnd.boundaries:
    grid.dnodes[[i - 1 for i in boundary.nodes], 3] = 1
    
for boundary in grid.landbnd.boundaries:
    grid.dnodes[[i - 1 for i in boundary.nodes], 3] = 1

# Writing the grid file
grid.writetofile(path=os.path.join(prjpath, 'wwmbnd.gr3'), overwrite=True, nodevalfmt='%4i')