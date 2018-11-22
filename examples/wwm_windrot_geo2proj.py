# -*- coding: utf-8 -*-
"""
Sample script to create the windrot_geo2proj file for SCHISM using schism toolbox.

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

# Writing the grid file
grid.writetofile(path=os.path.join(prjpath, 'windrot_geo2proj.gr3'), overwrite=True, writebounds=False, nodevalfmt='%4i')