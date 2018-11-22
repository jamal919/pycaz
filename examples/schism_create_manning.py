# -*- coding: utf-8 -*-
"""
Demo script to create mannings.gr3 file for SCHISM.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import sys
sys.path.append('/home/khan/MEGA/Codes/SCHISMMB')

import os
from schism import io as schismio

prjpath = '/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/'
grid = schismio.Gr3()
grid.readfromfile(path=os.path.join(prjpath, 'hgrid.gr3'))

# First setting all points to 0 (not on boundary)
grid.dnodes[:, 3] = 0.02

# Writing the grid file
grid.writetofile(path=os.path.join(prjpath, 'manning.gr3'), overwrite=True, writebounds=False, nodevalfmt='%4i')