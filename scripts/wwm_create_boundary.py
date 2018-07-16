# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 15:04:20 2018

@author: khan
"""

import schismio
grid = schismio.Gr3(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3', readflag=(True, True, True))
nodes = grid.readnodes()
elems= grid.readelems()
openbnd, landbnd = grid.readbounds()
