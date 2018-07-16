# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 15:04:20 2018

@author: khan
"""

from SCHISM import schismio
grid = schismio.Gr3()
grid.initfromfile(path='/run/media/khan/Workbench/Educations/SCHISM/Simulations/CycloneMora/Tide_Wind_Pressure/hgrid.gr3')
nodes = grid.readnodes()
elems= grid.readelems()
openbnd, landbnd = grid.readbounds()
