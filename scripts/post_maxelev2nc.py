#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reading maxelev files and putting them in a netcdf along with the bathymetry
mesh.

Dimensitons - 
    SCHISM_hgrid_node_x
    SCHISM_hgrid_node_y
    SCHISM_hgrid_element
    experiment

Variables - 
    SCHISM_hgrid_node_x
    SCHISM_hgrid_node_y
    SCHISM_hgrid_element
    experiment
    depth
    elev

@license: GPL3
@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
