#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create the bctides.in file for SCHISM with proper convention.

bctides.in is one of the must input files for SCHISM modelling system. It contains
the information for the tidal potential inside the domain of the grid as well as
the boundary condition and nudging at various open boundaries.

The open boundary can take many forms and documented as a table in the SCHISM 
manual. For practical reason, the following boundary conditions are currently 
implemented.

|--------------+------------------+-------------------------------+-------------------------|
| Variable     | eta              | Salt,Temp,Tracers             | u,v                     |
|--------------+------------------+-------------------------------+-------------------------|
| Type 1       | elev.th          | [MOD]_[ID].th                 | flux.th                 |
|              | (time history)   | (relax to time history)       | (discharge)             |
|              | (uniform at bnd) | (uniform at bnd for inflow)   | (<0 for inflow)         |
|--------------+------------------+-------------------------------+-------------------------|
| Type 2       | constant         | Relaxt to constant for inflow | discharge               |
|              |                  |                               | (<0 for inflow)         |
|--------------+------------------+-------------------------------+-------------------------|
| Type 3       | Tidal A/P        | Relax to i.c. for inflow      | Tides                   |
|              |                  |                               | (uniform A/P along bnd) |
|--------------+------------------+-------------------------------+-------------------------|
| Type -1      | Must=0           | N/A                           | Flather                 |
|--------------+------------------+-------------------------------+-------------------------|
| Nudging      | inu_elev=1       | inu_[MOD]=1 or 2              | inu_uv=1                |
| Sponge Layer |                  |                               |                         |
|--------------+------------------+-------------------------------+-------------------------|

The approach here is to set-up various types of boundary interms of boundary type,
i.e., eta, s/t, u/v etc.

The general format of the header line of the boundary is following - 
nnodes, elev, velocity, temperature, salinity ... and other modules if needed.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
class Node(object):
    pass

class Element(object):
    pass

class Mesh(object):
    pass

class Boundary(object):
    pass

class Hgrid(object):
    pass

class BC(Boundary):
    pass