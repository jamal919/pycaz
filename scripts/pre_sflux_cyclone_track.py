# -*- coding: utf-8 -*-
"""
Create wind and pressure field from cyclone track and write it into SCHISM
complient sflux file. 

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from netCDF4 import Dataset

class Track(object):
    def __init__(self, time, lon, lat, P=None, Vmax=None, Rmax=None):
        pass

class Grid(object):
    def __init__(self, x, y):
        pass

class WindInterpolator(object):
    def __init__(self, grid, method='E11'):
        pass

    def interpolate(self, at):
        pass

class PressureInterpolator(object):
    def __init__(self, grid, method='H80'):
        pass

    def interpolate(self, at):
        pass

class Sflux(object):
    def __init__(self, grid, inittime, nstep, path='./'):
        pass

    def __create_netcdf(self, filename):
        pass
    
    def putvalue(self, uwind, vwind, prmsl, stmp, spfh):
        pass