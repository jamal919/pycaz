#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Necessary conversion functions

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import numpy as np

def hpa2pa(hpa):
    '''
    Takes pressure value in hecta Pascal and return in Pascal. 
    '''
    return(hpa*100)

def knot2mps(knot):
    '''
    Takes velocity in knot and returns velocity in mps.
    '''
    return(knot*1.852/3.6)

def km2m(km):
    '''
    Takes distance in Km and converts it to meter.
    '''
    return(km*1000)

def ft2m(ft):
    '''
    Takes distance in Ft and converts it to meter
    '''
    return(ft*0.3048)

def ntm2m(ntm):
    '''
    Takes distance in Nautical mile and converts to meter.
    '''
    return(ntm*1852.001)

def lon180(lon360):
    '''
    Change lon range from 0-360 to -180-180
    '''
    lon360[lon360 > 180] = lon360[lon360 > 180] - 360
    return(lon360)

def gc_distance(of, origin, isradians=False):
    '''
    Calculates the great circle distance of 'of' from 'origin'
    
    of: list of lon lat of the point
    origin: list of lon lat of the origin point
    '''
    __dfac = 60*1.852*1000
    
    if isradians:
        __dtrans_x = __dfac*np.cos(origin[1])*(np.rad2deg(of[0])-np.rad2deg(origin[0]))
        __dtrans_y = __dfac*(np.rad2deg(of[1])-np.rad2deg(origin[1]))
    else:
        __dtrans_x = __dfac*np.cos(np.deg2rad(origin[1]))*(of[0]-origin[0])
        __dtrans_y = __dfac*(of[1]-origin[1])

    return((__dtrans_x, __dtrans_y))
