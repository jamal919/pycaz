#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Necessary conversion functions

@author: khan
"""
import numpy as np

def hpa2pa(hpa:float) -> float:
    '''
    Takes pressure value in hecta Pascal and return in Pascal. 
    '''
    return(hpa*100)

def pa2mb(pa:float) -> float:
    '''
    Takes pressure as Pa an return in milibar (hPa)
    '''
    return(pa/100.0)

def knot2mps(knot:float) -> float:
    '''
    Takes velocity in knot and returns velocity in mps.
    '''
    return(knot*1.852/3.6)

def km2m(km:float) -> float:
    '''
    Takes distance in Km and converts it to meter.
    '''
    return(km*1000)

def ft2m(ft:float) -> float:
    '''
    Takes distance in Ft and converts it to meter
    '''
    return(ft*0.3048)

def ntm2m(ntm:float) -> float:
    '''
    Takes distance in Nautical mile and converts to meter.
    '''
    return(ntm*1852.001)

def lon180(lon360:float) -> float:
    '''
    Change lon range from 0-360 to -180-180
    '''
    lon360[lon360 > 180] = lon360[lon360 > 180] - 360
    return(lon360)

def gc_distance(of_x:float, of_y:float, origin_x:float, origin_y:float, isradians:bool=False):
    '''
    Calculates the great circle distance of 'of' from 'origin'
    
    of: list of lon lat of the point
    origin: list of lon lat of the origin point
    '''
    dfac = 60*1.852*1000
    
    if isradians:
        dtrans_x = dfac*np.cos(origin_y)*(np.rad2deg(of_x)-np.rad2deg(origin_x))
        dtrans_y = dfac*(np.rad2deg(of_y)-np.rad2deg(origin_y))
    else:
        dtrans_x = dfac*np.cos(np.deg2rad(origin_y))*(of_x-origin_x)
        dtrans_y = dfac*(of_y-origin_y)

    return((dtrans_x, dtrans_y))
