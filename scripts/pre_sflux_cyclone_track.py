# -*- coding: utf-8 -*-
"""
Create wind and pressure field from cyclone track and write it into SCHISM
complient sflux file. 

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import sys
import numpy as np
from scipy import optimize
from datetime import datetime, timedelta
import calendar
import time
from netCDF4 import Dataset

class Converter(object):
    def __init__(self):
        pass

    @staticmethod
    def timestamp(d):
        return(calendar.timegm(d.timetuple()))

    @staticmethod
    def hpa2pa(hpa):
        return(hpa*100)
    
    @staticmethod
    def knot2mps(knot):
        return(knot*1.852/3.6)

    @staticmethod
    def km2m(km):
        return(km*1000)

    @staticmethod
    def lon180(lon360):
        lon360[lon360 > 180] = lon360[lon360 > 180] - 360
        return(lon360)

class Reader(object):
    '''
    Reader(readername) create a reader object for cyclone track file.
    '''
    def __init__(self, readername):
        self.readername = readername
        self.converter = Converter()

    def read(self, path):
        self.path = path
        self.track = getattr(self, self.readername)()
        return(self.track)

    def kerry(self):
        # Reading the track file
        __values = np.genfromtxt(fname=self.path, dtype=float, delimiter=',', skiprows=3)
        
        # Creating the datetime object
        __datetime = [datetime(int(__values[i, 0]),\
                            int(__values[i, 1]),\
                            int(__values[i, 2]),\
                            int(__values[i, 3]),\
                            int(__values[i, 4]),\
                            int(__values[i, 5])) \
                            for i in np.arange(len(__values))]
        # Creating the datetime index
        __time = np.array([self.converter.timestamp(__datetime[i]) \
                        for i in np.arange(len(__datetime))])
        # Processing other variables
        __lon = self.converter.lon180(lon360=__values[:, 6]) # -180 to 180
        __lat = __values[:, 7]
        __p = self.converter.hpa2pa(hpa=__values[:, 8]) # hPa to Pa
        __rm = self.converter.km2m(km=__values[:, 9]) # Km to m
        __vmax = self.converter.knot2mps(knot=__values[:, 10]) # knot to mps
        __track = np.array([dict(timeindex=__time[i],\
                    time=__datetime[i],\
                    lon=__lon[i],\
                    lat=__lat[i],\
                    p=__p[i],\
                    rm=__rm[i],\
                    vmax=__vmax[i]) for i in np.arange(len(__time))])
        return(__track)

class Track(object):
    def __init__(self, track, clipby=None):
        self.track = track
        self.starttime = track[0]['time']

        # Calculate translation speed
        self.__calc_translation()
        
        # Clipping the track
        if clipby is None:
            print('No clipping is done!')
        else:
            self.__clip(by=clipby)

    def __calc_translation(self):
        __dfac = 60*1.852*1000
        
        # Calculating __utrans and __vtrans for all timestep except the last
        for i in np.arange(0, len(self.track)-1):
            __dt = self.track[i+1]['time'] - self.track[i]['time']
            __dt = __dt.total_seconds()
            __lon1 = np.deg2rad(self.track[i]['lon'])
            __lon2 = np.deg2rad(self.track[i+1]['lon'])
            __lat1 = np.deg2rad(self.track[i]['lat'])
            __lat2 = np.deg2rad(self.track[i+1]['lat'])
            __dtrans_x = __dfac*np.cos(__lat1)*(__lon2-__lon1)
            __dtrans_y = __dfac*(__lat2-__lat1)
            __utstorm = __dtrans_x/__dt
            __vtstorm = __dtrans_y/__dt
            __tstorm = np.sqrt(__utstorm**2+__vtstorm**2)
            self.track[i]['utstorm'] = __utstorm
            self.track[i]['vtstorm'] = __vtstorm

        # For the last time step, we are keeping it to the same

    def __clip(self, by):
        '''
        Clip the track to a given boundary
        TODO: Check if the track time is OK
        '''
        if isinstance(by, list) and len(by)==4:
            # by lonlatbox cropping
            __lonmin = by[0]
            __lonmax = by[1]
            __latmin = by[2]
            __latmax = by[3]

            __lonlist = np.array([self.track[i]['lon'] for i in np.arange(len(self.track))])
            print(__lonlist)
            __lonselected = np.where((__lonlist >= __lonmin) & (__lonlist <= __lonmax))[0]
            print(__lonselected)
            __mod_track = self.track[__lonselected]
            __latlist = np.array([__mod_track[i]['lat'] for i in np.arange(len(__mod_track))])
            print(__latlist)
            __latselected = np.where((__latlist >= __latmin) & (__latlist <= __latmax))
            print(__latselected)
            __mod_track = __mod_track[__latselected]

            self.track = __mod_track

        elif isinstance(by, Grid):
            #TODO Implement clipping by grid selection
            print('Clipping by grid not yet implemented')
            sys.exit(1)
        else:
            print('Track clipping failed! Give a lonlatbox or Grid input. Aborting...')
            sys.exit(1)

    def interpolate(self, var, at):
        __var = var
        if isinstance(at, int):
            # already in timestamp format (seconds from epoch)
            __at = at
        elif isinstance(at, datetime):
            __at = calendar.timegm(at.timetuple())
        elif isinstance(at, timedelta):
            __at = self.starttime + at
            __at = calendar.timegm(__at.timetuple())
        else:
            print('Interpolation time must be in datetime/timedelta/timegm format')

        print(time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(__at)))
        
        if __var in self.track.keys():
            if __var=='lat' or __var=='lon':
                pass
            elif __var=='vtrans' or __var=='utrans':
                pass
            else:
                pass
        else:
            print('Variable {:s} not found in the track'.format(__var))

class Point(object):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def print(self):
        print(self.x, self.y)

    def __lt__(self, other):
        return(self.x < other.x and self.y < other.y)

class Grid(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = self.X.shape
        self.length = len(self.X.flat)

        self.points = np.array([Point(x=self.X.flat[i], y=self.Y.flat[i]) for i in np.arange(self.length)])
        self.points = np.reshape(self.points, self.shape)

    def distance_from(self, x, y):
        pass


class Interpolator(object):
    def __init__(self, grid):
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
    
    def __putvalue(self, stepi, uwind, vwind, prmsl, stmp, spfh):
        pass

if __name__=='__main__':
    reader = Reader(readername='kerry')
    trackfile = reader.read('/run/media/khan/Workbench/Projects/Surge Model/Emmanuel et al/tracks_csv/Track_0001.csv')
    track = Track(track=trackfile, clipby=[79, 99, 10.5, 24.5])
    print('clipped : ', len(track.track))
    track = Track(track=trackfile)
    print('unclipped : ', len(track.track))