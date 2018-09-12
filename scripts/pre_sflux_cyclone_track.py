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

class WindModel(object):
    def __init__(self, radialwind, isradians=False, SWRF=0.9):
        '''
        Wind model object contains the wind models and can be used to solve for
        radius, speed etc.

        args:
            radialwind (dict) : contains vmax, lat, rmax or other radial
                                information for a given time step
            isradians (bool)  : if the lat, long information is in radians
                                Default is False
            SWRF (float)      : Surface wind reduction factor to convert 10m
                                surface wind speed to boundary layer wind speed
                                where the wind models are valid.
                                Deafult is 0.9 (Chavaz and Lin 2015)
        methods:
            find_radius :   calculate the outer or inner radius at a given wind
                            speed
        '''

        self.isradians = isradians
        self.SWRF = 0.9
        self.converter = Converter()

        if isinstance(radialwind, dict):
            if 'vmax' in radialwind.keys():
                self.vmax = radialwind['vmax']
                self.vmaxb = self.__to_boundary(self.vmax)
            else:
                print('Must supply the maximum velocity (vmax)! Aborting...')
                sys.exit(1)

            if 'lat' in radialwind.keys():
                self.lat = radialwind['lat']
            else:
                print('Must supply latitude (lat)! Aborting...')
                sys.exit(1)

            if 'rmax' in radialwind.keys():
                self.rmax = radialwind['rmax']

            if 'r34' in radialwind.keys():
                self.r34 = radialwind['r34']
                self.v34 = self.converter.knot2mps(34)

            if 'r50' in radialwind.keys():
                self.r50 = radialwind['r50']

            if 'r64' in radialwind.keys():
                self.r64 = radialwind['r64']

            if 'f' in radialwind.keys():
                self.f = radialwind['f']
        else:
            print('Input a dictionary of radial wind information to Wind Model!')
            sys.exit(1)

    def __to_boundary(self, v10m):
        return(v10m/self.SWRF)

    def __calc_coriolis(self):
        ''' 
        Calculate the coriolis coefficient.
        Uses the simple function - 2*7.292e-5*sin(lat)
        '''
        if hasattr(self, 'lat'):
            # Calculate coriolis factor with latitude
            if self.isradians:
                self.f = 2*7.292e-5*np.sin(self.lat)
            else:
                self.f = 2*7.292e-5*np.sin(np.deg2rad(self.lat))
        else:
            print('Must provide lat or coriolis coefficient in the radialwind dict! Aborting...')
            sys.exit(1)

    def find_radius(self, v, model='E11', using='scan', limit=[500000, 0], step=-100, within='outer'):
        '''
        using: 'scan' or 'fsolve'
        limits: only available with using='scan'
        step: only used with using='scan' 
        r0: only used with using='fsolve'
        '''
        # Loading the essential variables
        __v = v
        __model = model

        # check if the coriolis is in the class
        if hasattr(self, 'f'):
            # f is spplied with the radial wind information
            pass
        else:
            self.__calc_coriolis()

        # Cheking the model options and set the function to solve
        if __model=='H80':
            pass
        elif __model=='E11':
            if hasattr(self, 'rmax'):
                __resfunc = lambda __r: __v - (2*__r*(self.vmax*self.rmax + 0.5*self.f*self.rmax**2)/\
                        (self.rmax**2+__r**2)-self.f*__r/2)
            else:
                print('Must provide lat or coriolis coefficient in the radialwind dict! Aborting...')
                sys.exit(1)
        else:
            print('Model {:s} not found'.format(__model))
            sys.exit(1)

        # solving for given v in the outer radius
        if within=='outer':
            if using=='vector':
                __rrange = np.arange(start=limit[0], stop=limit[1], step=step)
                __res = np.array([__resfunc(i) for i in __rrange])
                __loc = np.where(__res < 0)[0]
                __rsolved = __rrange[__loc[0]]
                return(__rsolved)
            elif using=='scan':
                __rrange = np.arange(start=limit[0], stop=limit[1], step=step)
                for i in __rrange:
                    __res = __resfunc(i)
                    if __res < 0:
                        __rsolved = i
                        break
                return(__rsolved)
            elif using=='fsolve':
                __rsolved = optimize.fsolve(__resfunc, x0=self.rmax*1.5)
                return(__rsolved[0])
            elif using=='bisect':
                __rsolved = optimize.bisect(__resfunc, a=limit[0], b=self.rmax)
                return(__rsolved)
            else:
                print('Method {:s} not found! Aborting...'.format(using))
        elif within=='inner':
            print('Radius calculation in the inner structure is not available yet!')
            sys.exit(1)


class PressureModel(object):
    def __init__(self):
        pass

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
        __lon = self.converter.lon180(lon360=__values[:, 6]) # -180 to 180 format
        __lat = __values[:, 7]
        __p = self.converter.hpa2pa(hpa=__values[:, 8]) # hPa to Pa conversion
        __rm = self.converter.km2m(km=__values[:, 9]) # Km to m conversion
        __vmax = self.converter.knot2mps(knot=__values[:, 10]) # knot to mps conversion
        __track = dict(timeindex=__time,\
                    time=__datetime,\
                    lon=__lon,\
                    lat=__lat,\
                    p=__p,\
                    rm=__rm,\
                    vmax=__vmax)
        return(__track)

class Track(object):
    def __init__(self, track, clipby=None):
        self.track = track
        self.starttime = track['time'][0]
        
        # calculate translation speed
        self.__translation_speed()
        print('Translation speed is calculated.')

        if clipby is None:
            # No need to clip the track
            pass
        else:
            self.clip(by=clipby)

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

    def clip(self, by):
        pass

    def __translation_speed(self):
        pass

    def __at_timedelta(self):
        pass
    
    def __at_datetime(self):
        pass

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

        __X, __Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = __X.shape
        self.length = len(__X.flat)

        self.points = np.array([Point(x=__X.flat[i], y=__Y.flat[i]) for i in np.arange(self.length)])
        self.points = np.reshape(self.points, self.shape)

    def distance_from(self, x, y):
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
    
    def __putvalue(self, stepi, uwind, vwind, prmsl, stmp, spfh):
        pass

def perf():
    windmodel = WindModel(radialwind=dict(rmax=20000, vmax=40, lat=24.5))
    windmodel.find_radius(v=50*1.852/3.6, using='bisect')


if __name__=='__main__':
    import timeit
    setup = 'from __main__ import perf'
    print(timeit.timeit('perf()', setup=setup))
    # reader = Reader(readername='kerry')
    # trackfile = reader.read('/run/media/khan/Workbench/Projects/Surge Model/Emmanuel et al/tracks_csv/Track_0001.csv')
    # track = Track(track=trackfile)
    # for i in range(1, 15):
    #     track.interpolate(var='lon', at=timedelta(minutes=15*i))
    # print([time.gmtime(i) for i in track['timeindex']])
    # print(track.keys())
    # print(isinstance(track['timeindex'][0], int))w


    # windmodel = WindModel(radialwind=dict(rmax=20000, vmax=40, lat=24.5))
    # print(windmodel.find_radius(v=50*1.852/3.6, using='scan'))
    # print(windmodel.find_radius(v=50*1.852/3.6, using='fsolve'))
    # print(windmodel.find_radius(v=50*1.852/3.6, using='bisect'))
