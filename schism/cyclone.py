# -*- coding: utf-8 -*-
"""
Cyclone wind field models. 

TODO:
    * Implement Holland 1980
    * Implement Emmanuel and Rotanno 2011

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
from scipy import optimize
import sys

class Holland1980(object):
    def __init__(self, radialwind, isradians=False, SWRF=0.9, rhoair=1.15, pn=101300):
        '''
        Wind model object contains the wind models and can be used to solve for
        radius, speed etc. Except otherwise stated, metric system is used as
        unit convention.

        args:
            radialwind (dict)   : contains vmax (@10m), lat, rmax or other radial
                                information for a given time step
            isradians (bool)    : if the lat, long information is in radians
                                Default is False
            SWRF (float)        : Surface wind reduction factor to convert 10m
                                surface wind speed to boundary layer wind speed
                                where the wind models are valid.
                                Deafult is 0.9 (Chavaz and Lin 2015)
            rhoair (float)      : Density of air. Default is 1.15 km/m**3
            pn (float)          : Ambient pressure. Default is 101300 Pa
        
        methods:
            find_radius :   calculate the outer or inner radius at a given wind
                            speed
        '''

        self.isradians = isradians
        self.SWRF = 0.9
        self.rhoair = rhoair
        self.pn = float(pn)

        if isinstance(radialwind, dict):
            if 'vmax' in radialwind.keys():
                self.vmax = radialwind['vmax']
                self.vmax = self.__to_boundary(self.vmax)
            else:
                print('Must supply the maximum velocity (vmax)! Aborting...')
                sys.exit(1)

            if 'lat' in radialwind.keys():
                self.lat = radialwind['lat']
            else:
                print('Must supply lat of the point! Aborting...')
                sys.exit(1)

            if 'rmax' in radialwind.keys():
                self.rmax = radialwind['rmax']

            if 'r34' in radialwind.keys():
                self.r34 = radialwind['r34']

            if 'r50' in radialwind.keys():
                self.r50 = radialwind['r50']

            if 'r64' in radialwind.keys():
                self.r64 = radialwind['r64']

            if 'p' in radialwind.keys():
                self.p = radialwind['p']

            if 'f' in radialwind.keys():
                self.f = radialwind['f']
        else:
            print('Input a dictionary of radial wind information to Wind Model!')
            sys.exit(1)

    def __holland_B(self, p, bmax=2.5, bmin=0.5):
        '''
        Calculate Holland B parameter for a given pressure p and ambient pressure
        defined in the class definition.
        '''
        __B = self.vmax**2*self.rhoair*np.exp(1)/(self.pn-p)

    def __to_boundary(self, v10m):
        '''
        Convert velocity at 10 to velocity at the boundary layer.
        '''
        return(v10m/self.SWRF)

    def __to_10m(self, vbound):
        '''
        Convert velocity at boundary layer to velocity at 10m.
        '''
        return(vbound*self.SWRF)

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

    

    def find_radius(self, v, at='boundary', model='E11', using='scan', limit=[500000, 0], step=-100, within='outer'):
        '''
        find_radius(v) returns the radius corresponding to the given v. It has
        several options to control the model, solving methods.
        
        args:
            v (float)   :   Velocity at which the radius is to be found. As the
                            models are generally defined at the boundary layer
                            it is expected this value is in the boundary. For
                            surface speed input change the at parameter.
                            In case v is higher than vmax, the function will
                            show an error and sys.exit(1) will be called.
            at (str)    :   Location of the specified velocity. If at='surface'
                            it will be factor up to estimate the boundary velocity
            model(str)  :   Cyclone model to be used. Current options are -
                                * 'E11' : Emanuel and Rotanno (2011)
                                * 'H80' : Holland (1980)
            using(str)  :   solver method to find the radius. Currently there are
                            for options -
                                * 'vector'  : find the minimum residual in vector
                                * 'scan'    : find the minimum residual in for loop
                                * 'fsolve'  : using scipy.optimize.fsolve
                                * 'bisect'  : using scipy.optimize.bisect
            limit(list) :   Two item list of [max, min] of search limit. Used for
                            'vector', 'scan', 'bisect'
        '''
        # Checking the velocity input
        if at=='boundary':
            __v = v
        elif at=='surface':
            __v = self.__to_boundary(v)
        else:
            print('The velocity can be at=boundary (default) or at=surface! Default is used.')
            __v = v

        if __v > self.vmax:
            print('Input velocity can not be more than vmax! Aborting...')
            sys.exit(1)
        
        # Loading the model name
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
                #TODO: calculate rmax from other radial wind info
                print('Must provide rmax in the radialwind dict! Aborting...')
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
            #TODO: implement inner structure
            print('Radius calculation in the inner structure is yet to be implemented!')
            sys.exit(1)

class Emnauel2011(object):
    def __init__(self):
        pass

if __name__=='__main__':
    import timeit
    setup = 'from __main__ import perf'
    print(timeit.timeit('perf()', setup=setup))