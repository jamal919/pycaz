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

class Eelement(object):
    pass

class Mesh(object):
    pass

class Boundary(object):
    pass

class Hgrid(object):
    pass

class Point(object):
    '''
    Point(x, y=0, a=None, p=None, isradians=False) is point object to hold point information
    including the amplitude and phase.

    args:
        x (float)       :   x position
        y (float)       :   y position
        a (float)       :   amplitude
        p (float)       :   phase in degree or radians (default: degrees)
        isradians (bool):   if input is in degrees or radians (default: false)

    returns:
        An instance of Point class.

    attributes:
        x (float)       : x position
        y (float)       : y position
        a (float)       : amplitude
        p (float)       : phase in radians

    methods:
        print() : prints the attributes

    TODO:
        * add typecheck and error handling
    '''
    def __init__(self, x, y=0, a=None, p=None, isradians=False):
        self.x = float(x)
        self.y = float(y)
        if a is not None:
            self.a = float(a)
        else:
            self.a = float(0)
        
        self.isradians = isradians
        
        if p is not None:
            if self.isradians:
                self.p = float(p)
            else:
                self.p = float(p)*np.pi/180.0
        else:
            self.p = float(0)

    def print(self):
        print(self.x, self.y, self.a, self.p)

    def __lt__(self, other):
        return(self.x < other.x and self.y < other.y)

class Grid(object):
    '''
    Grid(x, y, A=None, P=None, isradians=False) is the grid object to hold points in a meshgrid.

    args:
        x ([float])         : x positon array in the structured grid
        y ([float])         : y position array in the structured grid

        A ([[float]])       : 2D array of size (x, y) containing amplitude
        P ([[float]])       : 2D array of size (x, y) containing phase
        isradians (bool)    : if the phase in in radians

    '''
    def __init__(self, x, y, A=None, P=None, isradians=False):
        # Initiating variables
        self.x = x
        self.y = y
        __X, __Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = __X.shape
        self.length = len(__X.flat)
        self.isradians = isradians

        if A is None:
            __A = np.zeros(shape=self.shape)
        else:
            __A = A
        
        if P is None:
            __P = np.zeros(shape=self.shape)
        else:
            __P = P

        # Creating point meshgrid
        self.points = np.array([Point(x=__X.flat[i], y=__Y.flat[i], a=__A.flat[i], p=__P.flat[i], isradians=self.isradians) for i in np.arange(self.length)])
        self.points = np.reshape(self.points, self.shape)

    def getpoints(self, reshaped=True):
        __points = np.array([point for point in self.points.flat])
        if reshaped:
            return(np.reshape(self.points, self.shape))
        else:
            return(self.points)

    def getx(self, reshaped=True):
        __X = np.array([point.x for point in self.points.flat])
        if reshaped:
            __X = np.reshape(__X, self.shape)
        return(__X)

    def gety(self, reshaped=True):
        __Y = np.array([point.y for point in self.points.flat])
        if reshaped:
            __Y = np.reshape(__Y, self.shape)
        return(__Y)

    def getamplitude(self, reshaped=True):
        __A = np.array([point.a for point in self.points.flat])
        if reshaped:
            __A = np.reshape(__A, self.shape)
        return(__A)

    def getphase(self, reshaped=True, degrees=False):
        __P = np.array([point.p for point in self.points.flat])
        if degrees:
            __P = __P*180/np.pi

        if reshaped:
            __P = np.reshape(__P, self.shape)
        return(__P)

    def print(self, degrees=False):
        print('X =\n', self.getx())
        print('Y =\n', self.gety())
        print('A =\n', self.getamplitude())
        print('P =\n', self.getphase(degrees=degrees))


    def plot(self, degrees=False):
        __X = self.getx(reshaped=False)
        __Y = self.gety(reshaped=False)
        __A = self.getamplitude()
        if degrees:
            __P = self.getphase(degrees=True)
            __P[__P < 0] = 360 + __P[__P < 0]
        else:
            __P = self.getphase()

        __xy = np.array([(__X[i], __Y[i]) for i in np.arange(self.length)])
        
        __plot = plt.subplot(121)
        __plot.matshow(__A)
        plt.title('Amplitude')
        __s = [str(i) for i in __A.flat]
        for i in np.arange(len(__s)):
            __plot.annotate(s=__s[i], xy=__xy[i], ha='center', va='center')

        __plot = plt.subplot(122)
        __plot.matshow(__P)
        plt.title('Phase')
        __s = [str(i) for i in __P.flat]
        for i in np.arange(len(__s)):
            __plot.annotate(s=__s[i], xy=__xy[i], ha='center', va='center')

        plt.show()

class Interpolator1D(object):
    '''
        Interpolator1D(points, axis=1, sort=True) creates the 1-D interpolation
        object using the given points of amplitudes and phases.

        args:
            points ([Point]) : Array of given points
            axis ([1, 2]) : along which axis the interpolation will be done
                            1 : x axis
                            2 : y axis
            sort (boolean) : if sorting of the points is needed.
                            Set to True if the points are not structured
                            Set to False if the points are structured
    '''

    def __init__(self, points, axis=1, sort=True):
        self.points = points
        
        if sort:
            self.points.sort()
        
        self.axis = axis
        self.alpha = np.nan
        self.beta = np.nan

        if axis==1:
            # the line is along the x axis
            self.x = np.array([point.x for point in self.points])
        else:
            # the line is along the y axis
            self.x = np.array([point.y for point in self.points])

        self.a = np.array([point.a for point in self.points])
        self.p = np.array([point.p for point in self.points])
        
        # Check for length and unique values
        if len(self.x) < 2:
            print('Interpolant needs at least two input points.')
        else:
            if np.any(np.diff(self.x) == 0):
                print('The input must be unique elments')
            else:
                pass

    def findindices(self, point):
        __point = point

        if self.axis==1:
            __xi = __point.x
        else:
            __xi = __point.y
        
        __i1 = np.argmax([__xi <= __i for __i in self.x]) - 1
        __i2 = np.argmin([__xi >= __i for __i in self.x])

        return([__i1, __i2])

    def genweight(self, point, indices):
        __point = point

        if self.axis==1:
            __xi = __point.x
        else:
            __xi = __point.y

        __i1 = indices[0]
        __i2 = indices[1]

        __x1 = self.x[__i1]
        __x2 = self.x[__i2]
        __dx = __x2 - __x1

        __nantest = np.isnan(np.array([self.a[__i1], self.a[__i2]]))

        if np.all(__nantest):
            print('NaN value found around [x, y] = ', [__point.x, __point.y])
            print('Output will be nan values...')
            __alpha = np.nan
            __beta = np.nan

        elif np.any(__nantest):
            print('NaN value found for [x, y] = ', [__point.x, __point.y])
            print('Output will be nan values...')

            if np.array(np.where(__nantest)).flat[0]:
                # First element is available
                __alpha = 1
                __beta = 0
            else:
                # Second element is available
                __alpha = 0
                __beta = 1
        else:
            # Both elements are available
            if __dx == 0:
                __alpha = 1
            else:
                __alpha = (__x2 - __xi)/float(__dx)
            
            __beta = 1 - __alpha

        return([__alpha, __beta])

    def interpolate(self, point):
        __point = point
        __indices = self.findindices(point=__point)
        __alpha, __beta = self.genweight(point=__point, indices=__indices)

        __point1 = self.points[__indices[0]]
        __point2 = self.points[__indices[1]]

        __sinval = __alpha*__point1.a*np.sin(__point1.p)+__beta*__point2.a*np.sin(__point2.p)
        __cosval = __alpha*__point1.a*np.cos(__point1.p)+__beta*__point2.a*np.cos(__point2.p)

        __point.p = np.arctan2(__sinval, __cosval)
        __point.a = __sinval/np.sin(__point.p)

        return(__point)

class Interpolator2D(object):
    '''
    Interpolator2D(grid) create the 2D interpolator from the given grid.

    args:
        grid (Grid) :   Grid object containing amplitude and phase from which
                        the interpolation will be made
    '''
    def __init__(self, grid):
        self.sourcegrid = grid

    def interpolatepoint(self, point):
        __point = point
        __pointx = [Point(x=__point.x, y=self.sourcegrid.y[i]) for i in np.arange(self.sourcegrid.shape[0])]

        for i in np.arange(self.sourcegrid.shape[0]):
            # Finding all the interpolated points along y axis
            __points = self.sourcegrid.points[i, :]
            __interpolator = Interpolator1D(points=__points, axis=1, sort=False)
            __pointx[i] = __interpolator.interpolate(point=__pointx[i])

        __interpolator = Interpolator1D(points=__pointx, axis=2, sort=False)
        __point = __interpolator.interpolate(point=__point)

        return(__point)

    def interpolategrid(self, grid):
        __grid = grid
        __points = np.array([self.interpolatepoint(point) for point in __grid.points.flat])
        __grid.points = np.reshape(__points, __grid.shape)

        return(__grid)


if __name__=='__main__':
    print('Under development')