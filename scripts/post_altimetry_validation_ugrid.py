# -*- coding: utf-8 -*-
"""
Plot the validation with altimetry derived tidal constant. Using the previously
written pre_interp_amp_phase.py classes.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from netCDF4 import Dataset
import os
import sys

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
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = self.X.shape
        self.length = len(self.X.flatten())
        self.isradians = isradians

        if A is None:
            self.A = np.zeros(shape=self.shape)
        else:
            self.A = A
        
        if P is None:
            self.P = np.zeros(shape=self.shape)
        else:
            self.P = P

    def genpoints(self, reshaped=True):
        self.points = np.array(
            [
                Point(
                    x=self.X.flatten()[i],
                    y=self.Y.flatten()[i],
                    a=self.A.flatten()[i],
                    p=self.P.flatten()[i],
                    isradians=self.isradians
                ) for i in np.arange(self.length)
            ]
        )

        if reshaped:
            self.points = np.reshape(self.points, self.shape)
            return(self.points)
        else:
            return(self.points)

    def getx(self, reshaped=True):
        __X = self.X
        if reshaped:
            __X = np.reshape(__X, self.shape)
        return(__X)

    def gety(self, reshaped=True):
        __Y = self.Y
        if reshaped:
            __Y = np.reshape(__Y, self.shape)
        return(__Y)

    def getamplitude(self, reshaped=True):
        if hasattr(self, 'points'):
            __A = np.array([point.a for point in self.points.flatten()])
        else:
            __A = self.A

        if reshaped:
            __A = np.reshape(__A, self.shape)
        return(__A)

    def getphase(self, reshaped=True, degrees=False):
        if hasattr(self, 'points'):
            __P = np.array([point.p for point in self.points.flatten()])
        else:
            __P = self.P

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

class UGrid(object):
    '''
    UGrid (Unstructured Grid) is the class to store the Nodes, Elements and 
    triangulation information of unstructured grid and give similar functionality 
    of interpolation as the Grid class.

    The triangulate class of matplotlib has been used to work out the data
    structure.

    Args:
        nodex (array) : (N) x values of node
        nodey (array) : (N) y values of node
        elements (array) : (N, 3) array of element table, starting index at 1
        A (array) : (N) Value of amplitude
        P (array) : (N) Value of phase
        isradians (bool) : if the Phase is in radians
    '''
    def __init__(self, nodex=np.array([]), nodey=np.array([]), elements=np.array([]), A=np.array([]), P=np.array([]), isradians=False):
        self.nodex = nodex
        self.nodey = nodey
        self.elements = elements - 1
        self.A = A
        self.P = P
        self.isradians = isradians

        try:
            self.triang = mtri.Triangulation(x=self.nodex, y=self.nodey, triangles=self.elements)
            self.trifinder = self.triang.get_trifinder()
        except:
            print('Error! Element table must be counter clockwise with index starting at 1')
            sys.exit(1)

    def get_points(self):
        '''
        Create a Point object for each node with corresponding x, y, a, p value
        and return it. 

        return:
            numpy array of points
        '''
        points = np.array(
            [
                Point(
                    x=x, 
                    y=y, 
                    a=a, 
                    p=p, 
                    isradians=self.isradians
                ) for x, y, a, p in zip(self.nodex, self.nodey, self.A, self.P)
            ]
        )

        return(points)

    def sorrounding_points(self, of):
        '''
        arg:
            of: Point

        return:
            array of sorrounding Point of arg point
        '''
        element = self.trifinder(of.x, of.y)
        nodes = self.elements[element]
        points = np.array(
            [
                Point(
                    x=self.nodex[node],
                    y=self.nodey[node],
                    a=self.A[node],
                    p=self.P[node],
                    isradians=self.isradians
                ) for node in nodes
            ]
        )
        return(points)

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
        self.sourcegrid.genpoints(reshaped=True)

    def interpolatepoint(self, point):
        __point = point

        __pointx = np.array(
            [
                Point(
                    x=__point.x, 
                    y=self.sourcegrid.y[i]
                    ) for i in np.arange(self.sourcegrid.shape[0])
            ]
        )

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
        __grid.print()
        __grid.genpoints(reshaped=False)
        __grid.points = np.array([self.interpolatepoint(point) for point in __grid.points])

        return(__grid)


class FastInterpolator2D(object):
    '''
    Interpolator2D(grid) create the 2D interpolator from the given grid.

    It uses a lookup for individual x,y values, thus minimize the requirement
    for building a large grid.

    args:
        grid (Grid) :   Grid object containing amplitude and phase from which
                        the interpolation will be made
    '''
    def __init__(self, grid):
        self.sourcegrid = grid

    def interpolatepoint(self, point):
        __point = point
        __index_x1 = np.argmax([__point.x <= self.sourcegrid.x]) - 1
        __index_x2 = np.argmin([__point.x >= self.sourcegrid.x]) + 1
        __index_y1 = np.argmax([__point.y <= self.sourcegrid.y]) - 1
        __index_y2 = np.argmin([__point.y >= self.sourcegrid.y]) + 1

        # print(self.sourcegrid.x[__index_x1], __point.x, self.sourcegrid.x[__index_x2])
        # print(self.sourcegrid.y[__index_y1], __point.y, self.sourcegrid.y[__index_y2])

        self.trimgrid = Grid(
            x=self.sourcegrid.x[__index_x1: __index_x2],
            y=self.sourcegrid.y[__index_y1: __index_y2],
            A=self.sourcegrid.A[__index_y1:__index_y2, __index_x1:__index_x2],
            P=self.sourcegrid.P[__index_y1:__index_y2, __index_x1:__index_x2],
            isradians=self.sourcegrid.isradians
        )
        self.trimgrid.genpoints()
        __pointx = np.array(
            [
                Point(
                    x=__point.x, 
                    y=self.trimgrid.y[i]
                ) for i in np.arange(self.trimgrid.shape[0])
            ]
        )
        self.trimgrid.shape[0]
        for i in np.arange(self.trimgrid.shape[0]):
            # Finding all the interpolated points along y axis
            __points = self.trimgrid.points[i, :]
            __interpolator = Interpolator1D(points=__points, axis=1, sort=False)
            __pointx[i] = __interpolator.interpolate(point=__pointx[i])

        __interpolator = Interpolator1D(points=__pointx, axis=2, sort=False)
        __point = __interpolator.interpolate(point=__point)

        return(__point)

    def interpolategrid(self, grid):
        __grid = grid
        __points = __grid.genpoints(reshaped=False)
        __points = np.array([self.interpolatepoint(point) for point in __points])
        __grid.points = np.reshape(__points, __grid.shape)
        return(__grid) 

class UGridInterpolator(object):
    def __init__(self, ugrid):
        self.ugrid = ugrid
    
    def interpolatepoint(self, point):
        # Finding the triangle points that holds the point
        p1, p2, p3 = self.ugrid.sorrounding_points(of=point)

        # Finding linear weight applied to the three points
        denom = (p2.y-p3.y)*(p1.x-p3.x)+(p3.x-p2.x)*(p1.y-p3.y)
        alpha = ((p2.y-p3.y)*(point.x-p3.x)+(p3.x-p2.x)*(point.y-p3.y))/denom
        beta = ((p3.y-p1.y)*(point.x-p3.x)+(p1.x-p3.x)*(point.y-p3.y))/denom
        gamma = 1 - alpha - beta

        # Calculating amplitude and phase
        sin_value = alpha*p1.a*np.sin(p1.p) + beta*p2.a*np.sin(p2.p) + gamma*p3.a*np.sin(p3.p)
        cos_value = alpha*p1.a*np.cos(p1.p) + beta*p2.a*np.cos(p2.p) + gamma*p3.a*np.cos(p3.p)
        
        point.p = np.arctan2(sin_value, cos_value)
        point.a = sin_value/np.sin(point.p)

        return(point)

    def interpolategrid(self, grid):
        grid = grid
        points = grid.genpoints(reshaped=False)
        points = np.array([self.interpolatepoint(point) for point in points])
        grid.points = np.reshape(points, grid.shape)
        return(grid)


def bytes_to_string(chars):
    string = ''
    for char in chars:
        string = string+char.decode('UTF-8')
    return(string)


if __name__=='__main__':
    # Setting parameters
    atlas_loc = '/run/media/khan/Workbench/Projects/Tide/Experiments/EXP05_2014_Sa/Atlas'
    alt_loc = '/run/media/khan/Workbench/Data/Altimetry/Constituents/TP+J1+J2'
    save_loc = '/run/media/khan/Workbench/Projects/Tide/Experiments/EXP05_2014_Sa/AltimetryValidation'
    alt_tracks = np.array([40, 116, 192, 14, 90, 155, 231, 53, 129, 205, 27])
    wave_names = np.array(['M2', 'S2', 'O1', 'K1', 'Sa', 'Ssa'])
    
    for wave_name in wave_names:
        # Loading the netcdf grid
        atlas_file = '{:s}-elevation-atlas.nc'.format(wave_name)
        atlas_nc = Dataset(os.path.join(atlas_loc, atlas_file))
        atlas_lon = atlas_nc['lon'][:]
        atlas_lat = atlas_nc['lat'][:]
        atlas_elem = atlas_nc['element'][:]
        atlas_amp = atlas_nc['elevation_a'][:]
        atlas_pha = atlas_nc['elevation_G'][:]
        atlas_grid = UGrid(nodex=atlas_lon, nodey=atlas_lat, elements=atlas_elem, A=atlas_amp, P=atlas_pha, isradians=False)

        # Loading the altimetry grid
        for alt_track in alt_tracks:
            print('Track - {:d}, Wave - {:s}'.format(alt_track, wave_name))
            alt_file = 'ctoh.harmonics.ref.TP+J1+J2.nindian.{:03d}.nc'.format(alt_track)
            alt_nc = Dataset(os.path.join(alt_loc, alt_file), 'r')
            alt_consts = np.array([bytes_to_string(const) for const in alt_nc['constituentname'][:]])
            alt_wave_index = np.array(np.where(alt_consts==wave_name)).flatten()
            alt_lat = alt_nc['lat'][:]
            alt_lon = alt_nc['lon'][:]
            alt_amp = alt_nc['amplitude']
            alt_pha = alt_nc['phase_lag']
            alt_lat_criterion = (alt_lat >= np.min(atlas_lat)) & (alt_lat <= np.max(atlas_lat))
            alt_lon_criterion = (alt_lon >= np.min(atlas_lon)) & (alt_lon <= np.max(atlas_lon))
            indice = np.array(np.where(alt_lat_criterion & alt_lon_criterion)).flatten()
            alt_points = np.array(
                [
                    Point(
                        x=alt_lon[index],
                        y=alt_lat[index],
                        a=alt_amp[index, alt_wave_index],
                        p=alt_pha[index, alt_wave_index]
                    )
                    for index in indice
                ]
            )
            
            # Interpolation from the netcdf grid using Tide Interpolation algorithm
            atlas_points = np.array(
                [
                    Point(x=alt_lon[index], y=alt_lat[index]) for index in indice
                ]
            )
            atlas_interpolator = UGridInterpolator(ugrid=atlas_grid)
            atlas_interp_points = np.array(
                [atlas_interpolator.interpolatepoint(point=point) for point in atlas_points]
            )

            interpolated = np.array(
                [
                    [
                        i.x, i.y, i.a, i.p, j.a, j.p
                    ] for i,j in zip(alt_points, atlas_interp_points)
                ]
            )
            interpolated_file = os.path.join(save_loc, f'{alt_track:d}_{wave_name:s}.csv')
            np.savetxt(fname=interpolated_file, X=interpolated, delimiter=',')