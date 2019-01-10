#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interpolation of Amplitude and Phase using the interpolant derived by Zhigang Xu.
For details see - https://doi.org/10.1007/s10236-017-1122-8

This new update to algorithm tries to minimize the time of interpolation from a 
large grid. The previous approach was to build the grid at the first step and 
then use that to interpolate. However, for a very large grid it can get very time 
consuming to build such a large grid for a small number of points. 

The interpolation scheme is same as before, but in the class FastInterpolator2D
looks at the individual points and build the point grid afterwards.

Author: khan
Email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
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

    def get_amplitude(self, factor=1):
        return(self.a * factor)

    def get_phase(self, in_degrees=False, in_positive=True):
        if not in_degrees:
            return(self.p)
        else:
            p = np.degrees(self.p)
            if not in_positive:
                return(p)
            else:
                if p < 0:
                    return(360 + p)
                else:
                    return(p)

    def __repr__(self):
        return('x : {:4f}, y: {:4f}, a: {:4f}, p: {:4f}'.format(self.x, self.y, self.a, self.p))

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
        '''
        '''
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
        '''
        '''
        X = self.X
        if reshaped:
            X = np.reshape(X, self.shape)
        else:
            X = np.array([x for x in self.X.flatten()])
        
        return(X)

    def gety(self, reshaped=True):
        '''
        '''
        Y = self.Y
        if reshaped:
            Y = np.reshape(Y, self.shape)
        else:
            Y = np.array([y for y in self.Y.flatten()])
        
        return(Y)

    def getamplitude(self, reshaped=True):
        if hasattr(self, 'points'):
            A = np.array([point.a for point in self.points.flatten()])
        else:
            A = self.A

        if reshaped:
            A = np.reshape(A, self.shape)
        return(A)

    def getphase(self, reshaped=True, in_degrees=False, in_positive=False):
        if hasattr(self, 'points'):
            P = np.array([point.p for point in self.points.flatten()])
        else:
            P = self.P

        if in_degrees:
            P = np.degrees(P)
            if in_positive:
                P[P < 0] = 360 + P[P < 0]

        if reshaped:
            P = np.reshape(P, self.shape)
        return(P)

    def plot(self, in_degrees=False, in_positive=True):
        X = self.getx(reshaped=False)
        Y = self.gety(reshaped=False)
        A = self.getamplitude(reshaped=True)
        P = self.getphase(reshaped=True, in_degrees=in_degrees, in_positive=in_positive)

        xy = np.array([(X[i], Y[i]) for i in np.arange(self.length)])
        
        plot = plt.subplot(121)
        plot.matshow(A)
        plt.title('Amplitude')
        s = [str('{:4f}'.format(i)) for i in A.flat]
        for i in np.arange(len(s)):
            plot.annotate(s=s[i], xy=xy[i], ha='center', va='center')

        plot = plt.subplot(122)
        plot.matshow(P)
        plt.title('Phase')
        s = [str('{:4f}'.format(i)) for i in P.flat]
        for i in np.arange(len(s)):
            plot.annotate(s=s[i], xy=xy[i], ha='center', va='center')

        plt.show()

    def to_netcdf(self, fname):
        pass

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
        point = point

        if self.axis==1:
            xi = point.x
        else:
            xi = point.y
        
        i1 = np.argmax([xi <= i for i in self.x]) - 1
        i2 = np.argmin([xi >= i for i in self.x])

        return([i1, i2])

    def genweight(self, point, indices):
        point = point

        if self.axis==1:
            xi = point.x
        else:
            xi = point.y

        i1 = indices[0]
        i2 = indices[1]

        x1 = self.x[i1]
        x2 = self.x[i2]
        dx = x2 - x1

        nantest = np.isnan(np.array([self.a[i1], self.a[i2]]))

        if np.all(nantest):
            print('NaN value found around [x, y] = ', [point.x, point.y])
            print('Output will be nan values...')
            alpha = np.nan
            beta = np.nan

        elif np.any(nantest):
            print('NaN value found for [x, y] = ', [point.x, point.y])
            print('Output will be nan values...')

            if np.array(np.where(nantest)).flat[0]:
                # First element is available
                alpha = 1
                beta = 0
            else:
                # Second element is available
                alpha = 0
                beta = 1
        else:
            # Both elements are available
            if dx == 0:
                alpha = 1
            else:
                alpha = (x2 - xi)/float(dx)
            
            beta = 1 - alpha

        return([alpha, beta])

    def interpolate(self, point):
        point = point
        indices = self.findindices(point=point)
        alpha, beta = self.genweight(point=point, indices=indices)

        point1 = self.points[indices[0]]
        point2 = self.points[indices[1]]

        sinval = alpha*point1.a*np.sin(point1.p)+beta*point2.a*np.sin(point2.p)
        cosval = alpha*point1.a*np.cos(point1.p)+beta*point2.a*np.cos(point2.p)

        point.p = np.arctan2(sinval, cosval)
        point.a = sinval/np.sin(point.p)

        return(point)

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
        point = point

        pointx = np.array(
            [
                Point(
                    x=point.x, 
                    y=self.sourcegrid.y[i]
                    ) for i in np.arange(self.sourcegrid.shape[0])
            ]
        )

        for i in np.arange(self.sourcegrid.shape[0]):
            # Finding all the interpolated points along y axis
            points = self.sourcegrid.points[i, :]
            interpolator = Interpolator1D(points=points, axis=1, sort=False)
            pointx[i] = interpolator.interpolate(point=pointx[i])

        interpolator = Interpolator1D(points=pointx, axis=2, sort=False)
        point = interpolator.interpolate(point=point)

        return(point)

    def interpolategrid(self, grid):
        grid = grid
        grid.genpoints(reshaped=False)
        grid.points = np.array([self.interpolatepoint(point) for point in grid.points])

        return(grid)


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
        point = point
        index_x1 = np.argmax([point.x <= self.sourcegrid.x]) - 1
        index_x2 = np.argmin([point.x >= self.sourcegrid.x]) + 1
        index_y1 = np.argmax([point.y <= self.sourcegrid.y]) - 1
        index_y2 = np.argmin([point.y >= self.sourcegrid.y]) + 1

        # print(self.sourcegrid.x[__index_x1], __point.x, self.sourcegrid.x[__index_x2])
        # print(self.sourcegrid.y[__index_y1], __point.y, self.sourcegrid.y[__index_y2])

        self.trimgrid = Grid(
            x=self.sourcegrid.x[index_x1: index_x2],
            y=self.sourcegrid.y[index_y1: index_y2],
            A=self.sourcegrid.A[index_y1:index_y2, index_x1:index_x2],
            P=self.sourcegrid.P[index_y1:index_y2, index_x1:index_x2],
            isradians=self.sourcegrid.isradians
        )
        self.trimgrid.genpoints()
        pointx = np.array(
            [
                Point(
                    x=point.x, 
                    y=self.trimgrid.y[i]
                ) for i in np.arange(self.trimgrid.shape[0])
            ]
        )
        self.trimgrid.shape[0]
        for i in np.arange(self.trimgrid.shape[0]):
            # Finding all the interpolated points along y axis
            points = self.trimgrid.points[i, :]
            interpolator = Interpolator1D(points=points, axis=1, sort=False)
            pointx[i] = interpolator.interpolate(point=pointx[i])

        interpolator = Interpolator1D(points=pointx, axis=2, sort=False)
        point = interpolator.interpolate(point=point)

        return(point)

    def interpolategrid(self, grid):
        grid = grid
        points = grid.genpoints(reshaped=False)
        points = np.array([self.interpolatepoint(point) for point in points])
        grid.points = np.reshape(points, grid.shape)
        return(grid)

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


if __name__=='__main__':
    # Implementation of the test cases from the original paper
    # Linear Interpolation
    point1 = Point(x=0, a=1.5, p=45)
    point2 = Point(x=1, a=2, p=-45)
    pointi = Point(x=0.65)

    Interp1D = Interpolator1D(points=[point1, point2])
    pointi = Interp1D.interpolate(pointi)

    print('{:=>40}\n{: ^40}\n{:=>40}'.format('', '1D Interpolation', ''))
    print('|{:-<38}|'.format(''))
    print('|{:^12}|{:^12}|{:^12}|'.format('Value', 'Paper', 'Script'))
    print('|{:-<38}|'.format(''))
    print('|{:^12s}|{:^12.4f}|{:^12.4f}|'.format('Amplitude', 1.4020, pointi.a))
    print('|{:^12s}|{:^12.4f}|{:^12.4f}|'.format('Phase', -0.4016, pointi.p))
    print('|{:-<38}|\n'.format(''))

    # Grid Interpolation
    # Input grid and amplitude phase
    # The amplitude and phase are meant to be input as row major format
    inx = np.arange(0, 5)
    iny = np.arange(0, 5)
    inA = np.array(
        [[0.1526, 0.6104, 1.1075, 1.0474, 0.4655],
        [0.2135, 0.9765, 1.4736, 1.4135, 0.8315],
        [0.5000, 1.2630, 1.7601, 1.7000, 1.1180],
        [0.6634, 1.4264, 1.9234, 1.8634, 1.2814],
        [0.6787, 1.4417, 1.9387, 1.8787, 1.2967]]
    )
    inA = np.transpose(inA) # Chaning to row major format
    inP = np.array(
        [[6.5937, 17.4027, 122.1608, 259.5022, 348.9741],
        [342.2851, 353.0941, 97.8522, 235.1953, 324.6654],
        [294.6112, 305.4202, 50.1783, 187.5196, 276.9915,],
        [230.8300, 241.6390, 346.3971, 123.7385, 213.2103],
        [160.6516, 171.4606, 276.2187, 53.5601, 143.0320]]
    )
    inP = np.transpose(inP) # Changing to row major format

    ingrid = Grid(x=inx, y=iny, A=inA, P=inP, isradians=False)

    # Output grid
    outx = np.arange(0.5, 4.5)
    outy = np.arange(0.5, 4.5)
    outgrid = Grid(x=outx, y=outy)

    ipl = Interpolator2D(ingrid)
    out = ipl.interpolategrid(grid=outgrid)
    # print('{:=>40}\n{: ^40}\n{:=>40}'.format('', '2D Interpolation', ''))
    # out.plot(in_degrees=True, in_positive=True)

    # Extended from paper
    # Triangular test
    x = np.array([0, 0.5, 1])
    y = np.array([0, 1, 0])
    triang = np.array([[1, 3, 2]])
    a = np.array([1, 2, 3])
    p = np.array([10, 190, 150])

    intp = Point(x=0.5, y=0.5)

    ugrid = UGrid(nodex=x, nodey=y, elements=triang, A=a, P=p, isradians=False)
    upl = UGridInterpolator(ugrid=ugrid)
    intp = upl.interpolatepoint(point=intp)
    print(ugrid.sorrounding_points(of=Point(x=0.4, y=0.5)))
    print(intp)