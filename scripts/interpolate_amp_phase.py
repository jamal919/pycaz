"""
Interpolation of Amplitude and Phase using the interpolant derived by Zhigang Xu.
For details see - https://doi.org/10.1007/s10236-017-1122-8

Author: khan
Email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

class Point(object):
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
        __points = [point for point in self.points.flat]
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

        __xy = [(__X[i], __Y[i]) for i in np.arange(self.length)]
        
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
    def __init__(self, points, axis=1, sort=True):
        self.points = points
        
        if sort:
            self.points.sort()
        
        self.axis = axis
        self.alpha = np.nan
        self.beta = np.nan

        if axis==1:
            # the line is along the x axis
            self.x = [point.x for point in self.points]
        else:
            # the line is along the y axis
            self.x = [point.y for point in self.points]
        
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

        if __dx == 0:
            __alpha = 1
        else:
            __alpha = (__x2 - __xi)/__dx

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
    # Implementation of the test cases from the original paper
    # Linear Interpolation
    point1 = Point(x=0, a=1.5, p=45)
    point2 = Point(x=1, a=2, p=-45)
    pointi = Point(x=0.65)

    Interp1D = Interpolator1D(points=[point1, point2])
    pointi = Interp1D.interpolate(pointi)

    print('')
    print('             1D Interpolation            ')
    print('-----------------------------------------')
    print('Value      |', 'Paper   |', 'Result      ')
    print('-----------------------------------------')
    print('Amplitude  |', '1.4020  |', pointi.a)
    print('Phase      |', '-0.4016 |', pointi.p)
    print('-----------------------------------------')

    # Grid Interpolation
    # Input grid and amplitude phase
    # The amplitude and phase are meant to be input as row major format
    x = np.arange(0, 5)
    y = np.arange(0, 5)
    A = np.array([[0.1526, 0.6104, 1.1075, 1.0474, 0.4655],
        [0.2135, 0.9765, 1.4736, 1.4135, 0.8315],
        [0.5000, 1.2630, 1.7601, 1.7000, 1.1180],
        [0.6634, 1.4264, 1.9234, 1.8634, 1.2814],
        [0.6787, 1.4417, 1.9387, 1.8787, 1.2967]])
    A = np.transpose(A) # Chaning to row major format
    P = np.array([[6.5937, 17.4027, 122.1608, 259.5022, 348.9741],
        [342.2851, 353.0941, 97.8522, 235.1953, 324.6654],
        [294.6112, 305.4202, 50.1783, 187.5196, 276.9915,],
        [230.8300, 241.6390, 346.3971, 123.7385, 213.2103],
        [160.6516, 171.4606, 276.2187, 53.5601, 143.0320]])
    P = np.transpose(P) # Changing to row major format

    ingrid = Grid(x=x, y=y, A=A, P=P, isradians=False)

    # Output grid
    outx = np.arange(0.5, 4.5)
    outy = np.arange(0.5, 4.5)
    outgrid = Grid(x=outx, y=outy)

    ipl = Interpolator2D(ingrid)
    out = ipl.interpolategrid(grid=outgrid)
    out.print()
    out.plot(degrees=True)