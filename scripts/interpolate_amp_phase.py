"""
Interpolation of Amplitude and Phase using the interpolant derived by Zhigang Xu.
For details see - https://doi.org/10.1007/s10236-017-1122-8

Author: khan
Email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

class APInterpolator(object):
    def __init__(self):
        pass


    class Point1D(object):
        def __init__(self, x, a=None, p=None, radians=False):
            self.x = x
            self.a = a
            if p is not None:
                if radians:
                    self.p = p
                else:
                    self.p = p*np.pi/180.0

        def __lt__(self, other):
            return(self.x < other.x)


    class Point2D(object):
        def __init__(self, x, y, a=None, p=None, radians=False):
            self.x = x
            self.y = y
            self.a = a

            if p is not None:
                if radians:
                    self.p = p
                else:
                    self.p = p*np.pi/180.0

        def __lt__(self, other):
            return(self.x < other.x and self.y < other.y)

    class Grid(object):
        def __init__(self, x, y, A=None, P=None):
            self.x = x
            self.y = y
            self.X, self.Y = np.meshgrid(self.x, self.y, indexing='xy')
            
            if A is None:
                self.A = np.zeros(shape=self.X.shape)
            else:
                self.A = A
            
            if P is None:
                self.P = np.zeros(shape=self.X.shape)
            else:
                self.P = P

        def setamplitude(self, A):
            self.A = A

        def setphase(self, P):
            self.P = P

        def plotamplitude(self):
            plt.contourf(self.X, self.Y, self.A)
            plt.colorbar()
            plt.show()

        def plotphase(self):
            plt.contourf(self.X, self.Y, self.P)
            plt.colorbar()
            plt.show()


    class LinearInterpolator(object):
        def __init__(self, points):
            self.points = points
            self.points.sort()
            
            self.x = [point.x for point in points]
            if len(points) < 2:
                print('Interpolant needs at least two input points.')
            else:
                if np.any(np.diff(self.x) == 0):
                    print('The input must be unique in length')
                else:
                    pass

        def genweight(self, point=None):
            if point is None:
                __xi = self.point.x
            else:
                __xi = point.x
            
            __i1 = np.argmax([__xi <= __i for __i in self.x]) - 1
            __i2 = np.argmin([__xi >= __i for __i in self.x])

            self.point1 = self.points[__i1]
            self.point2 = self.points[__i2]

            __dx = self.point2.x - self.point1.x

            if __dx==0:
                __alpha = 1
            else:
                __alpha = (self.point2.x - self.point.x)/__dx

            __beta = 1 - __alpha

            return(__alpha, __beta)

        def interpolate(self, point):
            self.point = point
            self.alpha, self.beta = self.genweight()
            self.sinval = self.alpha*self.point1.a*np.sin(self.point1.p)+self.beta*self.point2.a*np.sin(self.point2.p)
            self.cosval = self.alpha*self.point1.a*np.cos(self.point1.p)+self.beta*self.point2.a*np.cos(self.point2.p)
            self.point.p = np.arctan2(self.sinval, self.cosval)
            self.point.a = self.sinval/np.sin(self.point.p)

            return(self.point)

    
    class GridInterpolator(object):
        def __init__(self, grid):
            self.grid = grid

        def interpolatepoint(self, point):
            self.point = point

            # Finding the position in the x direction
            self.px1 = np.argmax([self.point.x <= i for i in self.grid.x]) - 1
            self.px2 = np.argmin([self.point.x >= i for i in self.grid.x])

            self.x1 = self.grid.x[self.px1]
            self.x2 = self.grid.x[self.px2]

            # Finding the position in the y direction
            self.py1 = np.argmax([self.point.y <= i for i in self.grid.y]) - 1
            self.py2 = np.argmin([self.point.y >= i for i in self.grid.y])
            print(self.py1)
            self.y1 = self.grid.y[self.py1]
            self.y2 = self.grid.y[self.py2]

            # Linear interpolation along y1
            self.x1y1 = APInterpolator.Point1D(x=self.x1, a=self.grid.A[self.py1, self.px1], p=self.grid.P[self.py1, self.px1])
            self.x2y1 = APInterpolator.Point1D(x=self.x2, a=self.grid.A[self.py1, self.px2], p=self.grid.P[self.py1, self.px2])
            self.inty1 = APInterpolator.Point1D(x=self.point.x)

            self.interp1d = APInterpolator.LinearInterpolator(points=[self.x1y1, self.x2y1])
            self.inty1 = self.interp1d.interpolate(self.inty1)

            # Linear interpolation along y2
            self.x1y2 = APInterpolator.Point1D(x=self.x1, a=self.grid.A[self.py2, self.px1], p=self.grid.P[self.py2, self.px1])
            self.x2y1 = APInterpolator.Point1D(x=self.x2, a=self.grid.A[self.py2, self.px2], p=self.grid.P[self.py2, self.px2])
            self.inty2 = APInterpolator.Point1D(x=self.point.x)

            self.interp1d = APInterpolator.LinearInterpolator(points=[self.x1y1, self.x2y1])
            self.inty2 = self.interp1d.interpolate(self.inty2)

            # Linear interpolation along y1y2
            self.intx1 = self.inty1
            self.intx1.x = self.y1
            self.intx2 = self.inty2
            self.intx2.x = self.y2
            
            self.interp1d = APInterpolator.LinearInterpolator(points=[self.intx1, self.intx2])
            self.intpoint = self.interp1d.interpolate(APInterpolator.Point1D(x=self.point.y))
            self.point.a = self.intpoint.a
            self.point.p = self.intpoint.p

            return(self.point)

        def interpolategrid(self, grid, radians=False):
            self.outgrid = grid

            # Interpolate over all points of the grid
            for i in np.arange(len(self.outgrid.X.flat)):
                self.point = APInterpolator.Point2D(x=self.outgrid.X.flat[i], y=self.outgrid.Y.flat[i], radians=radians)
                self.interpolatepoint(self.point)
                self.outgrid.A.flat[i] = self.point.a
                self.outgrid.P.flat[i] = self.point.p

            return(self.outgrid)



if __name__=='__main__':
    # Implementation of the test cases from the original paper
    # Linear Interpolation
    point1 = APInterpolator.Point1D(x=0, a=1.5, p=45)
    point2 = APInterpolator.Point1D(x=1, a=2.0, p=-45)
    point = APInterpolator.Point1D(x=0.65)

    Interp1D = APInterpolator.LinearInterpolator(points=[point1, point2])
    intval = Interp1D.interpolate(point)

    print('')
    print('             1D Interpolation            ')
    print('-----------------------------------------')
    print('Value      |', 'Paper   |', 'Result      ')
    print('-----------------------------------------')
    print('Amplitude  |', '1.2040  |', intval.a)
    print('Phase      |', '0.4016  |', intval.p)
    print('-----------------------------------------')

    # Grid Interpolation
    inx = np.arange(0, 5)
    iny = np.arange(0, 5)
    INX,INY = np.meshgrid(inx, iny, indexing='xy')
    INA = np.array([[0.1526, 0.6104, 1.1075, 1.0474, 0.4655],
        [0.2135, 0.9765, 1.4736, 1.4135, 0.8315],
        [0.5000, 1.2630, 1.7601, 1.7000, 1.1180],
        [0.6634, 1.4264, 1.9234, 1.8634, 1.2814],
        [0.6787, 1.4417, 1.9387, 1.8787, 1.2967]])
    INA = np.transpose(INA)
    INP = np.array([[6.5937, 17.4027, 122.1608, 259.5022, 384.9741],
        [342.2851, 353.0941, 97.8522, 235.1953, 324.6654],
        [293.6112, 305.4202, 50.1783, 187.5196, 276.9915,],
        [230.8300, 241.6390, 246.3971, 123.7385, 213.2103],
        [160.6516, 171.4606, 276.2187, 53.5601, 143.0320]])
    INP = np.transpose(INP)

    ingrid = APInterpolator.Grid(x=inx, y=iny, A=INA, P=INP)

    outx = np.arange(0.5, 4.5)
    outy = np.arange(0.5, 4.5)
    outgrid = APInterpolator.Grid(x=outx, y=outy)

    interp = APInterpolator.GridInterpolator(grid=ingrid)
    interp.interpolategrid(grid=outgrid)

    plt.matshow(ingrid.A)
    xy = [(ingrid.X.flat[i], ingrid.Y.flat[i]) for i in np.arange(len(ingrid.X.flat))]
    s = [str(i) for i in ingrid.A.flat]
    for i in np.arange(len(s)):
        plt.annotate(s=s[i], xy=xy[i], ha='center', va='center')
    xy = [(outgrid.X.flat[i], outgrid.Y.flat[i]) for i in np.arange(len(outgrid.X.flat))]
    s = [str(i) for i in outgrid.A.flat]
    for i in np.arange(len(s)):
        plt.annotate(s=s[i], xy=xy[i], ha='center', va='center')
    plt.colorbar()
    plt.show()



    # plt.matshow(P)
    # xy = [(X.flat[i], Y.flat[i]) for i in np.arange(len(X.flat))]
    # s = [str(i) for i in A.flat]
    # for i in np.arange(len(s)):
    #     plt.annotate(s=s[i], xy=xy[i], ha='center', va='center')
    # plt.colorbar()
    # plt.show()
    