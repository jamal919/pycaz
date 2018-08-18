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
        def __init__(self):
            pass


if __name__=='__main__':
    # Implementation of the test cases from the original paper
    # Linear Interpolation
    point1 = APInterpolator.Point1D(x=0, a=1.5, p=45)
    point2 = APInterpolator.Point1D(x=1, a=2.0, p=-45)
    point = APInterpolator.Point1D(x=0.65)

    Interp1D = APInterpolator.LinearInterpolator(points=[point1, point2])
    intval = Interp1D.interpolate(point)

    print(intval.a, intval.p)

    # Grid Interpolation
    x = np.arange(90, 95)
    y = np.arange(20, 25)
    X,Y = np.meshgrid(x, y, indexing='ij')
    A = np.array([[0.1526, 0.6104, 1.1075, 1.0474, 0.4655],
        [0.2135, 0.9765, 1.4736, 1.4135, 0.8315],
        [0.5000, 1.2630, 1.7601, 1.7000, 1.1180],
        [0.6634, 1.4264, 1.9234, 1.8634, 1.2814],
        [0.6787, 1.4417, 1.9387, 1.8787, 1.2967]])
    P = np.array([[6.5937, 17.4027, 122.1608, 259.5022, 384.9741],
        [342.2851, 353.0941, 97.8522, 235.1953, 324.6654],
        [293.6112, 305.4202, 50.1783, 187.5196, 276.9915,],
        [230.8300, 241.6390, 246.3971, 123.7385, 213.2103],
        [160.6516, 171.4606, 276.2187, 53.5601, 143.0320]])

    print([i for i in X.flat])
    print([i for i in Y.flat])