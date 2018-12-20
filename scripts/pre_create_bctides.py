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

The approach here is to set-up various types of boundary in terms of boundary type,
i.e., eta, s/t, u/v etc.

The general format of the header line of the boundary is following - 
nnodes, elev, velocity, temperature, salinity ... and other modules if needed.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import glob
import sys
import re

class Node(object):
    def __init__(self, id, x, y, z):
        self.id = int(id)
        self.x = x
        self.y = y
        self.z = z

    def __lt__(self, other):
        if isinstance(other, Node):
            return(self.id < other.id)
        else:
            return(self.id < other)

class Element(object):
    def __init__(self, id, nnode, connectivity=[]):
        self.id = id
        self.nnode = nnode
        self.connectivity = connectivity

    def __lt__(self, other):
        if isinstance(other, Element):
            return(self.id < other.id)
        else:
            return(self.id < other)

class Mesh(object):
    def __init__(self, grname=None, nelem=0, nnode=0, nodes=[], elems=[]):
        self.grname = grname
        self.nelem = nelem
        self.nnode = nnode
        self.nodes = nodes
        self.elems = elems

    def read(self, fname, path='./', readnodes=True, readelements=True):
        __file = os.path.join(path, fname)
        with open(__file) as f:
            ds = f.readlines()
            __line = 0
            
            # Reading the gr3 name
            self.grname = ds[__line].strip()
            __line = __line + 1
            
            # Reading the number of nodes and elements
            __nelem, __nnode = np.fromstring(ds[__line].split('\n')[0], count=2, sep=' ')
            self.nelem = int(__nelem)
            self.nnode = int(__nnode)
            __line = __line + 1

            if readnodes:
                # Reading the nodes
                __nodes = np.genfromtxt(fname=ds[__line:__line+self.nnode])
                __line = __line + self.nnode
                self.nodes = np.array([Node(i[0], i[1], i[2], i[3]) for i in __nodes])

            if readelements:
                # Reading the elements
                __elems = np.genfromtxt(fname=ds[__line:__line+self.nelem], dtype=int)
                __line = __line + self.nelem
                self.elems = np.array([Element(i[0], i[1], i[2:len(i)]) for i in __elems])

    def write(self, fname, path='./'):
        with open(os.path.join(path, fname), 'w') as f:
            f.write('{:s}\n'.format(self.grname))
            f.write('{:d}\t{:d}\n'.format(self.nelem, self.nnode))
            
            for __node in self.nodes:
                f.write('{:d}\t{:.10f}\t{:.10f}\t{:.10f}\n'\
                        .format(int(__node.id), __node.x, __node.y, __node.z))

            for __elem in self.elems:
                f.write('{:d}\t{:d}'.format(__elem.id, __elem.nnode))
                [f.write('\t{:d}'.format(i)) for i in __elem.connectivity]
                f.write('\n')

class Boundary(object):
    def __init__(self, nnodes, nodes, landflag=None, bndname='', condition=None):
        """ SCHISM complient bounary
        
        Args:
            nnodes (int)     :   Serial number of the boundary
            nodes(int []):   Numpy array of nodes forming the boundary
            landflg(int)    :   Boundary type in case of land boundary.
                                For open boundary, no information needed.
            bndname(str)    :   Name of the boundary (optional)
        """
        self.number = nnodes
        self.nodes = nodes
        self.landflag = landflag
        self.name = bndname

    def countnodes(self):
        """Number of nodes in a boundary
        
        Returns:
            int : Number of nodes of the boundary
        """
        return(len(self.nodes))

class Boundaries(object):
    def __init__(self, openbnd=[], landbnd=[]):
        self.open = openbnd
        self.nopen = len(self.open)
        self.land = landbnd
        self.nland = len(self.land)

    def read(self, fname, mesh, hgrid=True):
        __mesh = mesh
        try:
            __file = fname
            with open(__file) as f:
                __ds = f.readlines()
        except:
            print('File - {:s} - does not exist!'.format(__file))
        else:
            if hgrid:
                # Reading the gr3 name
                __line = 0
                __grname = __ds[__line].strip()
                # Reading the number of nodes and elements
                __line = __line + 1
                __nelem, __nnode = np.fromstring(__ds[__line].split('\n')[0], count=2, sep=' ')
                __nelem = int(__nelem)
                __nnode = int(__nnode)

                # Boundary definition chunk in ds
                __line = __line + 1 + __nelem + __nnode
                __ds = __ds[__line:]
            
            # Reading open boundary sagments
            __line = 0
            self.nopen = int(__ds[__line].split()[0])
            
            __line = __line + 1
            self.nopennodes = int(__ds[__line].split()[0])

            __line = __line + 1
            for bnd in np.arange(self.nopen):
                __nnodes = int(__ds[__line].split()[0])
                
                __line = __line + 1
                __nodeids = np.genfromtxt(fname=__ds[__line:__line+__nnodes], dtype=int)
                __nodes = __mesh.nodes[__nodeids-1] # Python numbering scheme
                print(__nodeids[0], __nodes[0].id)

                __line = __line + __nnodes
                self.open.append(Boundary(nnodes=__nnodes, nodes=__nodes, bndname=str(bnd)))

            # Reading land boundary sagments
            self.nland = int(__ds[__line].split()[0])

            __line = __line + 1
            self.nlandnodes = int(__ds[__line].split()[0])

            __line = __line + 1
            for bnd in np.arange(self.nland):
                __nnodes, __flag = [int(i) for i in __ds[__line].split()[0:2]]
                
                __line = __line + 1
                __nodeids = np.genfromtxt(fname=__ds[__line:__line+__nnodes], dtype=int)
                __nodes = __mesh.nodes[__nodeids-1] # Python numbering scheme

                __line = __line + __nnodes
                self.land.append(Boundary(nnodes=__nnodes, nodes=__nodes, landflag=__flag, bndname=str(bnd)))

class Hgrid(object):
    def __init__(self, mesh=None, boundaries=None):
        self.mesh = mesh
        self.boundaries = boundaries

    def read(self, fname, path='./'):
        self.mesh = Mesh()
        self.mesh.read(fname=fname, path=path)

        self.boundaries = Boundaries()
        self.boundaries.read(fname=os.path.join(path, fname), mesh=self.mesh, hgrid=True)
    
    def coverage(self, padding=1):
        __padding = padding
        if self.mesh is not None:
            __xmin = int(np.floor(np.min(np.array([node.x for node in self.mesh.nodes]))) - __padding)
            __xmax = int(np.ceil(np.max(np.array([node.x for node in self.mesh.nodes]))) + __padding)
            __ymin = int(np.floor(np.min(np.array([node.y for node in self.mesh.nodes]))) - __padding)
            __ymax = int(np.ceil(np.max(np.array([node.y for node in self.mesh.nodes]))) + __padding)

        return(__xmin, __xmax, __ymin, __ymax)

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
        self.length = len(__X.flatten())
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
        self.points = np.array([Point(x=__X.flatten()[i], y=__Y.flatten()[i], a=__A.flatten()[i], p=__P.flatten()[i], isradians=self.isradians) for i in np.arange(self.length)])
        self.points = np.reshape(self.points, self.shape)

    def getpoints(self, reshaped=True):
        __points = np.array([point for point in self.points.flatten()])
        if reshaped:
            return(np.reshape(self.points, self.shape))
        else:
            return(self.points)

    def getx(self, reshaped=True):
        __X = np.array([point.x for point in self.points.flatten()])
        if reshaped:
            __X = np.reshape(__X, self.shape)
        return(__X)

    def gety(self, reshaped=True):
        __Y = np.array([point.y for point in self.points.flatten()])
        if reshaped:
            __Y = np.reshape(__Y, self.shape)
        return(__Y)

    def getamplitude(self, reshaped=True):
        __A = np.array([point.a for point in self.points.flatten()])
        if reshaped:
            __A = np.reshape(__A, self.shape)
        return(__A)

    def getphase(self, reshaped=True, degrees=False):
        __P = np.array([point.p for point in self.points.flatten()])
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

class Atlas(object):
    def __init__(self, waves={}):
        self.waves = waves

    def loadwaves(self, path, pattern='{wave}_FES2012_SLEV.nc', x='lon', y='lat', A='Ha', P='Hg'):
        __path = path
        __pattern = re.compile(pattern.format(wave='(\S+)'))
        fnames = glob.glob(os.path.join(__path, '*'))
        waves = [__pattern.findall(os.path.basename(fname))[0] for fname in fnames]
        
        for wave,fname in zip(waves, fnames):
            self.waves[wave.upper()] = dict(fname=fname, x=x, y=y, A=A, P=P)

    def wave(self, wave, coverage = None, dim='yx'):
        __wave = wave.upper()
        if __wave in self.waves:
            __waveinfo = self.waves[__wave]
            __nc = Dataset(__waveinfo['fname'], mode='r')
            __x = __nc.variables[__waveinfo['x']][:]
            __y = __nc.variables[__waveinfo['y']][:]
            if dim is 'yx':
                if coverage is None:
                    __A = __nc.variables[__waveinfo['A']][:]
                    __P = __nc.variables[__waveinfo['P']][:]
                    __grid = Grid(x=__x, y=__y)
                    __nc.close()
                else:
                    __xmin, __xmax, __ymin, __ymax = coverage
                    __x_left = np.argwhere(__x >= __xmin)[0]
                    __x_right = np.argwhere(__x <= __xmax)[-1]
                    __y_left = np.argwhere(__y >= __ymin)[0]
                    __y_right = np.argwhere(__y <= __ymax)[-1]
                    __x_selected = __nc.variables[__waveinfo['x']][__x_left:__x_right]
                    __y_selected = __nc.variables[__waveinfo['y']][__y_left:__y_right]
                    __A = __nc.variables[__waveinfo['A']][__y_left:__y_right,__x_left:__x_right]
                    __P = __nc.variables[__waveinfo['P']][__y_left:__y_right,__x_left:__x_right]
                    print(__x_selected.shape, __y_selected.shape, __A.shape, __P.shape)
                    __grid = Grid(x=__x_selected, y=__y_selected, A=__A, P=__P)
                return(__grid)
            elif dim is 'xy':
                print('The variables are set to be in xy format!')
                print('Under development')
                return(None)
        else:
            sys.exit('Failed! Constituent {wave} is not available!'.format(wave=__wave))

class TidalPotential(object):
    """
    The tidal potential enters the momentum equation as a body force term.
    The self.wave variable is a set of values taken from Reid 1990.
    """
    def __init__(self):
        self.waves = {
            '2N2':{'jspc':2, 'tamp':0.006141, 'tfreq':0.0001352404964640, 'tnf':0.96720, 'tear':251.59},
            'K1':{'jspc':1, 'tamp':0.141565, 'tfreq':0.0000729211583580, 'tnf':1.10338, 'tear':324.30},
            'K2':{'jspc':2, 'tamp':0.030684, 'tfreq':0.0001458423172010, 'tnf':1.28346, 'tear':109.01},
            'L2':{'jspc':2, 'tamp':0.006931, 'tfreq':0.0001431581055310, 'tnf':0.00000, 'tear':325.06},
            'M2':{'jspc':2, 'tamp':0.242334, 'tfreq':0.0001405189025090, 'tnf':0.96720, 'tear':313.79},
            'MU2':{'jspc':2, 'tamp':0.007408, 'tfreq':0.0001355937006840, 'tnf':0.96720, 'tear':266.58},
            'N2':{'jspc':2, 'tamp':0.046397, 'tfreq':0.0001378796994870, 'tnf':0.96720, 'tear':102.69},
            'NU2':{'jspc':2, 'tamp':0.008811, 'tfreq':0.0001382329037070, 'tnf':0.96720, 'tear':117.68},
            'O1':{'jspc':1, 'tamp':0.100661, 'tfreq':0.0000675977441510, 'tnf':1.16763, 'tear':348.06},
            'P1':{'jspc':1, 'tamp':0.046848, 'tfreq':0.0000725229459750, 'tnf':1.00000, 'tear':39.25},
            'Q1':{'jspc':1, 'tamp':0.019273, 'tfreq':0.0000649585411290, 'tnf':1.16763, 'tear':136.96},
            'S2':{'jspc':2, 'tamp':0.112743, 'tfreq':0.0001454441043330, 'tnf':1.00000, 'tear':0.00},
            'T2':{'jspc':2, 'tamp':0.006608, 'tfreq':0.0001452450073530, 'tnf':1.00000, 'tear':52.32}
            }

    def value(self, wavelist='default'):
        __values = {}
        if wavelist is 'default':
            __values = {wave:self.waves[wave] for wave in self.waves.keys()}
            return(__values)
        else:
            for wave in wavelist:
                if wave in self.waves.keys():
                    __values[wave] = self.waves[wave]
                else:
                    print('Wave {:s} - Not found!'.format(wave))
            return(__values)

class TidalForcing(object):
    def __init__(self):
        '''
        For each node where tidal forcing is expected, the common information 
        for each tidal constituent is described by tidal forcing. For each constituent
        they are the angular frequency (rad/sec), nodal factor, and earth equili-
        brium argument.According to SCHISM manual these variables are denoted 
        with amig, ff, and face respectively.
        '''
        self.waves = {
            'Q1':{'amig':6.4958541128674080E-05, 'ff':1.167380, 'face':136.93},
            'S2':{'amig':1.4544410433286079E-04, 'ff':1.000000, 'face':0.00},
            'S1':{'amig':7.2722052166430395E-05, 'ff':1.000000, 'face':180.00},
            'S4':{'amig':2.9088820866572158E-04, 'ff':1.000000, 'face':0.00},
            'M4':{'amig':2.8103780501728730E-04, 'ff':0.935590, 'face':267.58},
            'K1':{'amig':7.2921158357870547E-05, 'ff':1.103230, 'face':324.33},
            'M6':{'amig':4.2155670801074461E-04, 'ff':0.904960, 'face':221.38},
            'M3':{'amig':2.1077835376296549E-04, 'ff':0.951380, 'face':110.69},
            'M2':{'amig':1.4051890250864360E-04, 'ff':0.967260, 'face':313.79},
            'MSF':{'amig':4.9252018242171696E-06, 'ff':0.967260, 'face':47.20},
            'M8':{'amig':5.6207561051938818E-04, 'ff':0.875340, 'face':175.17},
            'O1':{'amig':6.7597744150773076E-05, 'ff':1.167380, 'face':348.03},
            'MN4':{'amig':2.7839860199518830E-04, 'ff':0.935590, 'face':56.49},
            '2N2':{'amig':1.3524049646444560E-04, 'ff':0.967260, 'face':251.60},
            'P1':{'amig':7.2522945974990243E-05, 'ff':1.000000, 'face':39.25},
            'R2':{'amig':1.4564320131284101E-04, 'ff':1.000000, 'face':127.68},
            'NU2':{'amig':1.3823290370652551E-04, 'ff':0.967260, 'face':117.69},
            'T2':{'amig':1.4524500735288060E-04, 'ff':1.000000, 'face':52.32},
            'N2':{'amig':1.3787969948654460E-04, 'ff':0.967260, 'face':102.70},
            'K2':{'amig':1.4584231720055481E-04, 'ff':1.282940, 'face':109.07},
            'MF':{'amig':5.3234146919111531E-06, 'ff':1.407460, 'face':157.74},
            'MU2':{'amig':1.3559370068442651E-04, 'ff':0.967260, 'face':266.59},
            'MM':{'amig':2.6392030220989930E-06, 'ff':0.885420, 'face':211.10},
            'MS4':{'amig':2.8596300684150439E-04, 'ff':0.967260, 'face':313.79},
            'J1':{'amig':7.5560361379969530E-05, 'ff':1.152360, 'face':176.86},
            'SSA':{'amig':3.9821286769398282E-07, 'ff':1.000000, 'face':101.50}
        }

class BoundaryConditon(object):
    pass

if __name__=='__main__':
    # path = '/home/khan/MEGA/Models/SCHISM/Storm Surge/Mesh/02_Variable_Polder'
    # fname = 'Mesh_WGS84.grd'

    # hgrid = Hgrid()
    # hgrid.read(fname=fname, path=path)
    # coverage = hgrid.coverage()
    # print(coverage)

    # atlas = Atlas()
    # atlas.loadwaves(path='/run/media/khan/Workbench/Data/FES2012', pattern='{wave}_FES2012_SLEV.nc')
    # M2 = atlas.wave(wave='M2', coverage=coverage)
    # M2.print()

    tip = TidalPotential()
    tip.value(wavelist=['M2', 'S2', 'TP2'])