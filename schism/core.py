#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Core data abstraction class.

The classes implemented here provides functionality to handle triangular mesh, 
rectangular mesh data. No timestamp is attached with the dataset.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
import os

class Node(object):
    def __init__(self, id, x, y, z):
        self.id = id
        self.x = x
        self.y = y
        self.z = z

        Warning('Node class will be removed in the future')

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

        Warning('Element class will be removed in the future')

    def __lt__(self, other):
        if isinstance(other, Element):
            return(self.id < other.id)
        else:
            return(self.id < other)

class Gr3(object):
    def __init__(self, nodes=[], elems=[], data=[], epsg=4326):
        '''
        Gr3 data object.
        '''
        self.nodes = np.array(nodes)
        self.elems = np.array(elems)
        self.data = np.array(data)
        self.epsg = epsg

        self.nelem = len(self.elems)
        self.nnode = len(self.nodes)

    def read(self, fname):
        '''
        Calling the gr3 object to read the data directly.
        '''
        __file = fname
        with open(__file) as f:
            ds = f.readlines()

            __line = 0
            
            # Reading the gr3 name
            self.grname = ds[__line].strip()
            __line = __line + 1
            
            # Reading the number of nodes and elements
            __nelem, __nnode = np.fromstring(string=ds[__line].split('\n')[0], count=2, sep=' ')
            self.nelem = int(__nelem)
            self.nnode = int(__nnode)
            __line = __line + 1

            # Reading the nodes
            try:
                __nodes = np.genfromtxt(fname=ds[__line:__line+self.nnode])
            except:
                raise Exception('Node data error')
            else:
                __line = __line + self.nnode
                self.nodes = np.array(__nodes[:, 1:3])
                self.data = np.array(__nodes[:, 3:])

            if len(ds) >= self.nelem + self.nnode + 2:
                try:
                    # Reading the elements
                    __elems = np.genfromtxt(fname=ds[__line:__line+self.nelem], dtype=int)
                except:
                    raise Exception('Element table error')
                else:
                    __line = __line + self.nelem
                    self.elems = np.array(__elems, dtype=int)
                    self.elems = self.elems[:, 2:]
            else:
                Warning('Element table does not exist in the file.')
        
        return(self)

    @property
    def x(self):
        return(self.nodes[:, 0])

    @property
    def y(self):
        return(self.nodes[:, 1])

    def __add__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data + other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data + other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data + other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __sub__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data - other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data - other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data - other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __mul__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data * other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data * other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data * other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __truediv__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data / other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data/other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data/other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')


    def __gt__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data>other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data>other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data>other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __ge__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data>=other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data>=other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data>=other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __lt__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data<other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data<other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data<other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __le__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data<=other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data<=other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data<=other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __pow__(self, other):
        if isinstance(other, (Gr3)):
            try:
                assert np.all(other.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data**other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Gr3(
                    nodes=self.nodes,
                    elems=self.elems,
                    data=self.data**other
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert len(other) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=self.data**other
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def __abs__(self):
        '''
        Absolute value of the gr3 object
        '''
        return(
            Gr3(
                nodes=self.nodes,
                elems=self.elems,
                data=np.abs(self.data)
            )
        )

    def where(self, cond, other=np.nan):
        if isinstance(cond, (Gr3)):
            try:
                assert np.all(cond.nodes.shape == self.nodes.shape)
            except:
                raise ValueError('Uneuqal gr3 object')
            else:
                data = np.zeros_like(self.data)*other
                data[cond.data] = self.data[cond.data]
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=data
                    )
                )
        elif isinstance(cond, (np.array)):
            try:
                assert len(cond) == len(self.data)
            except:
                raise ValueError('Unequal data shape')
            else:
                data = np.zeros_like(self.data)*other
                data[cond] = self.data[cond]
                return(
                    Gr3(
                        nodes=self.nodes,
                        elems=self.elems,
                        data=data
                    )
                )
        else:
            raise ValueError('Not a Gr3, array, or a number!')

    def plot(self, ax=None, clevels=None, cmap=None, colorbar=False, subplot_kw=None, gridspec_kw=None, **fig_kw):
        '''
        plot data
        '''
        if ax is None:
            _, ax = plt.subplots(nrows=1, ncols=1, subplot_kw=subplot_kw, gridspec_kw=gridspec_kw, **fig_kw)

        if clevels is None:
            clevels = np.linspace(np.min(self.data), np.max(self.data), num=10)

        if cmap is None:
            cmap = 'jet'
        try:
            tc = ax.tricontourf(
                self.nodes[:, 0],
                self.nodes[:, 1],
                self.elems - 1, # python 0-based index
                self.data.flatten(),
                levels=clevels,
                cmap=cmap
            )
        except:
            raise Exception('Check data')
        
        if colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes(
                position="right", 
                size="5%", 
                pad=0.05, 
                axes_class=plt.Axes # Important if you use cartopy to create ax
            )
            plt.colorbar(tc, cax=cax)
        
        return(ax)

    def interpolate(self, at, method='linear'):
        pass

    def to_xarray(self):
        '''
        Return an xarray dataset
        '''
        pass

    def to_netcdf(self, fname):
        '''
        Save Gr3 data as netcdf file.
        '''
        pass

    def __repr__(self):
        '''
        String representation of Gr3 object.
        '''
        repr_string = f'''Gr3 object with {self.nnode} nodes, {self.nelem} elemenets, {self.data.shape[1]}  data columns'''
        return(repr_string)
            

class Grid(object):
    def __init__(self, x, y, data=None, epsg=4326, indexing='xy'):
        '''
        Grid object to generate grid and provides function to find various
        values at grid points. 
        '''
        try:
            assert isinstance(x, (list, tuple, np.ndarray))
            assert isinstance(y, (list, tuple, np.ndarray))
        except:
            raise TypeError(f'x, y must be python or numpy array')
        else:
            self.x = x
            self.y = y
            self.shape = (len(self.x), len(self.y))

        try:
            assert indexing in ('xy', 'ij')
        except:
            raise ValueError(f'indexing argument can be only xy and ij')
        else:
            self.indexing = indexing

        self.epsg = epsg

        # Creating the meshgrid
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing=self.indexing)
        self.shape = self.X.shape
        self.length = len(self.X.flatten())

        # Data
        if data is None:
            self.data = np.ones(shape=(self,x, self.y))*np.nan
        elif isinstance(data, (np.array)):
            try:
                assert np.all(data.shape == self.shape)
            except:
                raise ValueError('Inconsistent data and x,y size')

    def __add__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data+other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data+other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data+other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __sub__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data-other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data-other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data-other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __mul__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data*other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data*other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data*other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __truediv__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data/other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data/other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data/other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __lt__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data<other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __le__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<=other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data<=other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<=other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __gt__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data>other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __ge__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>=other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data>=other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>=other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')

    def __pow__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data**other.data,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data**other,
                    epsg=self.epsg,
                    indexing=self.indexing
                )
            )
        elif isinstance(other, (np.array)):
            try:
                assert np.all(other.shape == self.shape)
            except:
                raise ValueError('Unequal data shape')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data**other,
                        epsg=self.epsg,
                        indexing=self.indexing
                    )
                )
        else:
            raise ValueError('Not a Grid, array, or a number!')
    
    def __repr__(self):
        '''
        String representation.
        '''
        repr_string = f'''Grid object with size {self.shape}, and epsg {self.epsg}'''
        return(repr_string)


    def radial_distance(self, originx, originy):
        '''
        Calculates distance from a given point in degrees(origionx, originy) 
        and returns the (radial_distance, x_distance, y_distance)
        '''
        try:
            assert self.epsg == 4326
        except:
            NotImplementedError('Projected coordinate not implemented')
        else:
            dfac = 60*1.852*1000
            dist_x = dfac*np.cos(np.deg2rad(self.Y))*(self.X-originx)
            dist_y = dfac*(self.Y-originy)
            
            radial_distance = np.sqrt(dist_x**2 + dist_y**2)
        
        return(
            Grid(
                x=self.x,
                y=self.y,
                data=radial_distance,
                epsg=self.epsg,
                indexing=self.indexing
            ),
            Grid(
                x=self.x,
                y=self.y,
                data=dist_x,
                epsg=self.epsg,
                indexing=self.indexing
            ),
            Grid(
                x=self.x,
                y=self.y,
                data=dist_y,
                epsg=self.epsg,
                indexing=self.indexing
            )
        )

    def radial_quadrant(self, originx, originy):
        '''
        Calculates the quadrant location of radial positions.
        '''
        pass

    def interpolate(self, at, method='linear'):
        '''
        Interpolate at another x,y point or grid.
        '''
        pass

    def meshgrid(self, type='xy'):
        '''
        Return meshgrid of the x,y grid at xy, ij, or schism sflux format
        '''
        pass

    def transform(self, to_crs):
        '''
        Transformation of projection.
        '''
        pass