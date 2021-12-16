#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Data classes and functions to handle SCHISM pre- and post-processing.

The classes implemented here provides functionality to handle triangular mesh, 
rectangular mesh data. No timestamp is attached with the dataset.

@author: khan
"""
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from datetime import datetime, timedelta
from netCDF4 import Dataset
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
    def __init__(self, x, y, data=None):
        '''
        Grid object to generate grid and provides function to find various
        values at grid points.

        x: number of rows
        y: number of columns
        data: must be of the length len(x)*len(y). It will be flatten to get
        i,j,n formation
        '''
        try:
            assert isinstance(x, (list, tuple, np.ndarray))
            assert isinstance(y, (list, tuple, np.ndarray))
        except:
            raise TypeError(f'x, y must be python or numpy array')
        else:
            self.x = x
            self.y = y
            self.size = (len(self.x), len(self.y))

        # Data
        if data is None:
            self.depth = 1
            self.data = np.zeros(shape=(self.size[0], self.size[1], self.depth))
        elif isinstance(data, (np.ndarray)):
            self.depth = np.int(len(data.flatten())/self.size[0]/self.size[1])
            try:
                self.data = np.array(data).reshape((self.size[0], self.size[1], self.depth))
            except:
                raise Exception('Size mismatch')

    @property
    def meshgrid(self):
        X, Y = np.meshgrid(self.x, self.y, indexing='ij')
        return(X, Y)

    def reshape(self):
        '''
        Reshape the data to conform data structure.
        '''
        self.depth = np.int(len(self.data.flatten())/self.size[0]/self.size[1])
        self.data = self.data.reshape((self.size[0], self.size[1], self.depth))

    def __add__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data+other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data+other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data + other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] + other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __sub__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data - other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data - other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data - other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] - other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __mul__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data*other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data*other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data * other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] * other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __truediv__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data/other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data/other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data/other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i]/other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __lt__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data<other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i]<other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __le__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<=other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data<=other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data<=other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i]<=other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __gt__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data>other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data > other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] > other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __ge__(self, other):
        if isinstance(other, (Grid)):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data>=other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data>=other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return(
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data >= other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] >= other
                    return(
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __pow__(self, other):
        if isinstance(other, (float, int)):
            return(
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data**other
                )
            )
        else:
            raise ValueError('Only float or int as power')
    
    def __repr__(self):
        '''
        String representation.
        '''
        return(self.data.__repr__())

    def __getitem__(self, key):
        return self.data[key]

    def __getattr__(self, name):
        return getattr(self.data, name)

    def __iter__(self):
        '''
        Return a list to iterate over - i in object
        '''
        return(iter(self.data.reshape((self.size[0]*self.size[1], self.depth))))

    def apply(self, func, **kwargs):
        f = lambda x : func(x, **kwargs)

        data = np.array([f(x) for x in self])
                
        return(
            Grid(
                x=self.x,
                y=self.y,
                data=data
            )
        )

    def polar_coordinate(self, origin):
        '''
        Calculate the polar distance from a given point of interest.

        For lon,lat values, the distance is calculated using great circle distance.
        '''
        try:
            originx, originy = origin
        except:
            raise Exception('Origin must be a list of lon, lat')
        
        X, Y = self.meshgrid
        dfac = 60*1.852*1000
        dist_x = dfac*np.cos(np.deg2rad(Y))*(X-originx)
        dist_y = dfac*(Y-originy)
            
        r = np.sqrt(dist_x**2 + dist_y**2)
        theta = np.arctan2(dist_y, dist_x)

        return(
            Grid(
                x=self.x,
                y=self.y,
                data=np.array([(rr, tt) for rr, tt in zip(r.flatten(), theta.flatten())])
            )
        )

    def interpolate(self, at, depth=0, method='linear', fill_value=np.nan, rescale=False):
        '''
        Interpolate at another x,y point or grid using scipy.interpolate.griddata

        at: {list, tuple, Grid} instance
        depth: depth of grid data to interpolate
        method: {'linear', 'nearest', 'cubic'}, optional
        fill_value: value used to fill in for requested point outside of convex hull
        rescale: rescale points to unit cube before preforming interpolation

        return Grid
        '''
        X, Y = self.meshgrid

        points = np.array([(x, y) for x, y in zip(X.flatten(), Y.flatten())])
        values = self[:, :, depth].flatten()

        if isinstance(at, (list, tuple)):
            # For x, y list or tuple
            return(
                griddata(
                    points, values, 
                    at, 
                    method=method, 
                    fill_value=fill_value, 
                    rescale=rescale
                )
            )
        
        if isinstance(at, Grid):
            return(
                Grid(
                    x=at.x,
                    y=at.y,
                    data=griddata(
                        points, values, at.meshgrid,
                        method=method, 
                        fill_value=fill_value,
                        rescale=rescale
                    )
                )
            )
class Boundary(object):
    """ SCHISM complient boundary class
    There are in general two types of boundary - open and land. This class
    contains the information regarding a single boundary definition.
    
    Attributes:
        number(int)  :  Serial number of the boundary
        nodes(int []):  Numpy array of nodes forming the bounary
        bndtype(int) :  Boundary type in case of land bounary.
                        For open boundary, no information needed.
        bndname(str) :  Name of the boundary
    """

    def __init__(self, bndno, bndnodes, landflag=None, bndname=''):
        """ SCHISM complient bounary
        
        Args:
            bndno (int)     :   Serial number of the boundary
            bndnodes(int []):   Numpy array of nodes forming the boundary
            bndtype(int)    :   Boundary type in case of land boundary.
                                For open boundary, no information needed.
            bndname(str)    :   Name of the boundary (optional)
        """
        self.number = bndno
        self.nodes = bndnodes
        self.landflag = landflag
        self.name = bndname

    def countnodes(self):
        """Number of nodes in a boundary
        
        Returns:
            int : Number of nodes of the boundary
        """
        return(len(self.nodes))

class Boundaries(object):
    """Collection of Boundary Objects 
    
    Args:
        bndtype(str)    :  type of boundary (open or land)
        totalnodes(int) :  number of total nodes in the open boundary
    
    Returns:
        Instance of an empty Boundaries object.
    
    Methods:
        addboundary(Boundary boundary) : add a Boundary object to the touple
        nopen() : returns the number of open boundary added to the object
        
    TODO:
        * Add checktotalnodes method
    """

    def __init__(self, bndtype="open", totalnodes=None):
        self.bndtype = bndtype
        self.totalnodes = totalnodes
        self.boundaries = []

    def addboundary(self, boundary):
        """Add a new Boundary object
        
        Args:
            boundary(Boundary) : object of Boundary class to be added
        """
        self.boundaries.append(boundary)
        
    def count(self):
        """Number of boundary
        
        Returns:
            int : number of Boundary
        """
        return(len(self.boundaries))
        
class Global2Local(object):
    """
    Global2Local(path)
    
    This is the core class to read and hold global2local.prop files generated from SCHISM model. 

    args:
        path : path to the outputs folder of the schism model output

    returns:
        class instance of Global2Local object

    TODO:
        - change the structure to read the file directly
    
    """
    def __init__(self, path):
        self.path = os.path.join(path, 'global_to_local.prop')
        
    def load_global2local(self):
        self.mapping = np.loadtxt(fname=self.path, dtype='int32')
        return(self.mapping)

class Local2Global(object):
    """
    Local2Global(path)
    
    This is the core class to read and hold local_to_global_* files generated from SCHISM model.
    These files contains the subdomains of the models and mapping to the whole domain, thus 
    essential to merge the results back to the whole domain.

    args:
        path : path to the target local_to_global file

    returns:
        instance of Local2Global object

    Contains:
        - node mapping
        - element mapping
        - time

    methods:
        read_local2global(self) : read the given file

    TODO:
        - Fix the issue with reading the vertical coordinate information
        - In all cases checking the gfortran and ifort line limit convention
    """
    def __init__(self, path=None):
        self.path = path
        
    def read_local2global(self):
        with open(self.path) as f:
            ds = f.readlines()
            init = ds[0].split()
            self.globalside = int(init[0])
            self.globalelem = int(init[1])
            self.globalnode = int(init[2])
            self.nvrt = int(init[3])
            self.nproc = int(init[4])
            self.elemcount = int(ds[2].split()[0])
            self.elems = np.loadtxt(fname=ds[3:self.elemcount+3], dtype='int32')
            self.nodecount = int(ds[self.elemcount+3].split()[0])
            self.nodes = np.loadtxt(fname=ds[self.elemcount+4:self.elemcount+self.nodecount+4], dtype='int32')
            self.sidecount = int(ds[self.elemcount+self.nodecount+4])
            self.sides = np.loadtxt(fname=ds[self.elemcount+self.nodecount+5:self.elemcount+self.nodecount+self.sidecount+5], dtype='int32')
            # There is a difference between gcc-fortran and intel fortran. In intel fortran the value
            # is saved till 72 character and in gcc-fortran version the value is saved as requested.
            # As the critical part of the variables (i.e., time) can be extracted safely we are not
            # bothering about the rest of the variables. However, for robustness, the reading function
            # should be rewritten.
            # One way to achieve this is by actively keep looking for values till
            # expected number of values are found. This is probably the most
            # reasonable solution after using scipy.FortranFile.
            # TODO: Test number of values vs fortran file idea
            timestring = ds[self.elemcount+self.nodecount+self.sidecount+6].split()
            self.year = int(timestring[0])
            self.month = int(timestring[1])
            self.day = int(timestring[2])
            self.hour = float(timestring[3])
            self.minute = divmod(self.hour*60, 60)[1]
            self.hour = int(divmod(self.hour*60, 60)[0])
            self.second = int(divmod(self.minute*60, 60)[1])
            self.minute = int(divmod(self.minute*60, 60)[0])
            # self.gmt = float(ds[self.elemcount+self.nodecount+self.sidecount+7].split()[0])
            # model = ds[self.elemcount+self.nodecount+self.sidecount+7].split()
            # self.nrec = int(model[0])
            # self.dtout = float(model[1])
            # self.nspool = int(model[2])
            # self.nvrt = int(model[3])
            # self.kz = int(model[4])
            # self.h0 = float(model[5])
            # model = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            # self.h_s = float(model[0])
            # self.h_c = float(model[1])
            # self.theta_b = float(model[2])
            # self.theta_f = float(model[3])
            # self.ics = int(model[4])
            self.elemtable = np.loadtxt(fname=ds[len(ds)-self.elemcount:len(ds)], dtype='int16')
            self.nodetable = np.loadtxt(fname=ds[len(ds)-self.elemcount-self.nodecount:len(ds)-self.elemcount], dtype='float32')

class Hgrid(object):
    """ SCHISM .gr3 type object. 
    
    This class contains the methods and classes to handle Gr3 like data from
    SCHISM. Gr3 format is a structured text format and similar to ADCIRC fort.11
    file. 

    A full gr3 a have several components in the file. The components are - 
        1. Nodal position and value information
        2. Nodal connectivitiy information
        3. Boundary and Boundary information
        
    Args:
        name(str)       :   Name of the gr3 file
        nelem(int)      :   Number of elements
        nnode(int)      :   Number of nodes
        dnodes(int[])   :   Numpy array of nodes. [no] [x] [y] [value] [opt value]
        delems(int[])   :   Numpy array of element tables. [no] [#nodes] [node 1] [node 2] [node 3] [opt node4]
        openbnd(Boundaries) : Object of Boundaries class for open boundary sagments
        landbnd(Boundaries) : Object of Boundaries class for land boundary sagments
        
    Returns: 
        Object of class Gr3
        
    Methods:
        readfromfile(path) :
        readnodes():
        readelems():
        gettriangulation():
        findboundaries():
        readbounds():
        
    TODO:
        * Implement option for vortex more than 3 (i.e., triangular element)
    """
    def __init__(self, name=None, nelem=None, nnode=None, dnodes=None,
                 delems=None, openbnd=None, landbnd=None):
                     self.name = name
                     self.nelem = nelem
                     self.nnode = nnode
                     self.dnodes = dnodes
                     self.delems = delems
                     self.openbnd = openbnd
                     self.landbnd = landbnd
                     
                     # File loading related initialization
                     self.fileinitialized = False
            
    def readfromfile(self, path=None):
        """
        Check the existance of a given gr3 file and find the available 
        chunk information.
        """
        # path options
        self.path = path
        
        if os.access(self.path, os.F_OK):
            print('File found @ ' + self.path)
            self.findchunk()
            
            if self.readflag[0]:
                self.readnodes()
            if self.readflag[1]:
                self.readelems()
            if self.readflag[2]:
                self.readbounds()
        else:
            print('No file is found!')
            
    def findchunk(self):
        """ Read and find different chunk of the Gr3 """
        with open(self.path) as f:
            self.ds = f.readlines()
            self.filelength = len(self.ds)
            
            self.cline = 0
            self.name = ' '.join(self.ds[self.cline].split())
            print('Grid file:\t' + self.name)
            
            self.cline = 1
            self.nelem, self.nnode = self.ds[self.cline].split()
            self.nelem = int(self.nelem)
            self.nnode = int(self.nnode)
            self.cline = 2
            print('No of nodes:\t' + str(self.nnode))
            print('No of elems:\t' + str(self.nelem))
            
            self.nodeds = self.ds[2:self.nnode+2]
            if (len(self.nodeds) + 2) == self.filelength:
                print('Gr3 file contains - \n\t Nodal points')
                self.readflag = (True, False, False)
            else:
                self.elemds = self.ds[self.nnode+2:self.nnode+self.nelem+2]
                if (len(self.nodeds) + len(self.elemds) + 2) == self.filelength:
                    print('Gr3 file contains - \n\t Nodal points\n\t Nodal connectivity')
                    self.readflag = (True, True, False)
                else:
                    self.boundds = self.ds[self.nnode+self.nelem+2:self.filelength]
                    print('Gr3 file contains - \n\t Nodal points\n\t Nodal connectivity\n\t Boundary')
                    self.readflag = (True, True, True)
                    
        # File is initialized and ready to be read
        self.fileinitialized = True
    
    def readnodes(self):
        """ Extract the node information """
        if self.fileinitialized:
            self.dnodes = np.genfromtxt(fname=self.nodeds)
            print('Node informations reading successful... showing first 5 rows')
            print(self.dnodes[0:5,:])
            return(self.dnodes)
        else:
            print('File could not be initialized. Check the formatting.')
        
    def readelems(self):
        """ Extract the element information """
        if self.fileinitialized:
            self.delems = np.genfromtxt(fname=self.elemds)
            print('Element informations reading successful... showing first 5 rows')
            print(self.delems[0:5,:])
        
            return(self.delems)
        else:
            print('File could not be initialized. Check the formatting.')
        
    def gettriangulation(self):
        """ Triangulation from the element table
        
        Triangulation is calculated as element table - 1 because of the python 
        numbering starting from zero.
        
        Returns:
            Matplotlib complient triangulation
        """
        self.triang = self.delems[:, 2:5]
        self.triang = self.triang - 1
        
        return(self.triang)
        
    def findboundaries(self):
        """Separate the open and land boundaries """
        if self.fileinitialized:
            if len(self.boundds) <= 2:
                # First two line is checked to figure out the boundary type range
                print('Probably no boundary exists! Please check the file.')
            else:
                # Open boundary
                self.nopen = int(self.boundds[0].split()[0])
                self.nopennodes = int(self.boundds[1].split()[0])
                
                print('Open boundary :\t' + str(self.nopen))
                print('Open boundary nodes :\t' + str(self.nopennodes))
                
                self.openboundds = self.boundds[2:self.nopen+self.nopennodes+2]
                
                # Land boundary
                self.landboundds = self.boundds[self.nopen+self.nopennodes+2:len(self.boundds)]
                self.nland = int(self.landboundds[0].split()[0])
                self.nlandnodes = int(self.landboundds[1].split()[0])
                
                print('Land boundary :\t' + str(self.nland))
                print('Land boundary nodes :\t' + str(self.nlandnodes))
                
                self.landboundds = self.landboundds[2:self.nland+self.nlandnodes+2]
        else:
            print('File could not be initialized. Check the formatting.')
        
    def readbounds(self):
        """Extract the boundary informaiton """
        # Extracting the boundary sagments        
        self.findboundaries()
        
        # Open boundary
        ds = self.openboundds
        self.openbnd = Boundaries(bndtype='open', totalnodes=self.nopennodes)
        
        for i in range(self.nopen):
            bndno = i + 1
            nlength = int(ds[0].split()[0])
            bndnodes = np.genfromtxt(fname=ds[1:nlength+1])
            self.openbnd.addboundary(Boundary(bndno=bndno, bndnodes=bndnodes))
            ds = ds[nlength+1:len(ds)]
            print('Reading open boundary '+ str(i+1) + '... done.')
            
        # Land boundary
        ds = self.landboundds
        self.landbnd = Boundaries(bndtype='land', totalnodes=self.nlandnodes)
        
        for i in range(self.nland):
            bndno = i + 1
            nlength = int(ds[0].split()[0])
            landflag = int(ds[0].split()[1])
            bndnodes = np.genfromtxt(fname=ds[1:nlength+1])
            self.landbnd.addboundary(boundary=Boundary(bndno=bndno, bndnodes=bndnodes, landflag=landflag))
            ds = ds[nlength+1:len(ds)]
            print('Reading land boundary '+ str(i+1) + '... done.')
            
        return(self.openbnd, self.landbnd)
    
    def writetofile(self, path=None, overwrite=False, writebounds=True,
                    nodevalfmt='%16.10f'):
        """Write the grid to file
        
        This methods writes the grid to a file specified by path. The gr3
        format is implemented as appear in SCHISM manual. 
        
        Args:
            path(str)       :   path to the file to be written
            overwrite(bool) :   overwrite flag if file exist or same as input file
            writebounds(bool):  to write the boundary or not
            nodefmt         :   formatting string for the value at nodes
            
        TODO:
            * Check the grid object before writing for probable missing info
            * Add check for nodefmt
        """
        nodefmt=['%10i', '%16.10f', '%16.10f', nodevalfmt]        
        
        writefile = False
        if os.access(path, os.F_OK):
            if overwrite != True:
                print('File exists! Set overwrite=True to overwrite.')
            else:
                writefile = True
                self.outpath = path
                print('File is going to be over-written @ ' + self.outpath)
        else:
            writefile = True
            self.outpath = path
            print('File is going to be written @ ' + self.outpath)
        
        if writefile:
            # Header option
            with open(self.outpath, 'wb') as f:
                f.write(self.name + '\n')
                
            # Elements and nodes number
            with open(self.outpath, 'ab') as f:
                f.write(str(self.nelem) 
                + '\t' 
                + str(self.nnode) 
                + ' ! # of elements and nodes in the horizontal grid\n')
                
            # Nodes
            with open(self.outpath, 'ab') as f:
                np.savetxt(fname=f, X=self.dnodes, fmt=nodefmt)
                
            # Element table
            with open(self.outpath, 'ab') as f:
                np.savetxt(fname=f, X=self.delems, fmt='%i', delimiter='\t')
                
            # Boundary
            if writebounds:
                # open boundary
                with open(self.outpath, 'ab') as f:
                    f.write(str(self.openbnd.count()) + ' = Number of open boundaries\n')
                    f.write(str(self.openbnd.totalnodes) + ' = Total number of open boundary nodes\n')
                
                for boundary in self.openbnd.boundaries:
                    with open(self.outpath, 'ab') as f:
                        f.write(str(boundary.countnodes()) + ' = Number of nodes for open boundary ' + str(boundary.number) + ' ' + boundary.name +'\n')
                        np.savetxt(fname=f, X=boundary.nodes, fmt='%i')
                        
                # land boundary
                with open(self.outpath, 'ab') as f:
                    f.write(str(self.landbnd.count()) + ' = Number of land boundaries\n')
                    f.write(str(self.landbnd.totalnodes) + ' = Total number of land boundary nodes\n')
                
                for boundary in self.landbnd.boundaries:
                    with open(self.outpath, 'ab') as f:
                        f.write(str(boundary.countnodes()) + ' ' + str(boundary.landflag) + ' = Number of nodes for land boundary ' + str(boundary.number) + '\n')
                        np.savetxt(fname=f, X=boundary.nodes, fmt='%i')


class Tidefacout(object):
    def __init__(self, year=0, month=0, day=0, hour=0, rnday=0, const={}):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.rnday = rnday
        self.const = const
    
    def read(self, filepath):
        # Reading date information
        with open(filepath, 'r') as f:
            # Reading the date section
            __ds = f.readline()
            __date = np.fromstring(__ds, dtype=float, count=4, sep=',')
            self.year = __date[0]
            self.month = int(__date[1])
            self.day = int(__date[2]) 
            self.hour = int(__date[3])
            
            # Reading the run length section
            __ds = f.readline()
            __rnday = np.fromstring(__ds, dtype=float, count=1, sep=',')
            self.rnday = __rnday[0]

        # Reading the constants, node factor and eq. argument ref. to GM in deg.
        __const = np.genfromtxt(fname=filepath, dtype=None, skip_header=6, \
                                delimiter=None, autostrip=True)
        __const = np.array([[i for i in j] for j in __const])
        __const = {i[0].upper():[float(j) for j in i[1:3]] for i in __const}
        self.const = __const

        # Tidefac header information
        self.info = '{:.2f} days - {:4.0f}/{:02.0f}/{:02.0f} {:02.2f} UTC'.format(self.rnday,\
                self.year, self.month, self.day, self.hour)

    def __str__(self):
        return(self.info)
class Bctides(object):
    def __init__(self, info='', ntip=0, tip_dp=0, tip=[], nbfr=0, bfr=[], nope=0, boundaries=[]):
        self.info = info
        self.nitp = ntip
        self.tip_dp = tip_dp
        self.tip = tip
        self.nbfr = nbfr
        self.bfr = bfr
        self.nope = nope
        self.boundaries = boundaries

    def read(self, filepath):
        with open(filepath) as f:
            ds = f.readlines()
            # First the dates
            self.info = ds[0].split('\n')[0]
            __lnproc = 0

            # Then the tidal potential information
            self.ntip, self.tip_dp = np.fromstring(ds[1].split('!')[0], count=2, sep=' ')
            self.ntip = int(self.ntip)
            __lnproc = 1
            
            print('|{:3s}|{:10s}|{:10s}|{:10s}|{:10s}|{:10s}|{:10s}|'.format('no', 'talpha', 'jspc', 'tamp', 'tfreq', 'tnf', 'tear'))
            for i in np.arange(self.ntip):
                __talpha = ds[__lnproc+1].split('\n')[0]
                __jspc, __tamp, __tfreq, __tnf, __tear = np.fromstring(ds[__lnproc+2].split('\n')[0], count=5, sep=' ')
                print('|{:3d}|{:10s}|{:10d}|{:10f}|{:10f}|{:10f}|{:10f}|'.format(i, __talpha, __jspc, __tamp, __tfreq, __tnf, __tear))
                __rec = dict(talpha=__talpha, jspc=__jspc, tamp=__tamp, tfreq=__tfreq, tnf=__tnf, tear=__tear)
                self.tip.append(__rec)
                __lnproc = __lnproc + 2
            
            # Reading the boundary frequencies
            self.nbfr = np.fromstring(ds[__lnproc+1], count=1, sep=' ')
            self.nbfr = int(self.nbfr)
            __lnproc = __lnproc + 1
            
            self.bfr = []
            for i in np.arange(self.nbfr):
                __alpha = ds[__lnproc+1].split('\n')[0]
                __amig, __ff, __face = np.fromstring(ds[__lnproc+2].split('\n')[0], count=3, sep=' ')
                __rec = dict(alpha=__alpha, amig=__amig, ff=__ff, face=__face)
                self.bfr.append(__rec)
                __lnproc = __lnproc + 2
            
            # Open boundary sagments
            self.nope = ds[__lnproc+1].split(' ')[0]
            self.nope = int(self.nope)
            __lnproc = __lnproc + 1

            # For each open boundary sagment
            self.boundaries = ds[__lnproc+1:len(ds)]

    def update(self, tidefac):
        # Update time
        self.info = tidefac.info
        # Updating the tidal potential nodal factor and equilibrium argument
        for __tip in self.tip:
            __talpha = __tip['talpha'].strip().upper()
            if __talpha in tidefac.const.keys():
                __tip['tnf'] = tidefac.const[__talpha][0]
                __tip['tear'] = tidefac.const[__talpha][1]

        # Updating the Boundary frequency nodal factors and equilibrium argument
        for __bfr in self.bfr:
            __alpha = __bfr['alpha'].strip().upper()
            if __alpha in tidefac.const.keys():
                __bfr['ff'] = tidefac.const[__alpha][0]
                __bfr['face'] = tidefac.const[__alpha][1]

    def write(self, filepath):
        with open(filepath, 'w') as f:
            # Header information
            f.write('{:s}\n'.format(self.info))

            # Tidal potential
            f.write('{:d} {:3.2f} !ntip, tip_dp\n'.format(int(self.ntip), float(self.tip_dp)))

            for __tip in self.tip:
                f.write('{:s}\n{:d}\t{:.6f}\t{:.16f}\t{:.5f}\t{:.2f}\n'\
                        .format(__tip['talpha'].strip().upper(),\
                                int(__tip['jspc']),\
                                __tip['tamp'],\
                                __tip['tfreq'],\
                                __tip['tnf'],\
                                __tip['tear']))

            # Boundary frequencies
            f.write('{:d} !nbfr\n'.format(int(self.nbfr)))

            for __bfr in self.bfr:
                f.write('{:s}\n{:.16E}\t{:.6f}\t{:.2f}\n'\
                        .format(__bfr['alpha'].strip().upper(),\
                                __bfr['amig'],\
                                __bfr['ff'],\
                                __bfr['face']))

            # Open boundaries
            f.write('{:d} !Number of Open Boundaries\n'.format(self.nope))
            
            for __line in self.boundaries:
                f.write(__line)


class Sflux(object):
    def __init__(self, grid, basedate, sflux_type='air', nstep=96, priority=1, syncstep=10, path='./sflux'):
        self.grid = grid
        self.nstep = nstep # No of step
        self.basedate = basedate
        self.sflux_type = sflux_type
        self.nfile = 0 # No file at the beginning
        self.priority = priority # sflux_air_1 or 2
        self.syncstep = syncstep # Sync the netCDF each syncstep
        self.path = path

        sflux_func_map = {
            'air' : {
                'create_netcdf' : self.create_netcdf_air,
                'put_value' : self.put_value_air
            },
            'prc' : {
                'create_netcdf' : self.create_netcdf_prc,
                'put_value' : self.put_value_prc
            },
            'rad' : {
                'create_netcdf' : self.create_netcdf_rad,
                'put_value' : self.put_value_rad
            }
        }

        if sflux_type in sflux_func_map:
            self.create_netcdf = sflux_func_map[sflux_type]['create_netcdf']
            self.put_value = sflux_func_map[sflux_type]['put_value']
        else:
            raise Exception(f"sflux_type {self.sflux_type} not correct, one of 'air', 'prc', and 'rad'")
        
        # Directory creation
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def create_netcdf_air(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.filename = f'sflux_air_{self.priority:1d}.{self.nfile:03d}.nc'
        self.filepath = os.path.join(self.path, self.filename)

        # Creating the file first
        self.nc = Dataset(self.filepath, 'w', format='NETCDF4_CLASSIC')

        
        # Creating the dimensions
        self.nc.createDimension(dimname='nx_grid', size=len(self.grid.x))
        self.nc.createDimension(dimname='ny_grid', size=len(self.grid.y))
        self.nc.createDimension(dimname='ntime', size=None)

        # Creating the variables
        # Time
        self.v_time = self.nc.createVariable(
            varname='time',
            datatype=np.float32,
            dimensions=('ntime')
        )
        strf_basedate = self.basedate.strftime('%Y-%m-%d %H:%M:%S')
        self.v_time.units = f'days since {strf_basedate:s}'
        self.v_time.long_name = 'Time'
        self.v_time.calendar = 'standard'
        self.v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.v_lon = self.nc.createVariable(
            varname='lon',
            datatype=np.float32,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lon.units = 'degrees_north'
        self.v_lon.long_name = 'Longitude'
        self.v_lon.standard_name = 'longitude'

        # Latitude
        self.v_lat = self.nc.createVariable(
            varname='lat',
            datatype=np.float32,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lat.units = 'degrees_east'
        self.v_lat.long_name = 'Latitude'
        self.v_lat.standard_name = 'latitude'

        # Uwind
        self.v_uwind = self.nc.createVariable(
            varname='uwind',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.v_uwind.units = 'm/s'
        self.v_uwind.long_name = 'Surface Eastward Air Velocity (10m AGL)'
        self.v_uwind.standard_name = 'eastward_wind'

        # Vwind
        self.v_vwind = self.nc.createVariable(
            varname='vwind',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.v_vwind.units = 'm/s'
        self.v_vwind.long_name = 'Surface Northward Air Velocity (10m AGL)'
        self.v_vwind.standard_name = 'northward_wind'

        # Prmsl
        self.v_prmsl = self.nc.createVariable(
            varname='prmsl',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.v_prmsl.units = 'Pa'
        self.v_prmsl.long_name = 'Pressure Reduced to MSL'
        self.v_prmsl.standard_name = 'air_pressure_at_mean_sea_level'

        # stmp
        self.v_stmp = self.nc.createVariable(
            varname='stmp',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.v_stmp.units = 'K'
        self.v_stmp.long_name = 'Surface Temperature (2m AGL)'
        self.v_stmp.standard_name = 'surface_temperature'

        # spfh
        self.v_spfh = self.nc.createVariable(
            varname='spfh',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.v_spfh.units = 1
        self.v_spfh.long_name = 'Specific Humidity (2m AGL)'
        self.v_spfh.standard_name = 'surface_specific_humidity'
        
        # Writing lon-lat once
        X, Y = self.grid.meshgrid
        self.v_lon[:] = X.T
        self.v_lat[:] = Y.T
        
    
    def put_value_air(self, stepi, at, flux):
        if isinstance(at, (datetime, pd.DatetimeIndex)):
            at = pd.to_datetime(at) - self.basedate
        elif isinstance(at, (timedelta, pd.Timedelta)):
            at = at
        else:
            raise Exception(f'at must be datetime or timedelta object')

        self.v_time[stepi] = at.days + at.seconds/float(86400)
        self.v_uwind[stepi, :, :] = flux['uwind']
        self.v_vwind[stepi, :, :] = flux['vwind']
        self.v_prmsl[stepi, :, :] = flux['prmsl']
        self.v_stmp[stepi, :, :] = flux['stmp']
        self.v_spfh[stepi, :, :] = flux['spfh']

        self.step = self.step + 1

        # Syncing each 10 step
        if self.step%self.syncstep:
            self.nc.sync()

    def create_netcdf_prc(self):
        raise NotImplementedError

    def put_value_prc(self):
        raise NotImplementedError

    def create_netcdf_rad(self):
        raise NotImplementedError

    def put_value_rad(self):
        raise NotImplementedError

    def close_netcdf(self):
        self.nc.close()

    def write(self, at, flux):
        # First check if self.nc is available
        if hasattr(self, 'nc'):
            if self.step < self.nstep:
                self.put_value(self.step, at, flux)
            else:
                self.close_netcdf()
                self.create_netcdf()
                self.put_value(self.step, at, flux)
        else:
            self.create_netcdf()
            self.put_value(self.step, at, flux)

    def finish(self):
        if hasattr(self, 'nc'):
            self.close_netcdf()

    def sfluxtxt(self, dt):
        dt = dt.total_seconds()
        max_window = self.nstep*dt/float(3600)
        filepath = os.path.join(self.path, 'sflux_inputs.txt')
        with open(filepath, mode='w') as f:
            f.write('&sflux_inputs\n')
            f.write('air_1_relative_weight=1.,	!air_[12]_relative_weight set the relative ratio between datasets 1 and 2\n')
            f.write('air_2_relative_weight=99., \n')
            f.write(f'air_1_max_window_hours={max_window:.1f},	!max. # of hours (offset from start time in each file) in each file of set 1\n')
            f.write('air_1_fail_if_missing=.true.,	!set 1 is mandatory\n')
            f.write('air_2_fail_if_missing=.false., 	!set 2 is optional\n')
            f.write("air_1_file='sflux_air_1', 	!file name for 1st set of 'air'\n")
            f.write("air_2_file='sflux_air_2'\n")
            f.write("uwind_name='uwind', 		!name of u-wind vel.\n")
            f.write("vwind_name='vwind', 		!name of v-wind vel.\n")
            f.write("prmsl_name='prmsl', 		!name of air pressure (@MSL) variable in .nc file\n")
            f.write("stmp_name='stmp',  		!name of surface air T\n")
            f.write("spfh_name='spfh',  		!name of specific humidity\n")
            f.write('/\n')
        