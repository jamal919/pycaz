#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script is a standalone script to set nodal values using shape or tiff files.
It is also possible to set rule based values, but it is probably simpler to 
approach this with generic mesh class itself. 

@license: GPL3
@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import numpy as np
import os


class Node(object):
    def __init__(self, id, x, y, z):
        self.id = id
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

class Coverage(object):
    """ Coverage object to hold the information regarding spatially varying
    values to be assign to the nodes.

    TODO: Setup the data structure
    """
    def __init__(self, kind=None, fname=None):
        self.kind = kind
        self.fname = fname

        if self.kind == 'shp':
            self.readshape(self.fname)
        elif self.kind == 'tif' or self.kind == 'tiff':
            self.readtiff(self.fname)
        else:
            print('Incompatible file format. Only shapefile and tiff is supported.')

    def readshape(self, fname):
        pass

    def readtiff(self, fname):
        pass