#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sorting gr3 output from SCHISM in pythonic way. The script is intended to sort
a lot of max elevation results and sorting them at node based on the inundation
value.

This script was developed to tackle the problem of sorting the storm track
developed in Kerry Hydro simulations. As maxelev is a gr3 formatted file, a
generalized version is very much possible and is in the long list of future TODO

@license: GPL3
@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
import os
import glob
import sys

class Node(object):
    def __init__(self, id, x, y, z):
        self.id = id
        self.x = x
        self.y = y
        self.z = z

    def __lt__(self, other):
        if isinstance(other, Node):
            return(self.z < other.z)
        else:
            return(self.z < other)

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

class Gr3(object):
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

# Derived class made with Node and an Experiment name
class Pixel(object):
    def __init__(self, node, exp):
        self.id = node.id
        self.x = node.x
        self.y = node.y
        self.z = node.z
        self.exp = exp

    def __lt__(self, other):
        if isinstance(other, Pixel):
            return(self.z < other.z)
        else:
            return(self.z < other)

# Functions to sort the fullset or a subset
def gr3sort(fnames, consider='all'):
    try:
        # Check if the considered number of the files makes sense
        print('Total files to be sorted = {:d}'.format(len(fnames)))

        if consider != 'all':
            if(len(fnames)) <= consider:
                consider = 'all'
                print('The output will consider all ({:d}) values'.format(len(fnames)))
            else:
                print('The output will consider maximum {:d} values'.format(consider))
        else:
            print('The output will consider maximum {:d} values'.format(consider))

        # Loading first file and creating placeholder
        print('{:04d} - Reading {:s}'.format(0, os.path.basename(fnames[0])))
        exp = int(os.path.basename(fnames[0]).split('.gr3')[0].split('_')[1])
        gr3 = Gr3()
        gr3.read(fnames[0])
        gr3stack = np.array([Pixel(i, exp) for i in gr3.nodes])
        gr3stack = np.reshape(gr3stack, newshape=(1, len(gr3stack)))
    except:
        print('Problem with loading the first file! Exiting...')
        sys.exit(1)
    else:
        if consider == 'all':
            # Loading the rest of the files
            for i in np.arange(len(fnames))[1:len(fnames)]:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3nodes = np.array([Pixel(i, exp) for i in gr3.nodes])
                gr3stack = np.append(gr3stack, [gr3nodes], axis=0)

            # Sorting
            stackshape = gr3stack.shape
            gr3stack = np.sort(gr3stack, axis=0)
        elif consider > 1:
            # Loding upto first sagment of the files
            for i in np.arange(len(fnames))[1:consider]:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3nodes = np.array([Pixel(i, exp) for i in gr3.nodes])
                gr3stack = np.append(gr3stack, [gr3nodes], axis=0)

            # Initial sorting and setting the output shape of the stack
            stackshape = gr3stack.shape
            gr3stack = np.sort(gr3stack, axis=0)

            # Continue sorting the rest of the files
            for i in np.arange(len(fnames))[consider:len(fnames)]:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3nodes = np.array([Pixel(i, exp) for i in gr3.nodes])
                gr3stack = np.append(gr3stack, [gr3nodes], axis=0)
                gr3stack = np.sort(gr3stack, axis=0)
                gr3stack = gr3stack[1:consider+1, :]
    finally:
        # Saving the results
        sortelev = np.reshape([pixel.z for pixel in gr3stack.flatten()], stackshape)
        sortelev = np.flipud(sortelev)
        sortexp = np.reshape([pixel.exp for pixel in gr3stack.flatten()], stackshape)
        sortexp = np.flipud(sortexp)

        return(sortelev, sortexp)


if __name__=='__main__':
    folder = '/run/media/khan/Workbench/Projects/Surge Model/Kerry_Hydro/Maxelev'
    fnames = glob.glob(os.path.join(folder, 'Track_*.gr3'))
    
    sortelev, sortexp = gr3sort(fnames, consider=3)

    np.savetxt(fname=os.path.join(folder, 'sorted_elev.csv'), X=sortelev, fmt='%.3f', delimiter=',')
    np.savetxt(fname=os.path.join(folder, 'sorted_exp.csv'), X=sortexp, fmt='%d', delimiter=',')