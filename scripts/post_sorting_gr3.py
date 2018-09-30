# -*- coding: utf-8 -*-
"""
Sorting gr3 output from SCHISM in pythonic way. The script is intended to sort
a lot of max elevation results and sorting them at node based on the inundation
value.

This script was developed to tacke the problem of sorting the storm track
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


if __name__=='__main__':
    folder = '/run/media/khan/Workbench/Projects/Surge Model/Kerry_Hydro/Maxelev'
    fnames = glob.glob(os.path.join(folder, 'Track_*.gr3'))
    gr3nodes = []
    gr3tnames = []
    for i in np.arange(len(fnames))[0:10]:
        print('Reading - {:s}'.format(os.path.basename(fnames[i])))
        trackname = os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1]
        gr3 = Gr3()
        gr3.read(fnames[i])
        tname = np.repeat(trackname, len(gr3.nodes))
        gr3nodes.append(gr3.nodes)
        gr3tnames.append(tname)

    gr3nodes = np.array(gr3nodes)
    gr3tnames = np.array(gr3tnames)

    # sortindex = np.argsort(gr3nodes, axis=0)

    gr3nodeval = np.reshape(np.array([i.z for i in gr3nodes.flat]), gr3nodes.shape)
    # np.savetxt(fname=os.path.join(folder, 'elev.txt'), X=gr3nodeval, delimiter=',')
    # np.savetxt(fname=os.path.join(folder, 'tracknumber.txt'), X=gr3tnames, delimiter=',', fmt='%s')
    
    
    gr3sorted = np.sort(gr3nodes, axis=0)
    gr3sortedval = np.reshape(np.array([i.z for i in gr3sorted.flat]), gr3nodes.shape)
    
    gr3final = Gr3(grname='maximum', nelem=gr3.nelem, nnode=gr3.nnode, nodes=gr3sorted[len(gr3sorted)-1, :], elems=gr3.elems)
    gr3final.write(fname='maximum.gr3', path=folder)