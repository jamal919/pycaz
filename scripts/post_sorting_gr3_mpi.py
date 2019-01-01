#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sorting gr3 output from SCHISM in pythonic way. The script is intended to sort
a lot of max elevation results and sorting them at node based on the inundation
value.

This script was developed to tackle the problem of sorting the storm track
developed in Kerry Hydro simulations.

Based on the experience of this scripts twin sister, who only use one core, we 
experienced a very slow performance in sorting. So Finally, this script is born
to use the power of MPI (message passing interface) to speedup the work.

MPI Implementation Note:
    The original post_sorting_gr3 was implemented on a premise that the objects 
    are during the reading of the files. Since we are only allowed to communiate
    basic datastructure with mpi, it was not possible to pass the complex data
    structure, like Node and Elements. So finally, they were streamlined to 
    use only the numeric values. The pixel objects were formed, essentially,
    afterwards in the mpi processes.

@license: GPL3
@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
from mpi4py import MPI
import os
import glob
import sys

# Gr3 Object to read and write the gr3 formatted files
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
            __nelem, __nnode = np.fromstring(string=ds[__line].split('\n')[0], count=2, sep=' ')
            self.nelem = int(__nelem)
            self.nnode = int(__nnode)
            __line = __line + 1

            if readnodes:
                # Reading the nodes
                __nodes = np.genfromtxt(fname=ds[__line:__line+self.nnode])
                __line = __line + self.nnode
                self.nodes = np.array(__nodes)

            if readelements:
                # Reading the elements
                __elems = np.genfromtxt(fname=ds[__line:__line+self.nelem], dtype=int)
                __line = __line + self.nelem
                self.elems = np.array(__elems, dtype=int)

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

# Class to hold pixel values and help sorting
class Pixel(object):
    def __init__(self, id, x, y, z, exp):
        self.id = id
        self.x = x
        self.y = y
        self.z = z
        self.exp = exp

    def __lt__(self, other):
        if isinstance(other, Pixel):
            return(self.z < other.z)
        else:
            return(self.z < other)


if __name__=='__main__':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # Setting up
    folder = '/run/media/khan/Workbench/Projects/Surge Model/Kerry_Hydro/Maxelev'
    fnames = glob.glob(os.path.join(folder, 'Track_*.gr3'))
    fnames = fnames[0:4] # For testing, comment out to consider all files
    consider = 3 # How many files to consider by all processes

    if rank == 0:
        # Check if the considered number of the files makes sense
        print('Total files to be sorted = {:d}'.format(len(fnames)))

    if consider != 'all':
        if(len(fnames)) <= consider:
            consider = 'all'
            if rank == 0:
                print('The output will consider all ({:d}) values'.format(len(fnames)))
        else:
            if rank == 0:
                print('The output will consider maximum {:d} values'.format(consider))
    else:
        if rank == 0:
            print('The output will consider maximum {:d} values'.format(consider))

    # Loading First file and distribute to other process
    if rank == 0:
        print('{:04d} - Reading {:s}'.format(0, os.path.basename(fnames[0])))
        exp = int(os.path.basename(fnames[0]).split('.gr3')[0].split('_')[1])
        gr3 = Gr3()
        gr3.read(fnames[0])
        gr3data = gr3.nodes
        gr3shape = gr3data.shape
    else:
        gr3shape = None
        exp = None

    gr3shape = comm.bcast(gr3shape, root=0)
    exp = comm.bcast(exp, root=0)
    if rank != 0:
        gr3data = np.empty(gr3shape)
    
    comm.Bcast(gr3data, root=0)

    # Range of points used by each rank and initial Point array
    # Similar to use a parmetis library
    chunksize = int(np.ceil(float(gr3shape[0])/size))
    gr3stack = np.array([Pixel(i[0], i[1], i[2], i[3], exp) for i in gr3data[rank*chunksize:(rank+1)*chunksize]])
    gr3stack = np.reshape(gr3stack, newshape=(1, len(gr3stack)))

    # Now running over the files
    if consider == 'all':
        # Load all the files
        for i in np.arange(len(fnames))[1:len(fnames)]:
            if rank == 0:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3data = gr3.nodes
                gr3shape = gr3data.shape
            else:
                gr3shape = None
                exp = None

            gr3shape = comm.bcast(gr3shape, root=0)
            exp = comm.bcast(exp, root=0)
            if rank != 0:
                gr3data = np.empty(gr3shape)

            comm.Bcast(gr3data, root=0)
            gr3append = np.array([Pixel(i[0], i[1], i[2], i[3], exp) for i in gr3data[rank*chunksize:(rank+1)*chunksize]])
            gr3stack = np.append(gr3stack, [gr3append], axis=0)

        stackshape = gr3stack.shape
        gr3stack = np.sort(gr3stack, axis=0)
    
    elif consider != 'all' and consider > 1:
        # Preparing matrix upto first sagment of the files
        for i in np.arange(len(fnames))[1:consider]:
            if rank == 0:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3data = gr3.nodes
                gr3shape = gr3data.shape
            else:
                gr3shape = None
                exp = None
            
            gr3shape = comm.bcast(gr3shape, root=0)
            exp = comm.bcast(exp, root=0)
            if rank != 0:
                gr3data = np.empty(gr3shape)

            comm.Bcast(gr3data, root=0)
            gr3append = np.array([Pixel(i[0], i[1], i[2], i[3], exp) for i in gr3data[rank*chunksize:(rank+1)*chunksize]])
            gr3stack = np.append(gr3stack, [gr3append], axis=0)

        stackshape = gr3stack.shape
        gr3stack = np.sort(gr3stack, axis=0)

        # Continue sorting the rest of the files
        for i in np.arange(len(fnames))[consider:len(fnames)]:
            if rank == 0:
                print('{:04d} - Reading {:s}'.format(i, os.path.basename(fnames[i])))
                exp = int(os.path.basename(fnames[i]).split('.gr3')[0].split('_')[1])
                gr3 = Gr3()
                gr3.read(fnames[i])
                gr3data = gr3.nodes
                gr3shape = gr3data.shape
            else:
                gr3shape = None
                exp = None

            gr3shape = comm.bcast(gr3shape, root=0)
            exp = comm.bcast(exp, root=0)
            if rank != 0:
                gr3data = np.empty(gr3shape)

            comm.Bcast(gr3data, root=0)
            gr3append = np.array([Pixel(i[0], i[1], i[2], i[3], exp) for i in gr3data[rank*chunksize:(rank+1)*chunksize]])
            gr3stack = np.append(gr3stack, [gr3append], axis=0)
            gr3stack = np.sort(gr3stack, axis=0)
            gr3stack = gr3stack[1:consider+1, :]

    # Saving files with flipping the axis
    sortelev = np.reshape([pixel.z for pixel in gr3stack.flatten()], stackshape)
    sortelev = np.flipud(sortelev)
    sortexp = np.reshape([pixel.exp for pixel in gr3stack.flatten()], stackshape)
    sortexp = np.flipud(sortexp)

    np.savetxt(fname=os.path.join(folder, 'sorted_elev_'+str(rank)+'.csv'), X=sortelev, fmt='%.3f', delimiter=',')
    np.savetxt(fname=os.path.join(folder, 'sorted_exp_'+str(rank)+'.csv'), X=sortexp, fmt='%d', delimiter=',')