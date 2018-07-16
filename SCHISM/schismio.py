#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SCHISM input output file library (schismio)

SCHISMIO is the container of the classes and methods to read and write SCHISM type files. 

Input/Output formats - 
1. BP format
2. gr3 format
3. netcdf format
4. station format

@author: Md. Jamal Uddin Khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np

class Boundary(object):
    """ Class for SCHISM complient boundary 
    There are in general two types of boundary - open and land. This class
    contains the information regarding a single boundary definition.
    
    
    """

    def __init__(self, bndno, bndnodes, bndtype=None, bndname=None):
        self.number = bndno
        self.nodes = bndnodes
        self.bndtype = bndtype
        self.name = bndname

    def nodecount(self):
        return(len(self.nodes))

class Boundaries(object):
    """ Collection of the object of class Boundary """

    def __init__(self, bndtype="open", totalnodes=None):
        self.bndtype = bndtype
        self.tatalnodes = totalnodes
        self.boundaries = []

    def addboundary(self, boundary):
        self.boundaries.append(boundary)
        
    def nopen(self):
        return(len(self.boundaries))

class Gr3(object):
    """ SCHISM .gr3 type object. This object can read and write gr3 like data. 
Gr3 file is used for many purposes in SCHISM environment. It can be used as Input file for mesh, nodal properties, boundaries etc. It is also used as output format for node centered value.

A full gr3 a have several components in the file. The components are - 
1. Nodal position and value information
2. Nodal connectivitiy information
3. Boundary and Boundary information

The selected information to read is controlled by passing flags to the functions.
"""
    def __init__(self, path=None, readflag=None):
        # path options
        self.path = path
        
        if os.access(self.path, os.F_OK):
            print('File found @ ' + self.path)
            self.findchunk()
        else:
            print('No file is found!')
            
    def findchunk(self):
        """ Read and find different chunk of the Gr3 """
        with open(self.path) as f:
            self.ds = f.readlines()
            self.filelength = len(self.ds)
            
            self.cline = 0
            self.grname = self.ds[self.cline].split()[0]
            print('Grid file:\t' + self.grname)
            
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
                    self.read.flag = (True, True, False)
                else:
                    self.boundds = self.ds[self.nnode+self.nelem+2:self.filelength]
                    print('Gr3 file contains - \n\t Nodal points\n\t Nodal connectivity\n\t Boundary')
                    self.readflag = (True, True, True)
    
    def readnodes(self):
        """ Extract the node information """
        self.dnodes = np.genfromtxt(fname=self.nodeds)
        print('Node informations reading successful... showing first 5 rows')
        print(self.dnodes[0:5,:])
        return(self.dnodes)
        
    def readelems(self):
        """ Extract the element information """
        self.delems = np.genfromtxt(fname=self.elemds)
        print('Element informations reading successful... showing first 5 rows')
        print(self.delems[0:5,:])
        
        return(self.delems)
        
    def gettriangulation(self):
        """ Return the triangulation from the element table
Triangulation is element table - 1 because of the python numbering starting from zero.        
        """
        self.triang = self.delems[:, 2:5]
        self.triang = self.triang - 1
        
        return(self.triang)
        
    def findboundaries(self):
        """ Separate the open and land boundaries """
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
        
    def readbounds(self):
        """ Extract the boundary informaiton """
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
            bndtype = int(ds[0].split()[1])
            bndnodes = np.genfromtxt(fname=ds[1:nlength+1])
            self.landbnd.addboundary(boundary=Boundary(bndno=bndno, bndnodes=bndnodes, bndtype=bndtype))
            ds = ds[nlength+1:len(ds)]
            print('Reading land boundary '+ str(i+1) + '... done.')
            
        return(self.openbnd, self.landbnd)
    
    def read(self):
        """ Read gr3 file format """
        with open(self.path) as f:
            self.ds = f.readlines()
            if self.readflag[0]:
                # Reading the nodal information
                # Saving the current line in cline variable
                self.dnodes = np.genfromtxt(fname=self.ds[self.cline:self.nnode+2])

                self.cline = self.cline + self.nnode

            if self.readflag[1]:
                # Reading the nodal connectivity
                self.delems = np.genfromtxt(fname=self.ds[self.cline:self.cline+self.nelem])

                self.triang = self.delems[:, 2:5]
                self.triang = self.triang - 1

                self.cline = self.cline + self.nelem

            if self.readflag[2]:
                # Reading the boundaries
                # Type 1 - Typically open boundary
                self.temp = self.ds[self.cline].split()

                self.nbound = self.temp[0]
                self.boundtype = self.temp[4]
                self.cline = self.cline + 1

                self.temp = self.ds[self.cline].split()
                self.cline = self.cline + 1
                self.openbnd = Boundaries(bndtype=self.boundtype, totalnodes=int(self.temp[0]))

                for i in range(int(self.nbound)):
                    self.temp = self.ds[self.cline].split()
                    self.bndnode = int(self.temp[0])
                    self.cline = self.cline + 1
                    self.bndnodes = np.genfromtxt(fname=self.ds[self.cline:self.cline+self.bndnode])
                    self.openbnd.addboundary(Boundary(i+1, self.bndnodes))
                    self.cline = self.cline + self.bndnode
                    print('Open boundary ' + str(i + 1))

                # Type 2 - Typically land boundary
                self.temp = self.ds[self.cline].split()
                
                self.nbound = self.temp[0]
                self.boundtype = self.temp[4]
                self.cline = self.cline + 1

                self.temp = self.ds[self.cline].split()
                self.cline = self.cline + 1
                self.landbnd = Boundaries(bndtype=self.boundtype, totalnodes=int(self.temp[0]))

                for i in range(int(self.nbound)):
                    print('Land boundary ' + str(i + 1) + 'starts!')
                    self.temp = self.ds[self.cline].split()
                    self.bndnode = int(self.temp[0])
                    self.cline = self.cline + 1
                    self.bndnodes = np.genfromtxt(fname=self.ds[self.cline:self.cline+self.bndnode])
                    self.landbnd.addboundary(Boundary(i+1, self.bndnodes))
                    self.cline = self.cline + self.bndnode
                    print('Land boundary ' + str(i + 1) + 'done!')


                        
                
