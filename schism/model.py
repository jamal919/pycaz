#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SCHISM Model Builder (SCHISMMB)

SCHISMMB is the container of the the classes and methods to build the inputs 
for SCHISM Modelling system. 

The model input is first build using the following objects - 
    Hgrid:      Horizontal grid
    Vgrid:      Vertical grid
    Bctide:     Boundary conditions
    Roughness:  Roughness conditions
    
These objects are passed to the Project object with a directory to create the 
files to run the model.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np
from core import Gr3

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

class Hgrid(Gr3):
    """
    SCHISM hgrid.gr3 object. This object loads the hgrid.gr3 file, breaks it in
    it's subcomponent, and store it for the model object to use.
    
    Grid or *.gr3 is essentially the fort.14 file format in ADCIRC model in SMS. 
    The intended way for this model is to provide the original file from SMS.
    """
    def __init__(self, path = None):
        self.path = path
        Gr3.__init__()
        if os.access(self.path, os.F_OK):
            self.read(fname=path)
        else:
            print('Check the hgrid.gr3 file location')

    def read_mesh(self):
        pass
        
    def read_boundaries(self):
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
        
    def plot(self):
        """ Plot the hgrid.gr3 mesh """
