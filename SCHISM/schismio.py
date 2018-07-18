#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SCHISM input output file library (schismio)

This module contains the input output classes for SCHISM model. The Input/Output 
in SCHISM can be classfied in following formats: 
    1. BP format
    2. gr3 format
    3. netcdf format
    4. station format

Classes:
    * Gr3
    
TODO:
    * add station file read write functionality
    * add bctides.in read write functionality

@author: Md. Jamal Uddin Khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np
from schismcore import Boundary, Boundaries

class Gr3(object):
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
        * Implement rectangular element
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
                    self.read.flag = (True, True, False)
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

if __name__=='__main__':
    print('import this library using\n>>> from SCHISM import schismio')