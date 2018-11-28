#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SCHISM file input/output routines

This module contains the input output classes for SCHISM model. The Input/Output 
in SCHISM can be classfied in following formats: 
    1. BP format
    2. gr3 format (implemented)
    3. netcdf format (implemented)
    4. station format

TODO:
    * add station file read write functionality
    * add bctides.in read write functionality
    * make the path used by various classes clean

@author: Md. Jamal Uddin Khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np
import glob
from netCDF4 import Dataset, date2num
from datetime import datetime, timedelta
from core import Boundary, Boundaries, Local2Global

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

class Hgrid(object):
    pass

class Local2Globals(object):
    def __init__(self, path):
        self.path = path
    
    def load_files(self, prefix='local_to_global*'):
        self.filelist = glob.glob(os.path.join(self.path, prefix))
        self.filelist = sorted(self.filelist)
        
        self.files = []
        
        print('Loading local_to_global files...')
        for f in self.filelist:
            local2global = Local2Global(path=f)
            local2global.read_local2global()
            self.files.append(local2global)
            
        if(len(self.files)) == self.files[0].nproc:
            print('All local_to_global files are loaded!')
        else:
            print('Mismatch between number of expected and obtained local_to_global files!')
            raise(Exception)

    def merge_nodes(self):
        self.globalnodetable = np.empty(shape=(self.files[0].globalnode, 4))
        self.globalnodetable[:, 0] = np.array(np.arange(1, self.files[0].globalnode+1), dtype=int)
        for f in self.files:
            self.globalnodetable[f.nodes[:, 1] - 1, 1] = f.nodetable[:, 0]
            self.globalnodetable[f.nodes[:, 1] - 1, 2] = f.nodetable[:, 1]
            self.globalnodetable[f.nodes[:, 1] - 1, 3] = f.nodetable[:, 2]

    def merge_elements(self, vortex=3):
        self.globalelemtable = np.empty(shape=(self.files[0].globalelem, vortex+2))
        self.globalelemtable[:, 0] = np.array(np.arange(1,self.files[0].globalelem+1), dtype=int)
        for f in self.files:
            self.globalelemtable[f.elems[:, 1]-1, 1] = f.elemtable[:, 0]
            for i in np.arange(vortex):
                self.globalelemtable[f.elems[:, 1]-1, i+2] = f.nodes[f.elemtable[:, i+1]-1, 1]


class Schout(object):
    def __init__(self, path, local2globals, outfile='schout.nc', outpath='./'):
        self.path = path
        self.output = os.path.join(outpath, outfile)
        self.info = local2globals
    
    def list_inputs(self, inprefix='schout_*_', start=1, end=1):
        self.filelist = glob.glob(os.path.join(self.path, inprefix + '['+str(start) + '-' + str(end) + '].nc'))
        self.filelist = sorted(self.filelist)
        self.procs = [int(os.path.basename(i).split('_')[1]) for i in self.filelist]
        
        nc = Dataset(self.filelist[0])
        self.records = nc.variables['time'][:]
        nc.close()        
        
    def create_file(self):
        # Create netcdf file
        nc = Dataset(self.output, 'w', format='NETCDF4', clobber=True)       
        
        # Creating dimensions
        nc.createDimension(dimname='nSCHISM_hgrid_node', size=self.info.files[0].globalnode)
        nc.createDimension(dimname='time', size=None)
        
        # Variables
        vtime = nc.createVariable(varname='time', datatype=np.float64, dimensions=('time'))
        vtime.longname = 'Time'
        vtime.units = 'hours since 1970-01-01 00:00:00'
        vtime.calendar = 'standard'
        vtime.standard_name = 'Time hours'
        
        timearray = [datetime(self.info.files[0].year, self.info.files[0].month, self.info.files[0].day, self.info.files[0].hour, self.info.files[0].minute, self.info.files[0].second) + timedelta(seconds=t) for t in [int(i) for i in self.records]]
        vtime[:] = date2num(timearray, units = vtime.units, calendar=vtime.calendar)        
        
        vx = nc.createVariable(varname='SCHISM_hgrid_node_x', datatype=np.float64, dimensions=('nSCHISM_hgrid_node'))
        vx.longname = 'Node x-coordinate'
        vx.standard_name = 'longitude'
        vx.units = 'degrees_east'
        vx[:] = self.info.globalnodex
        
        vy = nc.createVariable(varname='SCHISM_hgrid_node_y', datatype=np.float64, dimensions=('nSCHISM_hgrid_node'))
        vy.longname = 'Node y-coordinate'
        vy.standard_name = 'latitude'
        vy.units = 'degrees_north'
        vy[:] = self.info.globalnodey
        
        vdepth = nc.createVariable(varname='depth', datatype=np.float64, dimensions=('nSCHISM_hgrid_node'))
        vdepth.longname = 'Bathymetry'
        vdepth.units = 'meter'
        vdepth.positive = 'down'
        vdepth.location = 'node'
        vdepth[:] = self.info.globalnodetable[:, 3]
        
        velev = nc.createVariable(varname='elev', datatype=np.float64, dimensions=('time', 'nSCHISM_hgrid_node'), chunksizes=(len(self.records), 1))
        velev.units = 'meter'
        velev.data_horizontal_center = 'node'
        velev.vertical_center = 'full'
        
        # Global Attirbute
        nc.title = 'Merged SCHISM model output'
        nc.institution = 'LEGOS'
        nc.source = 'SCHISM'
        nc.history = 'created by python netcdf library'
        nc.author = 'Jamal Uddin Khan'
        nc.email = 'jamal.khan@legos.obs-mip.fr'
        
        # Pulling and saving variables from sagmented netCDF files
        for proc in self.procs:
            infile = Dataset(self.filelist[proc])
            invalue = infile.variables['elev'][:]
            
            outindex = self.info.files[proc].nodes - 1
            velev[:, outindex[:, 1]] = invalue[:, outindex[:, 0]] 
            print(os.path.basename(self.filelist[proc]))
                
            infile.close()
            nc.sync()
        
        # Closing the output file
        nc.close()

class MaxVariable(object):
    def __init__(self, l2g, varname='maxelev', varnum=1):
        self.varname = varname
        self.varnum = varnum

        self.l2g = l2g
        self.path = self.l2g.path
        self.nnode = self.l2g.files[0].globalnode
        self.nelem = self.l2g.files[0].globalelem
        self.nodes = self.l2g.globalnodetable
        self.elems = self.l2g.globalelemtable

        # Setting the output to nan values
        for i in np.arange(0, self.varnum):
            self.nodes[:, i+3] = np.nan

        self.prefix = self.varname + '_*'
        

    def list_files(self):
        self.filelist = glob.glob(os.path.join(self.path, self.prefix))
        self.filelist = sorted(self.filelist)

    def read_file(self, path):
        with open(path) as f:
            __ds = f.readlines()
            __ds = __ds[1:len(__ds)]
            __table = np.loadtxt(__ds)
            return(__table)

    def merge_files(self):
        for f in self.filelist:
            __ds = self.read_file(f)
            __index = __ds[:, 0] - 1
            __index = [int(i) for i in __index]
            for i in np.arange(0, self.varnum):
                self.nodes[__index, i+3] = __ds[:, i+3]

    def write_file(self, path, fmt='%16.10f'):
        nodefmt = ['%10i', '%16.10f', '%16.10f']
        valfmt = np.repeat(fmt, self.varnum)
        for f in valfmt:
            nodefmt.append(f)

        with open(path, mode='wb') as f:
            f.write(self.varname + '\n')
            f.write(str(self.nelem) 
                + '\t' 
                + str(self.nnode) 
                + ' ! # of elements and nodes in the horizontal grid\n')
            np.savetxt(fname=f, X=self.nodes, fmt=nodefmt)
            np.savetxt(fname=f, X=self.elems, fmt='%i', delimiter='\t')

if __name__=='__main__':
    print('import this library using\n>>> from schism import io as schismio')