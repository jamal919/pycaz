# -*- coding: utf-8 -*-
"""
Merge the netCDF output from SCHISM model. 

This is a standalone script to merge the netcdf output from SCHISM model for
elevation only. It is faster than the fortran script provided with the
source code.This source code is developed for use with single output and testing 
purpose.

It was developed as a part of SCHISM model toolbox module for python. For 
more information visit - github.com/schismmb

@license: GPL3
@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import glob
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num
import numpy as np
import os

class Global2Local(object):
    def __init__(self, path):
        self.path = os.path.join(path, 'outputs', 'global_to_local.prop')
        
    def load_global2local(self):
        self.mapping = np.loadtxt(fname=self.path, dtype='int32')
        return(self.mapping)

class Local2Global(object):
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
            timestring = ds[self.elemcount+self.nodecount+self.sidecount+6].split()
            self.year = int(timestring[0])
            self.month = int(timestring[1])
            self.day = int(timestring[2])
            self.hour = float(timestring[3])
            self.minute = divmod(self.hour*60, 60)[1]
            self.hour = int(divmod(self.hour*60, 60)[0])
            self.second = int(divmod(self.minute*60, 60)[1])
            self.minute = int(divmod(self.minute*60, 60)[0])
            self.gmt = float(ds[self.elemcount+self.nodecount+self.sidecount+7].split()[0])
            model = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            self.nrec = int(model[0])
            self.dtout = float(model[1])
            self.nspool = int(model[2])
            self.nvrt = int(model[3])
            self.kz = int(model[4])
            self.h0 = float(model[5])
            model = ds[self.elemcount+self.nodecount+self.sidecount+9].split()
            self.h_s = float(model[0])
            self.h_c = float(model[1])
            self.theta_b = float(model[2])
            self.theta_f = float(model[3])
            self.ics = int(model[4])
            self.elemtable = np.loadtxt(fname=ds[len(ds)-self.elemcount:len(ds)], dtype='int16')
            self.nodetable = np.loadtxt(fname=ds[len(ds)-self.elemcount-self.nodecount:len(ds)-self.elemcount], dtype='float32')
            
class Local2Globals(object):
    def __init__(self, path):
        self.path = path
    
    def load_files(self, prefix='local_to_global*'):
        self.filelist = glob.glob(os.path.join(self.path, 'outputs', prefix))
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
        nodenumber = np.array(np.arange(1, self.files[0].globalnode+1), dtype=int)
        self.globalnodex = np.empty(shape=(self.files[0].globalnode))
        self.globalnodey = np.empty(shape=(self.files[0].globalnode))
        self.globaldepth = np.empty(shape=(self.files[0].globalnode))
        for f in self.files:
            self.globalnodex[f.nodes[:, 1] - 1] = f.nodetable[:, 0]
            self.globalnodey[f.nodes[:, 1] - 1] = f.nodetable[:, 1]
            self.globaldepth[f.nodes[:, 1] - 1] = f.nodetable[:, 2]
        self.globalnodetable = np.column_stack((nodenumber, self.globalnodex, self.globalnodey, self.globaldepth))
        
            
class Schout(object):
    def __init__(self, path, local2globals, outfile='schout.nc'):
        self.path = path
        self.output = os.path.join(self.path, outfile)
        self.info = local2globals
    
    def list_inputs(self, inprefix='schout_*_', start=1, end=1):
        self.filelist = glob.glob(os.path.join(self.path, 'outputs', inprefix + '['+str(start) + '-' + str(end) + '].nc'))
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
        self.intertidal = False
        for proc in self.procs:
            infile = Dataset(self.filelist[proc])
            invalue = infile.variables['elev'][:]
            
            if self.intertidal:
                for i in np.arange(invalue.shape[1]):
                    if(np.sum(invalue[:, i]) >= 1):
                        invalue[:, i] = np.nan
                    else:
                        # do nothing
                        invalue[:, i] = invalue[:, i]

            outindex = self.info.files[proc].nodes - 1
            velev[:, outindex[:, 1]] = invalue[:, outindex[:, 0]]
            print(os.path.basename(self.filelist[proc]))
                
            infile.close()
            nc.sync()
        
        # Closing the output file
        nc.close()
                    

# Sample script
if __name__=='__main__':
    path = './'
    l2gs = Local2Globals(path)
    l2gs.load_files()
    l2gs.merge_nodes()
    nc = Schout(path=path, local2globals=l2gs)
    nc.list_inputs()
    nc.create_file()
