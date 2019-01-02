# -*- coding: utf-8 -*-
"""
Merge the netCDF output from SCHISM model. 

This is a standalone script to merge the netcdf output from SCHISM model for
elevation only. It is faster than the fortran script provided with the
source code. However, it can only combine water level (elev) in a comodo compatible
format, which allows the use of comodo tools for tidal analysis. 

This source code is developed for use with single output and testing 
purpose. It was developed as a part of SCHISM model toolbox module for python. For 
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
            # There is a difference between gcc-fortran and intel fortran. In intel fortran the value
            # is saved till 72 character and in gcc-fortran version the value is saved as requested.
            # As the critical part of the variables (i.e., time) can be extracted safely we are not
            # bothering about the rest of the variables. However, for robustness, the reading function
            # must be rewritten.
            # TODO Rewrite the module using scipy.io.FortranFile
            timestring = ds[self.elemcount+self.nodecount+self.sidecount+6].split()
            self.year = int(timestring[0])
            self.month = int(timestring[1])
            self.day = int(timestring[2])
            self.hour = float(timestring[3])
            self.minute = divmod(self.hour*60, 60)[1]
            self.hour = int(divmod(self.hour*60, 60)[0])
            self.second = int(divmod(self.minute*60, 60)[1])
            self.minute = int(divmod(self.minute*60, 60)[0])
            #self.gmt = float(ds[self.elemcount+self.nodecount+self.sidecount+7].split()[0])
            #modelstring = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            #self.nrec = int(modelstring[0])
            #self.dtout = float(modelstring[1])
            #self.nspool = int(modelstring[2])
            #self.nvrt = int(modelstring[3])
            #self.kz = int(modelstring[4])
            #self.h0 = float(modelstring[5])
            #model = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            #self.h_s = float(model[0])
            #self.h_c = float(model[1])
            #self.theta_b = float(model[2])
            #self.theta_f = float(model[3])
            #self.ics = int(model[4])
            self.elemtable = np.loadtxt(fname=ds[len(ds)-self.elemcount:len(ds)], dtype='int16')
            self.nodetable = np.loadtxt(fname=ds[len(ds)-self.elemcount-self.nodecount:len(ds)-self.elemcount], dtype='float32')
            
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
        nodenumber = np.array(np.arange(1, self.files[0].globalnode+1), dtype=int)
        self.globalnodex = np.empty(shape=(self.files[0].globalnode))
        self.globalnodey = np.empty(shape=(self.files[0].globalnode))
        self.globaldepth = np.empty(shape=(self.files[0].globalnode))
        for f in self.files:
            self.globalnodex[f.nodes[:, 1] - 1] = f.nodetable[:, 0]
            self.globalnodey[f.nodes[:, 1] - 1] = f.nodetable[:, 1]
            self.globaldepth[f.nodes[:, 1] - 1] = f.nodetable[:, 2]
        self.globalnodetable = np.column_stack((nodenumber, self.globalnodex, self.globalnodey, self.globaldepth))

    def merge_elements(self, vortex=3):
        self.globalelemtable = np.empty(shape=(self.files[0].globalelem, vortex+2))
        self.globalelemtable[:, 0] = np.array(np.arange(1,self.files[0].globalelem+1), dtype=int)
        for f in self.files:
            self.globalelemtable[f.elems[:, 1]-1, 1] = f.elemtable[:, 0]
            for i in np.arange(vortex):
                self.globalelemtable[f.elems[:, 1]-1, i+2] = f.nodes[f.elemtable[:, i+1]-1, 1]
        
            
class Schout(object):
    def __init__(self, path, local2globals):
        self.path = path
        self.info = local2globals
    
    def list_inputs(self, inprefix='schout_*_', start=1, end=1):
        self.filelist = glob.glob(os.path.join(self.path, inprefix + '['+str(start) + '-' + str(end) + '].nc'))
        self.filelist = sorted(self.filelist)
        self.procs = np.array([int(os.path.basename(i).split('_')[1]) for i in self.filelist])
        self.outputs = np.array([Dataset(fname) for fname in self.filelist])
        self.timesteps = self.outputs[0].variables['time'][:]
        
    def create_file(self, outfile='schout.nc', path='./'):
        self.output = os.path.join(path, outfile)

        # Create netcdf file
        nc = Dataset(self.output, 'w', format='NETCDF4', clobber=True)       
        
        # Creating dimensions
        nc.createDimension(dimname='M', size=self.info.files[0].globalelem)
        nc.createDimension(dimname='N', size=self.info.files[0].globalnode)
        nc.createDimension(dimname='P', size=3)
        nc.createDimension(dimname='T', size=None)
        
        # Variables
        ## Time variable
        vtime = nc.createVariable(varname='time', datatype=np.float64, dimensions=('T'))
        vtime.longname = 'Time elapsed in seconds since time_origin'
        vtime.units = 'seconds'
        vtime.calendar = 'standard'
        vtime.title = 'Time'
        vtime.time_origin = '1970-01-01 00:00:00'
        
        timearray = [datetime(self.info.files[0].year, self.info.files[0].month, self.info.files[0].day, self.info.files[0].hour, self.info.files[0].minute, self.info.files[0].second) + timedelta(seconds=t) for t in [int(i) for i in self.timesteps]]
        vtime[:] = date2num(timearray, units = '{unit} since {origin}'.format(unit=vtime.units, origin=vtime.time_origin), calendar=vtime.calendar)        
        
        ## Longitude
        vx = nc.createVariable(varname='lon', datatype=np.float64, dimensions=('N'))
        vx.longname = 'longitude'
        vx.standard_name = 'longitude'
        vx.units = 'degrees_east'
        vx.valid_min = -180.0
        vx.valid_max = 180.0
        vx.subgrid ='point'
        vx.content = 'N'
        vx[:] = self.info.globalnodex
        
        ## Latitude
        vy = nc.createVariable(varname='lat', datatype=np.float64, dimensions=('N'))
        vy.longname = 'latitude'
        vy.standard_name = 'latitude'
        vy.units = 'degrees_north'
        vy.valid_min = -90.0
        vy.valid_max = 90.0
        vx.subgrid = 'point'
        vy.content = 'N'
        vy[:] = self.info.globalnodey

        ## Element table
        velem = nc.createVariable(varname='element', datatype=int, dimensions=('M', 'P'))
        velem.long_name = 'element_connectivity'
        velem.standard_name = 'element'
        velem.subgrid = 'cell'
        velem.content = 'MP'
        velem[:] = self.info.globalelemtable[:, 2:5]
        
        ## Bathymetry
        vdepth = nc.createVariable(varname='bathymetry', datatype=np.float64, dimensions=('N'))
        vdepth.longname = 'model_positive_bathymetry'
        vdepth.short_name = 'bathym'
        vdepth.units = 'm'
        vdepth.associate = 'lon lat'
        vdepth.positive = 'down'
        vdepth.location = 'node'
        vdepth.subgrid = 'point'
        vdepth.content = 'N'
        vdepth[:] = self.info.globalnodetable[:, 3]
        
        ## Water level
        velev = nc.createVariable(varname='elevation', datatype=np.float64, dimensions=('T', 'N'), chunksizes=(1, self.info.files[0].globalnode))
        velev.units = 'm'
        velev.standard_name = 'sea_surface_height_above_mean_sea_level'
        velev.short_name = 'elevation'
        velev.data_horizontal_center = 'node'
        velev.vertical_center = 'full'
        velev.content = 'TN'
        velev.associate = 'time lat lon'

        ## Pulling and saving variables from sagmented netCDF files to velev
        for i, timestep in enumerate(self.timesteps):
            invalues = np.empty(self.info.files[0].globalnode, dtype=float)

            for f, output in enumerate(self.outputs):
                invalue = output.variables['elev'][i, :]
                outindex = self.info.files[f].nodes - 1
                invalues[outindex[:, 1]] = invalue
                
            velev[i, :] = invalues
            print(i, timestep)

            if i%200 == 0:
                nc.sync()
        
        ## Global Attirbute
        nc.title = 'SCHISM model output'
        nc.institution = 'LEGOS'
        nc.source = 'SCHISM'
        nc.Conventions = 'CF 1.0'
        nc.Extensions = 'SIROCCO'
        nc.history = 'created by python netcdf library'
        nc.author = 'Jamal Uddin Khan'
        nc.email = 'jamal.khan@legos.obs-mip.fr'

        # Closing the merged file
        nc.close()

        # Closing the output files
        for output in self.outputs:
            output.close()
                    

# Sample script
if __name__=='__main__':
    path = '/run/media/khan/Storehouse/Projects/Tidal Modelling/Experiments/EXP07_2010_Sa_Boundary/outputs/outputs'
    l2gs = Local2Globals(path)
    l2gs.load_files()
    l2gs.merge_nodes()
    l2gs.merge_elements()
    nc = Schout(path=path, local2globals=l2gs)
    nc.list_inputs()
    nc.create_file(outfile='schout.nc', path='/run/media/khan/Workbench/Projects/Tide/Experiments/EXP07')