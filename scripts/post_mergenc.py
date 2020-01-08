#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Merge the netCDF output from SCHISM model.

This is a standalone script to merge the netcdf output from SCHISM model for
elevation only. It is faster than the fortran script provided with the
source code.This source code is developed for use with single output and
testing purpose.

It was developed as a part of SCHISM model toolbox module for python.
For more information visit - github.com/jamal919/pyschism

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
        '''
        path: str, path where the global_to_local file exists, typically 'outputs'
        '''
        self.path = os.path.join(path, 'global_to_local.prop')

    def load_global2local(self):
        '''
        loads global to local file
        '''
        self.mapping = np.loadtxt(fname=self.path, dtype='int32')
        return(self.mapping)

    
class Local2Global(object):
    def __init__(self, path=None, compiler='intel'):
        '''
        path: str, path to a local_to_global file
        compiler: str, flag to identify output from which compiler gnu or intel
        '''
        self.path = path

    def read_local2global(self):
        '''
        reading the local to global file

        There is a difference between gcc-fortran and intel fortran. In intel 
        fortran the value is saved till 72 character and in gcc-fortran version 
        the value is saved as requested. As the critical part of the variables 
        (i.e., time) can be extracted safely we are not bothering about the rest
        of the variables. However, for completeness, the reading function must be 
        rewritten.

        Currently a compiler flag is used to circumvent this issue.
        '''
        with open(self.path) as f:
            ds = f.readlines()
            init = ds[0].split()
            self.globalside = int(init[0])
            self.globalelem = int(init[1])
            self.globalnode = int(init[2])
            self.nvrt = int(init[3])
            self.nproc = int(init[4])
            self.elemcount = int(ds[2].split()[0])
            self.elems = np.loadtxt(fname=ds[3:self.elemcount+3],
                                    dtype='int32')
            self.nodecount = int(ds[self.elemcount+3].split()[0])
            self.nodes = np.loadtxt(fname=ds[self.elemcount+4:self.elemcount+self.nodecount+4], dtype='int32')
            self.sidecount = int(ds[self.elemcount+self.nodecount+4])
            self.sides = np.loadtxt(fname=ds[self.elemcount+self.nodecount+5:self.elemcount+self.nodecount+self.sidecount+5], dtype='int32')
            timestring = ds[self.elemcount+self.nodecount+self.sidecount+6].split()
            self.year = int(timestring[0])
            self.month = int(timestring[1])
            self.day = int(timestring[2])
            self.hour_model = float(timestring[3])
            self.minute = divmod(self.hour*60, 60)[1]
            self.hour = int(divmod(self.hour_model*60, 60)[0])
            self.second = int(divmod(self.minute*60, 60)[1])
            self.minute = int(divmod(self.minute*60, 60)[0])
            self.utc = float(ds[self.elemcount+self.nodecount+self.sidecount+7].split()[0])
            modelstring = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            self.nrec = int(modelstring[0])
            self.dtout = float(modelstring[1])
            self.nspool = int(modelstring[2])
            self.nvrt = int(modelstring[3])
            self.kz = int(modelstring[4])
            self.h0 = float(modelstring[5])
            vrt_s = ds[self.elemcount+self.nodecount+self.sidecount+8].split()
            self.h_s = float(vrt_s[0])
            self.h_c = float(vrt_s[1])
            self.theta_b = float(vrt_s[2])
            self.theta_f = float(vrt_s[3])
            self.ics = int(vrt_s[4])
            # afterwards vertical level definition
            self.elemtable = np.loadtxt(fname=ds[len(ds)-self.elemcount:len(ds)],
                                        dtype='int16')
            self.nodetable = np.loadtxt(fname=ds[len(ds)-self.elemcount-self.nodecount:len(ds)-self.elemcount],
                                        dtype='float32')


class Local2Globals(object):
    def __init__(self, path, prefix='local_to_global*', compiler='intel'):
        '''
        path: str, path to a bunch of local_to_global files
        compiler: str, type of compiler used in model, gnu or intel
                    this is to solve the issue to 72 char line and contineous
                    line in reading local_to_global files
        '''
        self.path = path
        self.load_files(prefix=prefix, compiler=compiler)
        self.merge_nodes()
        self.merge_elements()
        self.merge_edges()
    
    @property
    def year(self):
        return(self.files[0].year)

    @property
    def month(self):
        return(self.files[0].month)

    @property
    def day(self):
        return(self.files[0].day)

    @property
    def hour(self):
        return(self.files[0].hour) # Processed time
    
    @property
    def hour_model(self):
        return(self.files[0].hour_model) # As appear in param.in

    @property
    def minute(self):
        return(self.files[0].minute)

    @property
    def second(self):
        return(self.files[0].second)

    @property
    def utc(self):
        return(self.files[0].utc)

    @property
    def nglobalnode(self):
        '''
        Global node count from local_to_global file header
        '''
        return(self.files[0].globalnode)

    @property
    def nglobalelem(self):
        '''
        Global element count from local_to_global file header
        '''
        return(self.files[0].globalelem)

    @property
    def nglobaledge(self):
        '''
        Global edge count from local_to_global file header
        '''
        return(self.files[0].globalside)

    @property
    def nvrt(self):
        '''
        Number of vertical layer defined in local_to_global file header
        '''
        return(self.files[0].nvrt)

    @property
    def ics(self):
        return(self.files[0].ics)

    @property
    def h0(self):
        return(self.files[0].h0)

    @property
    def h_s(self):
        return(self.files[0].h_s)

    @property
    def h_c(self):
        return(self.files[0].h_c)

    @property
    def theta_b(self):
        return(self.files[0].theta_b)

    @property
    def theta_f(self):
        return(self.files[0].theta_f)

    @property
    def kz(self):
        return(self.files[0].kz)

    def load_files(self, prefix='local_to_global*', compiler='intel'):
        '''
        prefix: str, file perfix to list and process, default local_to_global*
        compiler: str, output from gnu or intel compilers
        '''
        self.filelist = glob.glob(os.path.join(self.path, prefix))
        self.filelist = sorted(self.filelist)

        self.files = []

        for f in self.filelist:
            local2global = Local2Global(path=f, compiler=compiler)
            local2global.read_local2global()
            self.files.append(local2global)

        if(len(self.files)) != self.files[0].nproc:
            print('Mismatch! expected != obtained local_to_global files!')
            raise(Exception)


    def merge_nodes(self):
        '''
        method - merge the nodes in all local to global files to global nodes
        '''
        self.globalnodex = np.empty(shape=(self.nglobalnode))
        self.globalnodey = np.empty(shape=(self.nglobalnode))
        self.globaldepth = np.empty(shape=(self.nglobalnode))

        for f in self.files:
            self.globalnodex[f.nodes[:, 1] - 1] = f.nodetable[:, 0]
            self.globalnodey[f.nodes[:, 1] - 1] = f.nodetable[:, 1]
            self.globaldepth[f.nodes[:, 1] - 1] = f.nodetable[:, 2]

    def merge_elements(self):
        '''
        method to merge the elements from all local_to_global files
        '''
        self.globalfacenodes = np.empty(shape=(self.nglobalelem, 4))*np.nan
        self.globalfacenodex = np.empty(shape=(self.nglobalelem))
        self.globalfacenodey = np.empty(shape=(self.nglobalelem))
        
        for f in self.files:
            for ilocal, iglobal in zip(f.elems[:, 0]-1, f.elems[:, 1]-1):
                # -1 is important, 0-based indexing
                if f.elemtable[ilocal, 0] == 3:
                    # triangular element
                    self.globalfacenodes[iglobal, 0:3] = f.elemtable[ilocal, 1:]
                else:
                    # rectangular element
                    self.globalfacenodes[iglobal, 0:4] = f.elemtable[ilocal, 1:]

    def merge_edges(self):
        '''
        method to merge the edges from all local_to_global files
        '''
        pass



class Schout(object):
    def __init__(self, path, local2globals, inprefix='schout', ispool=1, outfile='schout.nc'):
        '''
        path: str, path to location of schout files, typically 'outputs'
        local2globals: Local2Globals, instance of object Local2Globals
        '''
        self.path = path
        self.info = local2globals
        self.list_inputs(prefix=inprefix, ispool=ispool)
        self.create_file(outfile=outfile)

    def list_inputs(self, prefix='schout', ispool=1):
        '''
        prefix: str, prefix for the file, default: 'schout'
        ispool: int, spool number to combine, default: 1st spool
        '''
        file_pattern = f"{prefix}_*_{ispool:d}.nc"
        self.filelist = glob.glob(os.path.join(self.path, file_pattern))
        self.filelist = sorted(self.filelist)
        self.procs = [ 
            int(os.path.basename(i).split('_')[1]) for i in self.filelist
        ]

        # TODO to be used nrec from local_to_global_* later
        nc = Dataset(self.filelist[0])
        self.records = nc.variables['time'][:]
        nc.close()

    def create_file(self, outfile='schout.nc'):
        '''
        outfile: str, filename of the output file
        path: str, path to store the output file
        '''
        self.fname = outfile

        # Create netcdf file
        self.nc = Dataset(self.fname, 'w', format='NETCDF4', clobber=True)
        self.file_closed = False

        # Creating dimensions
        self.nc.createDimension(dimname='nSCHISM_hgrid_node',
                                size=self.info.nglobalnode)
        self.nc.createDimension(dimname='nSCHISM_hgrid_face',
                                size=self.info.nglobalelem)
        self.nc.createDimension(dimname='nSCHISM_hgrid_edge',
                                size=self.info.nglobaledge)
        self.nc.createDimension(dimname='nMaxSCHISM_hgrid_face_nodes',
                                size=4)
        self.nc.createDimension(dimname='nSCHISM_vgrid_layers',
                                size=self.info.nvrt)
        self.nc.createDimension(dimname='one', size=1)
        self.nc.createDimension(dimname='two', size=2)
        self.nc.createDimension(dimname='sigma',
                                size=self.info.nvrt-self.info.kz+1)
        if self.info.kz != 1:
            self.nc.createDimension(dimname='nz', size=self.info.kz-1)

        self.nc.createDimension(dimname='time', size=None)

        # time: Time variable
        _year = self.info.year
        _month = self.info.month
        _day = self.info.day
        _hour = self.info.hour
        _hour_model = self.info.hour_model
        _minute = self.info.minute
        _second = self.info.second
        _utc = self.info.utc
        vtime = self.nc.createVariable(varname='time',
                                       datatype=np.float64,
                                       dimensions=('time'))
        vtime.long_name = 'Time'
        vtime.units = f'seconds since {_year:4d}-{_month:02d}-{_day:02d} {_hour:02d}:{_minute:02d}:{_second:02d} {_utc*100:+04d}'
        vtime.calendar = 'standard'
        vtime.base_date = [_year, _month, _day, _hour_model, _utc]
        # timearray = [datetime(self.info.year,
        #                       self.info.month,
        #                       self.info.day,
        #                       self.info.hour,
        #                       self.info.minute,
        #                       self.info.second) + timedelta(seconds=t)
        #              for t in [int(i) for i in self.records]]
        # vtime[:] = date2num(timearray,
        #                     units=vtime.units,
        #                     calendar=vtime.calendar)
        vtime[:] = self.records

        # SCHISM_hgrid
        vhgrid = self.nc.createVariable(varname='SCHISM_hgrid',
                                        datatype=np.int32,
                                        dimensions=('one'))
        vhgrid.long_name = 'Topology data of 2d unstructured mesh'
        vhgrid.topology_dimensions = 2
        vhgrid.cf_role = 'mesh_topology'
        vhgrid.node_coordinates = 'SCHISM_hgrid_node_x SCHISM_hgrid_node_y'
        vhgrid.face_node_connectivity = 'SCHISM_hgrid_face_nodes'
        vhgrid.edge_coordinates = 'SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y'
        vhgrid.face_coordinates = 'SCHISM_hgrid_face_x SCHISM_hgrid_face_y'
        vhgrid.edge_node_connectivity = 'SCHISM_hgrid_edge_nodes'

        # SCHISM_hgrid_face_nodes
        vfacenodes = self.nc.createVariable(varname='SCHISM_hgrid_face_nodes',
                                            datatype=np.int32,
                                            dimensions=('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'))
        vfacenodes.long_name = 'Horizontal Element Table'
        vfacenodes.cf_role = 'face_node_connectivity'
        vfacenodes.start_index = 1
        vfacenodes[:] = self.info.globalfacenodes

        # SCHISM_hgrid_edge_nodes
        vedgenodes = self.nc.createVariable(varname='SCHISM_hgrid_edge_nodes',
                                            datatype=np.int32,
                                            dimensions=('nSCHISM_hgrid_edge', 'two'))
        vedgenodes.long_name = 'Map every edge to the two nodes that it connects'
        vedgenodes.cf_roles = 'edge_node_connectivity'
        vedgenodes.start_index = 1

        # SCHISM_hgrid_node_x
        vx = self.nc.createVariable(varname='SCHISM_hgrid_node_x',
                                    datatype=np.float32,
                                    dimensions=('nSCHISM_hgrid_node'))
        vx.long_name = 'Node x-coordinate'
        vx.standard_name = 'longitude'
        vx.units = 'degrees_east'
        vx[:] = self.info.globalnodex

        # SCHISM_hgrid_node_y
        vy = self.nc.createVariable(varname='SCHISM_hgrid_node_y',
                                    datatype=np.float32,
                                    dimensions=('nSCHISM_hgrid_node'))
        vy.long_name = 'Node y-coordinate'
        vy.standard_name = 'latitude'
        vy.units = 'degrees_north'
        vy[:] = self.info.globalnodey

        # node_bottom_index
        vnodebottomindex = self.nc.createVariable(varname='node_bottom_index',
                                                  datatype=np.int32,
                                                  dimensions=('nSCHISM_hgrid_node'))
        vnodebottomindex.long_name = 'x_coordinate of 2D mesh face'
        vnodebottomindex.units = 'non-dimensional'
        vnodebottomindex.mesh = 'SCHISM_hgrid'
        vnodebottomindex.location = 'node'
        vnodebottomindex.start_index = 1

        # SCHISM_hgrid_face_x
        vfacex = self.nc.createVariable(varname='SCHISM_hgrid_face_x',
                                        datatype=np.float32,
                                        dimensions=('nSCHISM_hgrid_face'))
        vfacex.long_name = 'x_coordinate of 2D mesh face'
        vfacex.standard_name = 'longitude'
        vfacex.units = 'degrees_east'
        vfacex.mesh = 'SCHISM_hgrid'

        # SCHISM_hgrid_face_y
        vfacey = self.nc.createVariable(varname='SCHISM_hgrid_face_y',
                                        datatype=np.float32,
                                        dimensions=('nSCHISM_hgrid_face'))
        vfacey.long_name = 'y_coordinate of 2D mesh face'
        vfacey.standard_name = 'latitude'
        vfacey.units = 'degrees_north'
        vfacey.mesh = 'SCHISM_hgrid'

        # ele_bottom_index
        velebottomindex = self.nc.createVariable(varname='ele_bottom_index',
                                                 datatype=np.int32,
                                                 dimensions=('nSCHISM_hgrid_face'))
        velebottomindex.long_name = 'bottom level index at each element'
        velebottomindex.units = 'non-dimensional'
        velebottomindex.mesh = 'SCHISM_hgrid'
        velebottomindex.location = 'elem'
        velebottomindex.start_index = 1

        # SCHISM_hgrid_edge_x
        vedgex = self.nc.createVariable(varname='SCHISM_hgrid_edge_x',
                                        datatype=np.float32,
                                        dimensions=('nSCHISM_hgrid_edge'))
        vedgex.long_name = 'x_coordinate of 2D mesh edge'
        vedgex.standard_name = 'longitude'
        vedgex.units = 'degrees_east'
        vedgex.mesh = 'SCHISM_hgrid'

        # SCHISM_hgrid_edge_y
        vedgey = self.nc.createVariable(varname='SCHISM_hgrid_edge_y',
                                        datatype=np.float32,
                                        dimensions=('nSCHISM_hgrid_edge'))
        vedgey.long_name = 'y_coordinate of 2D mesh edge'
        vedgey.standard_name = 'latitude'
        vedgey.units = 'degrees_north'
        vedgey.mesh = 'SCHISM_hgrid'

        # edge_bottom_index
        vedgebottomindex = self.nc.createVariable(varname='edge_bottom_index',
                                                datatype=np.int32,
                                                dimensions=('nSCHISM_hgrid_edge'))
        vedgebottomindex.long_name = 'bottom level index at each edge'
        vedgebottomindex.units = 'non-dimensional'
        vedgebottomindex.mesh = 'SCHISM_hgrid'
        vedgebottomindex.location = 'edge'
        vedgebottomindex.start_index = 1

        # depth
        vdepth = self.nc.createVariable(varname='depth', 
                                        datatype=np.float32, 
                                        dimensions=('nSCHISM_hgrid_node'))
        vdepth.long_name = 'Bathymetry'
        vdepth.units = 'm'
        vdepth.positive = 'down'
        vdepth.mesh = 'SCHISM_hgrid'
        vdepth.location = 'node'
        vdepth[:] = self.info.globaldepth

        # sigma
        vsigma = self.nc.createVariable(varname='sigma',
                                        datatype=np.float32,
                                        dimensions=('sigma'))
        vsigma.long_name = 'S coordinates at whole levels'
        vsigma.units = 1
        vsigma.standard_name = 'ocean_s_coordinate'
        vsigma.positive = 'up'
        vsigma.h_s = self.info.h_s
        vsigma.h_c = self.info.h_c
        vsigma.theta_b = self.info.theta_b
        vsigma.theta_f = self.info.theta_b

        # dry_value_flag
        # coordinae_system_flag
        vcoordflag = self.nc.createVariable(varname='coordinate_system_flag',
                                            datatype=np.int32,
                                            dimensions=('one'))
        vcoordflag[:] = self.info.ics

        # minimum_depth
        vmindepth = self.nc.createVariable(varname='minimum_depth',
                                            datatype=np.float32,
                                            dimensions=('one'))
        vmindepth[:] = self.info.h0
        
        # sigma_h_c
        vsigmahc = self.nc.createVariable(varname='sigma_h_c',
                                        datatype=np.float32,
                                        dimensions=('one'))
        vsigmahc[:] = self.info.h_c

        # sigma_theta_b
        vsigmatb = self.nc.createVariable(varname='sigma_theta_b',
                                        datatype=np.float32,
                                        dimensions=('one'))
        vsigmatb[:] = self.info.theta_b

        # sigma_theta_f
        vsigmatf = self.nc.createVariable(varname='sigma_theta_f',
                                        datatype=np.float32,
                                        dimensions=('one'))
        vsigmatf[:] = self.info.theta_f

        # sigma_maxdepth
        # Cs
        # wetdry_elem
        # zcor

        # Global Attirbute
        self.nc.conventions = 'CF-1.0, UGRID-1.0'
        self.nc.title = 'Merged SCHISM model output'
        self.nc.institution = 'LEGOS'
        self.nc.source = 'SCHISM v5.6.1 WWM Hybrid'
        self.nc.history = 'created by python netcdf library'
        self.nc.author = 'Jamal Uddin Khan'
        self.nc.email = 'jamal.khan@legos.obs-mip.fr'
        self.nc.VisIT_plugin = 'https://schism.water.ca.gov/library/-/document_library/view/3476283'

        self.nc.sync()

    def combine(self, varname, rename=None, datatype=np.float32, long_name=None, units=None, chunksizes=None, **kwargs):
        '''
        varname: str, name of the variable to be merged
        rename: str, renamed variable
        long_name: str, long_name of the variable
        units: str, custom units, if None then will try to get original units
        chunksizes: dict, dimname:size, if None then no chunking, size='full'
        '''
        # determine variable name
        in_varname = varname
        if rename is not None:
            out_varname = rename
        else:
            out_varname = in_varname

        # for checking if there is 'full' dimension
        if chunksizes is not None:
            fulldims = {}
            for d in self.nc.dimensions:
                fulldims[d] = self.nc.dimensions[d].size

        # check if variable exists and get its attributes
        tempfile = Dataset(self.filelist[0])
        if in_varname in tempfile.variables.keys():
            print(f"{in_varname} : processing as - {out_varname}")

            dims = {}
            for v in tempfile.variables[in_varname].get_dims():
                # setting dimname and 1 as chunking
                dims[v.name] = 1

            if chunksizes is not None:
                # update chunking info
                for dim in dims:
                    if dim in chunksizes:
                        if chunksizes[dim] == 'full':
                            dims[dim] = fulldims[dim]
                        else:
                            dims[dim] = chunksizes[dim]
                    else:
                        print(f"{dim} not in chunksizes, defaulting to 1")

            attrs = {}
            for attr in tempfile.variables[in_varname].ncattrs():
                attrs[attr] = tempfile.variables[in_varname].getncattr(attr)

            tempfile.close()

            # defining the variable
            out_dims = tuple(d for d in dims)
            out_chunks = tuple(dims[d] for d in dims)
            save_var = self.nc.createVariable(varname=out_varname,
                                              datatype=datatype,
                                              dimensions=out_dims,
                                              chunksizes=out_chunks)
            if long_name is not None:
                save_var.long_name = long_name

            if units is not None:
                save_var.units = units

            for attr in attrs:
                self.nc.variables[out_varname].setncattr(attr, attrs[attr])
            
            for kwarg in kwargs:
                self.nc.variables[out_varname].setncattr(kwarg, kwargs[kwarg])

            # Pulling and saving variables from sagmented netCDF files
            for proc in self.procs:
                infile = Dataset(self.filelist[proc])
                invalue = infile.variables[in_varname][:]

                outindex = self.info.files[proc].nodes - 1
                if len(out_dims) == 2:
                    save_var[:, outindex[:, 1]] = invalue[:, outindex[:, 0]]
                elif len(out_dims) == 3:
                    save_var[:, outindex[:, 1], :] = invalue[:, outindex[:, 0], :]

                print(os.path.basename(self.filelist[proc]))
                infile.close()
                self.nc.sync()
        else:
            print(f"variable {in_varname} does not exist!")

    def close_file(self):
        # Closing the output file
        # Make it automatic in case of destruction
        self.nc.close()
        self.file_closed = True

    def __del__(self):
        # If the object is deleted then calling close_file automatically
        # Should be called automatically at the end of the program
        if not self.file_closed:
            self.nc.close()


# Sample script
if __name__=='__main__':
    path = './outputs'
    l2gs = Local2Globals(path, prefix='local_to_global*', compiler='intel')
    nc = Schout(path=path, local2globals=l2gs, inprefix='schout', ispool=1, outfile='schout.nc')
    nc.combine(varname='elev', long_name='water level elevation', units='m', chunksizes={'time':'full'})
    nc.combine(varname='dahv', long_name='depth average velocity', units='m/s', chunksizes={'time':'full'})
    nc.combine(varname='WWM_1', rename='hs', long_name='significant wave height', units='m', chunksizes={'time':'full'})
    nc.combine(varname='WWM_2', rename='tm01', long_name='mean average period', units='s', chunksizes={'time':'full'})
    nc.combine(varname='WWM_3', rename='tm02', long_name='zero down crossing period', units='s', chunksizes={'time':'full'})
    nc.combine(varname='WWM_11', rename='tp', long_name='discrete peak period', units='s', chunksizes={'time':'full'})
    nc.combine(varname='WWM_energy_dir', long_name='energy direction vector', chunksizes={'time':'full'})
    nc.close_file()