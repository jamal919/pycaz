#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Combines the maxelev.gr3 files into a single netcdf with corresponding experiment
information.

It produces a netCDF file, similar to the structure of SHCISM output. 
The mesh structure is derived from the input hgrid.gr3 file.

Author: Jamal Uddin Khan
Email: jamal.khan@legos.obs-mip.fr
"""

import numpy as np
import netCDF4
import os
from glob import glob

# Gr3 Object to read and write the gr3 formatted files
class Gr3(object):
    def __init__(self, grname=None, nelem=0, nnode=0, nodes=[], elems=[]):
        self.grname = grname
        self.nelem = nelem
        self.nnode = nnode
        self.nodes = nodes
        self.elems = elems

    def read(self, fname, readnodes=True, readelements=True):
        __file = fname
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

if __name__=='__main__':
    # Inputs
    in_folder = './maxelev'
    in_files = glob(os.path.join(in_folder, 'Track_*.gr3'))

    mesh_file = './hgrid.gr3'
    mesh_data = Gr3()
    mesh_data.read(fname=mesh_file, readnodes=True, readelements=True)

    n_experiment = len(in_files)
    print(f'No of experiment to merge: {n_experiment:04d}')
    print(f'Depth file : {mesh_file:s}')

    # netCDF file
    try:
        nc = netCDF4.Dataset('maxelev.nc', 'w', format='NETCDF4', clobber=True)
        
        nc.createDimension(dimname='nSCHISM_hgrid_node', size=mesh_data.nnode)
        nc.createDimension(dimname='nSCHISM_hgrid_face', size=mesh_data.nelem)
        nc.createDimension(dimname='nMaxSCHISM_hgrid_face_nodes', size=4)
        nc.createDimension(dimname='n_experiment', size=n_experiment)

        var_node_x = nc.createVariable(varname='SCHISM_hgrid_nodex', datatype=float, dimensions=('nSCHISM_hgrid_node'))
        var_node_x.long_name = 'Node x-coordinate'
        var_node_x.standard_name = 'Longitude'
        var_node_x.units = 'degrees east'
        var_node_x[:] = mesh_data.nodes[:, 1]

        var_node_y = nc.createVariable(varname='SCHISM_hgrid_nodey', datatype=float, dimensions=('nSCHISM_hgrid_node'))
        var_node_y.long_name = 'Node y-coordinate'
        var_node_y.standard_name = 'Latitude'
        var_node_y.units = 'degrees north'
        var_node_y[:] = mesh_data.nodes[:, 2]

        var_elem = nc.createVariable(varname='SCHISM_hgrid_face_nodes', datatype=int, dimensions=('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'))
        var_elem.long_name = 'Element Connectivity Table'
        var_elem[:, 0:3] = mesh_data.elems[:, 1:4]
        
        var_depth = nc.createVariable(varname='depth', datatype=float, dimensions=('nSCHISM_hgrid_node'))
        var_depth.long_name = 'Bathymetry at node'
        var_depth.standard_name = 'Depth'
        var_depth.convention = 'Positive downard'
        var_depth[:] = mesh_data.nodes[:, 3]

        var_exps = nc.createVariable(varname='experiment', datatype=int, dimensions=('n_experiment'))
        var_exps.long_name = 'Track no of the experiment'
        var_exps.ensemble = 'Kerry Emmanuel'
        var_exps[:] = np.arange(start=1, stop=n_experiment+1, step=1)

        var_maxelev = nc.createVariable(varname='maxelev', datatype=float, dimensions=('n_experiment', 'nSCHISM_hgrid_node'))
        var_maxelev.long_name = 'Maximum elevation'
        var_maxelev.standard_name = 'Maximum elevation'
        var_maxelev.units = 'm'
        var_maxelev.datum = 'MSL'

        nc.experiment_setup = 'Hydro + 35cm Sa'
        nc.contact = 'Jamal Khan'
        nc.email = 'jamal.khan@legos.obs-mip.fr'
    except:
        print('netCDF Error!')
    else:
        track_data = Gr3()
        for i, in_file in enumerate(in_files):
            print(in_file)
            track_no = int(os.path.basename(in_file).split('.')[0].split('_')[1])
            track_data.read(in_file, readnodes=True, readelements=False)
            var_maxelev[track_no-1, :] = track_data.nodes[:, 3]
            
            if i%100 == 0:
                nc.sync()
    finally:
        nc.close()
