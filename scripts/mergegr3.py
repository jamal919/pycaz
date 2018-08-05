# -*- coding: utf-8 -*-
"""
Merge the gr3 formatted maximum elev and maximum velocity outputs. 

This script is a standalone script to merge the output variable containing maximum 
values. It is developed as a part of the schism toolbox module and intended to
be included in the toolbox itself. 

The program is compatible with both python 2 and 3.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
import os
import glob
import numpy as np

# from schism.core import Local2Global
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
                self.globalelemtable[f.elems[:, 1]-1, i+2] = f.elems[f.elemtable[:, 1]-1, 1]

class MaxVariable(object):
    def __init__(self, path, varname='maxelev', varnum=1):
        self.path = path
        self.varname = varname
        self.varnum = varnum

        self.l2g = Local2Globals(self.path)
        self.l2g.load_files()
        self.l2g.merge_nodes()
        self.l2g.merge_elements()
        self.nodes = self.l2g.globalnodetable
        self.elems = self.l2g.globalelemtable
        

    def list_files(self, prefix):
        self.filelist = glob.glob(os.path.join(self.path, prefix))
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
            for i in np.arange(0, self.varnum):
                self.nodes[__ds[:, 0] - 1, i+3] = __ds[i+3]



if __name__=='__main__':
    path = '/run/media/khan/Storehouse/Projects/201803_Surge Level Under Current Climate/Experiments/Sensitivity/Test01/outputs'
    maxelev = MaxVariable(path=path, varname='maxelev', varnum=1)
    maxelev.list_files(prefix='maxelev')
    maxelev.merge_files()
    print(maxelev.nodes)
    print(maxelev.elems)
