#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SCHISM Model Builder (SCHISMMB)

SCHISMMB is the container of the the classes and methods to build the the 
inputs for SCHISM Modelling system. 

@author: Md. Jamal Uddin Khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np

class project(object):
    """
    SCHISM model project object. It only contains Path information regarding
    the model project and passed to the model object for further processing.
    """
    def __init__(self, path = None):
        """ Initialize the project folder
        
            path - location of project folder
            
        """
        self.path = path
        if os.access(path, os.F_OK):
            self.direxist = True
        else:
            os.mkdir(path)
            self.direxist = True
            
    def loadhgrid(self, hgridfile):
        """ Load the hgrid file into the project """
        self.hgrid = hgrid(path = hgridfile)
            
    def delete(self):
        """ Delete project """
        os.rmdir(self.path)


class hgrid(object):
    """
    SCHISM hgrid.gr3 object. This object loads the hgrid.gr3 file, breaks it in
    it's subcomponent, and store it for the model object to use. 
    """
    def __init__(self, path = None):
        self.path = path
        if os.access(self.path, os.F_OK):
            self.read()
        else:
            print('Check the hgrid.gr3 file location')
        
    def read(self): 
        """ Read hgrid.gr3 file """
        print("Reading hgrid.gr3 file")
        
    def plot(self):
        """ Plot the hgrid.gr3 mesh """
        print('Plotting hgrid.gr3 file')