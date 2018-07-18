#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SCHISM data abstraction class for core data types

This class contains the core data structure to handle basic elements in 
SCHISM input output file. This library is a requirement for schismio and 
imported there.

Classes:
    * Boundary
    * Boundaries
    
TODO:
    * Check consistancy of the classes

@author: Md. Jamal Uddin Khan
@email: jamal.khan@legos.obs-mip.fr
"""

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