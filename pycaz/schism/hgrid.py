#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implements Hgrid and Gr3 related functionalities.
"""

from copy import deepcopy
from typing_extensions import Self
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr

# Data classes
class Gr3(dict):
    def __init__(self, **kwargs):
        """
        A gr3 object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.update(kwargs)

    # headers
    @property
    def header(self) -> str:
        return self['header']

    @header.setter
    def header(self, header:str):
        _header = str(header)
        self['header'] = _header

    # nodes related functionalities
    @property
    def nodes(self) -> np.ndarray:
        return self['nodes']

    @property
    def nodeid(self):
        return self.nodes[:, 0]

    @property
    def x(self):
        return self.nodes[:, 1]

    @property
    def y(self):
        return self.nodes[:, 2]

    @property
    def data(self):
        return self.nodes[:, 3:]

    @property
    def ndata(self):
        return self.nodes.shape[1] - 3

    @data.setter
    def data(self, data):
        try:
            _data = np.atleast_1d(data)
            assert(np.shape(_data)[0] == self.data.shape[0])
        except AssertionError:
            raise Exception(f'Size mismatch! Both length must be {self.data.shape[0]}')
        else:
            self['nodes'] = np.hstack([self.nodes[:, :3], _data])

    # Element related functionalities
    @property
    def elemtype(self) -> np.ndarray:
        try:
            return self['elemtype']
        except KeyError:
            raise Exception('Element table is not present!')

    @property
    def meshtype(self) -> str:
        if np.all([i in [3, 4] for i in np.unique(self['elemtype'])]):
            return('i34')
        else:
            return('i3')

    def copy(self) -> Self:
        return deepcopy(self)

    def split_quads(self) -> None:
        """Convert to a fully triangular grid from hybrid quad-tri grid

        TODO: Implemented based on split_quads_wwm.f90
        """
        pass
    
    # I/O functionalities
    def write(self, fname, fmt='%16.10f', overwrite=False) -> None:
        """Write the grid to file
        
        This methods writes the grid to a file specified by path. The gr3
        format is implemented as appear in SCHISM manual.
        """
        nodefmt=['%10i', '%16.10f', '%16.10f', fmt]

        if os.path.exists(fname) & ~overwrite:
            raise Exception('File exists! Set overwrite=True to overwrite.')
        else:
            with open(fname, 'w') as f:
                f.write(f"{self['header']}\n")
                f.write(f"{self['nelem']}\t{self['nnode']} ")
                f.write("! # of elements and nodes in the horizontal grid\n")
                np.savetxt(fname=f, X=self.nodes, fmt=nodefmt)

        # write element table if exist
        if 'elems' in self:
            with open(fname, 'a') as f:
                for i in np.arange(self['nelem']):
                    inodes = self['elems'][i]
                    f.write(f'{i+1:10d}\t')
                    i34 = self['elemtype'][i] # 3 or 4 node elements
                    f.write(f'{i34:2d}\t')
                    for n in np.arange(i34):
                        f.write(f'{inodes[n]:9d}\t')
                    f.write('\n')

    def to_xarray(self, varname:str='depth') -> xr.Dataset:
        """Returns a xarray dataset of gr3/hgrid
        
            Dimensions:
            nSCHISM_hgrid_node : number of nodes;
            nSCHISM_hgrid_face : number of faces;
            nMaxSCHISM_hgrid_face_nodes: 4; # hardcoded
            nData: number of data columns; only appears for velocity
            
            Variables:
            x column: SCHISM_hgrid_node_x
            y column: SCHISM_hgrid_node_y
            data columns: the name is to be taken as input
            
            Global attributes:
            Convention = "CF-1.0"
            Information = From mesh header
        """
        try:
            _necessary_fields = ['nodes', 'elems', 'nodes']
            _necessary_fields_found = [i in self for i in _necessary_fields]
            assert(np.all(_necessary_fields_found))
        except AssertionError:
            Exception(f'Necessary fields {_necessary_fields} not found!')
        
        # All element is in present
        nSCHISM_hgrid_node = np.arange(self['nnode']) + 1
        nSCHISM_hgrid_face = np.arange(self['nelem']) + 1
        nMaxSCHISM_hgrid_face_nodes = np.arange(4)
        nData = self.ndata
        ELEM_FILL_VALUE = self['elem_FillValue']

        SCHISM_hgrid_node_x = xr.DataArray(
            data=self.x, # first column
            dims=['nSCHISM_hgrid_node'],
            coords={
                'nSCHISM_hgrid_node':nSCHISM_hgrid_node
            }
        )

        SCHISM_hgrid_node_y = xr.DataArray(
            data=self.y, # second column
            dims=['nSCHISM_hgrid_node'],
            coords={
                'nSCHISM_hgrid_node':nSCHISM_hgrid_node
            }
        )

        SCHISM_hgrid_face_nodes = xr.DataArray(
            data=self['elems'],
            dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'],
            coords={
                'nSCHISM_hgrid_face':nSCHISM_hgrid_face,
                'nMaxSCHISM_hgrid_face_nodes':nMaxSCHISM_hgrid_face_nodes
            }
        )

        if nData == 1:
            vardata = xr.DataArray(
                data=self.data.flatten(),
                dims=['nSCHISM_hgrid_node'],
                coords={
                    'nSCHISM_hgrid_node':nSCHISM_hgrid_node
                },
                attrs={
                    'header':self['header']
                }
            )
        elif nData > 1:
            varcoord = np.arange(nData)
            vardata = xr.DataArray(
                data=self.data,
                dims=['nSCHISM_hgrid_node', varname],
                coords={
                    'nSCHISM_hgrid_node':nSCHISM_hgrid_node,
                    varname:varcoord
                },
                attrs={
                    'header':self['header']
                }
            )
        else:
            Exception('Data in Gr3 is not correct! Data size mismatch.')


        ds = xr.Dataset(
            data_vars={
                'SCHISM_hgrid_node_x':SCHISM_hgrid_node_x,
                'SCHISM_hgrid_node_y':SCHISM_hgrid_node_y,
                'SCHISM_hgrid_face_nodes':SCHISM_hgrid_face_nodes,
                varname:vardata
            },
            attrs={
                'header':self['header']
            }
        )

        ds.SCHISM_hgrid_face_nodes.encoding['_FillValue'] = ELEM_FILL_VALUE

        return ds

class Boundary(dict):
    def __init__(self, **kwargs):
        """
        A gr3 object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.update(kwargs)

    @property
    def name(self):
        return(self['name'])

    @name.setter
    def name(self, new_name: str):
        self['name'] = new_name
    
    @property
    def nodes(self):
        return self['nodes']

class Hgrid(Gr3):
    def __init__(self, **kwargs):
        """
        A Hgrid object extended from Gr3 object, with an overrided write functionality.
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(**kwargs)

    def write(self, fname:str, overwrite=False):
        nodefmt=['%10i', '%16.10f', '%16.10f', '%16.10f']

        if os.path.exists(fname) & ~overwrite:
            raise Exception('File exists! Set overwrite=True to overwrite.')
        else:
            with open(fname, 'w') as f:
                f.write(f"{self['header']}\n")
                f.write(f"{self['nelem']}\t{self['nnode']} ")
                f.write("! # of elements and nodes in the horizontal grid\n")
                np.savetxt(fname=f, X=self.nodes, fmt=nodefmt)

        # write element table if exist
        if 'elems' in self:
            with open(fname, 'a') as f:
                for i in np.arange(self['nelem']):
                    inodes = self['elems'][i]
                    f.write(f'{i+1:10d}\t')
                    i34 = self['elemtype'][i] # 3 or 4 node elements
                    f.write(f'{i34:2d}\t')
                    for n in np.arange(i34):
                        f.write(f'{inodes[n]:9d}\t')
                    f.write('\n')

        # open boundaries
        nbnds = len(self['open_bnds'])
        nbndnodes = np.sum([len(self['open_bnds'][bnd]['nodes']) for bnd in self['open_bnds']])
        with open(fname, 'a') as f:
            f.write(f'{nbnds:d} = Number of open boundaries\n')
            f.write(f'{nbndnodes:d} = Total number of open boundary nodes\n')

        for bnd in self['open_bnds']:
            bndname = self['open_bnds'][bnd]['name']
            bndnodes = self['open_bnds'][bnd]['nodes']

            with open(fname, 'a') as f:
                f.write(f'{len(bndnodes):d} = Number of nodes for open boundary {bnd:d} - {bndname}\n')
                np.savetxt(fname=f, X=bndnodes, fmt='%i')
        
        # land boundaries
        nbnds = len(self['land_bnds'])
        nbndnodes = np.sum([len(self['land_bnds'][bnd]['nodes']) for bnd in self['land_bnds']])

        with open(fname, 'a') as f:
            f.write(f'{nbnds:d} = Number of land boundaries\n')
            f.write(f'{nbndnodes:d} = Total number of land boundary nodes\n')

        for bnd in self['land_bnds']:
            bndname = self['land_bnds'][bnd]['name']
            bndnodes = self['land_bnds'][bnd]['nodes']
            bndtype = self['land_bnds'][bnd]['bndtype']
            with open(fname, 'a') as f:
                f.write(f'{len(bndnodes):d}\t{bndtype:d} = Number of nodes for land boundary {bnd:d} - {bndname}\n')
                np.savetxt(fname=f, X=bndnodes, fmt='%i')
    
    def describe(self):
        header = self['header']
        nelem = self['nelem']
        nnode = self['nnode']
        elemtype = 'hybrid' if np.sum(self['elemtype'] == 4) > 0 else 'tri'
        nopen = len(self['open_bnds'])
        nland = len(self['land_bnds'])

        print(f'{header}\n{nnode} nodes\n{nelem} elements of {elemtype} type\n{nopen} open, {nland} land boundaries')



# Parsing related functions
def _hgrid_find_chunks(fname: str):
    """
    Find different chunk of the Gr3/Hgrid file
    
    returns: dict(header, elem, nodes, [elems, [boundaries]])

    """
    with open(fname) as f:
        _ds = f.readlines()
        _length = len(_ds)

        # first check if there is minimum 2 lines to give us basic information
        try:
            assert(_length > 2)
        except AssertionError:
            raise Exception('Invalid length of file!')
        
        # Get the grid name at the first line, could be empty
        _name = ' '.join(_ds[0].split())
        
        # Get the element and node counts from the second line
        _nelem, _nnode = _ds[1].split()
        _nelem = int(_nelem)
        _nnode = int(_nnode)

        _return_chunks = {
            'header':_name,
            'nelem':_nelem,
            'nnode':_nnode
        }
        
        # Try reading the nodes sagment
        try:
            _nodes = _ds[2:_nnode+2]
        except IndexError:
            raise IndexError(f'Could not read {_nnode} nodes.')
        else:
            _return_chunks['nodes'] = _nodes
        
        # If we do not hit empty line just after nodes, then try reading the element sagment
        if _ds[_nnode+2].strip():
            try:
                _elems = _ds[_nnode+2:_nnode+_nelem+2]
            except IndexError:
                raise IndexError(f'Could not read {_nelem} elements.')
            else:
                _return_chunks['elems'] = _elems
            
        # If we do not hit empty line just after elements, then try reading the boundary sagment
        if _ds[_nnode+_nelem+2].strip():
            _boundaries = _ds[_nnode+_nelem+2:_length]
            _readflag = (True, True, True)
            _return_chunks['boundaries'] = _boundaries
                
    return _return_chunks

def _hgrid_parse_nodes(chunk: list) -> dict:
    """
    Parse nodes from the chunk found from _find_hgrid_chunks()

    returns: dict(nodes, data)
    """
    try:
        _nodes = np.genfromtxt(chunk)
    except:
        raise Exception('Problem with parsing nodes. Check the output of _find_hgrid_chunks().')
    else:
        return {
            'nodes':_nodes
        }

def _hgrid_parse_elements(chunk: list) -> dict:
    """
    Parse elements from the chunk found from _find_hgrid_chunks()

    return: dict(elemtype, data)
    """
    # Element table can be of 3 or 4 nodes, needs parsing line by line
    _nelem = len(chunk)
    _elemtype = np.zeros(_nelem, dtype=int)
    _elem_FillValue = -99999
    _elems = np.ones((_nelem, 4), dtype=int)*_elem_FillValue

    for i, _line in enumerate(chunk):
        _columns = [int(i) for i in _line.split()]
        _elemtype[i] = _columns[1] # i3 or 4 indicator
        _elems[i, :_columns[1]] = _columns[2:]

    return {
        'elemtype':_elemtype,
        'elems':_elems,
        'elem_FillValue':_elem_FillValue
    }

def _hgrid_parse_boundaries(chunk: list) -> dict:
    """
    Parse boundaries from the chunk found from _find_hgrid_chunks()

    return: dict(open_bnds, land_bnds)
    """
    # First the open boundaries
    _nopen = int(chunk[0].split('=')[0])
    _nopennodes = int(chunk[1].split('=')[0])
    _lnum = 2 # move cursor to first line in the loop
    _open = {}
    for _n in np.arange(_nopen):
        _nnodes =  int(chunk[_lnum].split('=')[0])
        _lnum = _lnum + 1
        _nodes = np.genfromtxt(chunk[_lnum:_lnum+_nnodes], dtype=int)
        _open.update({_n+1:Boundary(name=_n+1, nodes=_nodes)})
        _lnum = _lnum+_nnodes # move cursor to next sagment

    # Then the land boundaries
    _nland = int(chunk[_lnum].split('=')[0])
    _lnum = _lnum + 1
    _nlandnodes = int(chunk[_lnum].split('=')[0])
    _lnum = _lnum + 1
    _land = {}
    for _n in np.arange(_nland):
        _nnodes, _bndtype = [int(i) for i in chunk[_lnum].split('=')[0].split()]
        _lnum = _lnum + 1
        _nodes = np.genfromtxt(chunk[_lnum:_lnum+_nnodes], dtype=int)
        _land.update({_n+1:Boundary(name=_n+1, bndtype=_bndtype, nodes=_nodes)})
        _lnum = _lnum+_nnodes # move cursor to next sagment

    _land.keys()

    return {
        'open_bnds':_open,
        'land_bnds':_land
    }

# Main exposed functions
def read_gr3(fname: str) -> Gr3:
    """
    Reads a gr3 file and return a Gr3 dict object.

    Does not through error if elems are missing.
    """
    chunks = _hgrid_find_chunks(fname)
    
    # Name, elem/node count, and nodes are by design available
    gr3 = Gr3(header=chunks['header'], nelem=chunks['nelem'], nnode=chunks['nnode'])
    gr3.update(_hgrid_parse_nodes(chunks['nodes']))

    # Add the element table only if element is available
    if 'elems' in chunks:
        gr3.update(_hgrid_parse_elements(chunks['elems']))

    return gr3

def read_hgrid(fname: str) -> Hgrid:
    """
    Reads a hgrid file
    """
    chunks = _hgrid_find_chunks(fname)

    try:
        _necessary_fields = ['header', 'nelem', 'nnode', 'nodes', 'elems']
        _necessary_fields_found = [i in chunks for i in _necessary_fields]
        assert(np.all(_necessary_fields_found))
    except AssertionError:
        Exception('hgrid file does not have all the  necessary fields. Currently {_necessary_fields} are mandatory.')
    else:
        hgrid = Hgrid(header=chunks['header'], nelem=chunks['nelem'], nnode=chunks['nnode'])
        hgrid.update(_hgrid_parse_nodes(chunks['nodes']))
        hgrid.update(_hgrid_parse_elements(chunks['elems']))

    try:
        assert('boundaries' in chunks)
    except AssertionError:
        warnings.warn('Boundaries are not found in the file! Will set to empty values')
        boundaries = {
            'open_bnds':{},
            'land_bnds':{}
        }
    else:
        boundaries = _hgrid_parse_boundaries(chunks['boundaries'])

    hgrid.update(boundaries)

    return(hgrid)

