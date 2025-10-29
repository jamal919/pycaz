# -*- coding: utf-8 -*-

"""
Implements Hgrid and Gr3 related functionalities.
"""

import os
import warnings

import numpy as np

from pycaz.typing import PathLike
from .openbnd import OpenBoundary
from .landbnd import LandBoundary
from .gr3 import Gr3


class Hgrid(Gr3):
    def __init__(self, **kwargs):
        """
        A Hgrid object extended from Gr3 object, with an overrided write functionality.
        
        **kwargs: {
            'header' -> str,
            'nodes' -> np.ndarray,
            'elems' -> np.ndarray,
            'open_bnds' -> OpenBoundary,
            'land_bnds' -> LandBoundary
        }
        """
        super().__init__(**kwargs)

    def write(self, fname: str, overwrite=False):
        nodefmt = ['%10i', '%16.10f', '%16.10f', '%16.10f']

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
                    f.write(f'{i + 1:10d}\t')
                    i34 = self['elemtype'][i]  # 3 or 4 node elements
                    f.write(f'{i34:2d}\t')
                    for n in np.arange(i34):
                        f.write(f'{inodes[n]:9d}\t')
                    f.write('\n')

        # open boundaries
        nbnds = len(self['open_bnds'])
        nbndnodes = np.sum([len(self['open_bnds'][bnd]['nodes']) for bnd in self['open_bnds']]).astype(int)
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
        nbndnodes = np.sum([len(self['land_bnds'][bnd]['nodes']) for bnd in self['land_bnds']]).astype(int)

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

    @property
    def gr3(self):
        """
        Only the Gr3 object form the Hgrid object. In this case, element table is always 
        available.
        """
        return (
            Gr3(
                header=self.header,
                nnode=self.nnode,
                nelem=self.nelem,
                elemtype=self.elemtype,
                nodes=self.nodes,
                elems=self.elems
            )
        )

    @property
    def elems(self):
        """
        The nodal connectivity table of the Hgrid object.
        """
        return self['elems']

    @property
    def open_bnds(self):
        """
        A dictionary containing the open boundaries.
        """
        return self['open_bnds']

    @property
    def land_bnds(self):
        """
        A dictionary containing the land boundaries.
        """
        return self['land_bnds']

    @property
    def hgrid_info(self) -> dict:
        """
        A dictionary containing information related to hgrid to be passed to other classes, particularly Bctides

        :return: Dict
        """
        info = {
            "name": self.header,
            "center": self.center
        }
        return info

    def get_bctides(self):
        """Return an empty Bctides object with the definitions of the Boundaries."""
        # First we load the Bctides class
        # It should not be loaded at the top to avoid cyclic import
        from pycaz.schism.bctides import Bctides
        return Bctides(open_bnds=self.open_bnds, hgrid_info=self.hgrid_info)

    def describe(self):
        """A human-redable description of the hgrid."""
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
            assert (_length > 2)
        except AssertionError:
            raise Exception('Invalid length of file!')

        # Get the grid name at the first line, could be empty
        _name = ' '.join(_ds[0].split())

        # Get the element and node counts from the second line
        _nelem, _nnode = _ds[1].split('!')[0].strip().split()
        _nelem = int(_nelem)
        _nnode = int(_nnode)

        _return_chunks = {
            'header': _name,
            'nelem': _nelem,
            'nnode': _nnode
        }

        # Try reading the nodes sagment
        try:
            _nodes = _ds[2:_nnode + 2]
        except IndexError:
            raise IndexError(f'Could not read {_nnode} nodes.')
        else:
            _return_chunks['nodes'] = _nodes

        # If we do not hit empty line just after nodes, then try reading the element sagment
        if _ds[_nnode + 2].strip():
            try:
                _elems = _ds[_nnode + 2:_nnode + _nelem + 2]
            except IndexError:
                raise IndexError(f'Could not read {_nelem} elements.')
            else:
                _return_chunks['elems'] = _elems

        # If we do not hit empty line just after elements, then try reading the boundary sagment
        try:
            if _ds[_nnode + _nelem + 2].strip():
                _boundaries = _ds[_nnode + _nelem + 2:_length]
                _readflag = (True, True, True)
                _return_chunks['boundaries'] = _boundaries
        except IndexError:
            pass

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
            'nodes': _nodes
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
    _elems = np.ones((_nelem, 4), dtype=int) * _elem_FillValue

    for i, _line in enumerate(chunk):
        _columns = [int(i) for i in _line.split()]
        _elemtype[i] = _columns[1]  # i3 or 4 indicator
        _elems[i, :_columns[1]] = _columns[2:]

    return {
        'elemtype': _elemtype,
        'elems': _elems,
        'elem_FillValue': _elem_FillValue
    }


def _hgrid_parse_boundaries(chunk: list) -> dict:
    """
    Parse boundaries from the chunk found from _find_hgrid_chunks()

    return: dict(open_bnds, land_bnds)
    """
    # First the open boundaries
    _nopen = int(chunk[0].split('=')[0])
    _nopennodes = int(chunk[1].split('=')[0])
    _lnum = 2  # move cursor to first line in the loop
    _open = {}
    for _n in np.arange(_nopen):
        _nnodes = int(chunk[_lnum].split('=')[0])
        _lnum = _lnum + 1
        _nodes = np.genfromtxt(chunk[_lnum:_lnum + _nnodes], dtype=int)
        _open.update({_n + 1: OpenBoundary(name=_n + 1, nodes=_nodes, neta=len(_nodes))})
        _lnum = _lnum + _nnodes  # move cursor to next sagment

    # Then the land boundaries
    _nland = int(chunk[_lnum].split('=')[0])
    _lnum = _lnum + 1
    _nlandnodes = int(chunk[_lnum].split('=')[0])
    _lnum = _lnum + 1
    _land = {}
    for _n in np.arange(_nland):
        _nnodes, _bndtype = [int(i) for i in chunk[_lnum].split('=')[0].split()]
        _lnum = _lnum + 1
        _nodes = np.genfromtxt(chunk[_lnum:_lnum + _nnodes], dtype=int)
        _land.update({_n + 1: LandBoundary(name=_n + 1, bndtype=_bndtype, nodes=_nodes)})
        _lnum = _lnum + _nnodes  # move cursor to next sagment

    _land.keys()

    return {
        'open_bnds': _open,
        'land_bnds': _land
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


def read_hgrid(fname: PathLike) -> Hgrid:
    """
    Reads a hgrid file
    """
    chunks = _hgrid_find_chunks(fname)

    try:
        _necessary_fields = ['header', 'nelem', 'nnode', 'nodes', 'elems']
        _necessary_fields_found = [i in chunks for i in _necessary_fields]
        assert (np.all(_necessary_fields_found))
    except AssertionError:
        Exception('hgrid file does not have all the  necessary fields. Currently {_necessary_fields} are mandatory.')
    else:
        hgrid = Hgrid(header=chunks['header'], nelem=chunks['nelem'], nnode=chunks['nnode'])
        hgrid.update(_hgrid_parse_nodes(chunks['nodes']))
        hgrid.update(_hgrid_parse_elements(chunks['elems']))

    try:
        assert ('boundaries' in chunks)
    except AssertionError:
        warnings.warn('Boundaries are not found in the file! Will set to empty values')
        boundaries = {
            'open_bnds': {},
            'land_bnds': {}
        }
    else:
        boundaries = _hgrid_parse_boundaries(chunks['boundaries'])

    # Add nodes x,y location in the boundaries
    for open_bnd in boundaries['open_bnds']:
        boundaries['open_bnds'][open_bnd]['xy'] = hgrid.subset_nodes(boundaries['open_bnds'][open_bnd]['nodes'])

    for land_bnd in boundaries['land_bnds']:
        boundaries['land_bnds'][land_bnd]['xy'] = hgrid.subset_nodes(boundaries['land_bnds'][land_bnd]['nodes'])

    # Update the hgrid structure
    hgrid.update(boundaries)

    return hgrid
