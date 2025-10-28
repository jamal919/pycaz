# -*- encoding: utf-8 -*-
import os
from copy import deepcopy

import numpy as np
import xarray as xr
from typing_extensions import Self

from pycaz.typing import PathLike


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
    def header(self, header: str):
        _header = str(header)
        self['header'] = _header

    # nodes related functionalities
    @property
    def nnode(self) -> int:
        return self['nnode']

    @property
    def nodes(self) -> np.ndarray:
        return self['nodes']

    @property
    def nodeid(self):
        return self.nodes[:, 0].astype(int)

    @property
    def x(self):
        return self.nodes[:, 1].astype(float)

    @property
    def y(self):
        return self.nodes[:, 2].astype(float)

    @property
    def center(self):
        _x = (np.min(self.x) + np.max(self.x)) / 2
        _y = (np.min(self.y) + np.max(self.y)) / 2
        return _x, _y

    @property
    def xy(self):
        return self.nodes[:, 1:3].astype(float)

    @property
    def data(self):
        return self.nodes[:, 3:].astype(float)

    @property
    def ndata(self):
        return self.nodes.shape[1] - 3

    @data.setter
    def data(self, data):
        try:
            _data = np.atleast_1d(data)
            assert (np.shape(_data)[0] == self.data.shape[0])
        except AssertionError:
            raise Exception(f'Size mismatch! Both length must be {self.data.shape[0]}')
        else:
            self['nodes'] = np.hstack([self.nodes[:, :3], _data])

    # Element related functionalities
    @property
    def elemtype(self) -> np.ndarray:
        return self['elemtype']

    @property
    def nelem(self) -> int:
        return self['nelem']

    @property
    def meshtype(self) -> str:
        if np.all([i in [3, 4] for i in np.unique(self.elemtype)]):
            return 'i34'
        else:
            return 'i3'

    def subset_nodes(self, nodeid: np.ndarray):
        # check the nodeid is 1-based index
        try:
            assert np.all(nodeid >= 1)
        except:
            raise AssertionError('Node ids must start from 1')

        return self.xy[nodeid - 1, :]

    def extent(self, buffer: float = 0):
        extent = np.array([np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)])

        if buffer != 0:
            extent = extent + np.array([buffer * -1, buffer, buffer * -1, buffer])

        return extent

    def copy(self) -> Self:
        return deepcopy(self)

    def split_quads(self) -> None:
        """Convert to a fully triangular grid from hybrid quad-tri grid

        TODO: Implement based on split_quads_wwm.f90
        """
        raise NotImplementedError

    # I/O functionalities
    def write(self, fname: PathLike, fmt='%16.10f', overwrite=False) -> None:
        """Write the grid to file

        This methods writes the grid to a file specified by path. The gr3
        format is implemented as appear in SCHISM manual.
        """
        nodefmt = ['%10i', '%16.10f', '%16.10f', fmt]

        if os.path.exists(fname) & ~overwrite:
            raise Exception('File exists! Set overwrite=True to overwrite.')
        else:
            with open(fname, 'w') as f:
                f.write(f"{self['header']}\n")
                f.write(f"{self['nelem']}\t{self['nnode']} ")
                f.write("! # of elements and nodes in the horizontal grid\n")
                np.savetxt(fname=f, X=self.nodes, fmt=nodefmt, delimiter='\t')

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

    def to_xarray(self, varname: str = 'depth') -> xr.Dataset:
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
            assert (np.all(_necessary_fields_found))
        except AssertionError:
            Exception(f'Necessary fields {_necessary_fields} not found!')

        # All element is in present
        nSCHISM_hgrid_node = np.arange(self['nnode']) + 1
        nSCHISM_hgrid_face = np.arange(self['nelem']) + 1
        nMaxSCHISM_hgrid_face_nodes = np.arange(4)
        nData = self.ndata
        ELEM_FILL_VALUE = self['elem_FillValue']

        SCHISM_hgrid_node_x = xr.DataArray(
            data=self.x,  # first column
            dims=['nSCHISM_hgrid_node'],
            coords={
                'nSCHISM_hgrid_node': nSCHISM_hgrid_node
            }
        )

        SCHISM_hgrid_node_y = xr.DataArray(
            data=self.y,  # second column
            dims=['nSCHISM_hgrid_node'],
            coords={
                'nSCHISM_hgrid_node': nSCHISM_hgrid_node
            }
        )

        SCHISM_hgrid_face_nodes = xr.DataArray(
            data=self['elems'],
            dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'],
            coords={
                'nSCHISM_hgrid_face': nSCHISM_hgrid_face,
                'nMaxSCHISM_hgrid_face_nodes': nMaxSCHISM_hgrid_face_nodes
            }
        )

        if nData == 1:
            vardata = xr.DataArray(
                data=self.data.flatten(),
                dims=['nSCHISM_hgrid_node'],
                coords={
                    'nSCHISM_hgrid_node': nSCHISM_hgrid_node
                },
                attrs={
                    'header': self['header']
                }
            )
        elif nData > 1:
            varcoord = np.arange(nData)
            vardata = xr.DataArray(
                data=self.data,
                dims=['nSCHISM_hgrid_node', varname],
                coords={
                    'nSCHISM_hgrid_node': nSCHISM_hgrid_node,
                    varname: varcoord
                },
                attrs={
                    'header': self['header']
                }
            )
        else:
            Exception('Data in Gr3 is not correct! Data size mismatch.')

        ds = xr.Dataset(
            data_vars={
                'SCHISM_hgrid_node_x': SCHISM_hgrid_node_x,
                'SCHISM_hgrid_node_y': SCHISM_hgrid_node_y,
                'SCHISM_hgrid_face_nodes': SCHISM_hgrid_face_nodes,
                varname: vardata
            },
            attrs={
                'header': self['header']
            }
        )

        ds.SCHISM_hgrid_face_nodes.encoding['_FillValue'] = ELEM_FILL_VALUE

        return ds
