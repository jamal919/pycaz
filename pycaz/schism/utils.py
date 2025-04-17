# -*- coding: utf-8 -*-

from scipy.spatial import KDTree
import xarray as xr
import numpy as np


def find_nearest_nodes(ds: xr.Dataset, xy: np.ndarray, deeper_than: float = None) -> np.ndarray:
    """
    Find the nearest node numbers for a given set of xy values from a SCHISM out2d output

    :param ds: A out2d SCHISM output dataset read with xr.open_dataset() or xr.open_mfdataset()
    :param xy: xy locations where the values are searched with shape (n, 2)
    :param deeper_than: Only search on the nodes that are deeper than deeper_than, downward positive.
    :return: Index of nSCHISM_hgrid_node to be selected

    TODO:
        - ds only takes xr.Dataset, needs to add schism.hgrid.Hgrid as an accepted format
    """
    # processing variables
    x = ds.SCHISM_hgrid_node_x
    y = ds.SCHISM_hgrid_node_y
    depth = ds.depth
    idx = ds.nSCHISM_hgrid_node

    # process deeper_than constraint
    if deeper_than is not None:
        is_deep = (depth >= deeper_than).compute()  # TODO: potential problem with compute()?
        _x = x[is_deep]
        _y = y[is_deep]
        _idx = idx[is_deep]
    else:
        _x = x
        _y = y
        _idx = idx

    kdtree = KDTree(np.vstack((_x, _y)).T)
    _nn_dist, _nn_idx = kdtree.query(xy)
    bnd_idx = _idx[_nn_idx]

    return bnd_idx
