# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
from scipy.spatial import KDTree


def find_nearest_nodes(ds: xr.Dataset, xy: np.ndarray, deeper_than: float = None) -> np.ndarray:
    """
    Find the nearest node numbers for a given set of xy values from a SCHISM out2d output

    Args:
        ds (xr.Dataset):  A out2d SCHISM output dataset read with xr.open_dataset() or xr.open_mfdataset()
        xy (np.ndarray): xy locations where the values are searched with shape (n, 2)
        deeper_than (float): Only search on the nodes that are deeper than deeper_than, downward positive.

    Returns: Index of nSCHISM_hgrid_node to be selected
    """
    # processing variables
    x = ds.SCHISM_hgrid_node_x
    y = ds.SCHISM_hgrid_node_y
    try:
        depth = ds.depth.isel(time=0).drop_vars("time")
    except ValueError:
        depth = ds.depth

    idx = ds.nSCHISM_hgrid_node

    # process deeper_than constraint
    if deeper_than is not None:
        is_deep = (depth >= deeper_than).compute()
        _x = x[is_deep]
        _y = y[is_deep]
        _idx = idx[is_deep]
    else:
        _x = x
        _y = y
        _idx = idx

    kdtree = KDTree(np.vstack((_x, _y)).T)
    _nn_dist, _nn_idx = kdtree.query(xy)
    bnd_idx = _idx[_nn_idx] # returns xr.DataArray
    bnd_idx = np.atleast_1d(bnd_idx.values) # extracting the array

    return bnd_idx


def compute_maxelev(ds: xr.Dataset) -> xr.Dataset:
    """
    Compute maximum elevation from the schism output dataset

    Args:
        ds (xr.Dataset): SCHISM output dataset loaded with xr.open_dataset() or xr.open_mfdataset()

    Returns: xr.Dataset of maximum elevation

    """
    maxelev = ds.elevation.max(dim="time")
    try:
        depth = ds.depth.isel(time=0).drop_vars("time")
    except ValueError:
        depth = ds.depth
    maxelev = (maxelev + depth) * (depth < 0) + (maxelev) * (depth >= 0)

    try:
        elements = ds.SCHISM_hgrid_face_nodes.isel(time=0).drop_vars("time")
    except ValueError:
        elements = ds.SCHISM_hgrid_face_nodes

    ds_out = xr.Dataset({
        "maxelev": maxelev,
        "depth": depth,
        "SCHISM_hgrid_face_nodes": elements,
    }).squeeze()

    return ds_out.compute()


def compute_depth_mask(ds: xr.Dataset, threshold: float = 0) -> np.array:
    """
    Compute mask based on depth for schism dataset, true: masked, false: not masked

    Args:
        ds (xr.Dataset): SCHISM output dataset, or SCHISM grid netcdf dataset
        threshold (float): Value below this threshold will be masked

    Returns:

    """
    try:
        triangles = ds.SCHISM_hgrid_face_nodes.isel(time=0)
    except ValueError:
        triangles = ds.SCHISM_hgrid_face_nodes
    finally:
        triangles = triangles[:, 0:3] - 1
        triangles = triangles.values.astype(int)

    try:
        depth = ds.depth.isel(time=0)
    except ValueError:
        depth = ds.depth
    finally:
        depth = depth.values

    depth_at_triangles = depth[triangles.flatten()].reshape(triangles.shape)
    depth_above_threshold = depth_at_triangles < threshold
    depth_mask = np.all(depth_above_threshold, axis=1)

    return depth_mask
