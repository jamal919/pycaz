#!/usr/bin/env python

import numpy as np
import matplotlib.tri as mtri

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from numpy.typing import ArrayLike


def grid_around(
        xy: ArrayLike,
        grid_x: ArrayLike,
        grid_y: ArrayLike,
        extrapolate: Literal['spherical', 'nearest'] = 'spherical') -> np.ndarray:
    """Returns the sorrounding grid nodes.

    Finds the sorrounding nodes from the structured grid defined by `grid_lon` and `grid_lat` at the location defined 
    by `point_lon`, `point_lat`. The returning results are a tuple [[lon_west, lon_east], [lat_south, lat_north]]

    If `extrapolate` is `spherical`, then the west and south rounds the world. For example, if a point (lon=0.01) is 
    chosen in a lon grid [0.025, 360], then the right side of the grid pixel will be 0.025, and the left pixel will be 360.

    If `extrapolate` is `nearest`, it will just select the nearest neighbour.

    Args:
        xy (ArrayLike): Point location [x, y]
        grid_x (ArrayLike): Grid x points list
        grid_y (ArrayLike): Grid y points list
        extrapolate (Literal['spherical', 'nearest'], optional): Extrapolate method. Defaults to 'spherical'.

    Returns:
        np.ndarray: [[lon_west, lon_east], [lat_south, lat_north]]
    """
    grid_x = np.atleast_1d(grid_x)
    grid_y = np.atleast_1d(grid_y)
    xy = np.atleast_1d(xy)

    ieast = np.argmax(xy[0] <= grid_x)
    if ieast == 0:
        iwest = (ieast - 1) * (extrapolate == 'spherical') + (ieast) * (extrapolate == 'nearest')
    else:
        iwest = ieast - 1
    sorrounding_lon = grid_x[[iwest, ieast]]

    inorth = np.argmax(xy[1] <= grid_y)
    if inorth == 0:
        isouth = (inorth - 1) * (extrapolate == 'spherical') + (inorth) * (extrapolate == 'nearest')
    else:
        isouth = inorth - 1

    sorrounding_lat = grid_y[[isouth, inorth]]

    # returning sorrounding nodes in counter-clockwise order
    sorrounding_nodes = np.array([
        [sorrounding_lon[0], sorrounding_lat[0]],  # [west, south]
        [sorrounding_lon[1], sorrounding_lat[0]],  # [east, south]
        [sorrounding_lon[1], sorrounding_lat[1]],  # [east, north]
        [sorrounding_lon[0], sorrounding_lat[1]]  # [west, north]
    ])

    return (sorrounding_nodes)


def tri_around(
        xy: np.ndarray,
        triang: mtri.Triangulation) -> np.ndarray:
    '''
    Finds the sorroundig nodes in a triangular mesh.
    '''
    raise NotImplementedError
