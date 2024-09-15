#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from typing import Dict, Tuple

import numpy as np
from pyproj import Transformer

from pycaz.utils.geometry import get_utm_code

logger = logging.getLogger(__name__)


def get_basic_config(extent, dx, dy=None, rotation=0, wgs_code=4326):
    """
    Create a basic SFINCS grid configuration.

    This function is not recommended to use, but was written for learning purpose. Use the dedicated
    hydromt_sfincs toolbox for actual working configuration.

    :param extent: Extent of the model in lon,lat as array of [west, east, south, north]
    :param dx: dx size, meters
    :param dy: dy size, meters
    :param rotation: rotation, degrees, default 0
    :param wgs_code: default 4326, do not change
    :return: (basic_config, wgs2utm, utm2wgs)
    """
    # Find the corresponding utm_code
    utm_code = get_utm_code(extent)
    logger.info(f'UTM code for the extent: {utm_code}')

    # Create a transfer to compute the corresponding transformation
    wgs2utm = Transformer.from_crs(crs_from=wgs_code, crs_to=utm_code, always_xy=True)
    utm2wgs = Transformer.from_crs(crs_from=utm_code, crs_to=wgs_code, always_xy=True)

    # Lower left corner
    utm_x1, utm_y1 = wgs2utm.transform(extent[0], extent[2])
    # Upper right corner
    utm_x2, utm_y2 = wgs2utm.transform(extent[1], extent[3])

    # Now generating model variable
    x0 = utm_x1
    y0 = utm_y1

    # Model resolution
    dx = dx
    if dy is None:
        dy = dx
    else:
        dy = dy

    # mmax, number of "grid cell" in x
    mmax = np.ceil((utm_x2 - utm_x1) / dx).astype(int)

    # nmax, number of "grid cell" in y
    nmax = np.ceil((utm_y2 - utm_y1) / dy).astype(int)

    # rotation is set to 0 for now
    if rotation != 0:
        rotation = 0
        logger.info(f'rotation feature is implemented, rotation value is set to 0')

    basic_config = {
        'x0': x0,
        'y0': y0,
        'mmax': mmax,
        'nmax': nmax,
        'dx': dx,
        'dy': dy,
        'rotation': rotation,
        'epsg': utm_code
    }

    return basic_config, wgs2utm, utm2wgs


def get_ordinate_grid(config: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns the ordinate grid of the model.

    :param config: Model configuration as a dictionary
    :return: (X_ord, Y_ord)
    """
    # Get values from config
    x0 = config['x0']
    y0 = config['y0']
    dx = config['dx']
    dy = config['dy']
    mmax = config['mmax']
    nmax = config['nmax']

    # model grid coordinates
    x_ord = x0 + np.arange(mmax + 1) * dx
    y_ord = y0 + np.arange(nmax + 1) * dy

    X_ord, Y_ord = np.meshgrid(x_ord, y_ord)

    return X_ord, Y_ord


def get_cell_grid(config: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns the cell grid of the model. Most of the variable in the model is based on the cell grid.

    :param config: Model configuration as a dictionary
    :return: (X_cell, Y_cell)
    """
    # Get values from config
    x0 = config['x0']
    y0 = config['y0']
    dx = config['dx']
    dy = config['dy']
    mmax = config['mmax']
    nmax = config['nmax']

    # model grid coordinates
    x_ord = x0 + np.arange(mmax + 1) * dx
    y_ord = y0 + np.arange(nmax + 1) * dy

    # model cell coordinates, 1 less than grid size
    x_cell = (x_ord[0:-1] + x_ord[1:]) / 2
    y_cell = (y_ord[0:-1] + y_ord[1:]) / 2

    X_cell, Y_cell = np.meshgrid(x_cell, y_cell)

    return X_cell, Y_cell
