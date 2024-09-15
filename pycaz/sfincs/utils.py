#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from hydromt_sfincs import SfincsModel, utils
from pyproj import Transformer
from typing import Literal


def get_bnd_points(model: SfincsModel, bnd_type: Literal['wl'] = 'wl', to_crs=None) -> np.ndarray:
    """
    Get the boundary points from the SfincsModel

    :param model: A SfincsModel instance with msk
    :param bnd_type: Type of boundary, currently only 'wl' is supported, default
    :param to_crs: A target epsg code, typically 4326 (lon-lat)
    :return:
    """
    if bnd_type == 'wl':
        msk_value = 2
    else:
        raise Exception(f'Not defined!')

    try:
        utm_code = model.config['epsg']
    except KeyError:
        raise KeyError(f'Missing epsg key in the model.config!')

    bnds = utils.get_bounds_vector(model.grid.msk)

    xy = np.empty((0, 2))

    for geometry in bnds[bnds.value == 2].geometry:
        x, y = geometry.xy
        xy = np.append(xy, np.vstack([x, y]).T, axis=0)

    if to_crs is not None:
        transformer = Transformer.from_crs(crs_from=utm_code, crs_to=to_crs, always_xy=True)
        xy = np.vstack(transformer.transform(xy[:, 0], xy[:, 1])).T

    return xy

