#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from hydromt_sfincs import SfincsModel, utils


def get_bnd_points(model: SfincsModel, bnd_type='wl') -> np.ndarray:
    """
    Get the boundary points from the SfincsModel

    :param model: A SfincsModel instance with msk
    :param bnd_type: Type of boundary, currently only 'wl' is supported, default
    :return:
    """
    if bnd_type == 'wl':
        msk_value = 2
    else:
        raise Exception(f'Not defined!')

    bnds = utils.get_bounds_vector(model.grid.msk)

    xy = np.empty((0, 2))

    for geometry in bnds[bnds.value == 2].geometry:
        x, y = geometry.xy
        xy = np.append(xy, np.vstack([x, y]).T, axis=0)

    return xy
