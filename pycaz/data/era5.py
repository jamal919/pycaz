# -*- coding: utf-8 -*-

import numpy as np
from typing import Tuple
from pycaz.typing import ArrayLike

ERA5_WIND_CORRECTIONS = {
    'io': np.array([
        [0, 1],
        [17, 1.03],
        [20, 1.08],
        [23, 1.14],
        [26, 1.20],
        [30, 1.27]
    ])
}


def apply_wind_correction(
        u: ArrayLike,
        v: ArrayLike,
        qfactors: ArrayLike = None,
        basin: str = None) -> Tuple[ArrayLike, ArrayLike]:
    """
    Apply wind correction to ERA5 data. A set of pre-defined correction is given in ERA5_WIND_CORRECTIONS.

    :param u: Input wind u10
    :param v: Input wind v10
    :param qfactors: A 2d array of size(nwindspeed, nfactor). Check ERA5_WIND_CORRECTIONS for sample.
    :param basin: String to indicate basin, as listed in ERA5_WIND_CORRECTIONS.
    :return: Corrected u and v
    """
    _qfactors = None
    if basin is not None:
        if basin.lower() in ERA5_WIND_CORRECTIONS:
            _qfactors = ERA5_WIND_CORRECTIONS[basin]

    if qfactors is not None:
        _qfactors = qfactors

    if _qfactors is None:
        raise ValueError(f'Either qfactors or basin name must be provided')

    u = np.asarray(u)
    v = np.asarray(v)

    uv_lows = _qfactors[:, 0]
    uv_highs = np.append(_qfactors[1:, 0], np.inf)
    factors = _qfactors[:, 1]

    norm = np.sqrt(u ** 2 + v ** 2)
    theta = np.arctan2(v, u)
    norm_corr = np.full(norm.shape, fill_value=np.nan)

    for uv_low, uv_high, factor in zip(uv_lows, uv_highs, factors):
        idx = np.logical_and(norm >= uv_low, norm < uv_high)
        norm_corr[idx] = norm[idx] * factor

    u_corr = np.cos(theta) * norm_corr
    v_corr = np.sin(theta) * norm_corr

    return u_corr, v_corr
