#!/usr/bin/env python

import numpy as np
from numpy.typing import ArrayLike


def interp_complex_1D(x: float, bnds: ArrayLike, amp_pha: ArrayLike, pha_unit='degrees') -> np.ndarray:
    '''
    Interpolate amplitude and phase in 1D between two sorrounding point.
    :x: x
    :bnds: [x1, x2]
    :amp_pha: [[amp_x1, pha_x1], [amp_x2, pha_x2]]
    :pha_unit: 'degrees', 'radians'

    :output: [amp_x, pha_x]
    '''
    # Variables
    bnds = np.atleast_1d(bnds).astype(float)
    amp_pha = np.atleast_2d(amp_pha).astype(float)

    if pha_unit is 'degrees':
        amp_pha[:, 1] = np.deg2rad(amp_pha[:, 1])

    # Test for nan values
    nantest = np.isnan(amp_pha[:, 0])

    if np.all(nantest):
        alpha = np.nan
        beta = np.nan
        # rest of the computation is not necessary
        # retuen a tuple of 
        return (np.array([np.nan, np.nan]))

    # Compute proper dx
    if bnds[0] < bnds[1]:
        # does not go around the world
        dx = bnds[1] - bnds[0]
    else:
        # goes around the world
        dx = (bnds[1] % 90) - (bnds[0] % 90)

    # Find the weight for x1 and x2 point based on dx   
    alpha = (bnds[1] - x) / float(dx)
    beta = 1 - alpha

    # Consider for cases of missing values
    if np.isnan(amp_pha[0, 0]):
        alpha = 1
        beta = 0

    if np.isnan(amp_pha[1, 0]):
        alpha = 0
        beta = 1

    # Interpolate
    sinval = alpha * amp_pha[0, 0] * np.sin(amp_pha[0, 1]) + beta * amp_pha[1, 0] * np.sin(amp_pha[1, 1])
    cosval = alpha * amp_pha[0, 0] * np.cos(amp_pha[0, 1]) + beta * amp_pha[1, 0] * np.cos(amp_pha[1, 1])

    pha = np.arctan2(sinval, cosval)
    if np.sin(pha) == 0:
        amp = cosval / np.cos(pha)
    else:
        amp = sinval / np.sin(pha)

    if pha_unit == 'degrees':
        pha = np.rad2deg(pha)

    return (np.array([amp, pha]))


def interp_complex_2D(xy: ArrayLike, bnds: ArrayLike, amp_pha: ArrayLike, pha_unit='degrees') -> np.ndarray:
    '''
    :xy: [x, y]
    :bnds: [[w, s], [e, s], [e, n], [w, n]] -> counter-clockwise
    :amp_pha: [[amp_ws, pha_ws], [amp_es, pha_es], [amp_en, pha_en], [amp_wn, pha_wn]]
    '''
    xy = np.atleast_1d(xy)
    bnds = np.atleast_2d(bnds)
    amp_pha = np.atleast_2d(amp_pha)

    # in x direction
    amp_pha_s = interp_complex_1D(
        x=xy[0],
        bnds=bnds[0:2, 0],
        amp_pha=amp_pha[0:2, :],
        pha_unit=pha_unit)
    amp_pha_n = interp_complex_1D(
        x=xy[0],
        bnds=bnds[0:2, 0],
        amp_pha=amp_pha[2:4, :],
        pha_unit=pha_unit)

    # in y direction
    amp, pha = interp_complex_1D(
        x=xy[1],
        bnds=bnds[1:3, 1],
        amp_pha=np.array([amp_pha_s, amp_pha_n]),
        pha_unit=pha_unit)

    return (np.array([amp, pha]))


def interp_complex_tri(xy: ArrayLike, bnds: ArrayLike, amp_pha: ArrayLike, pha_unit='degrees'):
    '''
    TODO: Implement the triangular interpolation schemes from pycaz/pyschism (old) scripts
    '''
    raise NotImplementedError


if __name__ == '__main__':
    pass
