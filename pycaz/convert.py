# -*- coding: utf-8 -*-

"""
Useful conversion and distance related function

@author: pycaz authors
"""

import numpy as np

from .typing import ArrayLike, Tuple


def hpa2pa(hpa: ArrayLike) -> ArrayLike:
    """
    Takes pressure value in hecta Pascal and return in Pascal

    :param hpa: Pressure in hPa
    :return: Pressure in Pascal (Pa)
    """
    return hpa * 100


def pa2mb(pa: ArrayLike) -> ArrayLike:
    """
    Takes pressure as Pa and return in mili bar (hPa)

    :param pa: Pressure in Pascal (Pa)
    :return: Pressure in mili bar
    """
    return pa / 100.0


def knot2mps(knot: ArrayLike) -> ArrayLike:
    """
    Takes velocity in knot and returns velocity in mps.

    :param knot: Velocity in knot
    :return: Velocity in meters per second
    """
    return knot * 1.852 / 3.6


def km2m(km: ArrayLike) -> ArrayLike:
    """
    Takes distance in Kilometers and converts it to meter.

    :param km: distance in kilometers
    :return: distance in meters
    """
    return km * 1000


def ft2m(ft: ArrayLike) -> ArrayLike:
    """
    Takes distance in Ft and converts it to meter

    :param ft: Distance in feet
    :return: Distance in meter
    """
    return ft * 0.3048


def ntm2m(ntm: ArrayLike) -> ArrayLike:
    """
    Takes distance in Nautical Mile and converts to meter.

    :param ntm: Distance in Nautical Mile
    :return: Distance in meter
    """
    return ntm * 1852.001


def lon180(lon360: ArrayLike, sort_array: bool = False) -> ArrayLike:
    """
    Change lon formatting range from 0-360 to -180-180.

    It does not change the order. Use sort_array=True for ordering.

    :param lon360: Input longitude value/array in degrees
    :param sort_array: If output needs to be sorted, default is False
    :return: Longitude value/array in degrees in -180 to 180 range
    """
    lon360[lon360 > 180] = lon360[lon360 > 180] - 360

    if sort_array:
        lon360 = np.sort(lon360)

    return lon360


def gc_distance(
        of_x: ArrayLike,
        of_y: ArrayLike,
        origin_x: float,
        origin_y: float,
        isradians: bool = False) -> Tuple[ArrayLike, ArrayLike]:
    """Calculates the great circle distance of 'of' from 'origin'

    of: list of lon lat of the point
    origin: list of lon lat of the origin point

    :param of_x: longitude(s)
    :param of_y: latitude(s)
    :param origin_x: origin longitude
    :param origin_y: origin latitude
    :param isradians: if input is in radians, default is False
    :return:
    """
    dfac = 60 * 1.852 * 1000

    if isradians:
        dtrans_x = dfac * np.cos(origin_y) * (np.rad2deg(of_x) - np.rad2deg(origin_x))
        dtrans_y = dfac * (np.rad2deg(of_y) - np.rad2deg(origin_y))
    else:
        dtrans_x = dfac * np.cos(np.deg2rad(origin_y)) * (of_x - origin_x)
        dtrans_y = dfac * (of_y - origin_y)

    return dtrans_x, dtrans_y
