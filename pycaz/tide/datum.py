# -*- coding: utf-8 -*-

import numpy as np
from typing import List
from numpy.typing import ArrayLike
from pycaz.tide.utilities import coutier_classification, dirunal_ineuality


def compute_distance_msl2cd(
        m2: List[ArrayLike, ArrayLike],
        s2: List[ArrayLike, ArrayLike],
        o1: List[ArrayLike, ArrayLike],
        k1: List[ArrayLike, ArrayLike],
        n2: List[ArrayLike, ArrayLike],
        k2: List[ArrayLike, ArrayLike],
        p1: List[ArrayLike, ArrayLike]) -> ArrayLike:
    """
    Compute the distance from MSL to chart datum (always positive distance) using the Baley, 1952.

    Ref: Balay, M. A. (1952). Determination of Plane of Reduction of Soundings in any Place.
    The International Hydrographic Review.

    :param m2: (amplitude, phase) of m2.
    :param s2: (amplitude, phase) of s2.
    :param o1: (amplitude, phase) of o1.
    :param k1: (amplitude, phase) of k1.
    :param n2: (amplitude, phase) of n2.
    :param k2: (amplitude, phase) of k2.
    :param p1: (amplitude, phase) of p1.
    :return: Distance from MSL to chart datum.
    """

    # splitting amplitude and phase data
    m2_a, m2_g = m2
    s2_a, s2_g = s2
    o1_a, o1_g = o1
    k1_a, k1_g = k1
    n2_a, n2_g = n2
    k2_a, k2_g = k2
    p1_a, p1_g = p1

    # ensuring np.array
    m2_a, m2_g = np.asarray(m2_a), np.asarray(m2_g)
    s2_a, s2_g = np.asarray(s2_a), np.asarray(s2_g)
    o1_a, o1_g = np.asarray(o1_a), np.asarray(o1_g)
    k1_a, k1_g = np.asarray(k1_a), np.asarray(k1_g)
    n2_a, n2_g = np.asarray(n2_a), np.asarray(n2_g)
    k2_a, k2_g = np.asarray(k2_a), np.asarray(k2_g)
    p1_a, p1_g = np.asarray(p1_a), np.asarray(p1_g)

    # computing coefficients
    c = coutier_classification(m2_a, s2_a, o1_a, k1_a)
    two_k = dirunal_ineuality(m2_g, o1_g, k1_g)

    msl2cd = np.full(shape=m2_a.shape, fill_value=np.nan)

    # eqn 3, 0 < C <= 0.25
    idx = np.logical_and(c > 0, c <= 0.25)
    msl2cd[idx] = m2_a[idx] + s2_a[idx] + n2_a[idx] + k2_a[idx]

    # 0.25 < C <= 1.5
    idx = np.logical_and(c > 0.25, c <= 1.5)
    msl2cd[idx] = m2_a[idx] + s2_a[idx] + k1_a[idx] + o1_a[idx] + p1_a[idx]
    # particular case: two_K == 0
    idx_two_k0 = np.logical_and(idx, two_k == 0)
    msl2cd[idx_two_k0] = m2_a[idx_two_k0] + s2_a[idx_two_k0] + n2_a[idx_two_k0]
    # particular case: two_K == 180
    idx_two_k180 = np.logical_and(idx, two_k == 180)
    msl2cd[idx_two_k180] = m2_a[idx_two_k180] + s2_a[idx_two_k180] + n2_a[idx_two_k180]

    # 1.5 < C <= 3
    idx = np.logical_and(c > 1.5, c <= 3)
    msl2cd[idx] = m2_a[idx] + s2_a[idx] + k1_a[idx] + o1_a[idx]

    # C > 3
    idx = c > 3
    msl2cd[idx] = m2_a[idx] + s2_a[idx] + k1_a[idx] + o1_a[idx] + p1_a[idx]

    return msl2cd
