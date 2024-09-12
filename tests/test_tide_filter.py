#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycaz.tide.filters import filters
import numpy as np


def test_filter():
    """
    Sum of filter should be 1 (or very close to it)

    :return:
    """
    for f in filters:
        filtertotal = np.sum(filters[f]['Filter'] / filters[f]['Denom'])
        print(f, filtertotal)
        np.testing.assert_almost_equal(filtertotal, 1, decimal=5)
