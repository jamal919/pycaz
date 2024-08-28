#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl


class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = np.max([0, 1 / 2 * (1 - np.abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax)))])
        normalized_max = np.min([1, 1 / 2 * (1 + np.abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin)))])
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))

