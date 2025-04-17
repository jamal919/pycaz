# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
from cartopy.io.img_tiles import GoogleTiles
from typing import Literal, get_args

ESRILayer = Literal[
    'NatGeo_World_Map',
    'World_Imagery',
    'World_Physical_Map',
    'World_Shaded_Relief',
    'World_Street_Map',
    'World_Terrain_Base',
    'World_Topo_Map'
]


def esri_tiles(layer: ESRILayer = 'World_Imagery', cache: bool = False) -> GoogleTiles:
    """
    Returns an ESRI Tiles object to be used with cartopy add_image()

    :param layer: Name of the ESRI tile layer, see ESRI_LAYERS list
    :param cache: If the files to be saved in cache, default is False
    """
    if layer not in get_args(ESRILayer):
        raise ValueError(f'Layer must be one of {ESRILayer}')

    url = f'https://server.arcgisonline.com/ArcGIS/rest/services/{layer}/' + 'MapServer/tile/{z}/{y}/{x}.jpg'

    return GoogleTiles(url=url, cache=cache)


class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        """
        A normalize class for color scaling with a defined midpoint

        :param vmin: minimum value
        :param vmax: maximum value
        :param midpoint: midpoint value
        :param clip: True/False passed to Normalize class
        """
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = np.max([0, 1 / 2 * (1 - np.abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax)))])
        normalized_max = np.min([1, 1 / 2 * (1 + np.abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin)))])
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))
