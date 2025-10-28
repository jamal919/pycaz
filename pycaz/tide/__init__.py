# -*- coding: utf-8 -*-

"""
# pycaz.tide
Tide modules provides related functionalities for reading/writing tide atlas data, interpolate it, computing tidal datum,
applying tidal filter. It also provides some custom functions to compute nodal corrections from utide and generate
waterlevel using utide in a vectorized fashion.
"""

from .atlas import read_atlas
