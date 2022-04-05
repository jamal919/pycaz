#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pycaz package __init__ script

@author: khan
"""

import pkg_resources

try:
    __version__ = pkg_resources.get_distribution('pycaz').version
except Exception:
    __version__ = 'unknown'