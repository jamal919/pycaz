#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pkg_resources

try:
    __version__ = pkg_resources.get_distribution('pycaz').version
except Exception:
    __version__ = 'unknown'