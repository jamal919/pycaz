#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Typical types used in pycaz codes.

@author: pycaz authors
"""

import os
from typing import Union, TypeAlias
from numpy.typing import ArrayLike

PathLike: TypeAlias = str | bytes | os.PathLike
ArrayLike: TypeAlias = ArrayLike

