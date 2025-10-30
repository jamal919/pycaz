# -*- coding: utf-8 -*-

"""
Typical types used in pycaz codes.

@author: pycaz authors
"""

import os
from typing import List, TypeAlias
from numpy.typing import ArrayLike
from pandas._typing import TimestampConvertibleTypes


try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

PathLike: TypeAlias = str | bytes | os.PathLike
ArrayLike: TypeAlias = ArrayLike
Literal: TypeAlias = Literal
TimestampConvertibleTypes: TypeAlias = TimestampConvertibleTypes
List:TypeAlias = List

