# -*- coding: utf-8 -*-

"""
Typical types used in pycaz codes.

@author: pycaz authors
"""

import os
from typing import Tuple, List, TypeAlias
from numpy.typing import ArrayLike
from pandas._typing import TimestampConvertibleTypes
from pandas._typing import TimedeltaConvertibleTypes
from pathlib import Path


try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

PathLike: TypeAlias = str | bytes | os.PathLike | Path
ArrayLike: TypeAlias = ArrayLike
Literal: TypeAlias = Literal
TimestampConvertibleTypes: TypeAlias = TimestampConvertibleTypes
TimedeltaConvertibleTypes: TypeAlias = TimedeltaConvertibleTypes
List:TypeAlias = List

