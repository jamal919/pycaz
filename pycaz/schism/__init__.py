# -*- coding: utf-8 -*-
"""
# pycaz.schism
The `schism` module of `pycaz` concerns handling of SCHISM model, paricularly the input
in SCHISM model. The model code is available at https://github.com/schism-dev/schism. An
official toolbox is available at https://github.com/schism-dev/pyschism.
"""

from .hgrid import read_gr3
from .hgrid import read_hgrid
from .staout import read_staout