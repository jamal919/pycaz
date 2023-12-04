#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
The `schism` module of `pycaz` concerns handling of SCHISM model, paricularly the input
in SCHISM model. The model code is available at https://github.com/schism-dev/schism. An
official toolbox is available at https://github.com/schism-dev/pyschism. The functionalities
of the current toolbox is not one-to-one, and likely much less than the official pyschism.
However, the philosophy of the current toolbox is comparatively simple, and based on the
idea that research activities requires more freedom.

Further documentation will be added here in the future.

'''

from pycaz.schism.hgrid import read_gr3
from pycaz.schism.hgrid import read_hgrid