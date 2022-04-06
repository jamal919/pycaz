#!/usr/bin/env python
# -*- coding: utf-8 -*-
from copy import deepcopy
import numpy as np

class Tidefac(dict):
    def __init__(self, **kwargs):
        """ A bctides object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.update(kwargs)

    def copy(self):
        return(deepcopy(self))
    
    @property
    def consts(self):
        return(self['const'])
    
    @property
    def rnday(self):
        return(self['rnday'])
    
    @property
    def start_date(self):
        year = self['year']
        month = self['month']
        day = self['day']
        hour = self['hour'] // 1
        minute = (self['hour'] % 1 * 60)
        second = np.round((minute % 1) * 60)
        minute = minute // 1

        return(f'{int(year):4d}-{int(month):02d}-{int(day):02d} {int(hour):02d}:{int(minute):02d}:{int(second):02d}')

    @property
    def info(self):
        rnday = self['rnday']
        start_date = self.start_date
        return(f'{rnday:.2f} days run starting from {start_date} UTC')

    def describe(self):
        print(self.info)

def read_tidefacout(fname: str) -> Tidefac:
    """
    A reader for the output from tide_fac.f program, used for generate tidefac information.
    """
    tidefac = Tidefac()
    # Reading date information
    with open(fname, 'r') as f:
        # Reading the date section
        _date = np.fromstring(f.readline(), dtype=float, count=4, sep=',')
        tidefac['year'] = _date[0]
        tidefac['month'] = int(_date[1])
        tidefac['day'] = int(_date[2]) 
        tidefac['hour'] = int(_date[3])
        
        # Reading the run length section
        tidefac['rnday'] = float(f.readline().strip())

    # Reading the constants, node factor and eq. argument ref. to GM in deg.
    _const = np.genfromtxt(fname=fname, dtype=None, skip_header=6, delimiter=None, autostrip=True, encoding='UTF8')
    _const = np.array([[i for i in j] for j in _const])
    _const = {i[0].upper():{'nf':float(i[1]), 'ear':float(i[2])} for i in _const}
    tidefac['const'] = _const

    return(tidefac)