#!/usr/bin/env python
# -*- coding: utf-8 -*-
from copy import deepcopy
import numpy as np
import pandas as pd
import logging
from pycaz.tide.utide import utide_names, nodal_factor

logger = logging.getLogger(__name__)


class Tidefac(dict):
    def __init__(self, **kwargs):
        """ A bctides object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.update(kwargs)

    def copy(self):
        return (deepcopy(self))

    @property
    def consts(self):
        return (self['const'])

    @property
    def rnday(self):
        return (self['rnday'])

    @property
    def start_date(self):
        year = self['year']
        month = self['month']
        day = self['day']
        hour = self['hour'] // 1
        minute = (self['hour'] % 1 * 60)
        second = np.round((minute % 1) * 60)
        minute = minute // 1

        return (f'{int(year):4d}-{int(month):02d}-{int(day):02d} {int(hour):02d}:{int(minute):02d}:{int(second):02d}')

    @property
    def info(self):
        rnday = self['rnday']
        start_date = self.start_date
        return (f'{rnday:.2f} days from {start_date} UTC')

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
    _const = {i[0].upper(): {'nf': float(i[1]), 'ear': float(i[2])} for i in _const}
    tidefac['const'] = _const

    return tidefac


def generate_tidefac_fortran(start_date, rnday=30, end_date=None) -> Tidefac:
    """
    Generate tidefac using the tide_fac.f fortran program.

    :param start_date: Start date of the simulation.
    :param rnday: Length of the simulation in days, default 30.
    :param end_date: End date of the simulation.
    :return:

    """
    raise NotImplementedError


def generate_tidefac_utide(consts, start_date, rnday=30, end_date=None, lat=None) -> Tidefac:
    """
    Generate tidefac from at the mid of the start and end day.

    :param consts: List of consts for which the nodal factors are to be computed.
    :param start_date: Start date of the simulation.
    :param rnday: Length of simulation in days, default 30.
    :param end_date: End date of the simulation.
    :param lat: Center latitude of the model grid.
    :return: Tidefac object
    """
    start_date = pd.to_datetime(start_date)
    if rnday <= 0:
        raise ValueError('rnday must be greater than 0!')
    if end_date is not None:
        end_date = pd.to_datetime(end_date)
    if end_date is None:
        end_date = start_date + pd.Timedelta(days=rnday)
    if end_date <= start_date:
        raise ValueError('end_date must be greater than start_date!')
    if lat is None:
        lat = 5
        logger.info('Lat is set to 5°N')

    available, missing = utide_names(consts)
    const_nodal_factor = nodal_factor(
        t=[start_date, end_date],
        consts=available,
        lat=lat,
        correct_phase=True)
    tidefac = Tidefac(
        year=start_date.year,
        month=start_date.month,
        day=start_date.day,
        hour=start_date.hour,
        rnday=rnday,
        const=const_nodal_factor
    )

    return tidefac
