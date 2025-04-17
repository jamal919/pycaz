# -*- coding: utf-8 -*-

"""
Implements the classes that provides the reading and writing functionalities of the boundary
conditions.
"""

consts = [
    'M2', 'M3', 'M4', 'M6', 'M8', 'MF', 'MM',
    'MN4', 'MS4', 'MSF', 'MU2', 'N2', 'NU2',
    'O1', 'P1', 'Q1', 'R2', 'S1', 'S2', 'S4',
    'SSA', 'T2', 'K2', 'K1', 'J1', '2N2']


class Iettype1:
    pass


class Iettype2:
    pass


class Iettype3:
    pass


class Iettype4:
    """
    netcdf elev2D.th {
        dimensions:
            time = UNLIMITED ; // (73 currently)
            nOpenBndNodes = 748 ;
            nLevels = 1 ;
            nComponents = 1 ;
            one = 1 ;
        variables:
            double time_series(time, nOpenBndNodes, nLevels, nComponents) ;
            float time_step(one) ;
            double time(time) ;
        }
    """
    pass


class Iettype5:
    pass


class Ifltype1:
    pass


class Ifltype2:
    pass


class Ifltype3:
    pass


class Ifltype4:
    """
    netcdf uv3D.th {
        dimensions:
            nOpenBndNodes = 748 ;
            one = 1 ;
            time = UNLIMITED ; // (73 currently)
            nLevels = 44 ;
            nComponents = 2 ;
        variables:
            float time_step(one) ;
            double time(time) ;
            float time_series(time, nOpenBndNodes, nLevels, nComponents) ;
        }
    """
    pass


class Ifltype5:
    pass


class Ifltype_1:
    pass
