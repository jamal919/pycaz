# -*- coding: utf-8 -*-
'''
tide_fac is a program to compute nodal factors and equilibrium arguments.

the intention of this script is to replace the tide_fac program distributed
with the SCHISM model and provide a flexible way to incorporate the 
functionality in other scripted program.
'''

import numpy as np

class Tidefac(object):

    def __init__(self, year, month, day, hour, runtime):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.runtime = runtime

        # Constituents as [Constituent, Speed, Number of Cycle]
        self.constituents = np.array([['M2', 28.984104252, 2],
                                    ['S2', 30.0000000, 2],
                                    ['N2', 28.439729568, 2],
                                    ['K1', 15.041068632, 1],
                                    ['M4', 57.968208468, 4],
                                    ['O1', 13.943035584, 1],
                                    ['M6', 86.952312720, 6],
                                    ['MK3', 44.025172884, 3],
                                    ['S4', 60.0, 4],
                                    ['MN4', 57.423833820, 4],
                                    ['nu2', 28.5125831, 2],
                                    ['S6', 90.0, 6],
                                    ['mu2', 27.9682084, 2],
                                    ['2N2', 27.8953548, 2],
                                    ['OO1', 16.139101680, 1],
                                    ['lambda2', 29.4556253, 2],
                                    ['S1', 15.0, 1],
                                    ['M1', 14.496693984, 1],
                                    ['J1', 15.5854433, 1],
                                    ['Mm', 0.5443747, 0],
                                    ['Ssa', 0.0821373, 0],
                                    ['Sa', 0.0410686, 0],
                                    ['MSf', 1.0158957720, 0],
                                    ['Mf', 1.0980331, 0],
                                    ['rho1', 13.4715145, 1],
                                    ['Q1', 13.3986609, 1],
                                    ['T2', 29.9589333, 2],
                                    ['R2', 30.0410667, 2],
                                    ['2Q1', 12.854286252, 1],
                                    ['P1', 14.9589314, 1],
                                    ['2SM2', 31.01589576, 2],
                                    ['M3', 43.476156360, 3],
                                    ['L2', 29.5284789, 2],
                                    ['2MK3', 42.927139836, 3],
                                    ['K2', 30.0821373, 2],
                                    ['M8', 115.936416972, 8],
                                    ['MS4', 58.984104240, 4]])

    @staticmethod
    def arctan(parameter_list):
        pass


    @staticmethod
    def julianday(year, month, day):
        dayt = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334], dtype=float)
        days = np.zeros(shape=12)
        days[[0, 1]] = [0, 31]
        dinc = 0
        leapyear = (year-1900)%4
        if not leapyear:
            dinc = 1

        for i in np.arange(2, 12):
            days[i] = dayt[i] + dinc

        return(days[int(month)-1]+day)

if __name__=='__main__':
    tidefac = Tidefac(2007, 12, 12, 0, 30)
    print(tidefac.julianday(2007, 12, 12))
    print(tidefac.constituents[1, 0])