# -*- coding: utf-8 -*-
'''
tide_fac is a program to compute nodal factors and equilibrium arguments.

based on tide_fac.f of the NOS/NOAA
'''

import numpy as np

class Tidefac(object):

    def __init__(self, year, month, day, hour, runtime):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.runtime = runtime

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
                                    ['S6', 90.0, 6]])

        self.constituents = np.array([, , , , , , , , 
                                , , , , 'mu2', '2N2', 'OO1', 'lambda2', 
                                'S1', 'M1', 'J1', 'Mm', 'Ssa', 'Sa', 'MSf', 'Mf', 
                                'rho1', 'Q1', 'T2', 'R2', '2Q1', 'P1', '2SM2', 'M3', 
                                'L2', '2MK3', 'K2', 'M8', 'MS4'])
        self.speed = np.array([,,,, ,,,,
                            ,,,, 27.9682084,27.8953548,16.139101680,29.4556253,
                            15.0,14.496693984,15.5854433,0.5443747, 0.0821373,0.0410686,1.0158957720,1.0980331,
                            13.4715145,13.3986609,29.9589333,30.0410667, 12.854286252,14.9589314,31.01589576,43.476156360,
                            29.5284789,42.927139836,30.0821373,115.936416972, 58.984104240], dtype=float)
        self.ncycle = np.array([ , , , , , , ,,
                                 , , , ,2 ,2 ,1 ,2 ,
                                1 ,1, 1, 0, 0, 0, 0, 0, 
                                1, 1, 2, 2, 1, 1, 2, 3, 
                                2, 3, 2, 8, 4], dtype=float)

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