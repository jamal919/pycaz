# -*- coding: utf-8 -*-
"""
Reads tidefac output and adapt given bctides file as required by the tidefac 
outputs.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
import sys
import re

class Bctides(object):
    def __init__(self, info='', ntip=0, tip_dp=0, tip=[], nbfr=0, bfr=[], nope=0, boundaries=[]):
        self.info = info
        self.nitp = ntip
        self.tip_dp = tip_dp
        self.tip = tip
        self.nbfr = nbfr
        self.bfr = bfr
        self.nope = nope
        self.boundaries = boundaries

    def read(self, filepath):
        with open(filepath) as f:
            ds = f.readlines()
            # First the dates
            self.info = ds[0].split('\n')[0]
            __lnproc = 0

            # Then the tidal potential information
            self.ntip, self.tip_dp = np.fromstring(ds[1].split('!')[0], count=2, sep=' ')
            self.ntip = int(self.ntip)
            __lnproc = 1
            
            for i in np.arange(self.ntip):
                __talpha = ds[__lnproc+1].split('\n')[0]
                __jspc, __tamp, __tfreq, __tnf, __tear = np.fromstring(ds[__lnproc+2].split('\n')[0], count=5, sep=' ')
                __rec = dict(talpha=__talpha, jspc=__jspc, tamp=__tamp, tfreq=__tfreq, tnf=__tnf, tear=__tear)
                self.tip.append(__rec)
                __lnproc = __lnproc + 2
            
            # Reading the boundary frequencies
            self.nbfr = np.fromstring(ds[__lnproc+1], count=1, sep=' ')
            self.nbfr = int(self.nbfr)
            __lnproc = __lnproc + 1
            
            self.bfr = []
            for i in np.arange(self.nbfr):
                __alpha = ds[__lnproc+1].split('\n')[0]
                __amig, __ff, __face = np.fromstring(ds[__lnproc+2].split('\n')[0], count=3, sep=' ')
                __rec = dict(alpha=__alpha, amig=__amig, ff=__ff, face=__face)
                self.bfr.append(__rec)
                __lnproc = __lnproc + 2
            
            # Open boundary sagments
            self.nope = ds[__lnproc+1].split(' ')[0]
            self.nope = int(self.nope)
            __lnproc = __lnproc + 1

            # For each open boundary sagment
            self.boundaries = ds[__lnproc+1:len(ds)]

    def update(self, tidefac):
        # Update time
        self.info = tidefac.info
        # Updating the tidal potential nodal factor and equilibrium argument
        for __tip in self.tip:
            __talpha = __tip['talpha'].strip().upper()
            if __talpha in tidefac.const.keys():
                __tip['tnf'] = tidefac.const[__talpha][0]
                __tip['tear'] = tidefac.const[__talpha][1]

        # Updating the Boundary frequency nodal factors and equilibrium argument
        for __bfr in self.bfr:
            __alpha = __bfr['alpha'].strip().upper()
            if __alpha in tidefac.const.keys():
                __bfr['ff'] = tidefac.const[__alpha][0]
                __bfr['face'] = tidefac.const[__alpha][1]

    def write(self, filepath):
        with open(filepath, 'w') as f:
            # Header information
            f.write('{:s}\n'.format(self.info))

            # Tidal potential
            f.write('{:d} {:3.2f} !ntip, tip_dp\n'.format(int(self.ntip), float(self.tip_dp)))

            for __tip in self.tip:
                f.write('{:s}\n{:d}\t{:.6f}\t{:.16f}\t{:.5f}\t{:.2f}\n'\
                        .format(__tip['talpha'].strip().upper(),\
                                int(__tip['jspc']),\
                                __tip['tamp'],\
                                __tip['tfreq'],\
                                __tip['tnf'],\
                                __tip['tear']))

            # Boundary frequencies
            f.write('{:d} !nbfr\n'.format(int(self.nbfr)))

            for __bfr in self.bfr:
                f.write('{:s}\n{:.16E}\t{:.6f}\t{:.2f}\n'\
                        .format(__bfr['alpha'].strip().upper(),\
                                __bfr['amig'],\
                                __bfr['ff'],\
                                __bfr['face']))

            # Open boundaries
            f.write('{:d} !Number of Open Boundaries\n'.format(self.nope))
            
            for __line in self.boundaries:
                f.write(__line)

class Tidefacout(object):
    def __init__(self, year=0, month=0, day=0, hour=0, rnday=0, const={}):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.rnday = rnday
        self.const = const
    
    def read(self, filepath):
        # Reading date information
        with open(filepath, 'r') as f:
            # Reading the date section
            __ds = f.readline()
            __date = np.fromstring(__ds, dtype=float, count=4, sep=',')
            self.year = __date[0]
            self.month = int(__date[1])
            self.day = int(__date[2]) 
            self.hour = int(__date[3])
            
            # Reading the run length section
            __ds = f.readline()
            __rnday = np.fromstring(__ds, dtype=float, count=1, sep=',')
            self.rnday = __rnday[0]

        # Reading the constants, node factor and eq. argument ref. to GM in deg.
        __const = np.genfromtxt(fname=filepath, dtype=None, skiprows=6, \
                                delimiter=None, autostrip=True)
        __const = np.array([[i for i in j] for j in __const])
        __const = {i[0].upper():[float(j) for j in i[1:3]] for i in __const}
        self.const = __const

        # Tidefac header information
        self.info = '{:.2f} days - {:4.0f}/{:02.0f}/{:02.0f} {:02.2f} UTC'.format(self.rnday,\
                self.year, self.month, self.day, self.hour)

    def __str__(self):
        return(self.info)

if __name__=='__main__':
    bctide_source = 'bctides.ini'
    bctide_update = 'bctides.in'
    tfacfile = 'tide_fac.out'
    bctides = Bctides()
    bctides.read(filepath=bctide_source)
    tfac = Tidefacout()
    tfac.read(filepath=tfacfile)
    bctides.update(tfac)
    bctides.write(filepath=bctide_update)
