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
import re

class Bctides(object):
    def __init__(self, filepath):
        self.filepath = filepath

    def read(self):
        with open(self.filepath) as f:
            ds = f.readlines()
            # First the dates
            self.info = ds[0].split('\n')[0]
            print(self.info)
            __lnproc = 0

            # Then the tidal potential information
            self.ntip, self.tip_dp = np.fromstring(ds[1].split('!')[0], count=2, sep=' ')
            self.ntip = int(self.ntip)
            __lnproc = 1
            
            self.tip = []
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
    def __init__(self, filepath):
        self.filepath = filepath
        self.r = re.compile(r'[0-9]?.[0-9]\d+')
        self.__read()
    
    def __read(self):
        # Reading date information
        with open(self.filepath, 'r') as f:
            # Reading the date section
            __ds = f.readline()
            __date = np.array(self.r.findall(__ds), dtype='float')
            self.hour = __date[0]
            self.day = __date[1]
            self.month = __date[2]
            self.year = __date[3]            
            
            # Reading the run length section
            __ds = f.readline()
            __rnday = self.r.findall(__ds)
            
            if(len(__rnday)==0):
                __ds = f.readline()
                
            __rnday = np.array(self.r.findall(__ds), dtype='float')
            self.rnday = __rnday[0]

        # Reading the constants, node factor and eq. argument ref. to GM in deg.
        __const = np.genfromtxt(fname=self.filepath, dtype=None, skiprows=8, \
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
    bctides = Bctides(filepath=bctide_source).read()
    tfac = Tidefacout(filepath=tfacfile)
    bctides.update(tfac)
    bctides.write(filepath=bctide_update)