#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reads tidefac output and adapt given bctides file as required by the tidefac 
outputs.

This script uses a dictionary to store the constituent values. The problem with
dictionary is it does not maintain the order of the dictionary key as they are 
added. So consequently when we used the distionary structure on tidal constituent
information and not on the boundary information, it produced a tidefac file which
is not consistent - the listing of constituent in the tidal potential portion
is not the same as the listing of the tidal amplitude and phases at the boundary
nodes. The soltuion from this is to read also the boundaries while reading the 
initial constituent and write them as a whole data structure or keep track of the
constituents as they are read from the file. Currently none is implemented and
it is better not to use this version until this problem is fixed. 

It is to be noted that, the reason behing changing the data structure is to make
use of more pythonic way writing codes using dictionaries, instead of numpy array.

TODO: Read boundary information

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
import os
from collections import OrderedDict

class Bctides(object):
    def __init__(
        self,
        info='', 
        ntip=0, 
        tip_dp=0, 
        tip=OrderedDict(), 
        nbfr=0, 
        bfr=OrderedDict(), 
        nope=0, 
        boundaries=[]
    ):
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
                print('Reading Wave {:d} - {:s}'.format(i, __talpha))
                __jspc, __tamp, __tfreq, __tnf, __tear = np.fromstring(ds[__lnproc+2].split('\n')[0], count=5, sep=' ')
                self.tip[__talpha.strip().upper()] = OrderedDict(jspc=int(__jspc), tamp=__tamp, tfreq=__tfreq, tnf=__tnf, tear=__tear)
                __lnproc = __lnproc + 2
            
            # Reading the boundary frequencies
            self.nbfr = np.fromstring(ds[__lnproc+1], count=1, sep=' ')
            self.nbfr = int(self.nbfr)
            __lnproc = __lnproc + 1
            
            for i in np.arange(self.nbfr):
                __alpha = ds[__lnproc+1].split('\n')[0]
                __amig, __ff, __face = np.fromstring(ds[__lnproc+2].split('\n')[0], count=3, sep=' ')
                self.bfr[__alpha.strip().upper()] = OrderedDict(amig=__amig, ff=__ff, face=__face)
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
        for talpha in self.tip.keys():
            if talpha in tidefac.const.keys():
                self.tip[talpha]['tnf'] = tidefac.const[talpha][0]
                self.tip[talpha]['tear'] = tidefac.const[talpha][1]

        # Updating the Boundary frequency nodal factors and equilibrium argument
        for alpha in self.bfr.keys():
            if alpha in tidefac.const.keys():
                self.bfr[alpha]['ff'] = tidefac.const[alpha][0]
                self.bfr[alpha]['face'] = tidefac.const[alpha][1]

    def write(self, filepath):
        with open(filepath, 'w') as f:
            # Header information
            f.write('{:s}\n'.format(self.info))

            # Tidal potential
            f.write('{:d} {:3.2f} !ntip, tip_dp\n'.format(int(self.ntip), float(self.tip_dp)))

            for alpha in self.tip.keys():
                f.write('{:s}\n{:d}\t{:.6f}\t{:.16f}\t{:.5f}\t{:.2f}\n'\
                        .format(alpha,\
                                int(self.tip[alpha]['jspc']),\
                                self.tip[alpha]['tamp'],\
                                self.tip[alpha]['tfreq'],\
                                self.tip[alpha]['tnf'],\
                                self.tip[alpha]['tear']))

            # Boundary frequencies
            f.write('{:d} !nbfr\n'.format(int(self.nbfr)))

            for alpha in self.bfr.keys():
                f.write('{:s}\n{:.16E}\t{:.6f}\t{:.2f}\n'\
                        .format(alpha,\
                                self.bfr[alpha]['amig'],\
                                self.bfr[alpha]['ff'],\
                                self.bfr[alpha]['face']))

            # Open boundaries
            f.write('{:d} !Number of Open Boundaries\n'.format(self.nope))
            
            for __line in self.boundaries:
                f.write(__line)

class Tidefacout(object):
    def __init__(self, year=0, month=0, day=0, hour=0, rnday=0, const=OrderedDict()):
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
        __const = np.genfromtxt(fname=filepath, dtype=None, skip_header=6, \
                                delimiter=None, autostrip=True)
        __const = np.array([[i for i in j] for j in __const])
        __const = OrderedDict({i[0].upper():[float(j) for j in i[1:3]] for i in __const})
        self.const = __const

        # Tidefac header information
        self.info = f'{self.rnday:.2f} days - {self.year:4.0f}/{self.month:02.0f}/{self.day:02.0f} {self.hour:02.2f} UTC'
    def __str__(self):
        return(self.info)

if __name__=='__main__':
    path = '/run/media/khan/Workbench/Projects/Surge Model/Bctides'
    bctide_source = os.path.join(path, 'bctides.ini')
    bctide_update = os.path.join(path, 'bctides.in')
    tfacfile = os.path.join(path, 'tide_fac.out')
    bctides = Bctides()
    bctides.read(filepath=bctide_source)
    tfac = Tidefacout()
    tfac.read(filepath=tfacfile)
    bctides.update(tfac)
    bctides.write(filepath=bctide_update)
