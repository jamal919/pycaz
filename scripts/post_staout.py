#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script reads a staout file, station.in file and generate named time series.
It can write a timeseries to REFMAR to be used with comodo-tools, or can directly
be used with utide.

author: Md Jamal Uddin Khan
email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
from subprocess import call
import os

class Param(object):
    def __init__(self, params=None):
        self.params = params

    def read(self, fname):
        __params = np.genfromtxt(fname, dtype=None, comments='!', delimiter='=', autostrip=True)
        __params = dict(__params)
        self.params = __params

    def get_param(self, param):
        if(self.params.has_key(param)):
            return(self.params[param])
        else:
            print('Key {:s} not found in the param dictionary!'.format(param))
            return(None)

    def get_starttime(self):
        __year = int(self.get_param('start_year'))
        __month = int(self.get_param('start_month'))
        __day = int(self.get_param('start_day'))
        __dhour = float(self.get_param('start_hour'))
        __hour = int(__dhour)
        __minutes = __dhour*60 % 60
        __seconds = __dhour*3600 % 60
        return(datetime(year=__year, month=__month, day=__day, hour=__hour))
        

class Timeseries(object):
    def __init__(self, timestamps, values):
        self.timestamps = timestamps
        self.values = values

    def to_num_array(self):
        __tmstr = self.timestamps
        __values = self.values
        __ts = np.transpose(np.array([__tmstr, __values]))
        return(__ts)

    def to_strftime_array(self, fmt='%Y/%m/%d %H:%M:%S'):
        __tmstr = np.array([i.strftime(fmt) for i in self.timestamps])
        __values = self.values
        __ts = np.transpose(np.array([__tmstr, __values]))
        return(__ts)

class Station(object):
    def __init__(self, id, x, y, z, name, ts=None):
        self.id = int(id)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.name = str(name)
        self.ts = ts

    def __lt__(self, other):
        return(self.id < other.id)

class Stations(object):
    def __init__(self):
        self.list = np.array([])
        self.nstation = 0
        self.flag = np.array([])

    def read(self, fname):
        try:
            with open(fname) as f:
                __ds = f.readlines()
        except:
            print('Error in station.in file!')
        else:
            self.flag = np.fromstring(__ds[0], dtype=int, count=9, sep=' ')
            self.nstation = int(__ds[1].split('\n')[0])
            __list = np.genfromtxt(fname, dtype=None, skip_header=2, autostrip=True)
            self.list = np.array([Station(id=entry[0], x=entry[1], y=entry[2], z=entry[3], name=entry[4]) for entry in __list])
        finally:
            print('Reading completed!')
            print('No of stations       = {:d}'.format(self.nstation))
            print('No of station entry  = {:d}'.format(len(self.list)))

class Staout(object):
    def __init__(self, stations, param, exp='EXP'):
        self.stations = stations
        self.param = param
        self.exp = exp

    def read(self, fname):
        __ds = np.genfromtxt(fname, dtype=float, autostrip=True)
        __time = np.array([self.param.get_starttime() + timedelta(seconds=s) for s in __ds[:, 0]])

        for i, station in enumerate(stations.list):
            station.ts = Timeseries(timestamps=__time, values=__ds[:, i+1])

    def write(self, path, format='REFMAR'):
        if not os.path.exists(path):
            os.mkdir(path)

        for i, station in enumerate(stations.list):
            print('Processing Station {:d} - {:s}'.format(i, station.name))
            fname = os.path.join(path, station.name)
            with open(fname, 'w') as f:
                f.write('# Station: {:s}\n'.format(station.name))
                f.write('# Longitude: {:f}\n'.format(station.x))
                f.write('# Latitude: {:f}\n'.format(station.y))
                f.write('# Unit: {:s}\n'.format('m'))
                f.write('# Experiment: {:s}\n'.format(self.exp))
                f.write('# rnday: {:.1f}\n'.format(float(self.param.params['rnday'])))
                f.write('# dt: {:.1f}\n'.format(float(self.param.params['dt'])))
                
                np.savetxt(fname=f, X=station.ts.to_strftime_array('%Y/%m/%d %H:%M:%S'), fmt='%s', delimiter=' ')

if __name__=='__main__':
    path = '/home/khan/MEGA/Models/SCHISM/Toy'
    param = Param()
    param.read(os.path.join(path, 'param.in'))
    
    stations = Stations()
    stations.read(os.path.join(path, 'station.in'))

    staout = Staout(stations=stations, param=param, exp='EXP07_Sa_35cm')
    staout.read(os.path.join(path, 'staout_1'))
    staout.write(path=os.path.join(path, 'timeseries'))