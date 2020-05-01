#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cyclone tracks object, wind field models

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
from .conversion import gc_distance
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
import sys

class Record(dict):
    def __init__(self, *args, **kwargs):
        '''
        A recrod object takes keyworded input as arguments by entending the python
        dictionary object. 

        The required fields are - 
            : timestamp:    datetime, pd.Datetimeindex
            : lon:          longitude, 0-359, float
            : lat:          latitude, -180, 180, float
            : mslp:         central pressure, Pa, float
            : vmax:         maximum velocity, m/s, float
        
        The optional (but recommended) fields are - 
            : rmax:         radius of maximum wind, m, float
            : vinfo:        list of velocity, m/s, list or numpy array, size n
            : radinfo:      2d list of radial info, m, list or numpy array, size nx4
            : ustorm:       speed of storm in x-direction, m/s
            : vstorm:       speed of strom in y-direction, m/s
        '''

        req_kw = ['timestamp', 'lon', 'lat', 'mslp', 'vmax']
        opt_kw = ['rmax', 'vinfo', 'radinfo', 'ustorm', 'vstorm']

        try:
            assert np.all([kw in kwargs for kw in req_kw])
        except:
            raise Exception(f'All required fields are not provided')

        for kw in opt_kw:
            if kw not in kwargs:
                kwargs[kw] = np.nan

        # Correcting datetime and convert to pandas datetime object
        allowed_datetime = (datetime, pd.DatetimeIndex)
        try:
            assert isinstance(kwargs['timestamp'], allowed_datetime)
        except:
            raise Exception(f'timestamp must be parsable by pd.to_datetime')
        else:
            kwargs['timestamp'] = pd.to_datetime(kwargs['timestamp'])

        # Initiate the dictionary
        super(Record, self).__init__(*args, **kwargs)
    
    def __setitem__(self, key, value):
        super(Record,self).__setitem__(key, value)

    def apply_translation_correction(self, fraction, angle):
        pass

    def to_boundary(self, swrf):
        pass

    def to_10m(self, swrf):
        pass

    def __str__(self):
        repr_str = f'\t'.join([str(self[key]) for key in self])
        return(repr_str)

class Track(object):
    def __init__(self, records):
        '''
        Take a record or an array of record and provide track related functionalities.
        '''
        try:
            assert isinstance(records, (Record, np.ndarray))
        except:
            raise Exception(f'Records must be a single or an array of records')

        self.records = np.atleast_1d(records)

        try:
            assert np.all([isinstance(i, Record) for i in records])
        except:
            raise Exception(f'The records must be an array of Record object')

        self.timeindex = pd.to_datetime([record['timestamp'] for record in self.records])

        self.sort() # To put all records in ascending order
        self.calc_translation() # Recalculate translation speed

    def __iter__(self):
        '''
        Called and return the timeindex for invoke like - for t in track
        '''
        return iter(self.timeindex)

    def __getitem__(self, key):
        '''
        Get an item at a given timeindex
        '''
        if isinstance(key, int):
            track = Track(self.records[key])
        elif isinstance(key, np.ndarray):
            try:
                np.all([isinstance(i, int) for i in key])
            except:
                raise Exception(f'Integer indexing array expected')
            
            track = Track(self.records[key])
        elif isinstance(key, slice):
            track = Track(self.records[key.start:key.stop:key.step])
        else:
            try:
                key = pd.to_datetime(key)
            except:
                raise Exception(f'Accessor must be parsable by pd.to_datetime')

            track = Track(self.records[self.timeindex == key])

        return track

    def __setitem__(self, key, value):
        '''
        Set an item at a given time index
        '''
        try:
            key = pd.to_datetime(key)
        except:
            raise Exception(f'Accessor must be parsable by pd.to_datetime')

        self.records[self.timeindex == key] = value

    def __contains__(self, key):
        '''
        Implements checking of a record
        '''
        try:
            key = pd.to_datetime(key)
        except:
            key = key
        
        return key in self.timeindex

    def __getattr__(self, name):
        '''
        Implements accessing the the dictionary objects as array.
        '''
        try:
            attr = np.array([record[name] for record in self.records])
        except:
            raise Exception(f'{name} not found in the records')
        
        return attr

    def sort(self):
        '''
        Apply time sorting and sort the records
        '''
        isort = np.argsort(self.timeindex)
        self.timeindex = self.timeindex[isort]
        self.records = self.records[isort]

    def calc_translation(self):
        '''
        Calculate and update the ustorm, vstorm fields.
        '''
        for i in np.arange(0, len(self.records)-1):
            dt = (self.timeindex[i+1] - self.timeindex[i]).total_seconds()
            origin = (self.records[i]['lon'], self.records[i]['lat'])
            of = (self.records[i+1]['lon'], self.records[i+1]['lat'])
            dtrans_x, dtrans_y = gc_distance(of=of, origin=origin, isradians=False)
            self.records[i]['ustorm'] = dtrans_x/dt
            self.records[i]['vstorm'] = dtrans_y/dt

            # For the last time step, we are keeping it to the same
            self.records[-1]['ustorm'] = self.records[-2]['ustorm']
            self.records[-1]['vstorm'] = self.records[-2]['vstorm']

    def append(self, records):
        try:
            assert isinstance(records, (Record, np.array))
        except:
            raise Exception(f'Records must be a single or an array of records')

        records = np.atleast_1d(records)

        try:
            assert np.all([isinstance(i, Record) for i in records])
        except:
            raise Exception(f'The records must be an array of Record object')

        self.records = np.append(self.records, records)
        self.timeindex = pd.to_datetime([record.timestamp for record in self.records])
        self.sort() # To put all record in ascending order
        self.calc_translation() # Recalculate translation speed

    def interpolate(self, at, **kwargs):
        try:
            at = pd.to_datetime(at, **kwargs)
        except:
            raise Exception(f'`at` must be parsable by pd.to_datetime')

        # Calculate the corresponding indices
        i_right = np.searchsorted(a=self.timeindex, v=at, side='left')

        if i_right == 0:
            i_left = 0
        else:
            i_left = i_right - 1

        # Calculate the corresponding weights
        v_left, v_right = self.timeindex[[i_left, i_right]]
        
        if(i_left < i_right):
            td_right = (v_right - at).total_seconds() # Timedelta to seconds
            td_total =  (v_right - v_left).total_seconds() # Timedelta to seconds
            w_left = td_right/float(td_total)
            w_right = 1-w_left
        else:
            w_left = float(1)
            w_right = float(0)

        # Create Record for the necessary fields with direct interpolation
        int_record = Record(
            timestamp=at,
            lon=self.records[i_left]['lon']*w_left + self.records[i_right]['lon']*w_right,
            lat=self.records[i_left]['lat']*w_left + self.records[i_right]['lat']*w_right,
            mslp=self.records[i_left]['mslp']*w_left + self.records[i_right]['mslp']*w_right,
            vmax=self.records[i_left]['vmax']*w_left + self.records[i_right]['vmax']*w_right,
            rmax=self.records[i_left]['rmax']*w_left + self.records[i_right]['rmax']*w_right,
            ustorm=self.records[i_left]['ustorm']*w_left + self.records[i_right]['ustorm']*w_right,
            vstorm=self.records[i_left]['vstorm']*w_left + self.records[i_right]['vstorm']*w_right
        )

        # now interpolate vinfo and radinfo
        # export np.nan if not available
        vinfo_left = self.records[i_left]['vinfo']
        lvinfo_left = len(vinfo_left)
        vinfo_right = self.records[i_right]['vinfo']
        lvinfo_right = len(vinfo_right)
        radinfo_left = np.atleast_2d(self.records[i_left]['radinfo'])
        sradinfo_left = radinfo_left.shape
        radinfo_right = np.atleast_2d(self.records[i_right]['radinfo'])
        sradinfo_right = radinfo_right.shape

        # Checking if the left and right are same size
        if lvinfo_left == lvinfo_right:
            vinfo_diff = False
        else:
            vinfo_diff = True

        # Interpolating
        if vinfo_diff:
            # vinfo size is different more check is needed
            min_vinfo_length = np.min([lvinfo_left, lvinfo_right])
            vinfo = vinfo_left[0:min_vinfo_length]

            # First find which one is the longer and set radinfo size
            if lvinfo_left > lvinfo_right:
                radinfo = np.zeros(shape=(min_vinfo_length, sradinfo_left[1]))
            elif lvinfo_left < lvinfo_right:
                radinfo = np.zeros(shape=(min_vinfo_length, sradinfo_right[1]))

            for i in np.arange(min_vinfo_length):
                    if np.any([np.isnan(vinfo_left[i]), np.isnan(vinfo_right[i])]):
                        radinfo = radinfo*np.nan
                    else:
                        vinfo[i] = vinfo_left[i]
                        for j in np.arange(radinfo.shape[1]):
                            radinfo[i, j] = radinfo_left[i, j]*w_left + radinfo_right[i, j]*w_right
        else:
            # vinfo size is same
            len_vinfo = lvinfo_left
            # If the both of them has length 1
            if len_vinfo == 1:
                # Return missing if vinfo itselft is missing
                if np.any([np.isnan(vinfo_left[0]), np.isnan(vinfo_right[0])]):
                    vinfo = np.nan
                    radinfo = np.nan
                # Or set to zero
                elif np.any([vinfo_left[0]==0, vinfo_right[0]==0]):
                    vinfo = np.nan
                    radinfo = np.nan
                # Otherwise calculate the linear interpolation
                else:
                    vinfo = vinfo_left[0]
                    radinfo = np.zeros_like(radinfo_left)
                    for j in np.arange(sradinfo_left[0]):
                        radinfo[j] = radinfo_left[j]*w_left + radinfo_right[j]*w_right
            elif len_vinfo > 1:
                # If both of them as length > 1
                vinfo = vinfo_left
                shape_radinfo = (lvinfo_left, sradinfo_left[1])
                radinfo = np.zeros(shape=shape_radinfo)
                for i, v in enumerate(vinfo):
                    if np.any([np.isnan(vinfo_left[i]), np.isnan(vinfo_right[i])]):
                        radinfo = radinfo*np.nan
                    else:
                        vinfo[i] = vinfo_left[i]
                        for j in np.arange(shape_radinfo[1]):
                            radinfo[i, j] = radinfo_left[i, j]*w_left + radinfo_right[i, j]*w_right
        
        # Providing the same as the dataset if interpolated on a given time
        if w_right == 1:
            vinfo = vinfo_right
            radinfo = radinfo_right
        elif w_left == 1:
            vinfo = vinfo_left
            radinfo = radinfo_left

        # Updating int_record with vinfo, radinfo
        int_record['vinfo'] = vinfo
        int_record['radinfo'] = radinfo

        # Return the record file for rest of the calculations
        return(int_record)

    def plot(self, ax=None, subplot_kw={'projection':ccrs.PlateCarree()}):
        if not isinstance(ax, plt.Axes):
            _, ax = plt.subplots(subplot_kw=subplot_kw)
        else:
            ax = ax
        
        ax.plot(self.lon, self.lat, '.-')
        
        return(ax)

    def __str__(self):
        '''
        String representation for print function.
        '''
        out_str = '\n'.join(str(record) for record in self.records)
        return(out_str)

if __name__=='__main__':
    print('Cyclone module of pyschism package.')