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

class Track(object):
    def __init__(self, **kwargs):
        '''
        Track object takes keyworded input as arguments. All keyworded arguments
        should contain an array of inputs. 

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
        '''
        # Variables
        try:
            self.timestamp = kwargs['timestamp'] # Timestamp, datetime, pd.Datetimeindex
        except:
            raise Exception(f'timestamp is a required data')

        try:
            self.lon = kwargs['lon'] # Longitude, float
        except:
            raise Exception(f'lon is a required data')

        try:
            self.lat = kwargs['lat'] # Latitude, float
        except:
            raise Exception(f'lat is a required data')
        
        try:
            self.mslp = kwargs['mslp'] # Central pressure in Pa, float
        except:
            raise Exception(f'mslp is a required data')

        try:
            self.vmax = kwargs['vmax'] # Maximum velocity m/s, float
        except:
            raise Exception(f'vmax is a required data')
        
        try:
            self.rmax = kwargs['rmax'] # Radius of maximum velocity m, float
        except:
            self.rmax = np.ones_like(self.lon)*np.nan
            raise Warning(f'Radius of max velocity not provided. Set to np.nan')

        try:
            self.vinfo = kwargs['vinfo']  # List of velocity where radial info available
        except:
            self.vinfo = np.ones_like(self.lon)*np.nan
            raise Warning(f'Radial velocity info is not provided. Set to np.nan')

        try:
            self.radinfo = kwargs['radinfo']  # List of radial info, one row for each vinfo
        except:
            self.radinfo = np.ones_like(self.lon)*np.nan
            raise Warning(f'Radial distance info not provided. Set to np.nan')

        # Check if the datetime object is timezone aware or not
        try:
            assert np.all([d.tzinfo is None for d in self.timestamp])
        except:
            raise Exception(f'Date should be timezone naive')
        
        # Try to convert to datetime
        allowed_datetime = (datetime, pd.DatetimeIndex)
        try:
            assert np.all([isinstance(i, allowed_datetime) for i in self.timestamp])
        except:
            raise Exception(f'Date must be one of datetime, pandas datetime, numpy datetime64')
        else:
            self.timeindex = pd.to_datetime(self.timestamp)

        # Create storm translation fields
        self.utstorm = np.empty_like(self.lon) # Storm translation speed u
        self.vtstorm = np.empty_like(self.lon) # Storm translation speed v

        # Calculating __utrans and __vtrans for all timestep except the last, i.e, -1
        for i in np.arange(0, len(self.timeindex)-1):
            dt = (self.timeindex[i+1] - self.timeindex[i]).total_seconds()
            origin = (self.lon[i], self.lat[i])
            of = (self.lon[i+1], self.lat[i+1])
            dtrans_x, dtrans_y = gc_distance(of=of, origin=origin, isradians=False)
            self.utstorm[i] = dtrans_x/dt
            self.vtstorm[i] = dtrans_y/dt

        # For the last time step, we are keeping it to the same
            self.utstorm[-1] = self.utstorm[-2]
            self.vtstorm[-1] = self.vtstorm[-2]

    def interpolate(self, at, **kwargs):
        try:
            at = pd.to_datetime(at, **kwargs)
        except:
            raise Exception(f'the value of at must be formatted parsable by pd.to_datetime')

        # Calculate the corresponding indices
        i_right = np.searchsorted(a=self.timeindex, v=at, side='left')

        if i_right == 0:
            i_left = 0
        else:
            i_left = i_right - 1

        # Calculate the corresponding weights
        v_left, v_right = self.timeindex[[i_left, i_right]]
        
        if(i_left < i_right):
            td_right = (v_right - at).total_seconds()
            td_total =  (v_right - v_left).total_seconds()
            w_left = td_right/float(td_total)
            w_right = 1-w_left
        else:
            w_left = float(1)
            w_right = float(0)

        # Lon
        lon = self.lon[i_left]*w_left + self.lon[i_right]*w_right

        # Lat
        lat = self.lat[i_left]*w_left + self.lat[i_right]*w_right

        # Vmax
        vmax = self.vmax[i_left]*w_left + self.vmax[i_right]*w_right

        # Mslp
        mslp = self.mslp[i_left]*w_left + self.mslp[i_right]*w_right

        # Rmax
        rmax = self.rmax[i_left]*w_left + self.rmax[i_right]*w_right

        # utstorm
        utstorm = self.utstorm[i_left]*w_left + self.utstorm[i_right]*w_right

        # vtrans
        vtstorm = self.vtstorm[i_left]*w_left + self.vtstorm[i_right]*w_right

        # vinfo and radinfo
        # export np.nan if not available
        vinfo_left = self.vinfo[i_left]
        lvinfo_left = len(vinfo_left)
        vinfo_right = self.vinfo[i_right]
        lvinfo_right = len(vinfo_right)
        radinfo_left = np.atleast_2d(self.radinfo[i_left])
        sradinfo_left = radinfo_left.shape
        radinfo_right = np.atleast_2d(self.radinfo[i_right])
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

        return(
            dict(
                timestamp = at,
                lon = lon,
                lat = lat,
                vmax = vmax,
                mslp = mslp,
                rmax = rmax,
                vinfo = vinfo,
                radinfo = radinfo,
                utstorm = utstorm,
                vtstorm = vtstorm
            )
        )
    
    def plot(self, ax=None, subplot_kw={'projection':ccrs.PlateCarree()}):
        if ax is None:
            _, ax = plt.subplots(subplot_kw=subplot_kw)
        else:
            ax.plot(self.lon, self.lat, '.-')
        return(ax)

    def __str__(self):
        msg = f'''
        Starttime : {self.timestamp[0]},
        Endtime : {self.timestamp[-1]},
        Minimum velocity : {np.min(self.vmax):.2f} m/s,
        Maximum velocity : {np.max(self.vmax):.2f} m/s,
        Minimum central pressure : {np.min(self.mslp):.2f} Pa
        '''
        return(msg)

    def __repr__(self):
        df = pd.DataFrame({
            'Datetime':self.timestamp,
            'Lon':self.lon,
            'Lat':self.lat,
            'Vmax':self.vmax,
            'Rmax':self.rmax,
            'Mslp':self.mslp
        })
        df = df.set_index('Datetime')
        return(df.__repr__())

if __name__=='__main__':
    print('Cyclone module of pyschism package.')