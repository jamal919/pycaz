#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cyclone tracks object, wind field models

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
from __future__ import print_function
from .conversion import gc_distance, pa2mb, km2m
import pandas as pd
import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
import sys

def coriolis(lat):
    '''
    Calculate the coriolis parameter
    '''
    f = 2*7.292e-5*np.sin(np.deg2rad(lat))
    return(f)

def calc_rmax_s02(mslp):
    '''
    Calculates rmax as a function of central pressure as described in Silva 2002.

    mslp: float, central pressure in Pascal 

    Ref: Silva,  R., G. Georges, S. Paulo, B. Gustavo and B. G. D. Gabriel (2002).
    Oceanographic vulnerability to hurricanes on the Mexican coast.
    Pro-ceedings of 28th Conference on Coastal Engineering, 39-51. 
    '''
    try:
        mslp = pa2mb(mslp)
    except:
        raise Exception(f'mslp is a required input for S02 method of rmax calculation')
    
    rmax = 0.4785 * mslp - 413 # In Km
    rmax = km2m(rmax)

    return rmax

def calc_rmax_w04(vmax, lat):
    '''
    Calculates rmax as a function of the maximum velocity as described in 
    Willoughby and Rahn (2004) paper.

    vmax: float, m/s
    lat: float, m/s

    Ref: 
    Willoughby and Rahn (2004) Parametric Representation of the Primary Hurricane
    Vortex. Part I: Observations andEvaluation of the Holland (1980) Model, Monthly
    Wearher Review
    '''
    try:
        rmax = 51.6*np.exp(-0.0223*vmax + 0.0281*lat) # Km
    except:
        raise Exception(f'vmax and lat are required input for W04 method of rmax calculation')

    rmax = km2m(rmax)

    return(rmax)

def calc_rmax_e11(v, r, vmax, f, solver='fsolve', limit=[5000, 500000], step=100):
    '''
    Calculate maximum radius using Emanuel 2011 model. 

    v: float, known velocity information at the boundary layer, m/s
    r: float, known radial information corresponding to v, m
    vmax: float, known maximum velocity, m/s
    f: coriolis coefficient, taken at the center as it is small
    solver: solver type -
        fsolve: use scipy.optimize.fsolve
        scan: use linear scanning
    limit: limit for scanning in scan solver and limit[0] is x0 for fsolve
    step: stepping for scan solver
    vlimit: (vmin, vmax) limit to return a result, otherwise exception thrown

    Ref: Emanuel 2011, Lin and Chavas 2012, Krien et al. 2017
    '''

    resfunc = lambda rmax: v - (2*r*(vmax*rmax+0.5*f*rmax**2)/(rmax**2+r**2)-f*r/2)
    try:
        if solver=='scan':
            rm_range = np.arange(start=limit[0], stop=limit[1]+step, step=step)
            for rm in rm_range:
                res = resfunc(rm)
                if res < 0:
                    break
            rmax_solved = rm
        elif solver=='fsolve':
            rmax_solved = optimize.fsolve(func=resfunc, x0=limit[0])
    except:
        raise Exception('Solver failed')
    
    return(rmax_solved)

def calc_holland_B(vmax, pc, pn, rhoair, bmax=2.5, bmin=0.5):
    '''
    Calculates the simplified version of Holland's B parameters excluding the
    term with coriolis.

    Arguments:
        vmax: Velocity of maximum wind at boundary layer, m/s
        pc: Central pressure, Pa
        rhoair: Density of air km/m**3
    '''
    B = vmax**2*rhoair*np.exp(1)/(pn-pc)

    if B<bmin:
        B = bmin
    elif B>bmax:
        B = bmax
    else:
        B = B
    
    return(B)

def calc_rmax_h80(v, r, pc, B, f, pn, rhoair, solver='fsolve', limit=[5000, 100000], step=100):
    '''
    Solve rmax with for a given r and v using Holland 1980 model.
    '''

    resfunc = lambda rmax: v - (np.sqrt(B/rhoair*(rmax/r)**B*(pn-pc)*np.exp(-(rmax/r)**B)+(r*f/2)**2) - r*f/2)

    try:
        if solver=='scan':
            rm_range = np.arange(start=limit[0], stop=limit[1]+step, step=step)
            for rm in rm_range:
                res = resfunc(rm)
                if res < 0:
                    break
            rmax_solved = rm
        elif solver=='fsolve':
            rmax_solved = optimize.fsolve(func=resfunc, x0=limit[0])
    except:
        raise Exception('Solver failed')
    
    return(rmax_solved)

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

    def gen_radial_fields(self, fraction=0.56, angle=19.2, swrf=0.9):
        '''
        Generate the interpolator for radial field for exisiting radial informations.
        To get correct radial field, atmospheric background is removed and
        the amount to be removed is controlled by fraction, and angle.  

        Arguments:
            fraction: float, fraction of translation speed to remove
            angle: float, angle of translation velocity change in degree
            swrf: float, surface wind reduction factor

        Due to translation a fraction of the wind-speed is embadded to the axysymmetric
        wind structure. To get the radial wind structure, we need to remove this
        translation velocity from the recorded wind velocity.
        
        The fraction is historically taken as unity (1). However, as shown in 
        lin and chavas (2012), more statistically consistent value of this fraction
        is 0.56 (default)

        Similarly, the translation is not fully forward (angle=0), rather the 
        background wind is slightly on southerly direction if we consider the
        storm is westerly. The angle is found around 19.2 degrees from statistical
        analysis. Thus the background wind is rotated anti-clockwise and needs 
        to be updated before applying as fraction to the wind values.

        In the advisories, radial fields are defined at 10m level, whereas 
        the analytical fields are defined at the boundary layer. Thus the velocity
        is needed to be divided by SWRF to shift it to boundary layer and multiplied
        by SWRF to get it back to 10m level.
        '''
        # Rotating ustorm vstorm by 19.2 degree anti-clockwise
        # Same operation as rotating a cartesian coordinate by 19.2 degree
        angle_rad = np.deg2rad(angle)
        ustorm = self['ustorm']*np.cos(angle_rad) - self['vstorm']*np.sin(angle_rad)
        vstorm = self['ustorm']*np.sin(angle_rad) + self['vstorm']*np.cos(angle_rad)

        # Apply correction to vmax for 4 quadrant clockwise from 0N
        # Careful about the sin/cos from north
        #            N
        #            | t /
        #            | /
        #    W ------*------- E
        #            |
        #            |
        #            S
        theta = np.linspace(0, 2*np.pi, 5)
        vmax_x = self['vmax'] * np.sin(theta)*(-1) - fraction*ustorm 
        vmax_y = self['vmax'] * np.cos(theta) - fraction*vstorm
        vmax_amp = np.sqrt(vmax_x**2 + vmax_y**2)/swrf
        self['fvmax'] = interp1d(theta, vmax_amp)

        # Apply correction to vinfo for 4 quadrant
        fvinfo = np.array([])
        fradinfo = np.array([])
        for i, vinfo in enumerate(np.atleast_1d(self['vinfo'])):
            if np.isnan(vinfo):
                self['vinfo'][i] = np.nan
                self['radinfo'][i] = np.nan
            else:
                vinfo_x = vinfo * np.sin(theta)*(-1) - fraction*ustorm
                vinfo_y = vinfo * np.cos(theta) - fraction*vstorm
                vinfo_amp = np.sqrt(vinfo_x**2 + vinfo_y**2)/swrf
                vinfo_func = interp1d(theta, vinfo_amp)
                fvinfo = np.append(fvinfo, vinfo_func)

                radinfo = np.append(self['radinfo'][i], self['radinfo'][i][0])
                radinfo_func = interp1d(theta, radinfo)
                fradinfo = np.append(fradinfo, radinfo_func)
        
        self['fvinfo'] = fvinfo
        self['fradinfo'] = fradinfo

    def calc_rmax(
        self, 
        methods=['E11', 'H80', 'S02', 'W04'], 
        discard=False, 
        vlimit_min=np.nan, 
        vlimit_max=np.nan,
        kw_atmos={'pn':101325, 'rhoair':1.15},
        kw_h80={'bmax':2.5, 'bmin':0.5, 'solver':'fsolve', 'limit':[5000, 500000], 'step':100},
        kw_e11={'solver':'fsolve', 'limit':[5000, 500000], 'step':100}
    ):
        '''
        Calculate Radius of maximum wind based on a given method on the selected
        radial informations.

        methods: list, list of methods to try one after another
        discard: bool, existing rmax to be used directly or not
        vlimit_min: list (4) or nan, list of minimum v for method to be used
        vlimit_max: list (4) or nan, list of maximum v for method to be used
        kw_atmos: dict, atmospheric information, generally needed for H80 model
        kw_h80: dict, solver options for H80 model
        kw_e11: dict, solver options for E11 model

        Following methods are implemented - 
            'E11' : Emanuel 2011 model as described in Lin and Chavas 2012
            'H80' : Holland 1980 model based on radial info
            'S02' : Silva et al. 2002 based on regression of central pressure
                    rmax = 0.4785 * Pc (mb) - 413 (km)
            'W04' : Willoughby and Rahn (2004) regression of Vmax
                    rmax = 51.6*exp(-0.0223*Vmax + 0.0281*lat) (km)

        Models to be implemented in the future - 
            'H80c' : Holland 1980c (f related values are removed)
            'S92' : Slosh model, by jeleniaski et al . 1992
            'E04' : Emanuel 2004 model
            'M16' : Murty et al 2016 model
            'W06' : Willoughby et al. 2006 model 
        
        '''
        methods = np.atleast_1d(methods)

        calc_rmax = {
            'E11': calc_rmax_e11,
            'H80': calc_rmax_h80,
            'S02': calc_rmax_s02,
            'W04': calc_rmax_w04
        }
        
        if np.isnan(vlimit_min):
            vlimit_min = np.ones(shape=methods)*np.nan

        if np.isnan(vlimit_max):
            vlimit_max = np.ones(shape=methods)*np.nan
        
        try:
            assert len(vlimit_min) == len(methods)
        except:
            vlimit_min = np.ones(shape=methods)*np.nan
            raise Warning(f'vmin must be the size of methods. Set to no limit')

        try:
            assert len(vlimit_max) == len(methods)
        except:
            vlimit_max = np.ones(shape=methods)*np.nan
            raise Warning(f'vmin must be the size of methods. Set to no limit')


        if not (np.isnan(self['rmax']) & ~discard):
            # Use the provided vmax
            theta = np.linspace(0, 2*np.pi, 5)
            rmax = np.ones_like(theta)*self['rmax']
            frmax = interp1d(theta, rmax)
            self['frmax'] = frmax
        elif np.isnan(self['vinfo']):
            theta = np.linspace(0, 2*np.pi, 5)
            if not np.isnan(self['mslp']):
                # Use S02 regression method for rmax
                rmax = np.ones_like(theta)*calc_rmax_s02(mslp=self['mslp'])
            elif not np.isnan(self['vmax']):
                # Use W04 regression method for rmax
                rmax = np.ones_like(theta)*calc_rmax_w04(vmax=self['vmax'], lat=self['lat'])
            
            frmax = interp1d(theta, rmax)
            self['frmax'] = frmax
        else:
            frmax = [] # length of fvinfo or atleast 1
            theta = np.linspace(0, 2*np.pi, 5)
            
            # Radial info is available and rmax is to be calculated from radinfo
            for fv, fr in zip(self['fvinfo'], self['fradinfo']):
                rmax_fv_theta = np.ones_like(theta) # Placeholder
                
                # Iterating over all the values interpolated by theta
                for i, thetai in enumerate(theta):
                    for method in methods:
                        try:
                            # Check the vlimit_min vlimit_max
                            assert fv(thetai) >= vlimit_min[i]
                            assert fv(thetai) <= vlimit_max[i]
                            
                            # Trying each methods
                            if method == 'E11':
                                kwargs = {
                                    'v':fv(thetai), 
                                    'r':fr(thetai), 
                                    'vmax':self['fvmax'](thetai), 
                                    'f':coriolis(self['lat']), 
                                    'solver':kw_e11['solver'], 
                                    'limit':kw_e11['limit'], 
                                    'step':kw_e11['step']
                                }
                                rmax_fv_theta[i] = calc_rmax[method](**kwargs)

                            if method == 'H80':
                                # First calculate holland B parameter
                                kwargs = {
                                    'vmax':self['fvmax'](thetai),
                                    'pc':self['mslp'],
                                    'rhoair':kw_atmos['rhoair'],
                                    'bmax':kw_h80['bmax'],
                                    'bmin':kw_h80['bmin']
                                }
                                B = calc_holland_B(**kwargs)

                                # Then calculate the rmax
                                kwargs = {
                                    'v':fv(thetai), 
                                    'r':fr(thetai), 
                                    'pc':self['mslp'], 
                                    'B':B, 
                                    'f':coriolis(self['lat']), 
                                    'pn':kw_atmos['pn'], 
                                    'rhoair':kw_atmos['rhoair'], 
                                    'solver':kw_h80['solver'], 
                                    'limit':kw_h80['limit'], 
                                    'step':kw_h80['step']
                                }
                                rmax_fv_theta[i] = calc_rmax[method](**kwargs)

                            if method == 'S02':
                                kwargs = {
                                    'mslp':self['mslp']
                                }
                                rmax_fv_theta[i] = calc_rmax[method](**kwargs)

                            if method == 'W04':
                                kwargs = {
                                    'vmax':self['fvmax'](thetai),
                                    'lat':self['lat']
                                }
                                rmax_fv_theta[i] = calc_rmax[method](**kwargs)

                            # When successful break the loop
                            break
                        except:
                            # Otherwise continue with the methods
                            continue
                
                # For a particular vinfo create the interpolation function
                # append to frmax
                frmax_fv = interp1d(theta, rmax_fv_theta)
                frmax = np.append(frmax, frmax_fv)
            
            # Save to Record dictionary
            self['frmax'] = frmax  

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
                for i, _ in enumerate(vinfo):
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