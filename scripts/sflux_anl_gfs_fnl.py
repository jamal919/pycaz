#!/usr/bin/env python
# coding: utf-8
import cartopy
import pandas as pd
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import sys
sys.path.append('/home/jkhan02/pyschism')
from schism.cyclone import Record, Track
from schism.core import Grid
from schism.flux import Sflux
from schism import read_jtwc
from schism.cyclone import coriolis, calc_holland_B
from schism.cyclone import calc_vcirc_e11, calc_vcirc_h80
from schism.cyclone import calc_mslp_h80
from schism.conversion import gc_distance

import os
import logging

# Logger
logging.basicConfig(filename=f'sflux_fnl_gfs_best.log', level=logging.INFO, filemode='w')


# Track
track_loc = './tracks'
track = read_jtwc(os.path.join(track_loc, f'track_best.txt'))

# GFS
model_loc = './gfs'
model = xr.open_dataset(os.path.join(model_loc, f'gfs_best.nc'))
model_u10 = model['u10']
model_v10 = model['v10']
model_prmsl = model['prmsl']
model_spfh = model['spfh']
model_stmp = model['tmp2m']

# FNL dataset
fnl_loc = './fnl'
fnl = xr.open_dataset(os.path.join(fnl_loc, 'fnl.nc'))
fnl_u10 = fnl['u10']
fnl_v10 = fnl['v10']
fnl_prmsl = fnl['prmsl']

# Parameters
fraction = 0.56
angle = 19.2
swrf = 0.9
one2ten = 0.88
pn = 100800
rhoair = 1.15
temp_K = 300
spfh = 0.0175873

# Geographical limit
x = np.arange(79, 99, 0.025)
y = np.arange(10, 24, 0.025)

X, Y = np.meshgrid(x, y, indexing='xy')
grid = Grid(x=x, y=y)
basedate = '2020-05-13 00:00:00'
sflux = Sflux(
    grid=grid,
    basedate=pd.to_datetime(basedate),
    sflux_type='air',
    nstep=96,
    path=f'./sflux_fnl_gfs_best'
)

# First FNL only
starttime = '2020-05-13 00:00:00'
endtime = '2020-05-16 05:45:00'
timestep = '15min'

timestamps = pd.date_range(start=starttime, end=endtime, freq=timestep)

for timestamp in timestamps:
    # Merging with model
    u = fnl_u10.interp(lon=x, lat=y, time=timestamp)
    v = fnl_v10.interp(lon=x, lat=y, time=timestamp)
    prmsl = fnl_prmsl.interp(lon=x, lat=y, time=timestamp)
    stmp = u*0+273+25 # some value, irrelevant for us
    spfh = u*0+0.1 # some value, irrelevant for us

    
    # Writing to file
    sflux.write(at = timestamp, flux={
        'uwind': u,
        'vwind': v,
        'prmsl': prmsl,
        'stmp' : stmp,
        'spfh' : spfh
    })

# Merged best and gfs
starttime = '2020-05-16 06:00:00'
endtime = '2020-05-21 06:00:00'
timestep = '15min'

timestamps = pd.date_range(start=starttime, end=endtime, freq=timestep)


for timestamp in timestamps:
    logging.info(f'Processing - {timestamp}')
    # Advisory fields
    try:
        record = track.interpolate(at=timestamp)
        record = record.gen_radial_fields(
            fraction=fraction, 
            angle=angle, 
            swrf=swrf
        )
        record = record.calc_rmax(
            methods=['E11'], 
            use_rmax_info=False, 
            vlimit=[np.inf],
            kw_atmos={'pn':pn, 'rhoair':rhoair}
        )
        # calculate cartesian position
        dx, dy = gc_distance(X, Y, record['lon'], record['lat'])
        
        # calculate the radial distance
        r = np.sqrt(dx**2 + dy**2)
        r = r*(r>0) + (r+1e-8)*(r<=0)

        # calculate the theta
        theta = np.arctan2(dy, dx)
        theta_N = -1*theta + np.pi/2
        theta_N = (theta_N+2*np.pi)*(theta_N<0)+(theta_N)*(theta_N>=0)

        # vmax
        vmax = record['fvmax'](theta_N)

        # Other parameters
        f = coriolis(Y)
        B = (vmax**2)*rhoair*np.exp(1)/(pn-record['mslp'])
        
        # marging circular winds
        frmax = np.atleast_1d(record['frmax'])
        fradinfo = np.atleast_1d(record['fradinfo'])

        # Defining rmax based on first radial info
        rmax = frmax[0](theta_N)

        # Circular wind for number of rmaxs
        if len(frmax)==1:
            # with already calculated rmax
            vcirc = calc_vcirc_h80(r, rmax, record['mslp'], B, pn, rhoair, f)
        elif len(frmax)==2:
            # outer to inner circle
            rad_0 = fradinfo[0](theta_N)
            rmax_0 = frmax[0](theta_N)
            vcirc_0 = calc_vcirc_h80(r, rmax_0, record['mslp'], B, pn, rhoair, f)
            
            rad_1 = fradinfo[1](theta_N)
            rmax_1 = frmax[1](theta_N)
            vcirc_1 = calc_vcirc_e11(r, rmax_1, vmax, f)
            
            vcirc = vcirc_0*(r>=rad_0) + ((vcirc_0-vcirc_1)/(rad_0-rad_1)*(r-rad_1)+vcirc_1)*(np.logical_and(r<rad_0, r>rad_1)) + vcirc_1*(r<=rad_1)
        else:
            # outer to inner circle
            rad_0 = fradinfo[0](theta_N)
            rmax_0 = frmax[0](theta_N)
            vcirc_0 = calc_vcirc_h80(r, rmax_0, record['mslp'], B, pn, rhoair, f)
            
            rad_1 = fradinfo[1](theta_N)
            rmax_1 = frmax[1](theta_N)
            vcirc_1 = calc_vcirc_e11(r, rmax_1, vmax, f)

            rad_2 = fradinfo[2](theta_N)
            rmax_2 = frmax[2](theta_N)
            vcirc_2 = calc_vcirc_e11(r, rmax_2, vmax, f)

            vcirc = vcirc_0*(r>=rad_0) + ((vcirc_0-vcirc_1)/(rad_0-rad_1)*(r-rad_1)+vcirc_1)*(np.logical_and(r<rad_0, r>rad_1)) + ((vcirc_1-vcirc_2)/(rad_1-rad_2)*(r-rad_2)+vcirc_2)*(np.logical_and(r<rad_1, r>rad_2)) + vcirc_2*(r<=rad_2)

        # correcting circular wind
        vcirc = vcirc * swrf * one2ten # SWRF, 10min
        vcirc = vcirc*(vcirc>=0) + 0*(vcirc<0)
        
        utrans = record['ustorm']*np.cos(np.deg2rad(angle)) - record['vstorm']*np.sin(np.deg2rad(angle))
        vtrans = record['ustorm']*np.sin(np.deg2rad(angle)) + record['vstorm']*np.cos(np.deg2rad(angle))
        
        ua = -vcirc*r*np.sin(theta)/r + fraction*utrans
        va = vcirc*r*np.cos(theta)/r + fraction*vtrans

        # pressure
        prmsla = calc_mslp_h80(r, rmax, record['mslp'], B, pn)
        logging.info(f'\tTrack: Record found.')
    except:
        ua = model_u10.interp(lon=x, lat=y, time=timestamp)
        va = model_v10.interp(lon=x, lat=y, time=timestamp)
        prmsla = model_prmsl.interp(lon=x, lat=y, time=timestamp)
        logging.info(f'\tTrack: Not possible to calculate. Default.')

    # Merging with model
    um = model_u10.interp(lon=x, lat=y, time=timestamp)
    vm = model_v10.interp(lon=x, lat=y, time=timestamp)
    prmslm = model_prmsl.interp(lon=x, lat=y, time=timestamp)
    stmpm = model_stmp.interp(lon=x, lat=y, time=timestamp)
    spfhm = model_spfh.interp(lon=x, lat=y, time=timestamp)

    lr_min = 3*rmax
    lr_max = 10*rmax
    weight = (r-lr_min)/(lr_max - lr_min)
    u = ua*(r<lr_min) + (ua*(1-weight) + um*weight)*np.logical_and(r>=lr_min, r<=lr_max) + um*(r>lr_max)
    v = va*(r<lr_min) + (va*(1-weight) + vm*weight)*np.logical_and(r>=lr_min, r<=lr_max) + vm*(r>lr_max)
    prmsl = prmsla*(r<lr_min) + (prmsla*(1-weight) + prmslm*weight)*np.logical_and(r>=lr_min, r<=lr_max) + prmslm*(r>lr_max)
    
    # Writing to file
    sflux.write(at = timestamp, flux={
        'uwind': u,
        'vwind': v,
        'prmsl': prmsl,
        'stmp' : stmpm,
        'spfh' : spfhm
    })
sflux.finish()
sflux.sfluxtxt(dt=pd.Timedelta(timestep))
