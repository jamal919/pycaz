import sys
sys.path.append('/home/khan/MEGA/Codes/pyschism')
from schism.cyclone import Record, Track
from schism.core import Grid
from schism.flux import Sflux
from schism import read_jtwc
from schism.cyclone import coriolis, calc_holland_B
from schism.cyclone import calc_vcirc_e11, calc_vcirc_h80
from schism.cyclone import calc_mslp_h80
from schism.conversion import gc_distance

import cartopy
import pandas as pd
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import logging

# Logger
logging.basicConfig(filename='sflux.log', level=logging.INFO, filemode='w')


# Track
track = read_jtwc('bio012020.dat')
ax = track.plot() # Geoaxes object
ax.set_xlim([83, 97])
ax.set_xticks(np.arange(84, 97, 2))
ax.set_xticklabels(np.arange(84, 97, 2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.set_ylim([9, 27])
ax.set_yticks(np.arange(10, 27, 2))
ax.set_yticklabels(np.arange(10, 27, 2))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS)
ax.stock_img()
plt.savefig('track.png', dpi=150, bbox_inches='tight')

# GFS
model = xr.open_dataset('gfs.nc')
model_u10 = model['u10']
model_v10 = model['v10']
model_prmsl = model['prmsl']
model_spfh = model['spfh']
model_stmp = model['tmp2m']

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

starttime = '2020-05-16 00:00:00'
endtime = '2020-05-21 06:00:00'
timestep = '15min'

X, Y = np.meshgrid(x, y, indexing='xy')
timestamps = pd.date_range(start=starttime, end=endtime, freq=timestep)


grid = Grid(x=x, y=y)
sflux = Sflux(
    grid=grid, 
    basedate=pd.to_datetime(starttime), 
    sflux_type='air', 
    nstep=96, 
    path='./sflux'
)

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
        if len(frmax)==1:
            rmax = frmax[0](theta_N)
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
