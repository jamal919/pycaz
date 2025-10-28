#!/usr/bin/env python
# coding: utf-8

from pycaz.core.grid import Grid
from pycaz.schism.sflux import Sflux

import pandas as pd
import numpy as np
import xarray as xr

# Geographical limit
x = np.arange(79, 99, 0.025)
y = np.arange(10, 24, 0.025)

starttime = '2020-05-16 00:00:00'
endtime = '2020-05-21 12:00:00'
timestep = '15min'

temp_K = 300
spfh = 0.0175873
pn = 100800

X, Y = np.meshgrid(x, y, indexing='xy')

model = xr.open_dataset('gfs.nc')
u10 = model['u10']
v10 = model['v10']
prmsl = model['prmsl']


grid = Grid(x=x, y=y)
sflux = Sflux(
    grid=grid, 
    basedate=pd.to_datetime(starttime), 
    sflux_type='air', 
    nstep=96, 
    path='./sflux'
)

timestamps = pd.date_range(start=starttime, end=endtime, freq=timestep)

for timestamp in timestamps:
    # Model fields
    try:
        # Model wind
        ui = u10.interp(lon=x, lat=y, time=timestamp)
        vi = v10.interp(lon=x, lat=y, time=timestamp)
        
        # Model pressure
        prmsli = prmsl.interp(lon=x, lat=y, time=timestamp)
        print(f'Timestamp {str(timestamp)} interpolated.')
    except:
        ui = X*0
        vi = X*0
        prmsli = X*pn
        print(f'Timestamp {str(timestamp)} not found. Default value.')

    # Writing to file
    sflux.write(at = timestamp, flux={
        'uwind': ui,
        'vwind': vi,
        'prmsl': prmsli,
        'stmp' : X*0+temp_K,
        'spfh' : X*0+spfh
    })
sflux.finish()
sflux.sfluxtxt(dt=pd.Timedelta(timestep))
