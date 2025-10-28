import sys

sys.path.append('/home/khan/Documents/pycaz')
import logging

import numpy as np
import pandas as pd
import xarray as xr

from pycaz.core.grid import Grid
from pycaz.schism.sflux import Sflux

# ERA5
# Load your own data here and change the variable name as necessary
model = xr.open_dataset('era5.nc')
model_u10 = model['u10']
model_v10 = model['v10']
model_prmsl = model['msl']

# Filler sflux variables
temp_K = 300
spfh = 0.0175873

# Geographical limit
# or Load it from era5 dataset
x = np.arange(30, 70.5, 0.25) # x = model.longitude
y = np.arange(-40, 0.5, 0.25) # y = model.latitude[::-1]

# Time
# Or you can load it from era5
starttime = '2023-02-10 00:00:00'
endtime = '2023-03-20 00:00:00'
timestep = '1h' # 15min

# Spatio-temporal
X, Y = np.meshgrid(x, y, indexing='xy')
timestamps = pd.date_range(start=starttime, end=endtime, freq=timestep) # timesteps = model.time


grid = Grid(x=x, y=y) # Create a grid
sflux = Sflux(
    grid=grid, 
    basedate=pd.to_datetime(starttime), 
    sflux_type='air', 
    nstep=96, # create new file after nstep
    path='./sflux'
)

for timestamp in timestamps:
    logging.info(f'Begun processing {timestamp}')

    # Interpolating the gridded data
    u = model_u10.interp(longitude=x, latitude=y, time=timestamp) # for era5 you can use sel
    v = model_v10.interp(longitude=x, latitude=y, time=timestamp)
    prmsl = model_prmsl.interp(longitude=x, latitude=y, time=timestamp)

    # Filler variables
    stmpm = X*0+temp_K
    spfhm = X*0+spfh
    
    # Writing to file
    sflux.write(at = timestamp, flux={
        'uwind': u,
        'vwind': v,
        'prmsl': prmsl,
        'stmp' : stmpm,
        'spfh' : spfhm
    })
    logging.info(f'Written {timestamp} to sflux')

sflux.finish()
sflux.sfluxtxt(dt=pd.Timedelta(timestep))
logging.info(f'Sflux generation done')
