#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interface to generate model input fluxes.

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""
import os
import numpy as np

from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date, date2num

class Sflux(object):
    def __init__(self, grid, basedate, sflux_type='air', nstep=96, priority=1, syncstep=10, path='./sflux'):
        self.grid = grid
        self.nstep = nstep # No of step
        self.basedate = basedate
        self.sflux_type = sflux_type
        self.nfile = 0 # No file at the beginning
        self.priority = priority # sflux_air_1 or 2
        self.syncstep = syncstep # Sync the netCDF each syncstep
        self.path = path

        sflux_func_map = {
            'air' : {
                'create_netcdf' : self.__create_netcdf_air,
                'put_value' : self.__put_value_air
            },
            'prc' : {
                'create_netcdf' : self.__create_netcdf_prc,
                'put_value' : self.__put_value_prc
            },
            'rad' : {
                'create_netcdf' : self.__create_netcdf_rad,
                'put_value' : self.__put_value_rad
            }
        }

        if sflux_type in sflux_func_map:
            self.__create_netcdf = sflux_func_map[sflux_type]['create_netcdf']
            self.__put_value = sflux_func_map[sflux_type]['put_value']
        else:
            raise Exception(f"sflux_type {self.sflux_type} not correct, one of 'air', 'prc', and 'rad'")
        
        # Directory creation
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __create_netcdf_air(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.__filename = f'sflux_air_{self.priority:1d}.{self.nfile:03d}.nc'
        self.__filepath = os.path.join(self.path, self.__filename)

        # Creating the file first
        self.nc = Dataset(self.__filepath, 'w', format='NETCDF4_CLASSIC')

        
        # Creating the dimensions
        self.__d_nx_grid = self.nc.createDimension(dimname='nx_grid', size=len(self.grid.x))
        self.__d_ny_grid = self.nc.createDimension(dimname='ny_grid', size=len(self.grid.y))
        self.__d_time = self.nc.createDimension(dimname='ntime', size=None)

        # Creating the variables
        # Time
        self.__v_time = self.nc.createVariable(
            varname='time',
            datatype=np.float32,
            dimensions=('ntime')
        )
        self.__v_time.units = f'days since {self.basedate.strftime('%Y-%m-%d %H:%M:%S'):s}'
        self.__v_time.long_name = 'Time'
        self.__v_time.calendar = 'standard'
        self.__v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.__v_lon = self.nc.createVariable(
            varname='lon',
            datatype=np.float32,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.__v_lon.units = 'degrees_north'
        self.__v_lon.long_name = 'Longitude'
        self.__v_lon.standard_name = 'longitude'

        # Latitude
        self.__v_lat = self.nc.createVariable(
            varname='lat',
            datatype=np.float32,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.__v_lat.units = 'degrees_east'
        self.__v_lat.long_name = 'Latitude'
        self.__v_lat.standard_name = 'latitude'

        # Uwind
        self.__v_uwind = self.nc.createVariable(
            varname='uwind',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.__v_uwind.units = 'm/s'
        self.__v_uwind.long_name = 'Surface Eastward Air Velocity (10m AGL)'
        self.__v_uwind.standard_name = 'eastward_wind'

        # Vwind
        self.__v_vwind = self.nc.createVariable(
            varname='vwind',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.__v_vwind.units = 'm/s'
        self.__v_vwind.long_name = 'Surface Northward Air Velocity (10m AGL)'
        self.__v_vwind.standard_name = 'northward_wind'

        # Prmsl
        self.__v_prmsl = self.nc.createVariable(
            varname='prmsl',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.__v_prmsl.units = 'Pa'
        self.__v_prmsl.long_name = 'Pressure Reduced to MSL'
        self.__v_prmsl.standard_name = 'air_pressure_at_mean_sea_level'

        # stmp
        self.__v_stmp = self.nc.createVariable(
            varname='stmp',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.__v_stmp.units = 'K'
        self.__v_stmp.long_name = 'Surface Temperature (2m AGL)'
        self.__v_stmp.standard_name = 'surface_temperature'

        # spfh
        self.__v_spfh = self.nc.createVariable(
            varname='spfh',
            datatype=np.float32,
            dimensions=('ntime', 'ny_grid', 'nx_grid')
        )
        self.__v_spfh.units = 1
        self.__v_spfh.long_name = 'Specific Humidity (2m AGL)'
        self.__v_spfh.standard_name = 'surface_specific_humidity'
        
        # Writing lon-lat once
        self.__v_lon[:] = self.grid.X
        self.__v_lat[:] = self.grid.Y
        
    
    def __put_value_air(self, stepi, at, flux):
        self.__v_time[stepi] = at.days + at.seconds/float(86400)
        self.__v_uwind[stepi, :, :] = flux['uwind']
        self.__v_vwind[stepi, :, :] = flux['vwind']
        self.__v_prmsl[stepi, :, :] = flux['prmsl']
        self.__v_stmp[stepi, :, :] = flux['stmp']
        self.__v_spfh[stepi, :, :] = flux['spfh']

        self.step = self.step + 1

        # Syncing each 10 step
        if self.step%self.syncstep:
            self.nc.sync()

    def __close_netcdf(self):
        self.nc.close()

    def write(self, at, flux):
        # First check if self.nc is available
        if hasattr(self, 'nc'):
            if self.step < self.nstep:
                self.__put_value(self.step, at, flux)
            else:
                self.__close_netcdf()
                self.__create_netcdf()
                self.__put_value(self.step, at, flux)
        else:
            self.__create_netcdf()
            self.__put_value(self.step, at, flux)

    def finish(self):
        if hasattr(self, 'nc'):
            self.__close_netcdf()

    def sfluxtxt(self, dt):
        __dt = dt.total_seconds()
        __max_window = self.nstep*__dt/float(3600)
        __filepath = os.path.join(self.path, 'sflux_inputs.txt')
        with open(__filepath, mode='w') as f:
            f.write('&sflux_inputs\n')
            f.write('air_1_relative_weight=1.,	!air_[12]_relative_weight set the relative ratio between datasets 1 and 2\n')
            f.write('air_2_relative_weight=99., \n')
            f.write(f'air_1_max_window_hours={__max_window:.1f},	!max. # of hours (offset from start time in each file) in each file of set 1\n')
            f.write('air_1_fail_if_missing=.true.,	!set 1 is mandatory\n')
            f.write('air_2_fail_if_missing=.false., 	!set 2 is optional\n')
            f.write("air_1_file='sflux_air_1', 	!file name for 1st set of 'air'\n")
            f.write("air_2_file='sflux_air_2'\n")
            f.write("uwind_name='uwind', 		!name of u-wind vel.\n")
            f.write("vwind_name='vwind', 		!name of v-wind vel.\n")
            f.write("prmsl_name='prmsl', 		!name of air pressure (@MSL) variable in .nc file\n")
            f.write("stmp_name='stmp',  		!name of surface air T\n")
            f.write("spfh_name='spfh',  		!name of specific humidity\n")
            f.write('/\n')