#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset
import os


def default_filename_formatter(sflux: Sflux) -> str:
    """
    Default name formatter for the Sflux file which follows the naming conventions from schism model.

    :param sflux: Sflux object
    :return: filename string
    """
    return f'sflux_{sflux.sflux_type}_{sflux.priority:1d}.{sflux.nfile:04d}.nc'


class Sflux(object):
    def __init__(self,
                 grid,
                 basedate,
                 sflux_type='air',
                 nstep=96,
                 priority=1,
                 syncstep=10,
                 nfile=0,
                 filename_formatter=default_filename_formatter,
                 path='./sflux'):
        """
        Generate SCHISM complient Sflux files

        :param grid:
        :param basedate:
        :param sflux_type:
        :param nstep:
        :param priority:
        :param syncstep:
        :param path:

        TODO:
            - Validate inputs
            - Update to a context interface
            - Separate into Sub-classes
            - Automate dt of the dataset
            - Allow numpy datetime64 objects
            - Separate netcdf file creation class, with data access, close function
        """
        self.grid = grid
        self.nstep = nstep  # No of step
        self.basedate = basedate
        self.sflux_type = sflux_type
        self.nfile = nfile  # No file at the beginning
        self.priority = priority  # sflux_air_1 or 2
        self.syncstep = syncstep  # Sync the netCDF each syncstep
        self.filename_formatter = filename_formatter
        self.path = path

        self.nc = None
        self.filepath = None
        self.filename = None
        self.step = None

        sflux_func_map = {
            'air': {
                'create_netcdf': self.create_netcdf_air,
                'put_value': self.put_value_air
            },
            'prc': {
                'create_netcdf': self.create_netcdf_prc,
                'put_value': self.put_value_prc
            },
            'rad': {
                'create_netcdf': self.create_netcdf_rad,
                'put_value': self.put_value_rad
            }
        }

        if sflux_type in sflux_func_map:
            self.create_netcdf = sflux_func_map[sflux_type]['create_netcdf']
            self.put_value = sflux_func_map[sflux_type]['put_value']
        else:
            raise Exception(f"sflux_type {self.sflux_type} not correct, one of 'air', 'prc', and 'rad'")

        # Directory creation
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def create_netcdf_air(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.filename = self.filename_formatter(self)
        self.filepath = os.path.join(self.path, self.filename)

        # Creating the file first
        self.nc = Dataset(self.filepath, 'w', format='NETCDF4_CLASSIC')

        # Creating the dimensions
        self.nc.createDimension(dimname='nx_grid', size=len(self.grid.x))
        self.nc.createDimension(dimname='ny_grid', size=len(self.grid.y))
        self.nc.createDimension(dimname='time', size=None) # unlimited

        # Creating the variables
        # Time
        self.v_time = self.nc.createVariable(
            varname='time',
            datatype=float,
            dimensions=('time')
        )
        strf_basedate = self.basedate.strftime('%Y-%m-%d %H:%M:%S')
        self.v_time.units = f'days since {strf_basedate:s}'
        self.v_time.long_name = 'Time'
        self.v_time.calendar = 'standard'
        self.v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.v_lon = self.nc.createVariable(
            varname='lon',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lon.units = 'degrees_north'
        self.v_lon.long_name = 'Longitude'
        self.v_lon.standard_name = 'longitude'

        # Latitude
        self.v_lat = self.nc.createVariable(
            varname='lat',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lat.units = 'degrees_east'
        self.v_lat.long_name = 'Latitude'
        self.v_lat.standard_name = 'latitude'

        # Uwind
        self.v_uwind = self.nc.createVariable(
            varname='uwind',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_uwind.units = 'm/s'
        self.v_uwind.long_name = 'Surface Eastward Air Velocity (10m AGL)'
        self.v_uwind.standard_name = 'eastward_wind'

        # Vwind
        self.v_vwind = self.nc.createVariable(
            varname='vwind',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_vwind.units = 'm/s'
        self.v_vwind.long_name = 'Surface Northward Air Velocity (10m AGL)'
        self.v_vwind.standard_name = 'northward_wind'

        # Prmsl
        self.v_prmsl = self.nc.createVariable(
            varname='prmsl',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_prmsl.units = 'Pa'
        self.v_prmsl.long_name = 'Pressure Reduced to MSL'
        self.v_prmsl.standard_name = 'air_pressure_at_mean_sea_level'

        # stmp
        self.v_stmp = self.nc.createVariable(
            varname='stmp',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_stmp.units = 'K'
        self.v_stmp.long_name = 'Surface Temperature (2m AGL)'
        self.v_stmp.standard_name = 'surface_temperature'

        # spfh
        self.v_spfh = self.nc.createVariable(
            varname='spfh',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_spfh.units = '1'
        self.v_spfh.long_name = 'Specific Humidity (2m AGL)'
        self.v_spfh.standard_name = 'surface_specific_humidity'

        # Writing lon-lat once
        X, Y = self.grid.meshgrid
        self.v_lon[:] = X.T
        self.v_lat[:] = Y.T

    def put_value_air(self, stepi, at, flux):
        if isinstance(at, (datetime, pd.DatetimeIndex)):
            at = pd.to_datetime(at) - self.basedate
        elif isinstance(at, (timedelta, pd.Timedelta)):
            at = at
        else:
            raise Exception(f'at must be datetime or timedelta object')

        self.v_time[stepi] = at.days + at.seconds / float(86400)
        self.v_uwind[stepi, :, :] = flux['uwind']
        self.v_vwind[stepi, :, :] = flux['vwind']
        self.v_prmsl[stepi, :, :] = flux['prmsl']
        self.v_stmp[stepi, :, :] = flux['stmp']
        self.v_spfh[stepi, :, :] = flux['spfh']

        self.step = self.step + 1

        self.sync()

    def create_netcdf_prc(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.filename = self.filename_formatter(self)
        self.filepath = os.path.join(self.path, self.filename)

        # Creating the file first
        self.nc = Dataset(self.filepath, 'w', format='NETCDF4_CLASSIC')

        # Creating the dimensions
        self.nc.createDimension(dimname='nx_grid', size=len(self.grid.x))
        self.nc.createDimension(dimname='ny_grid', size=len(self.grid.y))
        self.nc.createDimension(dimname='time', size=None) # unlimited

        # Creating the variables
        # Time
        self.v_time = self.nc.createVariable(
            varname='time',
            datatype=float,
            dimensions=('time')
        )
        strf_basedate = self.basedate.strftime('%Y-%m-%d %H:%M:%S')
        self.v_time.units = f'days since {strf_basedate:s}'
        self.v_time.long_name = 'Time'
        self.v_time.calendar = 'standard'
        self.v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.v_lon = self.nc.createVariable(
            varname='lon',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lon.units = 'degrees_north'
        self.v_lon.long_name = 'Longitude'
        self.v_lon.standard_name = 'longitude'

        # Latitude
        self.v_lat = self.nc.createVariable(
            varname='lat',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lat.units = 'degrees_east'
        self.v_lat.long_name = 'Latitude'
        self.v_lat.standard_name = 'latitude'

        # Prate
        self.v_prate = self.nc.createVariable(
            varname='prate',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_prate.units = 'kg/m^2/s'
        self.v_prate.long_name = 'Surface Precipitation Rate'
        self.v_prate.standard_name = 'precipitation_flux'

        # Writing lon-lat once
        X, Y = self.grid.meshgrid
        self.v_lon[:] = X.T
        self.v_lat[:] = Y.T

    def put_value_prc(self, stepi, at, flux):
        if isinstance(at, (datetime, pd.DatetimeIndex)):
            at = pd.to_datetime(at) - self.basedate
        elif isinstance(at, (timedelta, pd.Timedelta)):
            at = at
        else:
            raise Exception(f'at must be datetime or timedelta object')

        self.v_time[stepi] = at.days + at.seconds / float(86400)
        self.v_prate[stepi, :, :] = flux['prate']

        self.step = self.step + 1

        self.sync()

    def create_netcdf_rad(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.filename = self.filename_formatter(self)
        self.filepath = os.path.join(self.path, self.filename)

        # Creating the file first
        self.nc = Dataset(self.filepath, 'w', format='NETCDF4_CLASSIC')

        # Creating the dimensions
        self.nc.createDimension(dimname='nx_grid', size=len(self.grid.x))
        self.nc.createDimension(dimname='ny_grid', size=len(self.grid.y))
        self.nc.createDimension(dimname='time', size=None) # unlimited

        # Creating the variables
        # Time
        self.v_time = self.nc.createVariable(
            varname='time',
            datatype=float,
            dimensions=('time')
        )
        strf_basedate = self.basedate.strftime('%Y-%m-%d %H:%M:%S')
        self.v_time.units = f'days since {strf_basedate:s}'
        self.v_time.long_name = 'Time'
        self.v_time.calendar = 'standard'
        self.v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.v_lon = self.nc.createVariable(
            varname='lon',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lon.units = 'degrees_north'
        self.v_lon.long_name = 'Longitude'
        self.v_lon.standard_name = 'longitude'

        # Latitude
        self.v_lat = self.nc.createVariable(
            varname='lat',
            datatype=float,
            dimensions=('ny_grid', 'nx_grid')
        )
        self.v_lat.units = 'degrees_east'
        self.v_lat.long_name = 'Latitude'
        self.v_lat.standard_name = 'latitude'

        # Dlwrf
        self.v_dlwrf = self.nc.createVariable(
            varname='dlwrf',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_dlwrf.units = 'W/m^2'
        self.v_dlwrf.long_name = 'Downward Long Wave Radiation Flux'
        self.v_dlwrf.standard_name = 'surface_downwelling_longwave_flux_in_air'

        # Dswrf
        self.v_dswrf = self.nc.createVariable(
            varname='dswrf',
            datatype=float,
            dimensions=('time', 'ny_grid', 'nx_grid')
        )
        self.v_dswrf.units = 'W/m^2'
        self.v_dswrf.long_name = 'Downward Short Wave Radiation Flux'
        self.v_dswrf.standard_name = 'surface_downwelling_shortwave_flux_in_air'

        # Writing lon-lat once
        X, Y = self.grid.meshgrid
        self.v_lon[:] = X.T
        self.v_lat[:] = Y.T

    def put_value_rad(self, stepi, at, flux):
        if isinstance(at, (datetime, pd.DatetimeIndex)):
            at = pd.to_datetime(at) - self.basedate
        elif isinstance(at, (timedelta, pd.Timedelta)):
            at = at
        else:
            raise Exception(f'at must be datetime or timedelta object')

        self.v_time[stepi] = at.days + at.seconds / float(86400)
        self.v_dlwrf[stepi, :, :] = flux['dlwrf']
        self.v_dswrf[stepi, :, :] = flux['dswrf']

        self.step = self.step + 1

        self.sync()

    def close_netcdf(self):
        self.nc.close()

    def write(self, at, flux):
        # First check if self.nc is available
        if hasattr(self, 'nc'):
            if self.step < self.nstep:
                self.put_value(self.step, at, flux)
            else:
                self.close_netcdf()
                self.create_netcdf()
                self.put_value(self.step, at, flux)
        else:
            self.create_netcdf()
            self.put_value(self.step, at, flux)

    def sync(self):
        if self.step % self.syncstep:
            self.nc.sync()

    def finish(self):
        if hasattr(self, 'nc'):
            self.close_netcdf()

    def sfluxtxt(self, dt):
        dt = dt.total_seconds()
        max_window = self.nstep * dt / float(3600)
        filepath = os.path.join(self.path, 'sflux_inputs.txt')
        with open(filepath, mode='w') as f:
            f.write('&sflux_inputs\n')
            f.write(
                'air_1_relative_weight=1.,\t!air_[12]_relative_weight set the relative ratio between datasets 1 and 2\n')
            f.write('air_2_relative_weight=99.,\n')
            f.write(
                f'air_1_max_window_hours={max_window:.1f},\t!max. # of hours (offset from start time in each file) in each file of set 1\n')
            f.write('air_1_fail_if_missing=.true.,\t!set 1 is mandatory\n')
            f.write('air_2_fail_if_missing=.false.,\t!set 2 is optional\n')
            f.write("air_1_file='sflux_air_1',\t!file name for 1st set of 'air'\n")
            f.write("air_2_file='sflux_air_2'\n")
            f.write("uwind_name='uwind',\t!name of u-wind vel.\n")
            f.write("vwind_name='vwind',\t!name of v-wind vel.\n")
            f.write("prmsl_name='prmsl',\t!name of air pressure (@MSL) variable in .nc file\n")
            f.write("stmp_name='stmp',\t!name of surface air T\n")
            f.write("spfh_name='spfh',\t!name of specific humidity\n")
            f.write("dlwrf_name='dlwrf',\t!name of downward longwave radiation variable\n")
            f.write("dswrf_name='dswrf',\t!name of downward shortwave radiation variable (solar)\n")
            f.write("prate_name='prate',\t!name of precipitation rate variable\n")
            f.write('/\n')



