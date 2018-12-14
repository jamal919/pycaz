#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create SCHISM complient wind and pressure field using climate / atmospheric model
output.

Issues:
    The netCDF4 module in python 2.7 has a memory leak and thus will create a
    memory issue. Better to use python3

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
import sys
import os
import numpy as np
from scipy import interpolate
from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date, date2num, date2index
import glob

class Converter(object):
    def __init__(self):
        '''
        Collection of converter function to be used in the module. They are
        generally static methods and meant to be used directly with input variables.
        '''
        pass

    @staticmethod
    def hpa2pa(hpa):
        '''
        Takes pressure value in hecta Pascal and return in Pascal. 
        '''
        return(hpa*100)
    
    @staticmethod
    def knot2mps(knot):
        '''
        Takes velocity in knot and returns velocity in mps.
        '''
        return(knot*1.852/3.6)

    @staticmethod
    def km2m(km):
        '''
        Takes distance in Km and converts it to meter.
        '''
        return(km*1000)

    @staticmethod
    def lon180(lon360):
        '''
        Change lon range from 0-360 to -180-180
        '''
        lon360[lon360 > 180] = lon360[lon360 > 180] - 360
        return(lon360)

    @staticmethod
    def gc_distance(of, origin, isradians=False):
        '''
        Calculates the great circle distance of 'of' from 'origin'
        '''
        __dfac = 60*1.852*1000
        
        if isradians:
            __dtrans_x = __dfac*np.cos(origin[1])*(np.rad2deg(of[0])-np.rad2deg(origin[0]))
            __dtrans_y = __dfac*(np.rad2deg(of[1])-np.rad2deg(origin[1]))
        else:
            __dtrans_x = __dfac*np.cos(np.deg2rad(origin[1]))*(of[0]-origin[0])
            __dtrans_y = __dfac*(of[1]-origin[1])

        return((__dtrans_x, __dtrans_y))

class Grid(object):
    def __init__(self, x, y):
        '''
        Grid object to generate grid and provides function to find various
        values at grid points. 
        '''
        self.x = x
        self.y = y

        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = self.X.shape
        self.length = len(self.X.flatten())
        self.converter = Converter()

    def radial_distance(self, originx, originy, isradians=False):
        '''
        Calculates distance from a given point (origionx, originy) and returns
        the (radial_distance, x_distance, y_distance)
        '''
        __dfac = 60*1.852*1000
        __dist_x = __dfac*np.cos(np.deg2rad(self.Y))*(self.X-originx)
        __dist_y = __dfac*(self.Y-originy)
        
        __radial_distance = np.sqrt(__dist_x**2+__dist_y**2)
        
        return(__radial_distance, __dist_x, __dist_y)

class Sflux(object):
    def __init__(self, grid, basedate, nstep, path='./sflux'):
        self.grid = grid
        self.nstep = nstep # No of step
        self.basedate = basedate
        self.nfile = 0 # No file at the beginning
        self.path = path
        
        # Directory creation
        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __create_netcdf(self):
        self.step = 0
        self.nfile = self.nfile + 1
        self.__filename = 'sflux_air_1.{:03d}.nc'.format(self.nfile)
        self.__filepath = os.path.join(self.path, self.__filename)

        # Creating the file first
        self.nc = Dataset(self.__filepath, 'w', format='NETCDF4_CLASSIC')

        
        # Creating the dimensions
        self.__len_nx_grid = len(self.grid.x)
        self.__len_ny_grid = len(self.grid.y)
        self.__d_nx_grid = self.nc.createDimension(dimname='nx_grid', size=self.__len_nx_grid)
        self.__d_ny_grid = self.nc.createDimension(dimname='ny_grid', size=self.__len_ny_grid)
        self.__d_time = self.nc.createDimension(dimname='ntime', size=None)

        # Creating the variables
        # Time
        self.__v_time = self.nc.createVariable(varname='time', \
                                        datatype=np.float32, \
                                        dimensions=('ntime'))
        self.__v_time.units = 'days since {:s}'.format(self.basedate.strftime('%Y-%m-%d %H:%M:%S'))
        self.__v_time.long_name = 'Time'
        self.__v_time.calendar = 'standard'
        self.__v_time.base_date = self.basedate.timetuple()[0:4]

        # Longitude
        self.__v_lon = self.nc.createVariable(varname='lon', \
                                        datatype=np.float32, \
                                        dimensions=('ny_grid', 'nx_grid'))
        self.__v_lon.units = 'degrees_north'
        self.__v_lon.long_name = 'Longitude'
        self.__v_lon.standard_name = 'longitude'

        # Latitude
        self.__v_lat = self.nc.createVariable(varname='lat', \
                                        datatype=np.float32, \
                                        dimensions=('ny_grid', 'nx_grid'))
        self.__v_lat.units = 'degrees_east'
        self.__v_lat.long_name = 'Latitude'
        self.__v_lat.standard_name = 'latitude'

        # Uwind
        self.__v_uwind = self.nc.createVariable(varname='uwind', \
                                        datatype=np.float32, \
                                        dimensions=('ntime', 'ny_grid', 'nx_grid'))
        self.__v_uwind.units = 'm/s'
        self.__v_uwind.long_name = 'Surface Eastward Air Velocity (10m AGL)'
        self.__v_uwind.standard_name = 'eastward_wind'

        # Vwind
        self.__v_vwind = self.nc.createVariable(varname='vwind', \
                                        datatype=np.float32, \
                                        dimensions=('ntime', 'ny_grid', 'nx_grid'))
        self.__v_vwind.units = 'm/s'
        self.__v_vwind.long_name = 'Surface Northward Air Velocity (10m AGL)'
        self.__v_vwind.standard_name = 'northward_wind'

        # Prmsl
        self.__v_prmsl = self.nc.createVariable(varname='prmsl', \
                                        datatype=np.float32, \
                                        dimensions=('ntime', 'ny_grid', 'nx_grid'))
        self.__v_prmsl.units = 'Pa'
        self.__v_prmsl.long_name = 'Pressure Reduced to MSL'
        self.__v_prmsl.standard_name = 'air_pressure_at_mean_sea_level'

        # stmp
        self.__v_stmp = self.nc.createVariable(varname='stmp', \
                                        datatype=np.float32, \
                                        dimensions=('ntime', 'ny_grid', 'nx_grid'))
        self.__v_stmp.units = 'K'
        self.__v_stmp.long_name = 'Surface Temperature (2m AGL)'
        self.__v_stmp.standard_name = 'surface_temperature'

        # spfh
        self.__v_spfh = self.nc.createVariable(varname='spfh', \
                                        datatype=np.float32, \
                                        dimensions=('ntime', 'ny_grid', 'nx_grid'))
        self.__v_spfh.units = 1
        self.__v_spfh.long_name = 'Specific Humidity (2m AGL)'
        self.__v_spfh.standard_name = 'surface_specific_humidity'
        
        # Writing lon-lat once
        self.__v_lon[:] = self.grid.X
        self.__v_lat[:] = self.grid.Y
        
    
    def __putvalue(self, stepi, at, flux):
        self.__v_time[stepi] = date2num(dates=at, units=self.__v_time.units, calendar=self.__v_time.calendar)
        self.__v_uwind[stepi, :, :] = flux['uwind']
        self.__v_vwind[stepi, :, :] = flux['vwind']
        self.__v_prmsl[stepi, :, :] = flux['prmsl']
        self.__v_stmp[stepi, :, :] = flux['stmp']
        self.__v_spfh[stepi, :, :] = flux['spfh']

        self.step = self.step + 1

    def __close_netcdf(self):
        self.nc.sync()
        self.nc.close()

    def write(self, at, flux):
        # First check if self.nc is available
        if hasattr(self, 'nc'):
            if self.step < self.nstep:
                self.__putvalue(self.step, at, flux)
            else:
                self.__close_netcdf()
                self.__create_netcdf()
                self.__putvalue(self.step, at, flux)
        else:
            self.__create_netcdf()
            self.__putvalue(self.step, at, flux)

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
            f.write('air_1_max_window_hours={:.1f},	!max. # of hours (offset from start time in each file) in each file of set 1\n'.format(__max_window))
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

class DummyFlux(object):
    def __init__(self, grid):
        self.grid = grid

    def spfh(self):
        '''
        Generate surface specific humidity at a given timestamp. 
        Currently we are using a constant specific humidity which is taken as
        0.0175873 from SCHISM Sflux_nc scripts.
        '''
        __const = 0.0175873
        __spfh = np.ones(shape=self.grid.shape)*__const
        return(__spfh)

    def stmp(self):
        '''
        Generate surface temperature at a given timestep.
        Currently we are using a constant temperature at 300K according to 
        SCHISM Sflux_nc scripts.
        '''
        __const = 300
        __stmp = np.ones(shape=self.grid.shape)*__const
        return(__stmp)

class CFSR(object):
    def __init__(self, fname, var, to_grid):
        self.fname = fname
        self.var = var
        self.to_grid = to_grid

        self.__prep_file()

    def __prep_file(self):
        self.nc = Dataset(self.fname, mode='r')
        self.timeunits = self.nc.variables['time'].units
        self.time = num2date(times=self.nc['time'][:], units=self.timeunits, calendar='standard')
        self.lat = self.nc.variables['lat'][:]
        self.lon = self.nc.variables['lon'][:]
        self.in_grid = Grid(x=self.lon, y=self.lat)
    
    def at(self, timestamp):
        __tstamp = timestamp
        # Finding the proper or time interpolated data in original grid at timestamp
        if __tstamp in self.time:
            __tstamp = np.where(self.time == __tstamp)[0]
            __values = self.nc.variables[self.var][__tstamp, :, :]
            __values = np.reshape(__values, self.in_grid.shape)
        else:
            # Check if __timestamp outside of the input time bound
            if __tstamp > self.time[-1]:
                print('Error! Exceeds time bound, returning the nearest time value!')
                __index = np.searchsorted(a=self.time, v=__tstamp)
                __values = np.reshape(self.nc.variables[self.var][__index-1, :, :], self.in_grid.shape)
            elif __tstamp < self.time[0]:
                print('Error! Below time bound, returning the nearest time value!')
                __index = np.searchsorted(a=self.time, v=__tstamp)
                __values = np.reshape(self.nc.variables[self.var][__index, :, :], self.in_grid.shape)
            else:
                # Interpolate between timesteps
                __right = np.searchsorted(a=self.time, v=__tstamp)
                __left = __right - 1
                __right_value = np.reshape(self.nc.variables[self.var][__right, :, :], self.in_grid.shape)
                __left_value = np.reshape(self.nc.variables[self.var][__left, :, :], self.in_grid.shape)

                __w_left = (self.time[__right] - __tstamp).total_seconds()/float((self.time[__right] - self.time[__left]).total_seconds())
                __w_right = 1 - __w_left
                __values = __left_value*__w_left + __right_value*__w_right

        # Interpolation of __values
        # Interpolation makes the output consistent
        __points = np.array((self.in_grid.X.flatten(), self.in_grid.Y.flatten())).T
        __values = __values.flatten()
        __interp = interpolate.griddata(points=__points, values=__values, xi=(self.to_grid.X, self.to_grid.Y), method='linear')
        return(__interp)

if __name__=='__main__':
    # The grid definition and corpping of the track
    area = [79, 99, 10.5, 24.5]
    res = 0.2
    grid = Grid(x=np.arange(area[0], area[1]+res/2, res), y=np.arange(area[2], area[3]+res/2, res))

    # prmsl file
    fname_prmsl = '/run/media/khan/Workbench/Data/CFSRV2/Tide/prmsl/prmsl.nc'
    prmsl = CFSR(fname = fname_prmsl, var='PRMSL_L101', to_grid=grid)

    # uv file
    fname_uv = '/run/media/khan/Workbench/Data/CFSRV2/Tide/uv/uv.nc'
    uwind = CFSR(fname=fname_uv, var='U_GRD_L103', to_grid=grid)
    vwind = CFSR(fname=fname_uv, var='V_GRD_L103', to_grid=grid)

    # Dummy temperature and humidity flux
    dummy = DummyFlux(grid=grid)

    # Sflux definition
    sfluxpath = '/run/media/khan/Workbench/Data/CFSRV2/Tide/sflux'
    sfluxstart = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0)
    sflux = Sflux(grid=grid, basedate=sfluxstart, nstep=24*30, path=sfluxpath)
    dt = timedelta(hours=1)
    
    # Time loop
    starttime = datetime(year=2015, month=10, day=31, hour=0, minute=0, second=0)
    endtime = datetime(year=2018, month=11, day=1, hour=0, minute=0, second=0)
    at = timedelta()
    while starttime + at <= endtime:
        tstep = starttime + at
        print(tstep)
        flux = dict(uwind=uwind.at(tstep), vwind=vwind.at(tstep), prmsl=prmsl.at(tstep), stmp=dummy.stmp(), spfh=dummy.spfh())
        sflux.write(at=tstep, flux=flux)
        at = at + dt

    sflux.finish()
    sflux.sfluxtxt(dt=dt)