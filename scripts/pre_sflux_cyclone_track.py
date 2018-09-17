# -*- coding: utf-8 -*-
"""
Create wind and pressure field from cyclone track and write it into SCHISM
complient sflux file. 

@author: khan
@email: jamal.khan@legos.obs-mip.fr
"""

from __future__ import print_function
import sys
import os
import numpy as np
from scipy import optimize
from datetime import datetime, timedelta
import calendar
import time
from netCDF4 import Dataset

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
        self.length = len(self.X.flat)
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
        self.__v_time[stepi] = at.days + at.seconds/float(86400)
        self.__v_uwind[stepi, :, :] = flux['uwind']
        self.__v_vwind[stepi, :, :] = flux['vwind']
        self.__v_prmsl[stepi, :, :] = flux['prmsl']
        self.__v_stmp[stepi, :, :] = flux['stmp']
        self.__v_spfh[stepi, :, :] = flux['spfh']

        self.step = self.step + 1

    def __close_netcdf(self):
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

class Reader(object):
    '''
    Reader(readername) create a reader object for cyclone track file.
    '''
    def __init__(self, readername):
        self.readername = readername
        self.converter = Converter()

    def read(self, path):
        self.path = path
        self.track = getattr(self, self.readername)()
        return(self.track)

    def kerry(self):
        # Reading the track file
        __values = np.genfromtxt(fname=self.path, dtype=float, delimiter=',', skiprows=3)
        
        # Creating the datetime object
        __datetime = np.array([datetime(int(__values[i, 0]),\
                            int(__values[i, 1]),\
                            int(__values[i, 2]),\
                            int(__values[i, 3]),\
                            int(__values[i, 4]),\
                            int(__values[i, 5])) \
                            for i in np.arange(len(__values))])
        # Processing other variables
        __lon = np.array(self.converter.lon180(lon360=__values[:, 6])) # -180 to 180
        __lat = np.array(__values[:, 7])
        __p = np.array(self.converter.hpa2pa(hpa=__values[:, 8])) # hPa to Pa
        __rm = np.array(self.converter.km2m(km=__values[:, 9])) # Km to m
        __vmax = np.array(self.converter.knot2mps(knot=__values[:, 10])) # knot to mps
        __track = np.array([dict(time=__datetime[i],\
                    lon=__lon[i],\
                    lat=__lat[i],\
                    p=__p[i],\
                    rm=__rm[i],\
                    vmax=__vmax[i]) for i in np.arange(len(__datetime))])
        return(__track)

class Track(object):
    def __init__(self, track, clipby=None):
        self.track = track
        self.converter = Converter()

        # Calculate translation speed
        self.__calc_translation()
        
        # Clipping the track
        if clipby is None:
            self.__gen_basedate()
        else:
            self.__clip(by=clipby)
            self.__gen_basedate()

        # Calc time index
        self.__calc_timeindex()

    def __calc_translation(self):
        
        # Calculating __utrans and __vtrans for all timestep except the last
        for i in np.arange(0, len(self.track)-1):
            __dt = self.track[i+1]['time'] - self.track[i]['time']
            __dt = __dt.total_seconds()
            __origin = (self.track[i]['lon'], self.track[i]['lat'])
            __of = (self.track[i+1]['lon'], self.track[i+1]['lat'])
            __dtrans_x, __dtrans_y = self.converter.gc_distance(of=__of, origin=__origin, isradians=False)
            __utstorm = __dtrans_x/__dt
            __vtstorm = __dtrans_y/__dt
            __tstorm = np.sqrt(__utstorm**2+__vtstorm**2)
            self.track[i]['utstorm'] = __utstorm
            self.track[i]['vtstorm'] = __vtstorm

        # For the last time step, we are keeping it to the same
            self.track[len(self.track)-1]['utstorm'] = __utstorm
            self.track[len(self.track)-1]['vtstorm'] = __vtstorm

    def __clip(self, by):
        '''
        Clip the track to a given boundary
        TODO: Check if the track time is OK
        '''
        if isinstance(by, list) and len(by)==4:
            # by lonlatbox cropping
            __lonmin = by[0]
            __lonmax = by[1]
            __latmin = by[2]
            __latmax = by[3]

            __lonlist = np.array([self.track[i]['lon'] for i in np.arange(len(self.track))])
            __lonselected = np.where((__lonlist >= __lonmin) & (__lonlist <= __lonmax))[0]
            __mod_track = self.track[__lonselected]
            __latlist = np.array([__mod_track[i]['lat'] for i in np.arange(len(__mod_track))])
            __latselected = np.where((__latlist >= __latmin) & (__latlist <= __latmax))
            __mod_track = __mod_track[__latselected]

            self.track = __mod_track

        elif isinstance(by, Grid):
            #TODO Implement clipping by grid selection
            print('Clipping by grid not yet implemented')
            sys.exit(1)
        else:
            print('Track clipping failed! Give a lonlatbox or Grid input. Aborting...')
            sys.exit(1)

    def __gen_basedate(self):
        self.basedate = self.track[0]['time']
        self.lastdate = self.track[len(self.track)-1]['time']

    def __calc_timeindex(self):
        self.timeindex = np.array([(record['time']-self.basedate).total_seconds() for record in self.track])

    def __find_indice(self, at):
        if(isinstance(at, timedelta)):
            __at = at.total_seconds()
        else:
            __at = at
        
        __right = np.searchsorted(a=self.timeindex, v=__at, side='left')
        
        if __right == 0:
            __left = 0
        else:
            __left = __right - 1

        return([__left, __right])

    def __find_weight(self, at):
        if(isinstance(at, timedelta)):
            __at = at.total_seconds()
        else:
            __at = at

        __indices = dict(index = self.__find_indice(at))
        __i_left = __indices['index'][0]
        __i_right = __indices['index'][1]
        __v_right = self.timeindex[__i_right]
        __v_left = self.timeindex[__i_left]
        
        if(__i_left != __i_right):
            __w_left = (__v_right - __at)/float((__v_right - __v_left))
            __w_right = 1-__w_left
        else:
            __w_left = float(1)
            __w_right = float(0)

        __indices['weight'] = [__w_left, __w_right]
        return(__indices)


    def interpolate(self, var, at):
        __at = at
        __var = var
        __weights = self.__find_weight(__at)
        
        if __var in self.track[0].keys():
            if __var=='vtrans' or __var=='utrans':
                __index = __weights['index'][0]
                __interp_val = self.track[__index][__var]
                return(__interp_val)
            else:
                __v_left = self.track[__weights['index'][0]][__var]
                __v_right = self.track[__weights['index'][1]][__var]
                __w_left = __weights['weight'][0]
                __w_right = __weights['weight'][1]
                __interp_val = __v_left*__w_left + __v_right*__w_right
                return(__interp_val)
        else:
            print('Variable {:s} not found in the track'.format(__var))

    def trackinfo(self, filepath):
        with open(filepath, 'w') as f:
            __basedate = self.basedate.timetuple()
            __rnday = self.lastdate-self.basedate
            __rnday = __rnday.total_seconds()/timedelta(days=1).total_seconds()
            f.write('start_year={:d}\n'.format(__basedate[0]))
            f.write('start_month={:d}\n'.format(__basedate[1]))
            f.write('start_day={:d}\n'.format(__basedate[2]))
            f.write('start_hour={:d}\n'.format(__basedate[3]))
            f.write('rnday={:.2f}\n'.format(__rnday))

class Generator(object):
    def __init__(self, track, grid, pn=101300, rhoair=1.15, transfac=0.56, transangle=19.2):
        self.track = track
        self.grid = grid
        self.pn = pn
        self.rhoair = rhoair
        self.transfac = transfac # Lin and Chavaz
        self.transangle = np.deg2rad(transangle) # Lin and Chavaz
        self.converter = Converter()

    @staticmethod
    def __coriolis(lat, isradians=False):
        ''' 
        Calculate the coriolis coefficient.
        Uses the simple function - 2*7.292e-5*sin(lat)
        '''
        if isradians:
            __f = 2*7.292e-5*np.sin(lat)
        else:
            __f = 2*7.292e-5*np.sin(np.deg2rad(lat))

        return(__f)
    
    def __holland_B(self, vmax_bl, p, bmax=2.5, bmin=0.5):
        __B = vmax_bl**2*self.rhoair*np.exp(1)/float((self.pn-p))
        
        if __B<bmin:
            __B = bmin
        elif __B>bmax:
            __B = bmax
        else:
            __B = __B
        
        return(__B)

    @staticmethod
    def surf2bl(vsurf, swrf=0.9):
        return(vsurf/float(swrf))

    @staticmethod
    def bl2surf(vbl, swrf=0.9):
        return(vbl*swrf)

    def __find_r_e11(self, v, vmax, rmax, f, solver='bisect', limit=[500000, 0], step=-100):
        __resfunc = lambda __r: v - (2*__r*(vmax*rmax+0.5*f*rmax**2)/\
                    (rmax**2+__r**2)-f*__r/float(2))
        if solver=='scan':
            __rrange = np.arange(start=limit[0], stop=limit[1]+step, step=step)
            for __r in __rrange:
                __res = __resfunc(__r)
                if __res < 0:
                    break
            __rsolved = __r
        elif solver=='bisect':
            __rsolved = optimize.bisect(f=__resfunc, a=limit[0], b=rmax)
        return(__rsolved)

    def __find_rmax_h80(self, vx, rx, p, B, f, solver='scan', limit=[1000, 100000], step=100):
        __resfunc = lambda __R: vx - (np.sqrt(B/float(self.rhoair)*(__R/float(rx))**B*(self.pn-p)*np.exp(-(__R/float(rx))**B)+(rx*f/float(2))**2) - rx*f/float(2))

        __rmrange = np.arange(start=limit[0], stop=limit[1]+step, step=step)
        for __rm in __rmrange:
            __res = __resfunc(__rm)
            if __res < 0:
                break
        __rsolved = __rm
        return(__rsolved)

    def __calc_pressure(self, p, rmax, rgrid, B):
        '''
        Calculate the pressure using Holland 1980 model
        '''
        __pn = self.pn
        __h80_pressure = lambda r: (__pn-p)*np.exp(-(rmax/(r+1e-8))**B)+p
        __pressure = np.array([__h80_pressure(r) for r in rgrid])
        return(np.reshape(__pressure, newshape=self.grid.shape))
    
    def __calc_holland_wind(self, rgrid, p, rmax, B, f):
        '''
        Calculate the u- and v- wind component at a given timestep using Holland
        1980 Model. The returned wind speed is in boundary layer.
        '''
        __pn = self.pn
        __rhoair = self.rhoair
        # Applying Holland 80 formula
        __vcirc = np.sqrt(B/__rhoair*(rmax/rgrid)**B*(__pn-p)*np.exp(-(rmax/rgrid)**B)+(rgrid*f/2)**2)-rgrid*f/2

        # Transfer it to surface wind
        __vcirc = self.bl2surf(vbl=__vcirc)

        return(__vcirc)

    def __calc_marged_wind(self, rgrid, p, rmax_h80, vmax, rmax_e11, B, f, threshold):
        __pn = self.pn
        __rhoair = self.rhoair
        __vcirc = np.zeros(shape=rgrid.shape)

        __h80 = lambda r: np.sqrt(B/float(__rhoair)*(rmax_h80/r)**B*(__pn-p)*np.exp(-(rmax_h80/r)**B)+(r*f/float(2))**2)-r*f/float(2)
        __e11 = lambda r: 2*r*(vmax*rmax_e11+0.5*f*rmax_e11**2)/(rmax_e11**2+r**2)-f*r/float(2)

        __vcirc[rgrid > threshold] = __h80(rgrid[rgrid > threshold])
        __vcirc[rgrid <= threshold] = __e11(rgrid[rgrid <= threshold])

        # Make sure no negative velocity
        __vcirc[__vcirc < 0] = 0

        # Transfer it to surface wind
        __vcirc = self.bl2surf(vbl=__vcirc)

        return(__vcirc)

    
    def __apply_translation(self, vcirc, rgrid, dxgrid, dygrid, utstorm, vtstorm):
        '''
        Apply translation speed correction to circular wind speed and return
        uwind vwind.
        '''
        __uwind = -vcirc*dygrid/rgrid
        __vwind = vcirc*dxgrid/rgrid
        
        __ustorm_modif = utstorm*np.cos(self.transangle)-vtstorm*np.sin(self.transangle)
        __vstorm_modif = utstorm*np.sin(self.transangle)+vtstorm*np.cos(self.transangle)

        __uwind = __uwind+self.transfac*__ustorm_modif
        __vwind = __vwind+self.transfac*__vstorm_modif

        __uwind = np.reshape(__uwind, newshape=self.grid.shape)
        __vwind = np.reshape(__vwind, newshape=self.grid.shape)

        return(__uwind, __vwind)

    def __one2ten(self, wind_1min, factor=0.88):
        return(wind_1min*factor)
        

    def __generate_wind(self, vcirc, rgrid, dxgrid, dygrid, utstorm, vtstorm):
        '''
        Generate u- and v- wind from circular wind speed. 
        '''
        
        __uwind, __vwind = self.__apply_translation(vcirc=vcirc, rgrid=rgrid, \
                                                    dxgrid=dxgrid, dygrid=dygrid,\
                                                    utstorm=utstorm, vtstorm=vtstorm)
        
        __uwind = self.__one2ten(wind_1min=__uwind, factor=0.88)
        __vwind = self.__one2ten(wind_1min=__vwind, factor=0.88)

        return(__uwind, __vwind)
    
    def __generate_stmp(self):
        '''
        Generate surface temperature at a given timestep.
        Currently we are using a constant temperature at 300K according to 
        SCHISM Sflux_nc scripts.
        '''
        __const = 300
        __stmp = np.ones(shape=self.grid.shape)*__const
        return(__stmp)
    
    def __generate_spfh(self):
        '''
        Generate surface specific humidity at a given timestamp. 
        Currently we are using a constant specific humidity which is taken as
        0.0175873 from SCHISM Sflux_nc scripts.
        '''
        __const = 0.0175873
        __spfh = np.ones(shape=self.grid.shape)*__const
        return(__spfh)

    def generate(self, at):
        '''
        Generate the fields and return a dict of the field at a given time.
        '''
        # Loading the time
        __at = at

        # Interpolating various values
        __lon = self.track.interpolate(var='lon', at=__at)
        __lat = self.track.interpolate(var='lat', at=__at)
        __p = self.track.interpolate(var='p', at=__at)
        __rmax_e11 = self.track.interpolate(var='rm', at=__at)
        __vmax = self.track.interpolate(var='vmax', at=__at)
        __utstorm = self.track.interpolate(var='utstorm', at=__at)
        __vtstorm = self.track.interpolate(var='vtstorm', at=__at)
        
        # Calculating parameters
        __B = self.__holland_B(vmax_bl=self.surf2bl(__vmax), p=__p, bmax=2.5, bmin=0.5)
        __f = self.__coriolis(lat=__lat, isradians=False)

        # Calculating transitioning speed and radius
        # Transition of E11 and H80 at 50knots
        __v50 = self.converter.knot2mps(50)

        # Calculating the radial distance grid from __lon, __lat
        __r_grid, __dx_grid, __dy_grid = self.grid.radial_distance(originx=__lon, originy=__lat, isradians=False)
        
        # Calculating the pressure grid
        # Using Holland 1980 model
        __prmsl = self.__calc_pressure(p=__p, rmax=__rmax_e11, rgrid=__r_grid, B=__B)

        # Selection of model for wind speed
        if __vmax < __v50:
            # Calculate pressure and wind field both using Holland 1980 Model
            # And transfer it to the surface
            # See Krien et al. (unpublished) for justification
            __vcirc = self.__calc_holland_wind(rgrid=__r_grid, \
                                                p=__p, \
                                                rmax=__rmax_e11, \
                                                B=__B, f=__f)
        else:
            # Calculate upto r50 using Emmanuel 2011 Model
            # Calculate the rest using Holland 1980 Model
            # And finally transfer it to the surface
            __r50 = self.__find_r_e11(v=__v50, vmax=self.surf2bl(__vmax), \
                            rmax=__rmax_e11, f=__f, \
                            solver='bisect',\
                            limit=[500000, 0], step=-100)

            # Calculating rmax for H80 model
            __rmax_h80 = self.__find_rmax_h80(vx=__v50, rx=__r50, \
                                    p=__p, B=__B, f=__f, \
                                    solver='scan', \
                                    limit=[1000, 100000], step=100)
            
            # Calculating circular wind based on criteria
            __vcirc = self.__calc_marged_wind(rgrid=__r_grid, p=__p, \
                                            rmax_h80=__rmax_h80, \
                                            vmax=self.surf2bl(__vmax), \
                                            rmax_e11=__rmax_e11, \
                                            B=__B, f=__f,\
                                            threshold=__r50)
        
        # Calculating the wind grid
        __uwind, __vwind = self.__generate_wind(vcirc=__vcirc, rgrid=__r_grid, \
                                                dxgrid=__dx_grid, dygrid=__dy_grid,\
                                                utstorm=__utstorm, vtstorm=__vtstorm)
        
        __stmp = self.__generate_stmp()
        __spfh = self.__generate_spfh()

        # Returning values
        return(dict(uwind=__uwind, vwind=__vwind, prmsl=__prmsl, stmp=__stmp, spfh=__spfh))

        

if __name__=='__main__':
    # Track file
    trackpath = './Track.csv'
    sfluxpath = './sflux'
    
    # The grid definition and corpping of the track
    area = [79, 99, 10.5, 24.5]
    res = 0.025
    grid = Grid(x=np.arange(area[0], area[1]+res/2, res), y=np.arange(area[2], area[3]+res/2, res))

    # Track reading
    trackreader = Reader(readername='kerry')
    trackfile = trackreader.read(trackpath)
    track = Track(track=trackfile, clipby=area)

    # Writing trackinfo file with 30 minutes reduction in runtime to avoid
    # no-data issue
    track.trackinfo(filepath='./trackinfo')

    # Generator
    generator = Generator(track=track, grid=grid)

    # sflux object creation
    sflux = Sflux(grid=grid, basedate=track.basedate, nstep=96, path=sfluxpath)

    # Time loop
    at = timedelta() # Starts at 0
    dt = timedelta(minutes=15) # Time step
    while track.basedate + at <= track.lastdate:
        print(datetime.strftime(track.basedate+at, '%Y-%m-%d %H:%M:%S'))
        flux = generator.generate(at=at)
        sflux.write(at=at, flux=flux)
        at = at + dt
    # Adding an extra timestep and finishing the file
    sflux.write(at=at, flux=flux)
    sflux.finish()
    sflux.sfluxtxt(dt=dt)