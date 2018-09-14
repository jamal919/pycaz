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

class Grid(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='xy')
        self.shape = self.X.shape
        self.length = len(self.X.flat)

    def radial_distance(self, originx, originy):
        pass

class Sflux(object):
    def __init__(self, grid, basedate, nstep, path='./'):
        self.grid = grid
        self.path = path
        self.nstep = nstep # No of step
        self.basedate = basedate
        self.nfile = 0 # No file at the beginning

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

class Converter(object):
    def __init__(self):
        pass

    @staticmethod
    def timestamp(d):
        return(calendar.timegm(d.timetuple()))

    @staticmethod
    def hpa2pa(hpa):
        return(hpa*100)
    
    @staticmethod
    def knot2mps(knot):
        return(knot*1.852/3.6)

    @staticmethod
    def km2m(km):
        return(km*1000)

    @staticmethod
    def lon180(lon360):
        lon360[lon360 > 180] = lon360[lon360 > 180] - 360
        return(lon360)

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
        __dfac = 60*1.852*1000
        
        # Calculating __utrans and __vtrans for all timestep except the last
        for i in np.arange(0, len(self.track)-1):
            __dt = self.track[i+1]['time'] - self.track[i]['time']
            __dt = __dt.total_seconds()
            __lon1 = np.deg2rad(self.track[i]['lon'])
            __lon2 = np.deg2rad(self.track[i+1]['lon'])
            __lat1 = np.deg2rad(self.track[i]['lat'])
            __lat2 = np.deg2rad(self.track[i+1]['lat'])
            __dtrans_x = __dfac*np.cos(__lat1)*(__lon2-__lon1)
            __dtrans_y = __dfac*(__lat2-__lat1)
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


    def interpolate(self, var, at, printval=False):
        __at = at
        __var = var
        __weights = self.__find_weight(__at)
        
        if __var in self.track[0].keys():
            if __var=='lat' or __var=='lon':
                pass
            elif __var=='vtrans' or __var=='utrans':
                __index = __weights['index'][0]
                return(self.track[__index][__var])
            else:
                __v_left = self.track[__weights['index'][0]][__var]
                __v_right = self.track[__weights['index'][1]][__var]
                __w_left = __weights['weight'][0]
                __w_right = __weights['weight'][1]
                __interp_val = __v_left*__w_left + __v_right*__w_right
                if printval:
                    print(__at.total_seconds(), \
                        self.track[__weights['index'][0]][__var], \
                        self.track[__weights['index'][1]][__var], \
                        __interp_val)
                return(__interp_val)
        else:
            print('Variable {:s} not found in the track'.format(__var))


class Generator(object):
    def __init__(self, track, grid, pn=101300, rhoair=1.15, transfac=0.56):
        self.track = track
        self.grid = grid
        self.pn = pn
        self.rhoair = rhoair
        self.transfac = transfac
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
    
    @staticmethod
    def __holland_B(vmax, p, rhoair, pn, bmax=2.5, bmin=0.5):
        __B = vmax**2*rhoair*np.exp(1)/(pn-p)
        
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
    def __find_r_e11(v, vmax, rmax, f, solver='bisect', limit=[500000, 0], step=-100):
        __resfunc = lambda __r: v - (2*__r*(vmax*rmax + 0.5*f*rmax**2)/\
                    (rmax**2+__r**2)-f*__r/2)
        
        if solver=='bisect':
            __rsolved = optimize.bisect(__resfunc, a=limit[0], b=rmax)
            return(__rsolved)
        elif solver=='scan':
            __rrange = np.arange(start=limit[0], stop=limit[1], step=step)
            __res = np.array([__resfunc(i) for i in __rrange])
            __loc = np.where(__res < 0)[0]
            __rsolved = __rrange[__loc[0]]
            return(__rsolved)

    @staticmethod
    def __find_rmax_h80(vx, rx, p, pn, rhoair, B, f, solver='bisect', limit=[1000, 100000], step=100):
        __resfunc = lambda __R: vx - np.sqrt(B/rhoair*(__R/rx)**B*(pn-p)*np.exp(-(__R/rx)**B)+(rx*f/2)**2) - (rx*f/2)

        if solver=='bisect':
            __rsolved = optimize.bisect(__resfunc, a=limit[0], b=rx)
            return(__rsolved)
        elif solver=='scan':
            __rrange = np.arange(start=limit[0], stop=limit[1], step=step)
            __res = np.array([__resfunc(i) for i in __rrange])
            __loc = np.where(__res < 0)[0]
            __rsolved = __rrange[__loc[0]]
            return(__rsolved)
    
    def __generate_wind(self):
        '''
        Calculate the u- and v- wind component at a given timestep
        '''
        pass

    def __generate_pressure(self, at):
        '''
        Calculate the pressure at a given timestep
        '''
        __at = at
        __pc = self.track.interpolate(var='p', at=__at, printval=False)
    
    def __generate_stmp(self):
        '''
        Generate surface temperature at a given timestep.
        Currently we are using a constant temperature at 300K according to 
        SCHISM Sflux_nc scripts.
        '''
        __const = 300
        __stmp = np.ones(shape=self.grid.shape)*__const
        return(__stmp)
    
    def __generate_sphd(self):
        '''
        Generate surface specific humidity at a given timestamp. 
        Currently we are using a constant specific humidity which is taken as
        0.0175873 from SCHISM Sflux_nc scripts.
        '''
        __const = 0.0175873
        __sphd = np.ones(shape=self.grid.shape)*__const
        return(__sphd)

    def generate(self, at):
        '''
        Generate the fields and return a dict of the field at a given time.
        '''
        # Loading the time
        __at = at

        # Interpolating various values
        __lon = self.track.interpolate(var='lon', at=__at, printval=False)
        __lat = self.track.interpolate(var='lat', at=__at, printval=False)
        __p = self.track.interpolate(var='p', at=__at, printval=False)
        __rmax_e11 = self.track.interpolate(var='rm', at=__at, printval=True)
        __vmax = self.track.interpolate(var='vmax', at=__at, printval=False)
        __utstorm = self.track.interpolate(var='utstorm', at=__at, printval=False)
        __vtstorm = self.track.interpolate(var='vtstorm', at=__at, printval=False)

        # Calculating parameters
        __B = self.__holland_B(vmax=__vmax, p=__p, rhoair=self.rhoair, pn=self.pn, bmax=2.5, bmin=0.5)
        __f = self.__coriolis(lat=__lat, isradians=False)

        # Calculating transitioning speed and radius
        # Transition of E11 and H80 at 50knots
        __v50 = self.converter.knot2mps(50)
        __r50 = self.__find_r_e11(v=__v50, vmax=self.surf2bl(__vmax), \
                                rmax=__rmax_e11, f=__f, \
                                solver='bisect', \
                                limit=[500000, 0], step=-100)

        # Calculating rmax for H80 model
        __rmax_h80 = self.__find_rmax_h80(vx=__v50, rx=__r50, p=__p, pn=self.pn, \
                                        rhoair=self.rhoair, B=__B, f=__f, \
                                        solver='scan', \
                                        limit=[1000, 100000], step=100)

        # __uwind, __vwind = self.__generate_wind()
        # __prmsl = self.__generate_pressure()
        # __stmp = self.__generate_stmp()
        # __sphd = self.__generate_sphd()

        # # Returning values
        # return(dict(uwind=__uwind, vwind=__vwind, prmsl=__prmsl, stmp=__stmp, sphd=__sphd))

        

if __name__=='__main__':
    # The grid definition
    area = [79, 99, 10.5, 24.5]
    res = 0.025
    grid = Grid(x=np.arange(area[0], area[1]+res/2, res), y=np.arange(area[2], area[3]+res/2, res))
    
    # Track reading
    trackreader = Reader(readername='kerry')
    trackfile = trackreader.read('/run/media/khan/Workbench/Projects/Surge Model/Emmanuel et al/tracks_csv/Track_0001.csv')
    track = Track(track=trackfile, clipby=area)

    # Generator
    generator = Generator(track=track, grid=grid)

    # sflux file creation
    sflux = Sflux(grid=grid, basedate=track.basedate, nstep=96, path='/run/media/khan/Workbench/Projects/Surge Model/Emmanuel et al/Sflux SCHISM/')
    at = timedelta() # Starts at 0
    dt = timedelta(minutes=15) # Time step
    while track.basedate + at <= track.lastdate:
        flux = generator.generate(at=at)
        # sflux.write(at=at, flux=flux)
        at = at + dt
    sflux.finish()