#!/usr/bin/env python

import numpy as np
import pandas as pd
import xarray as xr

from utide.astronomy import ut_astron
from utide.harmonics import FUV
from utide._ut_constants import ut_constants, constit_index_dict
from utide._time_conversion import _normalize_time as date2num

from netCDF4 import Dataset

from typing import Union, Tuple, List
from numpy.typing import ArrayLike
from pandas._typing import TimestampConvertibleTypes

from .atlas import Atlas
import os

import logging
logging.getLogger(__name__)

PathLike = Union[str, os.PathLike]

# required information from utide package
sat = ut_constants.sat
shallow = ut_constants.shallow

nshallow = np.ma.masked_invalid(ut_constants.const.nshallow).astype(int)
ishallow = np.ma.masked_invalid(ut_constants.const.ishallow).astype(int) - 1
not_shallow = ishallow.mask  # True where it was masked.
nshallow = nshallow.compressed()
ishallow = ishallow.compressed()
kshallow = np.nonzero(~not_shallow)[0]


# constituents from utide
mapping_utide = {
    wave:{'utide_name':wave, 'mapping_origin':'utide'} for wave in constit_index_dict
}

# constituents mapped from comodo and FES
comodo_utide = {
    'RO1':{'utide_name':'RHO1', 'mapping_origin':'FES'},
    'E2':{'utide_name':'EPS2', 'mapping_origin':'FES'},
    'LA2':{'utide_name':'LDA2', 'mapping_origin':'FES'},
    'MQ3':{'utide_name':'MO3', 'mapping_origin':'FES'},
    'M1':{'utide_name':'NO1', 'mapping_origin':'FES'},
    'KI1':{'utide_name':'CHI1', 'mapping_origin':'FES'},
    'TTA1':{'utide_name':'THE1', 'mapping_origin':'FES'}
}
mapping_utide.update(comodo_utide)

def utide_names(consts:ArrayLike):
    """Returns mapped utide constituents and missing constituents
    
    Returns the constituent names as found in utide as a {const:utide_const} dictionary. A second list is returned
    with names that are not found.

    Args:
        consts (ArrayLike): A constituent name, or a list of constituents names
    """
    consts = np.atleast_1d(consts)
    available = {}
    missing = []
    for const in consts:
        if const in mapping_utide:
            available[const] = mapping_utide[const]['utide_name']
        else:
            missing.append(const)
            logging.info(f'{const} is not found')
        
    return(available, missing)

def utide_freqs(consts:ArrayLike):
    """Returns the frequency of the constituents in cycle/hour

    Args:
        consts (ArrayLike): A constituent name, or a list of constituents name
    """
    available_consts, missing_consts = utide_names(consts)

    freq_dict = {
        name: freq for name, freq in zip(ut_constants.const.name, ut_constants.const.freq)
        }
    
    available_freq_dict = {
        const: freq_dict[available_consts[const]] for const in available_consts
    }

    return(available_freq_dict, missing_consts)


def utide_lind(consts:Union[dict, ArrayLike]) -> list:
    """Return the internal indices of utide constituents. 
    
    This is useful for getting the right output from various utide routines. The consts list should only contain utide 
    names. use `utide_names` to find available names.

    Args:
        consts (Union[dict, ArrayLike]): A constituent name, or a list of constituents names

    Raises:
        TypeError: `consts` should be list of dict type

    Returns:
        list: utide internal indices for the consts
    """
    if isinstance(consts, dict):
        consts_list = [consts[const] for const in consts]
    elif isinstance(consts, (np.ndarray, list)):
        consts_list = np.atleast_1d(consts)
    else:
        raise TypeError('consts should be one of list or dict type')
    
    lind = [constit_index_dict[const] for const in consts_list]
    
    return lind


def grid_FUV(
        timestamp:TimestampConvertibleTypes, 
        lat:Union[float, ArrayLike], 
        consts:Union[None, ArrayLike]=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the Nodal correction on grid

    Compute the Nodal correction to amplitude (F), phase (U), and the equilibrium argument(V) for the whole list of 
    constituent from utide at a given time t for a list of lat values. Adopted from utide.harmonics.FUV()

    Args:
        timestamp (TimestampConvertibleTypes): **single timestamp**
            Any formatted string 'yyyy-mm-dd HH:MM:SS', numpy.datetime64, pandas.datetime
        lat (Union[float, ArrayLike]): Latitude
        consts (Union[None, ArrayLike], optional): Constituent list.. Defaults to None.

    Raises:
        Exception: If len(t) > 1

    Returns:
        typing.Tuple[ArrayLike, ArrayLike, ArrayLike]: (F, U, V)
            F (nwaves, nlat): Nodal correction to amplitude (unitless multiplier), 
            U (nwaves, nlat): Nodal correction to phase (cycle)
            V (nwaves, nlat): Equilibrium argument (cycle)
    """
    const = ut_constants.const
    t = np.atleast_1d(date2num(np.array(timestamp)))
    if len(t) > 1:
        raise Exception('Only single value of t is acceptable in grid_FUV')
    lat = np.atleast_1d(lat)
    
    astro, ader = ut_astron(t)
    
    slat = np.sin(np.deg2rad(lat))
    rr = sat.amprat.copy()
    SLAT, RR = np.meshgrid(slat, rr)

    # Latitudinal correction for amplitude rate
    j = sat.ilatfac==1
    RR[j, :] *= 0.36309 * (1.0 - SLAT[j, :]**2 * 5.0) / SLAT[j, :]
    j = sat.ilatfac==2
    RR[j, :] *= SLAT[j, :] * 2.59808

    # 
    uu = np.dot(sat.deldood, astro[3:6, :]) + sat.phcorr[:, None]
    uu = np.fmod(uu, 1)
    mat = RR * np.exp(1j * 2 * np.pi * uu) # amprat x lat
    
    # 
    nfreq = len(const.isat)
    F = np.ones((nfreq, len(lat)), dtype=complex) # nfreq * lat
    iconst = sat.iconst - 1
    ind = np.unique(iconst)
    for i in ind:
        F[i, :] = 1 + np.sum(mat[iconst==i, :], axis=0)

    U = np.angle(F) / (2*np.pi) # cycles
    F = np.abs(F)

    # 
    V = np.dot(const.doodson, astro) + const.semi[:, None]
    V = np.fmod(V, 1)

    for i0, nshal, k in zip(ishallow, nshallow, kshallow):
        ik = i0 + np.arange(nshal)
        j = shallow.iname[ik] - 1
        exp1 = shallow.coef[ik, None]
        V[k, :] = np.sum(V[j, :]*exp1, axis=0)

    if consts is not None:
        lind = utide_lind(consts)
        F = F[lind]
        U = U[lind]
        V = V[lind]
        

    return (F, U, V)

def nodal_factor(t, consts, lat, correct_phase=True):
    """
    Compute the nodal correction, and astronomical argument at time t using the ut_FUV() in utide.

    Returns: Dictionary of Nodal correction, Astronomical argument for each constituent
    """
    t =np.atleast_1d(date2num(t))
    
    try:
        tref = (t[0] + t[-1])/2
    except:
        tref = t[0]

    lind = utide_lind(consts)

    nf, U, V = FUV(t=t, tref=tref, lind=lind, lat=lat, ngflgs=[0, 0, 0, 0])

    if correct_phase:
        ear = (U + V) * 360 # cycles to degree
    else:
        ear = V * 360
        
    ear %= 360 # always positive from 0-360

    nf_dict = {
        const:{'nf':cnf, 'ear':cear} for const, cnf, cear in zip(consts, nf, ear)
    }

    return nf_dict

def reconstruct_waterlevel(
        fname:PathLike, 
        atlas:Atlas, 
        timestamps:ArrayLike, 
        epoch:Union[None, TimestampConvertibleTypes]=None) -> xr.Dataset:
    """Reconstruct waterlevel from an atlas

    Args:
        fname (PathLike): File path for saving the reconstructed water levels
        atlas (Atlas): A tidal atlas
        timestamps (ArrayLike): A series of times when the water levels are to be reconstructed.
        epoch (Union[None, TimestampConvertibleTypes], optional): Epoch to be used in the output nc file. Defaults to None.

    Raises:
        Exception: `epoch` should be a TimestampConvertibleTypes

    Returns:
        xr.Dataset: The reconstructed dataset

    TODO: Reorganize the netcdf dataset creation process
    """
    timestamps = pd.to_datetime(np.array(timestamps))
    if epoch is not None:
        try:
            init_time = pd.Timestamp(epoch)
        except:
            raise Exception('epoch should be an datetime convertible type, e.g., `yyyy-mm-dd HH:MM:SS`')
    else:
        init_time = timestamps[0]

    const_available, const_missing = utide_names(atlas.waves)
    atlas_utide = atlas.select(const_available)
    if len(const_missing):
        print(f'{len(const_missing)} constituents - {" ".join(const_missing)} - are missing in utide.')
        print(f'Reconstructing from {" ".join(list(const_available.keys()))}')

    lon = atlas_utide.lon
    lat = atlas_utide.lat

    amp, pha = atlas_utide.values(consts=atlas_utide.waves)
    lind = utide_lind(atlas_utide.waves)

    elev_unit = atlas_utide.units['amp']

    with Dataset(fname, 'w', format='NETCDF4_CLASSIC') as nc:
        # Dimensions
        nc.createDimension(dimname='lon', size=len(lon))
        nc.createDimension(dimname='lat', size=len(lat))
        nc.createDimension(dimname='time', size=None)
        if atlas.grid_type=='points':
            nc.createDimension(dimname='points', size=len(lon))

        # Dimensions variables
        if atlas.grid_type=='points':
            var_lon = nc.createVariable(
                varname='lon',
                datatype=np.float32,
                dimensions=('points')
            )
        else:
            var_lon = nc.createVariable(
                varname='lon',
                datatype=np.float32,
                dimensions=('lon')
            )
        var_lon.units = 'degrees_east'
        var_lon.long_name = 'Longitude'
        var_lon.standard_name = 'longitude'
        var_lon[:] = lon

        if atlas.grid_type=='points':
            var_lat = nc.createVariable(
                varname='lat',
                datatype=np.float32,
                dimensions=('points')
            )
        else:
            var_lat = nc.createVariable(
            varname='lat',
            datatype=np.float32,
            dimensions=('lat')
        )
        var_lat.units = 'degrees_north'
        var_lat.long_name = 'Latitude'
        var_lat.standard_name = 'latitude'
        var_lat[:] = lat


        var_time = nc.createVariable(
            varname='time',
            datatype=np.float32,
            dimensions=('time'),
            
        )
        var_time.units = f'hours since {init_time.strftime("%Y-%m-%d %H:%M:%S")}'
        var_time.long_name = 'Time'
        var_time.calendar = 'standard'
        var_time[:] = (timestamps - init_time).total_seconds()/3600

        # Elevations
        if atlas.grid_type=='points':
            var_elev = nc.createVariable(
                varname='elev',
                datatype=np.float32,
                dimensions=('time', 'points')
            )
        else:
            var_elev = nc.createVariable(
                varname='elev',
                datatype=np.float32,
                dimensions=('time', 'lat', 'lon')
            )
        
        var_elev.units = elev_unit
        var_elev.long_name = 'Tidal water level elevation'

        nc.description = f'Water level created from tidal constituents in python'
        nc.constituents = ','.join(list(const_available.keys()))

        for i, t in enumerate(timestamps):
            F, U, V = grid_FUV(timestamp=t, lat=lat)
            F = F[lind, :, None]
            U = U[lind, :, None]*360
            V = V[lind, :, None]*360

            if atlas.grid_type=='points':
                var_elev[i, :] = np.sum(amp * F * np.cos(np.deg2rad((U + V - pha) % 360)), axis=0)
            else:
                var_elev[i, :, :] = np.sum(amp * F * np.cos(np.deg2rad((U + V - pha) % 360)), axis=0)

            if i%100 == 0:
                nc.sync()
    
    return xr.open_dataset(fname)
