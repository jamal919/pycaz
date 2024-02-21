#!/usr/bin/env python
from __future__ import annotations
from abc import ABC, abstractmethod
import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os
from typing import Union
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
from numpy.typing import ArrayLike

from .utilities import grid_around
from .interpolate import interp_complex_2D

PathLike = Union[str, os.PathLike]

# Abstract Dataset Class
# Each dataset derived from this class must contain the defined properties and methods at minimum
class Dataset(ABC):
    @property
    @abstractmethod
    def lon(self):
        pass
    
    @property
    @abstractmethod
    def lat(self):
        pass

    @property
    @abstractmethod
    def amp(self):
        pass

    @property
    @abstractmethod
    def pha(self):
        pass

    @property
    @abstractmethod
    def values(self):
        pass

    @abstractmethod
    def interp(self):
        pass
    
    @abstractmethod
    def __repr__(self):
        pass
    
    @abstractmethod
    def _repr_html_(self):
        pass

    @abstractmethod
    def to_netcdf(self):
        pass

# Class definitions
class PointsDataset(Dataset):
    def __init__(self, dataset:xr.Dataset, units:dict):
        """A class to contain the `amp` and `pha` on a list of `points`.

        Args:
            dataset (xr.Dataset): xarray dataset containing the `amp` and `pha`.
            units (dict): units dictionary, e.g., {'amp':'m', 'pha':'degrees'}
        """
        self.ds = dataset
        self.units = units

    @property
    def lat(self):
        return self.ds['lat']
    
    @property
    def lon(self):
        return self.ds['lon']
    
    @property
    def amp(self):
        return self.ds['amp']
    
    @property
    def pha(self):
        return self.ds['pha']
    
    @property
    def values(self):
        return np.vstack([self.amp, self.pha]).T
    
    def interp(self):
        raise Exception('Interpolation is not yet implemented for PointsDataset')
    
    def to_netcdf(self, fname):
        self.ds.to_netcdf(fname)

    def __repr__(self):
        return(self.ds.__repr__())

    def _repr_html_(self):
        return(self.ds._repr_html_())

class GriddedDataset(Dataset):
    def __init__(self, dataset:xr.Dataset, units:dict):
        """A class to contain the `amp` and `pha` on a regular `lat-lon` grid.

        Args:
            dataset (xr.Dataset): xarray dataset containing the `amp` and `pha`.
            units (dict): units dictionary, e.g., {'amp':'m', 'pha':'degrees'}
        """
        self.ds = dataset
        self.units = units

    @property
    def lon(self):
        return self.ds['lon']

    @property
    def lat(self):
        return self.ds['lat']

    @property
    def amp(self):
        return self.ds['amp']

    @property
    def pha(self):
        return self.ds['pha']
    
    @property
    def values(self):
        return np.stack([self.amp, self.pha]).T # ()
    
    def _interp_xy(
            self, 
            xy:ArrayLike, 
            method:Literal['complex', 'linear', 'nearest']='complex', 
            extrapolate:Literal['spherical', 'nearest']='spherical') -> PointsDataset:
        xy = np.atleast_2d(xy)
        amp_pha = np.empty(shape=xy.shape)
        if method=='complex':
            for i, ixy in enumerate(xy):
                bnds = grid_around(
                    xy=ixy,
                    grid_x=self.lon,
                    grid_y=self.lat,
                    extrapolate=extrapolate
                )
                amp_pha_grid = self.sel(xy=bnds).values # returns npoints x 2 (amp, pha) values
                iamp, ipha = interp_complex_2D(
                    xy=ixy,
                    bnds=bnds,
                    amp_pha=amp_pha_grid,
                    pha_unit=self.units['pha']
                )

                amp_pha[i, :] = np.array([iamp, ipha])
        else:
            xy_df = pd.DataFrame(xy, columns=['lon', 'lat'])
            famp = RegularGridInterpolator(
                points=(self.amp[self.amp.dims[0]], self.amp[self.amp.dims[1]]), 
                values=self.amp.values, method=method)
            xy_amp = xy_df.loc[:, self.amp.dims].values
            amp = famp(xy_amp)

            fpha = RegularGridInterpolator(
                points=(self.pha[self.pha.dims[0]], self.pha[self.pha.dims[1]]), 
                values=self.pha.values, method=method)
            xy_pha = xy_df.loc[:, self.pha.dims].values
            pha = fpha(xy_pha)

            amp_pha = np.vstack([amp, pha]).T
        
        # Create Points dataset
        points = xr.DataArray(data=np.arange(len(xy)), dims='points')
        lon = xr.DataArray(data=xy[:, 0], dims=('points'), coords={'points':points})
        lat = xr.DataArray(data=xy[:, 1], dims=('points'), coords={'points':points})
        amp = xr.DataArray(
            data=amp_pha[:, 0], 
            dims='points', 
            coords={'points':points})
        pha = xr.DataArray(
            data=amp_pha[:, 1], 
            dims='points', 
            coords={'points':points})
        ds = xr.Dataset({'lon':lon, 'lat':lat, 'amp':amp, 'pha':pha})
        
        return PointsDataset(dataset=ds, units=self.units)

    def _interp_ds(
            self, 
            method:Literal['complex', 'linear', 'nearest']='linear', 
            extrapolate:Literal['spherical', 'nearest']='spherical',
            **kwargs) -> GriddedDataset:
        """_summary_

        Args:
            method (Literal['complex', 'linear', 'nearest'], optional): _description_. Defaults to 'linear'.
            extrapolate (Literal['spherical', 'nearest'], optional): _description_. Defaults to 'spherical'.

        Returns:
            GriddedDataset: _description_
        """
        if method=='complex':
            # Get a consistent dataset using linear interpolation
            ds_interp = self.ds.interp(method='linear', **kwargs)
            X, Y = np.meshgrid(ds_interp.lon.values, ds_interp.lat.values)
            xy = np.vstack([X.flatten(), Y.flatten()]).T
            ds_xy = self._interp_xy(xy, method=method, extrapolate=extrapolate)
            ds_interp.amp.values = np.reshape(ds_xy.amp.values, X.shape)
            ds_interp.pha.values = np.reshape(ds_xy.pha.values, X.shape)
        else:
            ds_interp = self.ds.interp(method=method, **kwargs)

        return(GriddedDataset(dataset=ds_interp, units=self.units))
        
    def interp(
            self, 
            xy:Union[None, ArrayLike]=None,
            method:Literal['complex', 'linear', 'nearest']='linear', 
            extrapolate:Literal['spherical', 'nearest']='spherical', 
            **kwargs) -> Dataset:
        """Interpolate to a grid or a list of points

        Args:
            xy (Union[None, ArrayLike], optional): List of points. Defaults to None.
            method (Literal['complex', 'linear', 'nearest'], optional): Method to use for interpolation. 
                Use 'complex' for interpolation using coupled interpolation method of Xu et al. 
                Defaults to 'linear'.
            extrapolate (Literal['spherical', 'nearest'], optional): Controls how out of the bound interpolation is done. 
                'spherical' extrapolation goes around the grid for interpolation. Correct choice for global grid.
                'nearest' extrapolation takes the nearest value for extrapolation.
                Defaults to 'spherical'.
            **kwargs: 'lon', 'lat' values for interpolate in a grid defined by provided values. xy must be None in this case.

        Raises:
            Exception: if xy is None, **kwargs must contain one of both of `lon` and `lat`

        Returns:
            Dataset: _description_
        """
        if xy is not None:
            # we are doing an xy interpolate
            ds_interp = self._interp_xy(xy=xy, method=method, extrapolate=extrapolate)
        else:
            if len(kwargs) == 0:
                raise Exception('Dataset intepolation requires keyworded lon and/or lat values')
            else:
                ds_interp = self._interp_ds(method=method, extrapolate=extrapolate, **kwargs)

        return ds_interp
    
    # routing xarray sel function
    def _sel_ds(self, **kwargs) -> GriddedDataset:
        dataset = self.ds.sel(**kwargs)
        return GriddedDataset(dataset=dataset, units=self.units)
    
    # returning a list of values
    def _sel_xy(self, xy:np.ndarray) -> PointsDataset:
        '''
        return a list of values a the given xy points.

        xy needs to be a [2, n], or [n, 2] (preferred) [lon, lat] pairs, otherwise exception is raised.
        '''
        xy = np.atleast_2d(xy)

        points = xr.DataArray(data=np.arange(len(xy)), dims='points')
        lon = xr.DataArray(data=xy[:, 0], dims=('points'), coords={'points':points})
        lat = xr.DataArray(data=xy[:, 1], dims=('points'), coords={'points':points})
        
        amp = xr.DataArray(
            data=[ self.amp.sel(lon=x[0], lat=x[1]).values for x in xy ], 
            dims='points', 
            coords={'points':points})
        pha = xr.DataArray(
            data=[ self.pha.sel(lon=x[0], lat=x[1]).values for x in xy ], 
            dims='points', 
            coords={'points':points})
        
        ds = xr.Dataset({'lon':lon, 'lat':lat, 'amp':amp, 'pha':pha})

        return PointsDataset(dataset=ds, units=self.units)
    
    def sel(self, xy:Union[None, np.ndarray]=None, **kwargs):
        '''
        '''
        if xy is not None:
            sel = self._sel_xy(xy=xy)
        else:
            if len(kwargs) == 0:
                raise Exception('Dataset selection requires keyworded lon and/or lat values')
            else:
                sel = self._sel_ds(**kwargs)

        return sel

    # routing xarray isel function
    def _isel_ds(self, **kwargs):
        dataset = self.ds.isel(**kwargs)
        return GriddedDataset(dataset=dataset, units=self.units)
    
    # returning a list of values
    def _isel_ij(self, ij:np.ndarray):
        '''
        return a list of values a the given ij points.

        ij needs to be a [2, n], or [n, 2] (preferred) [lon, lat] pairs, otherwise exception is raised.

        TODO: Check ij shape
        '''
        ij = np.atleast_2d(ij)

        points = xr.DataArray(data=np.arange(len(ij)), dims='points')
        lon = xr.DataArray(data=self.lon[ij[:, 0]], dims=('points'), coords={'points':points})
        lat = xr.DataArray(data=self.lat[ij[:, 1]], dims=('points'), coords={'points':points})
        
        amp = xr.DataArray(
            data=[ self.amp.isel(lon=x[0], lat=x[1]).values for x in ij ], 
            dims='points', 
            coords={'points':points})
        pha = xr.DataArray(
            data=[ self.pha.isel(lon=x[0], lat=x[1]).values for x in ij ], 
            dims='points', 
            coords={'points':points})
        
        ds = xr.Dataset({'lon':lon, 'lat':lat, 'amp':amp, 'pha':pha})

        return PointsDataset(dataset=ds, units=self.units)
    
    def isel(self, ij=None, **kwargs):
        if ij is not None:
            sel = self._isel_ij(ij=ij)
        else:
            if len(kwargs) == 0:
                raise Exception('Dataset isel requires keyworded lon and/or lat indices')
            else:
                sel = self._isel_ds(**kwargs)

        return sel

    def __repr__(self):
        return(self.ds.__repr__())

    def _repr_html_(self):
        return(self.ds._repr_html_())

    # save to a file
    def to_netcdf(self, fname:PathLike, **kwargs):
        self.ds.to_netcdf(fname, **kwargs)

# load a netcdf dataset and set the longitude convention to [-180, 180] if lon180 is true
def read_gridded_dataset(
        fname:PathLike, 
        lon180:bool=False, 
        units:Union[str, dict]='auto', 
        variables:Union[str, dict]='auto'):
    """Read a gridded netcdf tidal dataset

    Args:
        fname (PathLike): File path.
        lon180 (bool, optional): Convert to longitude convention -180/180 if True. Defaults to False.
        units (Union[str, dict], optional): Dict with units. 
            Example - {'amp':'m', 'pha':'degrees'}. 
            Defaults to 'auto'.
        variables (Union[str, dict], optional): Dict with variable mapping. 
            Example - {'amp':'amplitude', 'pha':'phase', 'lon':'lon', 'lat':'lat'}
            Defaults to 'auto'.

    Raises:
        Exception: if variable={...}, it sould contain 'amp', 'pha', 'lon', 'lat'
        Exception: allowed pha units are - 'degrees', 'radians'
    """    
    ds = xr.open_dataset(fname)

    # mapping lon, lat, amp, pha
    lon_varients = ['lon', 'longitude']
    lat_varients = ['lat', 'latitude']
    amp_varients = ['amp', 'amplitude', 'A']
    pha_varients = ['pha', 'phase', 'g']
    var_fields = ['lon', 'lat', 'amp', 'pha']
    unit_fields = ['amp', 'pha']
    pha_units = ['degrees', 'radians']

    # Variable mapping
    if isinstance(variables, str):
        ds_kw = {}
        
        # find lon name, can be in both coords and data_vars
        for lon in lon_varients:
            if lon in ds.coords.keys():
                ds_kw['lon'] = lon
        
        # find lat name, can be in both coords, and data_vars
        for lat in lat_varients:
            if lat in ds.coords.keys():
                ds_kw['lat'] = lat

        # find amp name, should be in data_vars
        for amp in amp_varients:
            if amp in ds.data_vars.keys():
                ds_kw['amp'] = amp

        # find pha name, should be in data_vars
        for pha in pha_varients:
            if pha in ds.data_vars.keys():
                ds_kw['pha'] = pha
    elif isinstance(variables, dict):
        for field in var_fields:
            if field not in variables:
                raise Exception(f'Variable mapping should contain variable name for {field}')
        
        ds_kw = variables
    else:
        raise Exception(f'Variable mapping can either be auto, or a dictionary.')
    
    # Rename variables to lon, lat, amp, pha
    ds = ds.rename_vars({
            ds_kw['lon']: 'lon',
            ds_kw['lat']: 'lat',
            ds_kw['amp']: 'amp',
            ds_kw['pha']: 'pha'
        })

    # Convertion to lon180
    if lon180:
        lon = ds['lon']
        ds['_lon_180'] = (lon - 360)*(lon > 180) + (lon)*(lon <=180)
        ds = (
            ds
            .swap_dims({ds_kw['lon']:'_lon_180'})
            .sel(**{'_lon_180':sorted(ds._lon_180)})
            .drop('lon')
        ).swap_dims({
            '_lon_180': 'lon'
        })
        
    # Find units
    if isinstance(units, str):
        units = {}
        units['amp'] = ds.amp.units
        units['pha'] = ds.pha.units
    elif isinstance(units, dict):
        for unit in unit_fields:
            if unit not in units:
                raise Exception(f'Unit mapping should contain variable name for {unit}')
            
        units = dict(units)
    else:
        raise Exception(f'Unit mapping can either be auto, or a dictionary.')
    
    if units['pha'] not in pha_units:
        raise Exception(f"{units['pha']} not in allowed pha units - {pha_units}")
    
    return(GriddedDataset(dataset=ds, units=units))

def read_points_dataset(self):
    '''
    TODO: Implement a points reader, not urgent
    '''
    raise NotImplementedError()
