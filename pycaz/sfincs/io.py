#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

import pandas as pd
import xarray as xr
from pyproj import Transformer
from hydromt_sfincs import SfincsModel


def write_schism_waterlevel_forcing(ds: xr.Dataset, model: SfincsModel, out_fn: str = 'waterlevel.nc') -> None:
    """
    Write SCHISM output dataset (out2d_*.nc) as SFINCS complient boundary condition file.

    :param ds: Input SCHISM dataset, typically only selected nodes (see pycaz.schism.utils.find_nearest_node)
    :param model: An instance of SfincsModel
    :param out_fn: Output file name, default 'waterlevel.nc' saved to model.root location
    :return: None
    """
    # geographic transformation
    utm_code = model.config['epsg']
    wgs2utm = Transformer.from_crs(crs_from=4326, crs_to=utm_code, always_xy=True)

    # time selection
    ref_time = pd.to_datetime(model.config['tref'], format='%Y%m%d %H%M%S')
    ref_time_str = ref_time.strftime('%Y-%m-%d %H:%M:%S')

    # dataset
    ds_ = ds.copy()
    ds_ = ds_.rename({
        'time': 'time',
        'nSCHISM_hgrid_node': 'points',
        'elevation': 'zs'
    })

    # add x, y
    x_bnd_utm, y_bnd_utm = wgs2utm.transform(ds_.SCHISM_hgrid_node_x, ds_.SCHISM_hgrid_node_y)
    ds_ = ds_.assign(
        x=xr.DataArray(
            data=x_bnd_utm,
            dims=['points'],
            name='x',
            attrs={
                'standard_name': 'projection_x_coordinate',
                'long_name': 'x coordinate according to UTM',
                'axis': 'X',
                'units': 'm'
            }),
        y=xr.DataArray(
            data=y_bnd_utm,
            dims=['points'],
            name='y',
            attrs={
                'standard_name': 'projection_y_coordinate',
                'long_name': 'y coordinate according to UTM',
                'axis': 'Y',
                'units': 'm'
            }),
    )

    # add attr to zs
    ds_['zs'] = ds_.zs.assign_attrs(
        standard_name='water_level',
        long_name='sea_surface_height_above_mean_sea_level',
        units='m',
        coordinates='x y time',
        grid_mapping='crs'
    )

    # create crs
    ds_ = ds_.assign(crs=xr.DataArray(
        data=[1],
        dims=['crs'],
        coords={'crs': [0]},
        attrs={
            'standard_name': 'coordinate_reference_system',
            'long_name': 'coordinate_reference_system_in_UTM38S',
            'epsg_code': f'EPSG:{utm_code}'
        }
    ))

    # update encoding
    encoding = {
        'time': {
            'units': f'minutes since {ref_time_str}'
        }
    }

    # writing
    ds_.to_netcdf(Path(model.root) / out_fn, encoding=encoding)
