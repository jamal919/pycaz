
# -*- encoding: utf-8 -*-
import numpy as np
import pandas as pd
import xarray as xr

from pycaz.typing import PathLike
from pycaz.typing import TimestampConvertibleTypes
from .bctides import Bctides


def generate_elev2d(ds: xr.Dataset, bctides: Bctides, start_time: TimestampConvertibleTypes, rnday: float,
                    bufday: float = 1, fn_out: PathLike = "elev2d.th.nc", mapping: dict = None):
    """Generate elev2d.th.nc for SCHISM from a xarray dataset

    Args:
        ds (xr.Dataset): Input xarray dataset containing lon, lat, time, elev variables. Provide mapping if the names are different.
        bctides (Bctides): pycaz Bctides object, get it from read_bctides(), or Hgrid.get_bctides()
        start_time (TimestampConvertible): Start time of the simulation
        rnday (float): Runtime of the simulation in days
        bufday (float, optional): Number of extra days at the end as buffer. Defaults to 1.
        fn_out (PathLike, optional): Path to save the generated file. Defaults to "out2d.th.nc".
        mapping (dict, optional): Mapping of variable name from ds to expected variables. Defaults to None.

    Raises:
        ValueError: If the expected variables are not found in the dataset.
        ValueError: If no iettype 4 or 5 is found in the boundary.
    """

    if mapping is None:
        mapping = dict()

    ds_mapped = ds.rename(mapping)

    var_expected = {"time", "lon", "lat", "elev"}
    var_available = set(ds_mapped.variables.keys())
    var_unavailable = list(var_expected.difference(var_available))
    if len(var_unavailable):
        raise ValueError(f"Variables: {var_unavailable} not found, provide variable mapping with `mapping` keyword")

    start_time = pd.to_datetime(start_time)
    end_time = start_time + pd.Timedelta(days=rnday) + pd.Timedelta(days=bufday)

    ds_target = ds_mapped.sel(time=ds_mapped.time.loc[start_time:end_time])
    df_time = pd.to_datetime(ds_target.time)
    da_timestep = xr.DataArray([(df_time[1] - df_time[0]).total_seconds()], dims="one")

    # Now we need to get the xy locations of the bnd
    xy = None
    for bnd in bctides.open_bnds:
        iettype = bctides.open_bnds[bnd]["iettype"]
        if iettype == 4 or iettype == 5:
            bnd_xy = bctides.open_bnds[bnd]["xy"]
            if xy is None:
                xy = bnd_xy
            else:
                xy = np.vstack([xy, bnd_xy])

    if len(xy) == 0:
        raise ValueError("No boundary found with iettype=[3,4]")

    # Creating x, y coordinate to interpolate
    da_x = xr.DataArray(data=xy[:, 0], dims="nOpenBndNodes", attrs={"units": "degrees east"})
    da_y = xr.DataArray(data=xy[:, 1], dims="nOpenBndNodes", attrs={"units": "degrees north"})

    # Interpolating and expanding dims
    da_timeseries = ds_target.interp(lon=da_x, lat=da_y)
    da_timeseries = da_timeseries.expand_dims("nLevels").expand_dims("nComponents")
    da_timeseries = da_timeseries.transpose("time", "nOpenBndNodes", "nLevels", "nComponents")

    # Output dataset
    ds_elev2d = xr.Dataset({
        "time_series": da_timeseries.elev,
        "time_step": da_timestep,
        "time": da_timeseries.time
    })

    # Encoding conforming the SCHISM manual
    encoding = {
        "time_series": {"dtype": np.float32},
        "time_step": {"dtype": np.float32},
        "lon": {"dtype": np.float32},
        "lat": {"dtype": np.float32},
        "time": {"dtype": np.float64}
    }

    ds_elev2d.to_netcdf(fn_out, encoding=encoding, unlimited_dims="time")

