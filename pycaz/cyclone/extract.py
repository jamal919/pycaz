# -*- encoding: utf-8 -*-
import logging

import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
from scipy.ndimage import center_of_mass
from shapely import Point

from .track import Record, Track


def normalize_minmax(arr: np.ndarray|xr.DataArray) -> np.ndarray|xr.DataArray:
    """Normalize an array using its min max range

    Args:
        arr (np.ndarray): Array to normalize

    Returns:
        np.ndarray: Normalized array [0, 1]
    """
    arr_normalized = (arr - np.nanmin(arr)) / (np.nanmax(arr) - np.nanmin(arr))
    return arr_normalized


def extract_center_min_prmsl(prmsl: xr.DataArray) -> tuple[float, float]:
    """Extract center location using min prmsl value

    Args:
        prmsl (xr.DataArray): Mean sea level pressure data array

    Returns:
        tuple[float, float]: Longitude, Latitude of min prmsl value
    """
    id_ravel = np.nanargmin(prmsl.values)
    idy, idx = np.unravel_index(id_ravel, prmsl.values.shape)
    lon = prmsl.lon[idx]
    lat = prmsl.lat[idy]
    return lon, lat


def compute_relative_vorticity(u10: xr.DataArray, v10: xr.DataArray) -> xr.DataArray:
    """Compute relative vorticity

    Relative vorticity is positive if anti-clockwise, and negative if clockwise

    Args:
        u10 (xr.DataArray): west-east, zonal, wind component
        v10 (xr.DataArray): south-north, meridional, wind component

    Returns:
        xr.DataArray: relative vorticity
    """
    dv_dx = v10.differentiate("lon")
    du_dy = u10.differentiate("lat")

    rvor = dv_dx - du_dy
    return rvor


def extract_center_hybrid(
        u10: xr.DataArray,
        v10: xr.DataArray,
        prmsl: xr.DataArray,
        threshold: float = 0.8,
) -> tuple[float, float]:
    """Extract the storm center from hybrid relative vorticity and pressure field

    Args:
        u10 (xr.DataArray): west-east, zonal, wind component
        v10 (xr.DataArray): south-north, meridional, wind componet
        prmsl (xr.DataArray): mean sea level pressure component
        threshold (float, optional): Threshold to select data in the hybrid field. Defaults to 0.8.

    Returns:
        tuple[float, float]: _description_
    """
    rvor = np.abs(compute_relative_vorticity(u10=u10, v10=v10))

    msl_norm = 1.0 - normalize_minmax(prmsl)
    vor_norm = normalize_minmax(rvor)

    hyb_norm = msl_norm * vor_norm

    # Find the location of the maximum 'Hybrid' signal
    y_idx, x_idx = center_of_mass(hyb_norm > threshold * np.max(hyb_norm))
    lon = np.interp(x_idx, np.arange(len(hyb_norm.lon)), hyb_norm.lon)
    lat = np.interp(y_idx, np.arange(len(hyb_norm.lat)), hyb_norm.lat)
    return lon, lat


def model_extract_track(
        model_ds: xr.Dataset,
        track: Track = None,
        search_radius: float = 2,
        mapping: dict = None,
) -> Track:
    """Extract track from modeled data

    Args:
        model_ds (xr.Dataset): Dataset from which the track is extracted, requires time, lon, lat, prmsl, u10, v10. Can be mapped to other names using mapping.
        track (pycaz.cyclone.track.Track, optional): JTWC or similar track to constrain the extraction. Defaults to None.
        search_radius (float): The search radius to take if track is present. 2 degree ~ 200 Km is generally sufficient.
        mapping (dict): Variable mapping var_ds:var_target, defaults to {"time":"time", "longitude":"lon", "latitude":"lat", "u10":"u10", "v10":"v10", "prmsl":"prmsl"}

    Returns:
        Track: Extracted track from the dataset
    """
    if mapping is None:
        mapping = {
            "time": "time",
            "lon": "lon",
            "lat": "lat",
            "u10": "u10",
            "v10": "v10",
            "prmsl": "prmsl",
        }

    ds = model_ds.rename(mapping)

    # Selecting timerange
    timestamps = ds.time
    logging.info(f"Model data: From {timestamps[0].values} to {timestamps[-1].values}")
    if track:
        logging.info(f"Track data: From {track.timeindex[0]} to {track.timeindex[-1]}")
        timeindex_min = np.max([timestamps[0], track.timeindex[0]])
        timeindex_max = np.min([timestamps[-1], track.timeindex[-1]])
        timestamps = timestamps.loc[timeindex_min:timeindex_max].values

    logging.info(f"Extracting: From {timestamps[0]} to {timestamps[-1]}")

    # Creating track
    records = list()

    for timestamp in timestamps:
        timestamp = pd.to_datetime(timestamp)
        logging.info(f"Processing {timestamp}")
        prmsl = ds.prmsl.sel(time=timestamp).compute()
        u10 = ds.u10.sel(time=timestamp).compute()
        v10 = ds.v10.sel(time=timestamp).compute()

        if track:
            track_record = track.interpolate(at=timestamp)
            geom = [
                Point(track_record["lon"], track_record["lat"]).buffer(search_radius)
            ]
            try:
                prmsl = prmsl.rio.write_crs(4326).rio.clip(geom)
                u10 = u10.rio.write_crs(4326).rio.clip(geom)
                v10 = v10.rio.write_crs(4326).rio.clip(geom)
            except rioxarray.exceptions.NoDataInBounds:
                logging.info(
                    f"No data found within the cyclone footprint at {timestamp}"
                )
                continue

        vcirc = np.sqrt(u10 ** 2 + v10 ** 2)

        center_lon, center_lat = extract_center_hybrid(
            u10=u10, v10=v10, prmsl=prmsl, threshold=0.8
        )

        mslp = np.nanmin(prmsl.values)
        vmax = np.nanmax(vcirc)

        records.append(
            Record(
                timestamp=timestamp,
                lon=center_lon,
                lat=center_lat,
                mslp=mslp,
                vmax=vmax,
            )
        )

    logging.info(f"Track extraction is complete")
    track = Track(np.array(records))

    return track
