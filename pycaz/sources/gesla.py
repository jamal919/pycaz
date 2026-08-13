# -*- coding: utf-8 -*-

import pandas as pd
import requests
import logging
from pathlib import Path

def query_stations() -> pd.DataFrame:
    """Query available stations in GESLA dataset

    Returns:
        DataFrame: Dataframe containing the GESLA station information
    """
    erddap_url = (
        "https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_gesla.csv?"
        "latitude,longitude,record_id,agency_id,station_name,station_code,station_country_code"
        "&distinct()"
    )

    df = pd.read_csv(erddap_url)
    df = df.drop(0).reindex()
    return df

def download_erddap(record_id:str, outdir=None) -> None:
    """Download GESLA timeseries given a record_id.
    
    The output file is named as record_id.csv and saved to outdir.

    Args:
        record_id (str): record_id from the gesla_query_stations
        outdir (path, optional): output directory for saving the file. Defaults to None.
    """
    if outdir is None:
        outdir = "./"
    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir(parents=True)
    url = (
        'https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_gesla.csv?' 
        'time,sea_level,flag1,flag2' 
        '&record_id="{record_id}"')
    url = url.format(record_id=record_id)
    r = requests.get(url)
    r.raise_for_status()
    fn = outdir / f"{record_id}.csv"
    with open(fn, "w") as f:
        f.writelines(r.text)

def read_erddap(fn) -> pd.DataFrame:
    """ Read GESLA data from a file downloaded with gesla_download_erddap

    Args:
        fn (str): filepath to read

    Returns:
        pd.DataFrame: Data in a pandas DataFrame
    """
    df = pd.read_csv(fn, skiprows=2, parse_dates=[0], header=None, names=["Datetime","Elev", "Flag1", "Flag2"])
    df["Datetime"] = df.Datetime.dt.tz_convert(None)
    df = df.set_index("Datetime")
    return df

def read(fn, source="erddap") -> pd.DataFrame:
    """Read GELSA data

    Args:
        fn (str): file path to read
        source (str, optional): Defaults to "erddap".

    Returns:
        pd.DataFrame: Timeseries read from the file
    """
    fn = Path(fn)
    if source=="erddap":
        read_erddap(fn)
    else:
        raise(f"Currently only 'erddap' source is available for GESLA")