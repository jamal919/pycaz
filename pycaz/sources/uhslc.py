# -*- coding: utf-8 -*-

import pandas as pd
import requests
import logging
from pathlib import Path

def format_stn_name(name:str) -> str:
    """Format the station name to remove space

    Args:
        name (str): name string

    Returns:
        str: formatted name string
    """
    name = str(name)
    name = name.replace(" ", "_")
    name = name.replace("-", "_")
    return name

def query_stations(delivery="research", frequency="hourly") -> pd.DataFrame:
    """Query the UHSLC station list

    Args:
        delivery (str, optional): "research" or "fast". Defaults to "research".
        frequency (str, optional): "hourly" or "daily". Defaults to "hourly".

    Returns:
        DataFrame: Station information in a DataFrame
    """
    if delivery == "fast":
        ds_id = f"global_{frequency}_fast"
        erddap_url = (
            f"https://uhslc.soest.hawaii.edu/erddap/tabledap/{ds_id}.csv?"
            "station_name,station_country,uhslc_id,latitude,longitude"
            "&distinct()"
        )
    else:
        ds_id = f"global_{frequency}_rqds"
        erddap_url = (
            f"https://uhslc.soest.hawaii.edu/erddap/tabledap/{ds_id}.csv?"
            "station_name,station_country,version,uhslc_id,latitude,longitude"
            "&distinct()"
        )
    
    
    df = pd.read_csv(erddap_url)
    df = df.drop(0).reindex()
    df["uhslc_id"] = df["uhslc_id"].astype(int)
    return df

def download_erddap(station_name, uhslc_id, version=None, delivery="research", frequency="hourly", outdir=None):
    """Download UHSLC data from erddap query

    Args:
        station_name (str): Station name
        uhslc_id (int): UHSLC id
        version (str, optional): Version of the data 'A', 'B', 'C' etc. Defaults to None.
        delivery (str, optional): "research" or "fast" delivery. Defaults to "research".
        frequency (str, optional): "hourly" or "daily". Defaults to "hourly".
        outdir (str, optional): Path to output directory. Defaults to None.
    """
    if outdir is None:
        outdir = "./"
    outdir = Path("./")
    if not outdir.exists():
        outdir.mkdir(parents=True)
    
    if delivery == "fast":
        ds_id = f"global_{frequency}_fast"
        url = (
            'https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_rqds.csv?' 
            'time,sea_level' 
            '&station_name="{station_name}"' 
            '&uhslc_id={uhslc_id}')
    else:
        ds_id = f"global_{frequency}_rqds"
        url = (
            'https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_rqds.csv?' 
            'time,sea_level' 
            '&station_name="{station_name}"' 
            '&uhslc_id={uhslc_id}' 
            '&version="{version}"')

    url = url.format(station_name=station_name, uhslc_id=uhslc_id, version=version)
    print(url)
    r = requests.get(url)
    r.raise_for_status()
    fn = outdir / f"{format_stn_name(station_name)}.csv"
    with open(fn, "w") as f:
        f.writelines(r.text)

def read_webdata(fn) -> pd.DataFrame:
    """Read UHSLC data from a file and return a pandas dataframe.

    This function is to be used for the data downloaded directly from the website.

    Args:
        fn (str): filepath to read in.

    Returns:
        pd.DataFrame: Data in a pandas DataFrame
    """
    df = pd.read_csv(
        fn, 
        header=None, 
        names=["Year", "Month", "Day", "Hour", "Elev"], 
        na_values=[-9999, 9999, -32767]
        )
    datetime = pd.to_datetime(
        df.apply(
            lambda r: f"{int(r.Year):04d}-{int(r.Month):02d}-{int(r.Day):02d} {int(r.Hour):02d}:00:00", 
            axis=1)
        )
    df = df.loc[:, ["Elev"]].assign(Datetime=datetime).set_index("Datetime")
    df = df/1000 # mm -> m
    return df

def read_erddap(fn) -> pd.DataFrame:
    """Read UHSLC data from a file downloaded with uhslc_download_erddap

    Args:
        fn (str): filepath to read in.

    Returns:
        pd.DataFrame: Data in a pandas DataFrame
    """
    df = pd.read_csv(fn, skiprows=2, parse_dates=[0], header=None, names=["Datetime","Elev"])
    df["Datetime"] = df.Datetime.dt.tz_convert(None)
    df = df.set_index("Datetime")
    df["Elev"]  = df["Elev"]/1000 # mm -> m
    return df

def read(fn, source="erddap") -> pd.DataFrame:
    """Read UHSLC data

    Args:
        fn (str): file path to read
        source (str, optional): One of "erddap" or "webdata". Defaults to "erddap".

    Returns:
        pd.DataFrame: Timeseries read from the file
    """
    fn = Path(fn)
    if source=="erddap":
        read_erddap(fn)
    elif source=="webdata":
        read_webdata(fn)
    else:
        raise(f"Source must be one of 'erddap' or 'webdata'")
