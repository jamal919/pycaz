# -*- coding: utf-8 -*-
import logging
from pathlib import Path
from typing import List, Literal

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

DATA_URL = "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs"
DEFAULT_VERSION = "v04r01"
BASINS = Literal["ALL", "EP", "NA", "NI", "SA", "SI", "SP", "WP", "ACTIVE", "last3years", "since1980"]

def list_versions(data_url=DATA_URL) -> List[str]:
    """
    Returns a list of available IBTRACS versions.

    :param data_url: Source data url
    :return: List of versions
    """
    logger.info(f"Retrieving IBTRACS versions from {data_url}")
    versions = pd.read_html(data_url)[0]\
        .Name\
        .dropna()\
        .to_list()
    versions = filter(lambda x: x!= "Parent Directory", versions)
    versions = [version.replace("/", "") for version in versions]
    logger.info(f"Found {len(versions)} IBTRACS versions")
    return versions

def last_version(data_url=DATA_URL) -> str:
    """
    Returns the last version of ibtracs.

    :param data_url: Source data url
    :return: Version string, return DEFAULT_VERSION if list_version fails
    """
    versions = list_versions(data_url)
    if len(versions) == 0:
        logger.info(f"No IBTRACS versions found for {data_url}, returning default {DEFAULT_VERSION}")
        return DEFAULT_VERSION
    else:
        logger.info(f"Last version of IBTRACS is {versions[-1]}")
        return versions[-1]

def get_csv_url(data_url:str=DATA_URL, version:str=None, basin:BASINS=None) -> str:
    """
    Returns csv url for ibtracs.

    :param data_url: Source data url, defaults to DATA_URL
    :param version: Version string, defaults to DEFAULT_VERSION
    :param basin: Oceanic basin, defaults to "ALL"
    :return: Link to csv data url
    """
    if version is None:
        version = last_version(data_url)

    if basin is None:
        basin = "ALL"

    csv_url = f"{data_url}/{version}/access/csv/ibtracs.{basin}.list.{version}.csv"

    return csv_url

def get_data(url, saveto="./", fname=None, chunk_size=1024, overwrite=False) -> Path:
    """
    Downloads data from url and saves it to saveto directory

    :param url: Download url
    :param saveto: Directory to save to, defaults to current directory "./"
    :param fname: Custom name of the file to save
    :param chunk_size: chunk size of the data stream, defaults to 1024
    :param overwrite: if the data to be overwritten if exists
    """
    saveto = Path(saveto)
    if not saveto.exists():
        saveto.mkdir(parents=True, exist_ok=True)

    if fname is None:
        fname = saveto / url.split("/")[-1]

    try:
        res = requests.get(url, stream=True)
        res.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading {url}: {e}")
    else:
        total_size = int(res.headers.get("content-length", 0))

        file_exists = fname.exists() and fname.stat().st_size == total_size
        to_download = not file_exists or overwrite

        if to_download:
            with open(fname, "wb") as f, tqdm(
                    desc=fname.name,
                    total=total_size,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024) as pbar:
                for chunk in res.iter_content(chunk_size=chunk_size):
                    size = f.write(chunk)
                    pbar.update(size)

            logger.info(f"Download completed: {fname.as_posix()}")
        else:
            logger.info(f"File exists: {fname.as_posix()}")

    return fname


