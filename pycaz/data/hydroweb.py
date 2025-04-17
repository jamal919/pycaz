# -*- coding: utf-8 -*-

import pandas as pd
import geopandas as gpd
from shapely import Point
from pathlib import Path
from typing import Dict
from pycaz.typing import PathLike

HEADER_MAP = {
    'date': 'DATE(YYYY-MM-DD)',
    'time': 'TIME(HH:MM)',
    'waterlevel': 'ORTHOMETRIC HEIGHT (M) OF WATER SURFACE AT REFERENCE POSITION',
    'std': 'ASSOCIATED UNCERTAINTY(M)',
    'lon': 'LONGITUDE OF ALTIMETRY MEASUREMENT (deg)',
    'lat': 'LATITUDE OF ALTIMETRY MEASUREMENT (deg)',
    'height': 'ELLIPSOIDAL HEIGHT OF ALTIMETRY MEASUREMENT (M)',
    'geoid': 'GEOIDAL ONDULATION (M) at location [5,6]',
    'dist_to_ref': 'DISTANCE OF ALTIMETRY MEASUREMENT TO REFERENCE POSITION(KM)',
    'satellite': 'SATELLITE',
    'orbit': 'ORBIT / MISSION',
    'groundtrack': 'GROUND-TRACK NUMBER',
    'cycle': 'CYCLE NUMBER',
    'retracker': 'RETRACKING ALGORITHM',
    'gdr': 'GDR VERSION'
}


def read_hydroweb(fname: PathLike, header_map: Dict = None, metadata_only: bool = False):
    """
    Read water level timeseries downloaded from hydrweb.next.

    :param fname: Path of the file, typically a .txt file
    :param header_map: The mapping of the header in the data to the dataframe
    :param metadata_only: If only the metadata is to be parsed
    :return: DataFrame, Dict | Dict
    """
    if header_map is None:
        header_map = HEADER_MAP
    with open(fname) as f:
        txt = f.readlines()

    exclude = ['#PRODUCT CONTENT::', '#FIELD SEPARATOR :', '#SOURCES::']

    metadata = {}
    columns = {}
    data = []

    for line in txt:
        if line[0] == '#' and line[1] != '#' and line.strip() not in exclude:
            # header lines
            if 'COL' in line:
                # column information, splitted by ' : ', ':' will cause error with 'HH:MM'
                key, value = line.strip().replace('#', '').split(' : ')
                _, colnum = key.strip().split(' ')
                columns[int(colnum)] = value
            else:
                # general metadata, splitted by '::'
                key, value = line.strip().replace('#', '').split('::')
                metadata[key.strip()] = value.strip()

    if metadata_only:
        return metadata
    else:
        data = pd.read_csv(fname, sep=r'\s+', header=None, comment='#', na_values=[9999.999])
        data = data.drop(columns=[4])  # empty column
        data.columns = columns.keys()
        data['datetime'] = data.apply(lambda r: pd.to_datetime(f'{r[1]} {r[2]}:00'), axis=1)
        data = data.set_index('datetime')

        try:
            data.columns = list(list(header_map.keys()))
        except ValueError as e:
            raise Exception(f'{e}. Check header_map with {fname}')

        return data, metadata


def tabulate_stations(project_dir: PathLike, as_gdf: bool = True) -> pd.DataFrame | gpd.GeoDataFrame:
    """
    Create a dataframe (or geodataframe if as_gdf=True) for all the stations found in the project_dir.

    :param project_dir: Path to data downloaded from hydrweb.next
    :param as_gdf: If a geodataframe to be returned
    :return: DataFrame | GeoDataFrame
    """
    project_dir = Path(project_dir)
    stations = list(project_dir.glob('*'))
    database = {}
    for station in stations:
        fname = list(station.glob('*.txt'))[0]
        metadata = read_hydroweb(fname, header_map=HEADER_MAP, metadata_only=True)
        database[station.name] = metadata
        database[station.name].update({
            'fname': fname.resolve().as_posix(),
            'longitude': float(metadata['REFERENCE LONGITUDE']),
            'latitude': float(metadata['REFERENCE LATITUDE'])
        })

    # Create dataframe
    df = pd.DataFrame(database).T

    # Create geodataframe if gdf is True
    if as_gdf:
        geometry = [Point(x, y) for x, y in zip(df.longitude, df.latitude)]
        df = gpd.GeoDataFrame(
            data=df,
            geometry=geometry,
            crs="EPSG:4326"
        )

    return df
