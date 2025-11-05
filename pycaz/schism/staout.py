# -*- encoding: utf-8 -*-
import numpy as np
import pandas as pd
import f90nml

from pycaz.typing import PathLike, TimestampConvertibleTypes


def get_start_time(fn_param: PathLike) -> pd.Timestamp:
    """
    Get start time from the param.nml file

    :param fn_param: Path to param.nml file
    :return:
    """
    param = f90nml.read(fn_param)
    start_year = param["opt"].get("start_year")
    start_month = param["opt"].get("start_month")
    start_day = param["opt"].get("start_day")
    start_hour = param["opt"].get("start_hour")
    start_time = pd.to_datetime(f"{start_year:4d}-{start_month:02d}-{start_day:02d} {start_hour:02d}:00:00")

    return start_time


def read_station_in(fn_stn_in: PathLike) -> dict:
    """
    Read station information from the station.in file

    :param fn_stn_in: Path to station.in file
    :return: {header, nstation, stations}
    """
    # Data structure
    dict_stations = dict(
        header=dict(),
        nstation=0,
        stations=list()
    )

    with open(fn_stn_in, "r") as f:
        # Parse header
        header_flags = f.readline().split("!")[0].strip().split()
        header_flags = [int(i) for i in header_flags]
        header_vars = ["elev", "psmsl", "windx", "windy", "T", "S", "u", "v", "w"]
        header = {header_var: header_flag for header_var, header_flag in zip(header_vars, header_flags)}
        dict_stations.update(header=header)

        # Number of stations
        nstation = f.readline().split("!")[0].strip()
        nstation = int(nstation)
        dict_stations.update(nstation=nstation)

        # Now will read nstation time
        stations = list()
        for _ in np.arange(nstation):
            entry = f.readline().split("!")
            [stn_num, stn_lon, stn_lat, depth] = entry[0].strip().split()
            stn_name = entry[1].strip()
            stn_dict = {
                "num": int(stn_num),
                "name": str(stn_name),
                "lon": float(stn_lon),
                "lat": float(stn_lat),
                "depth": float(depth)
            }
            stations.append(stn_dict)

        dict_stations.update(stations=stations)

    return dict_stations


def read_staout(fn_stn_out: PathLike, fn_stn_in: PathLike, fn_param: PathLike = None,
                start_time: TimestampConvertibleTypes = None) -> pd.DataFrame:
    """
    Read staout file to a dataframe. Either fn_param or start_time is required to set the timestamps

    :param fn_stn_out: Path to staout_* file
    :param fn_stn_in: Path to station.in file
    :param fn_param: Path to param.nml file
    :param start_time: Timestamp to start reading
    :return: Dataframe of station output
    """
    if fn_param is None and start_time is None:
        raise Exception("Atleast one of fn_param and start_time must be given")

    if start_time is not None:
        start_time = pd.to_datetime(start_time)

    if fn_param is not None:
        start_time = get_start_time(fn_param)

    stn_df = pd.DataFrame(read_station_in(fn_stn_in)["stations"])
    stnout_colnames = np.append(["Timestamp"], stn_df.name.values)
    df_stn_out = pd.read_csv(fn_stn_out, header=None, sep='\s+', index_col=None, names=stnout_colnames)
    df_stn_out["Timestamp"] = start_time + pd.to_timedelta(df_stn_out.Timestamp, unit='s')
    df_stn_out = df_stn_out.set_index("Timestamp")
    return df_stn_out
