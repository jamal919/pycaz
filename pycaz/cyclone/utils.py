# -*- encoding: utf-8 -*-
import numpy as np
import pandas as pd

from .track import Record, Track


def track2df(track: Track) -> pd.DataFrame:
    """Convert a track to a simple dataframe

    Args:
        track (Track): A pycaz.cyclone.track.Track object

    Returns:
        pd.DataFrame: DataFrame of the track
    """
    df = dict(
        datetime=track.timeindex,
        lon=track.lon,
        lat=track.lat,
        mslp=track.mslp,
        vmax=track.vmax,
    )
    df = pd.DataFrame(df).set_index("datetime")
    return df


def df2track(df: pd.DataFrame) -> Track:
    """Convert a simple dataframe with necessary fields to Track

    Args:
        df (pd.DataFrame): DataFrame to conver to Track

    Returns:
        Track: resulting pycaz.cyclone.track.Track object
    """
    records = list()
    for timestamp, record in df.iterrows():
        record_dict = record.to_dict()
        record_dict.update(timestamp=timestamp)
        records.append(Record(**record_dict))

    records = np.array(records)
    track = Track(records)
    return track
