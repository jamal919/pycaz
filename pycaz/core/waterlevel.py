# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import dataclasses
import utide
import logging

@dataclasses.dataclass(eq=False)
class WaterLevelSeries:
    timeseries: pd.Series = dataclasses.field(repr=False)
    name: str = dataclasses.field(default="Unassigned")
    source: str = dataclasses.field(default="Unassigned")
    lon: float = dataclasses.field(default="Unassigned")
    lat: float = dataclasses.field(default="Unassigned")
    datum: str = dataclasses.field(default="Unassigned")

    def __post_init__(self):

        # Convert DataFrame to Series
        if isinstance(self.timeseries, pd.DataFrame):
            try:
                self.timeseries = pd.Series(self.timeseries.loc[:, 0])
                logging.warning(f"Taking the first column of the DataFrame")
            except Exception as e:
                logging.fatal(f"DataFrame could not be converted to Series")
                raise e

        # Check if the timeseries is a Series
        try:
            assert isinstance(self.timeseries, pd.Series)
        except Exception as e:
            logging.fatal(f"Only DataFrame/Series are accepted")
            raise e

        # Check if the index is DateTimeIndex
        try:
            assert isinstance(self.timeseries.index, pd.DatetimeIndex)
        except Exception as e:
            logging.fatal(f"Index is not DatetimeIndex")
            raise e

    def subset(self, start=None, end=None):
        start = pd.to_datetime(start)
        end = pd.to_datetime(end)
        timeseries = self.timeseries.loc[start:end]
        return dataclasses.replace(self, timeseries=timeseries)

    def interpolate(self, at):
        xp = pd.to_datetime(at)
        yp = np.interp(xp, self.timeseries.index, self.timeseries)
        series = pd.Series(yp, index=xp)

        obj = WaterLevelSeries(
            timeseries=series,
            name=self.name,
            source=f"{self.source}-Interpolated",
            lon=self.lon,
            lat=self.lat,
            datum=self.datum)

        return obj

    def utide_harmonics(self, **kwargs):
        coef = utide.solve(t=self.timeseries.index, u=self.timeseries, lat=self.lat, **kwargs)
        return coef

    def to_tuple(self):
        return dataclasses.astuple(self)

    def to_dict(self):
        return dataclasses.asdict(self)

    @property
    def stn_info(self):
        obj = self.to_dict()
        info_keys = ["name", "source", "lon", "lat", "datum"]
        info_dict = {key: obj["key"] for key in info_keys}
        return info_dict

    def plot(self, **kwargs):
        ax = self.timeseries.plot(**kwargs)
        return ax