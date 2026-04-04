# -*- encoding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from typing import List, Callable, Union, Optional
from pycaz.typing import PathLike, TimestampConvertibleTypes, TimedeltaConvertibleTypes

class DateOutOfBoundsError(Exception):
    """
    Raised when a requested date is outside the range of the managed files
    """
    pass

DURATION = {
    "daily":pd.DateOffset(days=1),
    "monthly":pd.DateOffset(months=1),
    "yearly":pd.DateOffset(years=1)
}


class NetCDFRegistry:
    def __init__(self,
                 fdir: PathLike,
                 date_parser: str | Callable[[PathLike], TimestampConvertibleTypes],
                 duration: Optional[str | TimedeltaConvertibleTypes] = None):
        """
        File registry to catalog and query netcdf files
        Args:
            fdir (PathLike): Folder containing the netcdf files, prefer single types of data only
            date_parser (str | Callable[[PathLike], TimestampConvertibleTypes]): A predefined function name or a function that takes file path as input and returns the datetime of the file as output
            duration (str | TimedeltaConvertibleTypes): Duration of a file. "daily" or "monthly" for typical online files, or None to infer from the first file
        """

        self.fdir = Path(fdir)
        self.date_parser = date_parser

        # Building initial file registry
        self._registry = self._build_initial_registry()

        # Figuring out duration
        if duration is not None and isinstance(duration, str):
            self.duration = DURATION.get(duration)
            if self.duration is None:
                self.duration = pd.to_timedelta(duration)

        if duration is None and not self._registry.empty:
            self.duration = self._sample_duration(Path(self._registry.iloc[0]['file_path']))

        # Apply the duration to calculate end times for all files
        self._finalize_registry()

    def _build_initial_registry(self) -> pd.DataFrame:
        """
        Fast scan: filename parsing only
        """
        records = []
        for fn in sorted(self.fdir.glob("*.nc")):
            try:
                start_time = self.date_parser(fn)
                records.append({
                    "file_path": str(fn.absolute()),
                    "start_time": pd.to_datetime(start_time)
                })
            except Exception as e:
                print(f"Skipping {fn.name}: {e}")
        return pd.DataFrame(records)

    @staticmethod
    def _sample_duration(sample_file: Path) -> TimedeltaConvertibleTypes:
        """
        The 'One-Shot' Read: Extracts duration from metadata
        """
        with xr.open_dataset(sample_file) as ds:
            # Assuming 'time' coordinate exists
            t_min = pd.to_datetime(ds.time.values.min())
            t_max = pd.to_datetime(ds.time.values.max())

            # Identify the frequency (e.g., hourly) to add one last step
            # so the end_time is the 'boundary' of the file.
            if len(ds.time) > 1:
                delta = pd.to_timedelta(np.diff(ds.time.values)[0])
                return (t_max - t_min) + delta
            return t_max - t_min

    def _finalize_registry(self):
        """
        Vectorized calculation of end times
        """
        self._registry['end_time'] = self._registry['start_time'] + self.duration - pd.to_timedelta("1s")
        self._registry = self._registry.set_index('start_time', drop=False).sort_index()

    def get_files(self, start_time: TimestampConvertibleTypes=None, end_time: TimestampConvertibleTypes=None) -> List[Path]:
        """
        Get the list of files, sliced by start and end time necessary with greedy selection

        Args:
            start_time (TimestampConvertibleTypes): Start time of files to return
            end_time (TimestampConvertibleTypes): End time of files to return

        Returns: List
        """
        if start_time is None:
            start_time = self.start_time

        if end_time is None:
            end_time = self.end_time

        s_query = pd.to_datetime(start_time)
        if s_query < self.start_time:
            raise DateOutOfBoundsError("The requested start date is out of bound")

        e_query = pd.to_datetime(end_time)
        if e_query > self.end_time:
            raise DateOutOfBoundsError("The requested end date is out of bound")

        # THE GREEDY MASK containing
        # 1. File ending after the query starts (captures the first partial file)
        # 2. File starts before our query ends (captures the last partial file)

        mask = (self._registry['end_time'] >= s_query) & (self._registry['start_time'] <= e_query)
        selected = self._registry[mask]

        if selected.empty:
            raise DateOutOfBoundsError("The requested dates are outside the range of the managed files")

        return [Path(p) for p in selected['file_path'].tolist()]

    @property
    def start_time(self) -> TimestampConvertibleTypes:
        t_start = self._registry['start_time'].min()
        return t_start

    @property
    def end_time(self) -> TimestampConvertibleTypes:
        t_end = self._registry['end_time'].max()
        return t_end


    def __repr__(self) -> str:
        """
        Returns a high-level summary of the Managed Spooled Registry.
        """
        if self._registry.empty:
            return f"<{self.__class__.__name__}(Empty - Folder: {self.fdir})>"

        n_files = len(self._registry)
        t_start = self.start_time.strftime('%Y-%m-%d %H:%M:%S')
        t_end = self.end_time.strftime('%Y-%m-%d %H:%M:%S')

        # Include the duration to show the sampling strategy used
        duration_str = str(self.duration)

        return (f"<{self.__class__.__name__}(\n"
                f"  Files: {n_files},\n"
                f"  Span:  [{t_start}] to [{t_end}],\n"
                f"  File Duration: {duration_str},\n"
                f"  Source: '{self.fdir}'\n"
                f")>")