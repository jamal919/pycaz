# -*- coding: utf-8 -*-

"""
Contains the class description of Atlas and relevant reader function.
"""

import os
import re
from glob import glob
from pathlib import Path
from typing import Union, List, Dict, Callable

import numpy as np
import xarray as xr

from .dataset import GriddedDataset, read_gridded_dataset
from .dataset import PointsDataset, read_points_dataset

from pycaz.typing import ArrayLike, PathLike


# A name parser for getting the first item, strip(), and return in uppercase
def generic_name_parser(fn):
    return re.findall(r"[0-9A-Za-z]+\d*", fn)[0].strip().upper()


DATASET = {
    "structured": {"class": GriddedDataset, "reader": read_gridded_dataset},
    "points": {"class": PointsDataset, "reader": read_points_dataset}
}


# An atlas which contains a list of waves from a atlas dir and provide convenient access each waves
class Atlas:
    def __init__(
            self,
            database: Dict,
            grid_type: str,
            lon180: bool,
            units: Union[str, Dict] = "auto",
            variables: Union[str, Dict] = "auto"):
        """
        A class to contain an atlas of constituents.

        :param database: Dictionary with constituents as key. The content should contain a `dataset.Dataset` object,
        if not a dataset object then it should contain a `fn` pointing to the dataset.
        :param grid_type: Type of the data grid.
            - "structured": contains a lon, lat dataset with amp(lat, lon), and pha(lat, lon)
            - "points": contains a points dataset with amp(points), and pha(points)
        :param lon180: preprocess if needed to load the dataset in -180/180 format
        :param units: dict of {"amp":"m", "pha":"degrees"}. Defaults to "auto".
        :param variables: dict of {"amp":"amplitude", "pha":"phase", "lon":"lon", "lat":"lat"}. Defaults to "auto".
        """
        self.database = database
        self.grid_type = grid_type

        if self.grid_type not in DATASET:
            raise NotImplementedError(
                f"Dataset type {self.grid_type} not implemented. Available types - {list(DATASET.keys())}")
        else:
            self.dataset = DATASET[self.grid_type]

        self.lon180 = lon180
        self.units = units
        self.variables = variables

        # Update units from the first wave, requires reading atleast one file
        self._update_units()

    @property
    def waves(self) -> list:
        """
        List of waves in the Atlas dataset

        :return: List of waves in the Atlas dataset
        """
        return [wave for wave in self.database]

    @property
    def lon(self) -> xr.DataArray:
        """
        Returns the longitude DataArray from the first dataset

        :return: Longitude data array
        """
        return self[self.waves[0]].lon

    @property
    def lat(self) -> xr.DataArray:
        """
        Returns the longitude DataArray from the first dataset

        :return: Latitude data array
        """
        return self[self.waves[0]].lat

    def _update_units(self) -> None:
        """
        Update the unit from the first dataset in the database
        """
        self.units = self[self.waves[0]].units
        print(f"units are set to {self.units}")

    def interp(self, xy: Union[ArrayLike, None] = None, method="linear", **kwargs) -> "Atlas":
        """

        :param xy: xy points to be interpolated. If not provided, interpolates all points.
        :param method: Interpolation method among `linear`, `nearest`, `complex`. Defaults to `linear`.
        :param kwargs: Keyword arguments passed to interpolation methods.
        :return: Interpolated Atlas object
        """
        database = {}
        for wave in self.waves:
            database[wave] = {"dataset": self[wave].interp(xy=xy, method=method, **kwargs)}

        if xy is not None:
            grid_type = "points"
        else:
            grid_type = self.grid_type

        return Atlas(
            database=database,
            grid_type=grid_type,
            lon180=self.lon180,
            units=self.units,
            variables=self.variables)

    def select(self, consts: Union[List, Dict]):
        """
        Selects a list of constituents from the atlas from the consts list

        :param consts: List of constituents to be selected
        :return: An Atlas with the selected consts
        """
        database = {}
        if isinstance(consts, list):
            for const in consts:
                database[const] = self.database[const]
        elif isinstance(consts, dict):
            for const in consts:
                database[consts[const]] = self.database[const]

        return Atlas(
            database=database,
            grid_type=self.grid_type,
            lon180=self.lon180,
            units=self.units,
            variables=self.variables)

    def __contains__(self, wave):
        if wave in self.waves:
            return True

        return False

    def __iter__(self):
        self.n = 0

        return self

    def __next__(self):
        try:
            next_item = self.database[self.waves[self.n]]
        except IndexError:
            raise StopIteration

        self.n += 1

        return next_item

    def __getitem__(self, wave):
        """
        Gives access to the Dataset for the given name of the wave

        :param wave: Wave name
        :return:
        """
        try:
            assert wave in self.waves
        except AssertionError:
            raise Exception("{wave} wave us not found in atlas {self.atlas_dir}")

        try:
            ds = self.database[wave]["dataset"]
        except KeyError:
            fn = self.database[wave]["fn"]
            self.database[wave]["dataset"] = self.dataset["reader"](
                fn,
                self.lon180,
                units=self.units,
                variables=self.variables)
            ds = self.database[wave]["dataset"]

        return ds

    def values(self, consts: Union[list, str] = "auto"):
        """
        Return a stacked dataset of the amplitudes for the given set of constituents

        :param consts: List of constituents to be selected
        :return: (Amplitude, Phase)
        """
        if consts == "auto":
            consts = self.waves
        else:
            consts = consts

        if self.dataset["class"] is GriddedDataset:
            amp = np.stack([np.atleast_2d(self[const].amp.values) for const in consts])
            pha = np.stack([np.atleast_2d(self[const].pha.values) for const in consts])
        elif self.dataset["class"] is PointsDataset:
            amp = np.atleast_3d(np.stack([np.atleast_1d(self[const].amp.values) for const in consts]))
            pha = np.atleast_3d(np.stack([np.atleast_1d(self[const].pha.values) for const in consts]))
        else:
            amp, pha = np.atleast_3d([]), np.atleast_3d([])

        return amp, pha

    def to_netcdf(self, fdir: PathLike, suffix="", prefix=""):
        """
        Save the current atlas to a directory dir with filenames formatted as `suffix_wave_prefix.nc`

        :param fdir: Directory where the atlas to be saved
        :param suffix: Suffix to the wave name
        :param prefix: Prefix to the wave name
        :return: None
        """
        fdir = Path(fdir)
        if not fdir.exists():
            fdir.mkdir(parents=True)

        for wave in self.waves:
            fn = fdir / f"{suffix}{wave}{prefix}.nc"
            self[wave].to_netcdf(fn)


def read_atlas(
        atlas_dir: PathLike,
        name_parser: Callable = generic_name_parser,
        ext: str = ".nc",
        grid_type: str = "structured",
        lon180: bool = True,
        units: Union[str, dict] = "auto",
        variables: Union[str, dict] = "auto") -> Atlas:
    """
    Read a tidal atlas from `atlas_dir`

    :param atlas_dir: Directory containing the tidal atlas
    :param name_parser: Name parser function which extract the constituent name. Defaults to `generic_name_parser`
    :param ext: Extensions of the file. Defaults to ".nc"
    :param grid_type: Name of the grid type, "structured", "points". Defaults to "structured"
    :param lon180: If lon to be changed from 0/360 to -180/180. Defaults to True.
    :param units: dict of {"amp":"m", "pha":"degrees"}. Defaults to "auto".
    :param variables: dict of {"amp":"amplitude", "pha":"phase", "lon":"lon", "lat":"lat"}. Defaults to "auto".

    :return: Atlas
    """
    database = {}
    file_paths = glob(os.path.join(atlas_dir, f"*{ext}"))
    for file_path in file_paths:
        wave = name_parser(os.path.basename(file_path))
        database[wave] = dict(fn=file_path)

    return Atlas(database=database, grid_type=grid_type, lon180=lon180, units=units, variables=variables)
