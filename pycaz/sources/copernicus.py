#!/usr/bin/env python

import json
import logging
import shutil
import time
import zipfile
from pathlib import Path

import cdsapi
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from pycaz.typing import TimestampConvertibleTypes, ArrayLike, PathLike
from tqdm.autonotebook import tqdm

VAR_SETS = {
    "reanalysis-era5-single-levels": {
        "full": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_sea_level_pressure",
            "2m_dewpoint_temperature",
            "2m_temperature",
            "mean_evaporation_rate",
            "mean_total_precipitation_rate",
            "total_precipitation",
            "mean_surface_downward_long_wave_radiation_flux",
            "mean_surface_downward_short_wave_radiation_flux",
        ],
        "essential": ["10m_u_component_of_wind", "10m_v_component_of_wind", "mean_sea_level_pressure"],
    }
}


class ERA5ReanalysisSingleLevel:
    def __init__(self, prefix: str = "era5", outdir: PathLike = None):
        self.prefix = prefix
        if outdir is None:
            outdir = "./"
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir()
        self.dataset = "reanalysis-era5-single-levels"
        self.request = {
            "product_type": ["reanalysis"],
            "variable": [],
            "year": [],
            "month": [],
            "day": [],
            "time": [],
            "data_format": "netcdf",
            "download_format": "zip",
            "area": [90, -180, -90, 180],
        }
        self.requests = list()
        self.client = cdsapi.Client()
        self.logger = self.client.logger

        logging.basicConfig(filename=self.outdir.joinpath("log.txt"), level=logging.INFO)

    def generate_requests(self,
                          start: TimestampConvertibleTypes,
                          end: TimestampConvertibleTypes,
                          extent: ArrayLike = None,
                          chunk_by: str = None,
                          dt: int = 1,
                          var_set: str | list = "essential"):
        requests = list()
        start = pd.to_datetime(start)
        end = pd.to_datetime(end)
        if extent is None:
            extent = [-180, 180, -90, 90]

        w, e, s, n = extent
        area = [n, w, s, e]

        if chunk_by is None:
            chunk_by = "none"
        chunk_types = ["none", "year", "month"]
        if chunk_by not in chunk_types:
            raise Exception(f"Invalid chunking {chunk_by}, valid values {chunk_types}")

        dt = int(dt)

        var_sets = VAR_SETS[self.dataset]
        if var_set in var_sets:
            variables = var_sets[var_set]
        else:
            variables = var_set

        years = [f"{y:4d}" for y in range(start.year, end.year + 1, 1)]
        months = [f"{m:02d}" for m in pd.date_range(start, end, freq="ME").month.unique().values]
        days = [f"{d:02d}" for d in range(1, 32, 1)]
        hours = [f"{h:02d}" for h in range(0, 24, dt)]

        request = self.request.copy()
        request.update(
            {"variable": variables, "year": years, "month": months, "day": days, "time": hours, "area": area})

        if chunk_by == "none":
            fn = f"{self.prefix}.zip"
            requests.append((self.dataset, request.copy(), fn))

        if chunk_by == "year":
            for year in years:
                fn = f"{self.prefix}_{year}.zip"
                request.update({"year": [year]})
                requests.append((self.dataset, request.copy(), fn))

        if chunk_by == "month":
            for year in years:
                for month in months:
                    fn = f"{self.prefix}_{year}_{month}.zip"
                    request.update({"year": [year], "month": [month]})
                    requests.append((self.dataset, request.copy(), fn))

        self.requests = requests

    def save_requests(self, fname: PathLike, summary: bool = False):
        with open(fname, 'w') as f:
            res_str = json.dumps(self.requests, indent=2)
            f.write(res_str)

    def load_requests(self, fname: PathLike):
        with open(fname, 'r') as f:
            res_data = json.load(f)

        self.requests = res_data

    def download(self):
        for dataset, request, fn in tqdm(self.requests):
            fpath = self.outdir / fn
            fpath = fpath.absolute().as_posix()

            download(dataset=dataset, request=request, target=fpath, client=self.client)
            self.logger.info(f"Downloaded {fn} to {self.outdir}")

    def extract(self, delete: bool = False):
        fns = self.outdir.glob("*.zip")
        rename = lambda fn: fn.name.replace(".zip", ".nc")

        for fn in tqdm(fns):
            dir_tmp = Path(f"./tmp_{str(hash(fn))[-4:]}")
            fn_out = fn.parent / rename(fn)
            with zipfile.ZipFile(fn) as fnz:
                fnz.extractall(dir_tmp)
            merge(dir_in=dir_tmp, fn_out=fn_out)
            shutil.rmtree(dir_tmp)

            if delete:
                Path(fn).unlink()

    def filter_unavailable(self):
        fns = self.outdir.glob("*.nc")
        rename = lambda fn: fn.name.replace(".nc", ".zip")
        fns = [rename(fn) for fn in fns]

        requests = list()
        for dataset, request, fn in self.requests:
            if fn not in fns:
                requests.append((dataset, request, fn))

        self.requests = requests


def download(dataset: str, request: dict, target: PathLike, client: cdsapi.Client = None):
    if client is None:
        client = cdsapi.Client()
    logger = client.logger

    func_download = lambda: client.retrieve(dataset, request, target)
    retry(func_download, logger=logger)


def retry(func: callable, retries: int = 3, delay: int = 1, exceptions: Exception = (Exception,),
          logger: logging.Logger = None):
    if logger is None:
        logger = logging.getLogger("retry")
    else:
        logger = logger

    for attempt in range(1, retries + 1):
        try:
            return func()
        except exceptions as e:
            logger.info(f"Attempt {attempt} failed: {e}")
            if attempt == retries:
                raise
            time.sleep(delay)


def merge(dir_in:PathLike, fn_out:PathLike):
    with xr.open_mfdataset(dir_in.glob("*.nc")) as ds:
        ds = (
            ds.rename_dims({"valid_time": "time"})
            .rename_vars({"valid_time": "time"})
            .drop_vars(["number", "expver"])
        )

        with ProgressBar():
            ds.to_netcdf(fn_out)
