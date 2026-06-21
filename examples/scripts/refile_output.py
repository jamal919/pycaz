#!/usr/bin/env python3

import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
import shutil
import sys
import logging

logging.basicConfig(
        format="[%(asctime)s | %(levelname)s] %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S", 
        stream=sys.stdout, 
        level=logging.INFO)

logging.info("Import successful")

def map_to_refile(fdir_in, fdir_out, skip_fn=0, prefix="out2d"):
    logging.info(f"Mapping outputs with prefix={prefix} in {fdir_in}")
    fdir_outputs = Path(fdir_in) / "outputs"
    fns = list(fdir_outputs.glob(f"{prefix}_*.nc"))
    n_fns = len(fns)
    logging.info(f"Found {n_fns} files")
    fns_num = [int(fn.name.split(".")[0].split("_")[1]) for fn in fns]
    fns_num.sort()
    fns_refile_name = list()
    logging.info(f"Starting file mapping to refile name skipping {skip_fn} file")
    for ifn in fns_num[skip_fn:]:
        fn_in = fdir_outputs / f"{prefix}_{ifn}.nc"
        with xr.open_dataset(fn_in) as ds:
            timestamp = pd.to_datetime(ds.time.values[0])
            timestr = timestamp.strftime("%Y%m%d%H")
            fn_out = fdir_out / f"{prefix}_{timestr}.nc"
            fns_refile_name.append((fn_in, fn_out))
            logging.info(f"Mapped: {fn_in} -> {fn_out}")

    return fns_refile_name

def check_run_complete(fdir_in):
    status = False
    fn = Path(fdir_in) / "outputs" / "mirror.out"
    logging.info(f"Checking {fdir_in} for outputs/mirror.out") 
    if not fn.exists():
        logging.warning(f"Folder {fdir_in} does not exist... not yet simulated?")
        return status

    with open(fn, "r") as f:
        statusline = f.readlines()[-1]
        status = statusline.split()[2] == "successfully"

    if status:
        logging.info(f"Run completed in {fdir_in}")
    else:
        logging.warning(f"Simulation is still running or got a crash in {fdir_in}")
    
    return status

def apply_rename(fn_in, fn_out):
    logging.info(f"Moving {fn_in} -> {fn_out}")
    if fn_out.exists():
        logging.warning(f"File exists: {fn_out}")
        raise ValueError(f"File exists: {fn_out}")
    shutil.move(fn_in, fn_out)
    logging.info(f"Moved {fn_in} -> {fn_out}")
    return fn_out


if __name__=="__main__":
    fdir_simulation = Path("/users/p19002/khan/scratch/bengal_era5_hindcast/simulations/")
    fdir_years = list(fdir_simulation.glob("*/"))
    fdir_years_completed = list(filter(check_run_complete, fdir_years))
    skip_fn=1
    for fdir_in in fdir_years_completed:
        fdir_out = Path("./")
        fns_map = map_to_refile(fdir_in=fdir_in, fdir_out=fdir_out, skip_fn=1)
        fns_renamed = [apply_rename(*fn_map) for fn_map in fns_map]
    logging.info("Script completed successfully")

