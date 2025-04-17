# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from io import StringIO
from pathlib import Path
import copy

from pycaz.wave.spectra import compute_bulk_params, correct_shoaling


def read_sp1d(fn: str | Path) -> dict:
    with open(fn, 'r') as f:
        data = f.readlines()

    # Header section
    _msc = int(data[0])  # spectral bins
    _mdc = int(data[1])  # directional bins

    ln = 2  # pointer/cursor to line

    # TODO: can be replaced with a generalized function
    _omega = np.full(_msc, fill_value=np.nan)
    m = 0  # data count cursor
    while m < _msc:
        tmps_ = np.genfromtxt(StringIO(data[ln]))
        for v in tmps_:
            _omega[m] = v
            m = m + 1
        ln = ln + 1

    _dir = np.full(_mdc, fill_value=np.nan)
    m = 0
    while m < _mdc:
        tmps_ = np.genfromtxt(StringIO(data[ln]))
        for v in tmps_:
            _dir[m] = v
            m = m + 1
        ln = ln + 1

    _isum = int(data[ln])
    ln = ln + 1

    _header = {
        'msc': _msc,
        'mdc': _mdc,
        'omega': _omega,
        'dir': _dir
    }

    # Data section
    # how many lines of data
    # TODO:potential bug due to unprinted wkloc_stations[i, :] in wwm_output.f90
    dl = _msc + 3

    _spec = {}  # spec data container
    i = 0
    while True:
        try:
            # Data header section
            date_ = pd.to_datetime(data[ln + dl * i + 0].strip(), format='%Y%m%d.%H%M%S')
            dep_ = float(data[ln + dl * i + 1])
            u_, v_ = np.genfromtxt(StringIO(data[ln + dl * i + 2]))

            # Data section for date
            spectbl_ = np.genfromtxt(data[ln + dl * i + 3: ln + dl * (i + 1)])
            spectbl_[spectbl_ == -999] = np.nan  # set nan value
            specdf_ = pd.DataFrame(
                spectbl_,
                columns=['f', 'efth', 'dir', 'spr'],
            )
            _spec[date_] = {
                'depth': dep_,
                'u': u_,
                'v': v_,
                'spec': specdf_
            }
            i = i + 1
        except:
            break

    return {
        'header': _header,
        'data': _spec
    }


def correct_sp1d_shoaling(depth):
    """
    Correct for the shoaling for `depth` from the sp1d file

    :param depth: Depth where the spectra needs to be converted
    :return: Corrected sp1d
    """
    pass


def compute_sp1d_bulk_params(sp1d: dict, freq_low=None, freq_high=None) -> pd.DataFrame:
    """
    Compute bulk parameters from a sp1d file read from WWM station output.

    :param sp1d: 1D spectra from WWM station output read by read_sp1d()
    :param freq_low: Low frequency cutoff, consistent with buoy.
    :param freq_high: High frequency cutoff, consistent with buoy.
    :return: DataFrame of bulk parameters computed by pycaz.wave.spectra.compute_bulk_params()
    """
    if sp1d.get('data') is None:
        raise ValueError('sp1d must have data. Use sp1d to read wwm sp1d output.')

    bulk_params = {}
    for timestamp in sp1d['data']:
        _freq = sp1d['data'][timestamp]['spec']['f'].values
        _efth = sp1d['data'][timestamp]['spec']['efth'].values

        bulk_params[timestamp] = compute_bulk_params(_freq, _efth, freq_low, freq_high)

    bulk_params = pd.DataFrame(bulk_params).T

    return bulk_params


def correct_shoaling_sp1d(sp1d: dict, target_depth: float) -> dict:
    """
    Correct the spectrum in sp1d for shoaling in target_depth.

    Currently the velocity in return sp1d dictionary is set to np.nan.
    TODO: Reconsider which velocity to keep.

    :param sp1d: 1D spectra from WWM station output read by read_sp1d().
    :param target_depth: Depth to which the spectrum is corrected. Typically buoy depth.
    :return: Corrected sp1d dictionary.
    """
    sp1d_corrected = copy.deepcopy(sp1d)
    for timestamp in sp1d['data']:
        spec = sp1d['data'][timestamp]['spec']
        depth = sp1d['data'][timestamp]['depth']
        efth_corrected = correct_shoaling(
            freq=spec.loc[:, 'f'],
            efth=spec.loc[:, 'efth'],
            depth=depth,
            target_depth=target_depth
        )

        sp1d_corrected['data'][timestamp]['spec'].loc[:, 'efth'] = efth_corrected
        sp1d_corrected['data'][timestamp]['depth'] = target_depth
        sp1d_corrected['data'][timestamp]['u'] = np.nan
        sp1d_corrected['data'][timestamp]['v'] = np.nan

    return sp1d_corrected
