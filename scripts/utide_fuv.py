#!/usr/bin/env python

import numpy as np
import pandas as pd

from utide.harmonics import FUV
from utide._ut_constants import ut_constants, constit_index_dict
from utide._time_conversion import _normalize_time as date2num

from typing import Union, Tuple, List
from numpy.typing import ArrayLike
from pandas._typing import TimestampConvertibleTypes

import os

import logging
logging.getLogger(__name__)

PathLike = Union[str, os.PathLike]

# required information from utide package
sat = ut_constants.sat
const = ut_constants.const
shallow = ut_constants.shallow

nshallow = np.ma.masked_invalid(const.nshallow).astype(int)
ishallow = np.ma.masked_invalid(const.ishallow).astype(int) - 1
not_shallow = ishallow.mask  # True where it was masked.
nshallow = nshallow.compressed()
ishallow = ishallow.compressed()
kshallow = np.nonzero(~not_shallow)[0]


# constituents from utide
mapping_utide = {
    wave:{'utide_name':wave, 'mapping_origin':'utide'} for wave in constit_index_dict
}

# constituents mapped from comodo and FES
comodo_utide = {
    'RO1':{'utide_name':'RHO1', 'mapping_origin':'FES'},
    'E2':{'utide_name':'EPS2', 'mapping_origin':'FES'},
    'LA2':{'utide_name':'LDA2', 'mapping_origin':'FES'},
    'MQ3':{'utide_name':'MO3', 'mapping_origin':'FES'},
    'M1':{'utide_name':'NO1', 'mapping_origin':'FES'},
    'KI1':{'utide_name':'CHI1', 'mapping_origin':'FES'},
    'TTA1':{'utide_name':'THE1', 'mapping_origin':'FES'}
}
mapping_utide.update(comodo_utide)

def utide_names(consts:ArrayLike):
    """Returns mapped utide constituents and missing constituents
    
    Returns the constituent names as found in utide as a {const:utide_const} dictionary. A second list is returned
    with names that are not found.

    Args:
        consts (ArrayLike): A constituent name, or a list of constituents names
    """
    consts = np.atleast_1d(consts)
    available = {}
    missing = []
    for const in consts:
        if const in mapping_utide:
            available[const] = mapping_utide[const]['utide_name']
        else:
            missing.append(const)
            logging.info(f'{const} is not found')
        
    return(available, missing)

def utide_lind(consts:Union[dict, ArrayLike]) -> list:
    """Return the internal indices of utide constituents. 
    
    This is useful for getting the right output from various utide routines. The consts list should only contain utide 
    names. use `utide_names` to find available names.

    Args:
        consts (Union[dict, ArrayLike]): A constituent name, or a list of constituents names

    Raises:
        TypeError: `consts` should be list of dict type

    Returns:
        list: utide internal indices for the consts
    """
    if isinstance(consts, dict):
        consts_list = [constit_index_dict[consts[const]] for const in consts]
    elif isinstance(consts, (np.ndarray, list)):
        consts_list = np.atleast_1d(consts)
    else:
        raise TypeError('consts should be one of list or dict type')
    
    lind = [constit_index_dict[const] for const in consts_list]
    
    return lind

def nodal_factor(t, consts, lat, correct_phase=True):
    """
    Compute the nodal correction, and astronomical argument at time t using the ut_FUV() in utide.

    Returns: Dictionary of Nodal correction, Astronomical argument for each constituent
    """
    t = date2num(t)
    
    #FIXME tref is only meaningful if ngflgs are changed, so the following calculation is, for now, useless
    if len(t) > 1:
        tref = (t[0] + t[-1])/2
    else:
        tref = t

    lind = utide_lind(consts)

    nf, U, V = FUV(t=t, tref=tref, lind=lind, lat=lat, ngflgs=[0, 0, 0, 0])

    if correct_phase:
        ear = (U + V) * 360 # cycles to degree
    else:
        ear = V * 360
        
    ear %= 360 # always positive from 0-360

    nf_dict = {
        const:{'nf':nf, 'ear':ear}
    }

    return nf_dict