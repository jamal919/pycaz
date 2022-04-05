#!/usr/bin/env python
# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
from .hgrid import Boundary
from .tidefac import Tidefac
import warnings

class Bctides(dict):
    def __init__(self, **kwargs):
        """ A bctides object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.update(kwargs)

    @property
    def header(self):
        return(self['header'])
    
    @header.setter
    def header(self, header_text: str):
        self['header'] = header_text

    @property
    def potential(self):
        return(self['potential'])

    @property
    def tidefr(self):
        return(self['tidefr'])

    @property
    def boundaries(self):
        return(self['boundaries'])

    def update_nodal(self, tidefac: Tidefac):
        update_bctide(bctides=self, tidefac=tidefac, inplace=True)


def read_bctides(fname: str) -> Bctides:
    bctides = Bctides()
    with open(fname) as f:
        txt = f.readlines()

    txt = [t.split('!')[0].strip() for t in txt]

    # Header
    bctides['header'] = txt[0]

    # Tidal potential
    ntip, tip_dp = np.fromstring(txt[1], count=2, sep=' ')
    bctides['potential'] = {
        'ntip': int(ntip),
        'tip_dp': tip_dp,
        'const': {}
    }

    ln = 1

    for k in np.arange(bctides['potential']['ntip']):
        ln += 1
        talpha = txt[ln].strip()
        ln += 1
        jspc, tamp, tfreq, tnf, tear = np.fromstring(txt[ln], count=5, sep=' ')
        bctides['potential']['const'][talpha] = {
            'spc':int(jspc),
            'amp':tamp,
            'freq':tfreq,
            'nf':tnf,
            'ear':tear
        }

    # Tidal frequencies
    ln += 1
    nbfr = int(txt[ln])
    bctides['tidefr'] = {
        'nbfr': nbfr,
        'const': {}
    }

    for k in np.arange(bctides['tidefr']['nbfr']):
        ln += 1
        alpha = txt[ln].strip()
        ln += 1
        amig, ff, face = np.fromstring(txt[ln], count=3, sep=' ')
        bctides['tidefr']['const'][alpha] = {
            'amig': amig,
            'ff': ff,
            'face': face
        }

    # Boundaries
    ln += 1
    nopen = int(txt[ln])
    bctides['open_bnds'] = {}

    for j in np.arange(nopen):
        boundary = Boundary(name=f'{j+1}', et={}, fl={}, te={}, sa={})
        ln += 1
        neta, iettype, ifltype, itetype, isatype = np.fromstring(txt[ln], dtype=int, count=5, sep=' ')
        
        if iettype == 0 and ifltype == 0:
            warnings.warn(f'Boundary {j} : Both elevation and flow are set to 0! One of them must be active.')
        
        boundary.update(neta=neta, iettype=iettype, ifltype=ifltype, isatype=isatype)
        
        # Elevation boundary conditions
        if iettype == 0 or iettype == 1 or iettype == 4 :
            # 0 : elevations are not specified for this boundary (in this case the velocity must be specified)
            # 1 : no input in bctides.in; time history of elevation is read in from elev.th (ASCII)
            # 4 : no input in this file; time history of elevation is read in from elev2D.th.nc (netcdf)
            pass
        elif iettype == 2:
            # constant elevation value for this segment
            ln += 1
            ethconst = float(txt[ln])
            boundary['et'][iettype] = ethconst
        elif iettype == 3 or iettype == 5:
            # 3: tidal forcing
            # 5: combination of 3 and 4
            values = {}
            for k in np.arange(bctides['tidefr']['nbfr']):
                ln += 1
                alpha = txt[ln].strip()
                ln += 1
                emo_efa = np.genfromtxt(txt[ln:ln+neta], encoding='UTF8')
                values[alpha] = emo_efa
                ln += neta - 1 # removes 1 for 0-based indexing
            boundary['et'][iettype] = values

        # Velocity boundary conditions
        if ifltype == 0 or ifltype == 1 or ifltype == 4:
            # 0: no boundary specified, no input needed. Elev boundary must be specified.
            # 1: no input in this file; time history of discharge is read in from flux.th (ASCII)
            # 4: time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
            pass
        elif ifltype == 2:
            # constant discharge (note that a negative number means inflow)
            ln += 1
            vthconst = float(txt[ln])
            boundary['fl'][ifltype] = vthconst
        elif ifltype == 3 or ifltype == 5:
            # vel. (not discharge!) is forced in frequency domain
            # 3: tidal forcing
            # 5: combination of 3 and 4
            values = {}
            for k in np.arange(bctides['tidefr']['nbfr']):
                ln += 1
                alpha = txt[ln].strip()
                ln += 1
                emo_efa = np.genfromtxt(txt[ln:ln+neta], encoding='UTF8')
                values[alpha] = emo_efa
                ln += neta - 1 # removes 1 for 0-based indexing
            boundary['fl'][ifltype] = values
        elif ifltype == -4:
            # time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
            # rel1, rel2: relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)
            ln += 1
            rel1, rel2 = np.fromstring(txt[ln], count=2, sep=' ')
            boundary['fl'][ifltype] = {'rel1': rel1, 'rel2':rel2}
        elif ifltype == -1:
            # flather type boundary condition, iettype must be 0
            ln += 1 # should give a text value 'eta_mean'
            ln += 1 # starts eta_mean values
            eta_mean = np.genfromtxt(txt[ln:ln+neta])
            ln += neta - 1
            ln += 1 # should give a text value 'vn_mean'
            ln += 1 # starts vn_mean values
            vn_mean = np.genfromtxt(txt[ln:ln+neta])
            ln += neta -1
            boundary['fl'][ifltype] = {
                'eta_mean': eta_mean,
                'vn_mean': vn_mean
            }

        bctides['open_bnds'][j+1] = boundary
    
    return(bctides)

def update_bctide(bctides: Bctides, tidefac: Tidefac, inplace:bool = False):
    """
    Update nodal information from a tidefac object.
    """
    if inplace:
        bctides_new = bctides
    else:
        bctides_new = deepcopy(bctides)

    # update potential
    if 'const' in bctides['potential']:
        print('Updating tidal potential...')
        for const in bctides['potential']['const']:
            if const in tidefac['const']:
                bctides['potential']['const'][const].update(tidefac['const'][const])
                print(f'\t{const} -> Updated')
            else:
                print(f'\t{const} -> Not updated')

    # update potential
    if 'const' in bctides['tidefr']:
        print('Updating tidal constituents...')
        for const in bctides['tidefr']['const']:
            if const in tidefac['const']:
                bctides['tidefr']['const'][const].update(
                    ff=tidefac['const'][const]['nf'],
                    face=tidefac['const'][const]['ear']
                    )
                print(f'\t{const} -> Updated')
            else:
                print(f'\t{const} -> Not updated')
    
    # update header
    bctides_new.header = tidefac.info
    
    if not inplace:
        return bctides_new
