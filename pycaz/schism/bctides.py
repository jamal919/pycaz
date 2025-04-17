# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
from pycaz.schism.hgrid import OpenBoundary
from pycaz.schism.tidefac import Tidefac
import warnings
import os
import logging

from typing import List, Dict, Union
from pycaz.typing import PathLike

logger = logging.getLogger(__name__)


class Bctides(dict):
    def __init__(self, **kwargs):
        """ A bctides object extended from dictonaries
        
        Additional key-value pairs can be added using keyworded arguments.

        TODO: To update with a more dedicated class, based on dataclass or pydantic.
        """
        super().__init__(self)
        # initialize everything with empty stuff
        self.update(
            header='bctides',
            potential={
                'ntip': 0,
                'tip_dp': 0,
                'const': {}
            },
            tidefr={
                'nbfr': 0,
                'const': {}
            },
            open_bnds={}
        )
        # then update with the kwargs
        self.update(kwargs)

    def copy(self):
        return deepcopy(self)

    @property
    def header(self) -> str:
        return self['header']

    @header.setter
    def header(self, header_text: str):
        self['header'] = header_text

    @property
    def potential(self) -> Dict:
        return self['potential']

    @property
    def tidefr(self) -> Dict:
        return self['tidefr']

    @property
    def tidefr_consts(self) -> List[str]:
        if self['tidefr']['nbfr']:
            return list(self['tidefr']['const'].keys())
        else:
            return []

    @property
    def potential_consts(self) -> List[str]:
        if self['potential']['const']:
            return list(self['potential']['const'].keys())
        else:
            return []

    @property
    def consts(self) -> List[str]:
        """
        Returns a unique list of all the constituents (potential and forcing) used in the bctides.

        For individual constituents, use either tidefr_consts or potential_consts property.

        :return List: A list of constituents name
        """
        _all_consts = list()
        for const in self.tidefr_consts:
            _all_consts.append(const)

        for const in self.potential_consts:
            _all_consts.append(const)

        return list(set(_all_consts))

    @property
    def open_bnds(self) -> Dict:
        """
        Returns the OpenBoundary dictionary.

        :return: Dictionary of the OpenBoundary
        """
        return self['open_bnds']

    def describe(self) -> None:
        """
        Prints a description of the bctides in plain text.
        :return: None

        TODO: Override with __repr__?
        """
        print(self['header'])
        ntip = self['potential']['ntip']
        tip_dp = self['potential']['tip_dp']
        if self['potential']['ntip']:
            print('For depth >', tip_dp, ',', ntip, 'const for tidal potential', list(self['potential']['const']))

        nbfr = self['tidefr']['nbfr']
        if nbfr:
            print(nbfr, 'tidal const for BC:', list(self['tidefr']['const'].keys()))
        else:
            print('No tidal boundary constituents defined')

        print(len(self['open_bnds']), 'open boundaries')
        for bnd in self['open_bnds']:
            name = self['open_bnds'][bnd]['name']
            iettype = self['open_bnds'][bnd]['iettype']
            ifltype = self['open_bnds'][bnd]['ifltype']
            itetype = self['open_bnds'][bnd]['itetype']
            isatype = self['open_bnds'][bnd]['isatype']
            print(
                f'Boundary {bnd} [{name}] - iettype: {iettype}, ifltype: {ifltype}, itetype: {itetype}, isatype: {isatype}')

    def add_tidefr(self, tidefr: Dict) -> None:
        """
        Adds the tidefr into the Bctides.

        :param tidefr:
        :return: None
        """
        self.tidefr['const'].update(tidefr)
        self.tidefr['nbfr'] = len(self.tidefr['const'])

    def update_nodal(self, tidefac: Tidefac, nodal: bool = True, eq_arg: bool = True) -> None:
        """
        Update the nodal factor in the bctides with the given tidefac content.

        Tidefac can be generated using 'tidefac.f' and then read using read pycaz.schism.tidefac.read_tidefacout() or
        using the pycaz.schis.tidefac.generate_tidefac_utide() function.

        :param tidefac: Tidefac object
        :param nodal: If the nodal factor be updated, default True.
        :param eq_arg: If the equilibrium arguments to be updated, default True.
        :return: None

        TODO: Adds inplace functionality and return updated bctides otherwise. Currently defaults to inplace.
        """
        update_bctide(bctides=self, tidefac=tidefac, nodal=nodal, eq_arg=eq_arg, inplace=True)

    def write(self, fname: PathLike = 'bctides.in', replace: bool = False) -> None:
        """
        Write bctide to a file, typically after updates are applied on it.

        :param fname: Path of the file to be written, default bctides.in
        :param replace: Existing file, if exists, to be replaced or not, default False.
        :return: None
        """
        write_bctides(bctides=self, fname=fname, replace=replace)


def read_bctides(fname: PathLike) -> Bctides:
    """
    Read a given bctide from `fname`.

    :param fname: Path to bctide file
    :return: Bctides
    """
    bctides = Bctides()
    with open(fname) as f:
        txt = f.readlines()

    txt = [t.split('!')[0].strip() for t in txt]

    # Header
    bctides.update(header=txt[0])

    # Tidal potential
    ntip, tip_dp = np.fromstring(txt[1], count=2, sep=' ')
    bctides.update(potential={
        'ntip': int(ntip),
        'tip_dp': tip_dp,
        'const': {}
    })

    ln = 1

    for k in np.arange(bctides['potential']['ntip']):
        ln += 1
        talpha = txt[ln].strip()
        ln += 1
        jspc, tamp, tfreq, tnf, tear = np.fromstring(txt[ln], count=5, sep=' ')
        bctides['potential']['const'][talpha] = {
            'spc': int(jspc),
            'amp': tamp,
            'freq': tfreq,
            'nf': tnf,
            'ear': tear
        }

    # Tidal frequencies
    ln += 1
    nbfr = int(txt[ln])
    bctides.update(tidefr={
        'nbfr': nbfr,
        'const': {}
    })

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

    # we want keep the sequence of the boundary intact
    # so boundaries will be identified using j+1 dictionary key
    for j in np.arange(nopen):
        boundary = OpenBoundary(name=f'{j + 1}')
        ln += 1
        neta, iettype, ifltype, itetype, isatype = np.fromstring(txt[ln], dtype=int, count=5, sep=' ')

        if iettype == 0 and ifltype == 0:
            warnings.warn(f'Boundary {j} : Both elevation and flow are set to 0! One of them must be active.')

        boundary.update(neta=neta, iettype=iettype, ifltype=ifltype, itetpye=itetype, isatype=isatype)

        # Elevation boundary conditions, iettype
        if iettype == 0 or iettype == 1 or iettype == 4:
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
                emo_efa = np.genfromtxt(txt[ln:ln + neta], encoding='UTF8')
                values[alpha] = emo_efa
                ln += neta - 1  # removes 1 for 0-based indexing
            boundary['et'][iettype] = values

        # Velocity boundary conditions, ifltype
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
                emo_efa = np.genfromtxt(txt[ln:ln + neta], encoding='UTF8')
                values[alpha] = emo_efa
                ln += neta - 1  # removes 1 for 0-based indexing
            boundary['fl'][ifltype] = values
        elif ifltype == -4:
            # time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
            # rel1, rel2: relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)
            ln += 1
            rel1, rel2 = np.fromstring(txt[ln], count=2, sep=' ')
            boundary['fl'][ifltype] = {'rel1': rel1, 'rel2': rel2}
        elif ifltype == -1:
            # flather type boundary condition, iettype must be 0
            ln += 1  # should give a text value 'eta_mean'
            ln += 1  # starts eta_mean values
            eta_mean = np.genfromtxt(txt[ln:ln + neta])
            ln += neta - 1
            ln += 1  # should give a text value 'vn_mean'
            ln += 1  # starts vn_mean values
            vn_mean = np.genfromtxt(txt[ln:ln + neta])
            ln += neta - 1
            boundary['fl'][ifltype] = {
                'eta_mean': eta_mean,
                'vn_mean': vn_mean
            }

        # Temperature boundary condition
        if itetype == 0:
            # 0: Temperature not specified
            pass
        elif itetype == 1:
            # time history of temperature on this boundary
            # here only nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; 
            # time history of temperature will be read in from TEM_1.th (ASCII)
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][itetype] = tobc
        elif itetype == 2:
            # this boundary is forced by a constant temperature
            # constant temperature on this segment
            ln += 1
            tthconst = float(txt[ln])
            # nudging factor (between 0 and 1) for inflow
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][itetype] = {'tthconst': tthconst, 'tobc': tobc}
        elif itetype == 3:
            # initial temperature profile for inflow
            # nudging factor (between 0 and 1) for inflow
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][itetype] = tobc
        elif itetype == 4:
            # 3D input
            # time history of temperature is read in from TEM_3D.th.nc (netcdf)
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][itetype] = tobc

        # Salinity boundary condition
        if isatype == 0:
            # 0: Temperature not specified
            pass
        elif isatype == 1:
            # time history of salinity on this boundary
            # here only nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; 
            # time history of temperature will be read in from SAL_1.th (ASCII)
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][isatype] = tobc
        elif isatype == 2:
            # this boundary is forced by a constant salinity
            # constant salinity on this segment
            ln += 1
            tthconst = float(txt[ln])
            # nudging factor (between 0 and 1) for inflow
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][isatype] = {'tthconst': tthconst, 'tobc': tobc}
        elif isatype == 3:
            # initial salinity profile for inflow
            # nudging factor (between 0 and 1) for inflow
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][isatype] = tobc
        elif isatype == 4:
            # 3D input
            # time history of salinity is read in from SAL_3D.th.nc (netcdf)
            ln += 1
            tobc = float(txt[ln])
            boundary['te'][isatype] = tobc

        bctides['open_bnds'][j + 1] = boundary

    return (bctides)


def update_bctide(bctides: Bctides, tidefac: Tidefac, nodal: bool = True, eq_arg: bool = True, inplace: bool = False):
    """
    Update `bctides` with `tidefac`. If `inplace=True`, then `bctides` will be copied first.

    :param bctides: Bctides object, either read or generated.
    :param tidefac: Tidefac object, either read or generated.
    :param nodal: If the nodal factor should be updated.
    :param eq_arg: If the equilibrium argument to be updated.
    :param inplace: To replace the original in memory or not, ddefault False.
    :return:
    """
    if inplace:
        bctides_new = bctides  # refer to the original
    else:
        bctides_new = bctides.copy()

    if nodal:
        logger.debug('Nodal factor update requested')

    if eq_arg:
        logger.debug('Equilibrium argument update requested')

    # update potential
    _updated = []
    _not_updated = []
    if 'const' in bctides['potential']:
        for const in bctides['potential']['const']:
            if const in tidefac.consts:
                if nodal:
                    bctides['potential']['const'][const].update(nf=tidefac['const'][const]['nf'])
                if eq_arg:
                    bctides['potential']['const'][const].update(ear=tidefac['const'][const]['ear'])
                _updated.append(const)
            else:
                _not_updated.append(const)

    logger.info('Potential: updated - ' + ', '.join(_updated))
    logger.info('Potential: not updated - ' + ', '.join(_not_updated))

    # update forced tidal harmonics
    _updated = []
    _not_updated = []
    if 'const' in bctides['tidefr']:
        for const in bctides['tidefr']['const']:
            if const in tidefac.consts:
                if nodal:
                    bctides['tidefr']['const'][const].update(ff=tidefac['const'][const]['nf'])
                if eq_arg:
                    bctides['tidefr']['const'][const].update(face=tidefac['const'][const]['ear'])
                _updated.append(const)
            else:
                _not_updated.append(const)

    logger.info('Tidefr: updated - ' + ', '.join(_updated))
    logger.info('Tidefr: not updated - ' + ', '.join(_not_updated))

    # update header
    bctides_new.header = tidefac.info

    if not inplace:
        return bctides_new


def write_bctides(bctides: Bctides, fname: PathLike, replace: bool = False) -> None:
    """
    Write the `bctides` to file given by `fname`. If `replace=True`, then the existing file will be overwritten.

    :param bctides: Bctides object to be written.
    :param fname: Path to the file to be wrrien.
    :param replace: If the file exists, overwrite if True, defaults to False.
    :return: None
    """
    if os.path.exists(fname) and not replace:
        raise Exception(f'{fname} already exists! Set replace=True if you want to replace.')

    with open(fname, 'w') as f:
        # Header
        f.write('{header}\n'.format(**bctides))

        # Tidal potential
        f.write('{ntip}\t{tip_dp} !ntip, tip_dp\n'.format(**bctides['potential']))
        for const in bctides['potential']['const']:
            f.write(f'{const}\n')
            f.write('{spc:1d}\t{amp:.6f}\t{freq:.15f}\t{nf:.6f}\t{ear:.2f}\n'.format(
                **bctides['potential']['const'][const]))

        # Tidal frequencies
        f.write('{nbfr} !nbfr\n'.format(**bctides['tidefr']))
        for const in bctides['tidefr']['const']:
            f.write(f'{const}\n')
            f.write('{amig:.15f}\t{ff:.6f}\t{face:.2f}\n'.format(**bctides['tidefr']['const'][const]))

        # Open boundaries
        nopen = len(bctides['open_bnds'])
        f.write(f'{nopen} !Number of Open Boundaries\n')
        for bnd in np.arange(nopen) + 1:
            boundary = bctides['open_bnds'][bnd]
            name = boundary['name']
            neta = boundary['neta']
            iettype = boundary['iettype']
            ifltype = boundary['ifltype']
            itetype = boundary['itetype']
            isatype = boundary['isatype']
            f.write(f'{neta} {iettype} {ifltype} {itetype} {isatype} !Boundary {bnd} [{name}]\n')

            # Elevation boundary conditions
            if iettype == 0 or iettype == 1 or iettype == 4:
                # 0 : elevations are not specified for this boundary (in this case the velocity must be specified)
                # 1 : no input in bctides.in; time history of elevation is read in from elev.th (ASCII)
                # 4 : no input in this file; time history of elevation is read in from elev2D.th.nc (netcdf)
                pass
            elif iettype == 2:
                # constant elevation value for this segment
                ethconst = boundary['et'][iettype]
                f.write(f'{ethconst}\n')
            elif iettype == 3 or iettype == 5:
                # 3: tidal forcing
                # 5: combination of 3 and 4
                for const in bctides['tidefr']['const']:
                    f.write(f'{const}\n')
                    np.savetxt(
                        fname=f,
                        X=boundary['et'][iettype][const],
                        fmt=['%.15f', '%.15f'],
                        delimiter='\t'
                    )

            # Check if both of ifltype and iettype is set to 0, and raise exception
            if iettype == 0 and ifltype == 0:
                raise Exception(
                    f'Bad bctides! Both iettype and ifltype set to 0 for Boundary {bnd}, atleast one BC needed.')

            # Velocity boundary conditions
            if ifltype == 0 or ifltype == 1 or ifltype == 4:
                # 0: no boundary specified, no input needed. Elev boundary must be specified.
                # 1: no input in this file; time history of discharge is read in from flux.th (ASCII)
                # 4: time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
                pass
            elif ifltype == 2:
                # constant discharge (note that a negative number means inflow)
                vthconst = boundary['fl'][ifltype]
                f.write(f'{vthconst}\n')
            elif ifltype == 3 or ifltype == 5:
                # vel. (not discharge!) is forced in frequency domain
                # 3: tidal forcing
                # 5: combination of 3 and 4
                for const in bctides['tidefr']['const']:
                    f.write(f'{const}\n')
                    np.savetxt(
                        fname=f,
                        X=boundary['fl'][ifltype][const],
                        fmt=['%.15f', '%.15f'],
                        delimiter='\t'
                    )
            elif ifltype == -4:
                # time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
                # rel1, rel2: relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)
                f.write('{rel1} {rel2}\n').format(**boundary['fl'][ifltype])
            elif ifltype == -1:
                # flather type boundary condition, iettype must be 0
                f.write('eta_mean !mean elevation below\n')
                np.savetxt(
                    fname=f,
                    X=boundary['fl'][ifltype]['eta_mean'],
                    fmt=['%.2f']
                )
                f.write('vn_mean !mean normal velocity\n')
                np.savetxt(
                    fname=f,
                    X=boundary['fl'][ifltype]['vn_mean'],
                    fmt=['%.3f', '%0.3f']
                )

            # Temperature boundary condition
            if itetype == 0:
                # 0: Temperature not specified
                pass
            elif itetype == 1:
                # time history of temperature on this boundary
                # here only nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; 
                # time history of temperature will be read in from TEM_1.th (ASCII)
                tobc = boundary['te'][itetype]
                f.write(f'{tobc}\n')
            elif itetype == 2:
                # this boundary is forced by a constant temperature
                # constant temperature on this segment
                f.write('{tthconst}\n{tobc}\n'.format(**boundary['te'][itetype]))
            elif itetype == 3:
                # initial temperature profile for inflow
                # nudging factor (between 0 and 1) for inflow
                tobc = boundary['te'][itetype]
                f.write(f'{tobc}\n')

            # Salinity boundary condition
            if itetype == 0:
                # 0: Salinity not specified
                pass
            elif itetype == 1:
                # time history of salinity on this boundary
                # here only nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; 
                # time history of salinity will be read in from SAL_1.th (ASCII)
                tobc = boundary['te'][itetype]
                f.write(f'{tobc}\n')
            elif itetype == 2:
                # this boundary is forced by a constant salinity
                # constant salinity on this segment
                f.write('{tthconst}\n{tobc}\n'.format(**boundary['te'][itetype]))
            elif itetype == 3:
                # initial salinity profile for inflow
                # nudging factor (between 0 and 1) for inflow
                tobc = boundary['te'][itetype]
                f.write(f'{tobc}\n')
