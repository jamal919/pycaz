# -*- coding: utf-8 -*-

"""
Implements a static list of tidal potential that can be applied to a schism model.
"""
from pycaz.tide.utide import nodal_factor

PREDEFINED_WAVES = {
    'raid1990': {
        '2N2': {'spc': 2, 'amp': 0.006141, 'freq': 0.0001352404964640, 'nf': 0.96720, 'ear': 251.59},
        'K1': {'spc': 1, 'amp': 0.141565, 'freq': 0.0000729211583580, 'nf': 1.10338, 'ear': 324.30},
        'K2': {'spc': 2, 'amp': 0.030684, 'freq': 0.0001458423172010, 'nf': 1.28346, 'ear': 109.01},
        'L2': {'spc': 2, 'amp': 0.006931, 'freq': 0.0001431581055310, 'nf': 0.00000, 'ear': 325.06},
        'M2': {'spc': 2, 'amp': 0.242334, 'freq': 0.0001405189025090, 'nf': 0.96720, 'ear': 313.79},
        'MU2': {'spc': 2, 'amp': 0.007408, 'freq': 0.0001355937006840, 'nf': 0.96720, 'ear': 266.58},
        'N2': {'spc': 2, 'amp': 0.046397, 'freq': 0.0001378796994870, 'nf': 0.96720, 'ear': 102.69},
        'NU2': {'spc': 2, 'amp': 0.008811, 'freq': 0.0001382329037070, 'nf': 0.96720, 'ear': 117.68},
        'O1': {'spc': 1, 'amp': 0.100661, 'freq': 0.0000675977441510, 'nf': 1.16763, 'ear': 348.06},
        'P1': {'spc': 1, 'amp': 0.046848, 'freq': 0.0000725229459750, 'nf': 1.00000, 'ear': 39.25},
        'Q1': {'spc': 1, 'amp': 0.019273, 'freq': 0.0000649585411290, 'nf': 1.16763, 'ear': 136.96},
        'S2': {'spc': 2, 'amp': 0.112743, 'freq': 0.0001454441043330, 'nf': 1.00000, 'ear': 0.00},
        'T2': {'spc': 2, 'amp': 0.006608, 'freq': 0.0001452450073530, 'nf': 1.00000, 'ear': 52.32}
    }
}


class TidalPotential:
    """
    The tidal potential enters the momentum equation as a body force term.
    The self.wave variable is a set of values taken from Reid 1990.
    """

    def __init__(self):
        """
        Initialize a tidal potential object.

        The tidal potential enters the momentum equation as a body force term.
        The self.wave variable is a set of values taken from Reid 1990.
        """
        self.waves = {}

    def get_predefined_waves(self, name="default"):
        if name not in PREDEFINED_WAVES.keys() or name == "default":
            name = "raid1990"

        self.waves = PREDEFINED_WAVES[name]

    @property
    def consts(self):
        consts_list = [const for const in self.waves.keys()]
        return consts_list

    def update_nodal(self, at, lat, correct_phase=False):
        nf_consts = nodal_factor(t=at, consts=self.consts, lat=lat, correct_phase=correct_phase)

        for const in self.waves.keys():
            if const in nf_consts:
                self.waves[const].update(**nf_consts[const])

    def get_dict(self, wavelist='default') -> dict:
        dict_waves = {}
        if wavelist == 'default':
            dict_waves = {wave: self.waves[wave] for wave in self.waves.keys()}
            return dict_waves
        else:
            for wave in wavelist:
                if wave in self.waves.keys():
                    dict_waves[wave] = self.waves[wave]
                else:
                    print('Wave {:s} - Not found!'.format(wave))
            return dict_waves


def get_tidal_potential(at, consts="default", lat=0, correct_phase=False) -> dict:
    """
    Get the default set of tidal potentials with nodal correction

    :param at: DatetimeLike, time at which the nodal factor information will be updated
    :param consts: str|List, list of constituents or "default" for default set
    :param lat: float, latitude for nodal factor update
    :param correct_phase: bool, whether to correct the phase of the nodal factor
    :return:
    """
    potential = TidalPotential()
    potential.get_predefined_waves("default")
    potential.update_nodal(at=at, lat=lat, correct_phase=correct_phase)

    return potential.get_dict(wavelist=consts)
