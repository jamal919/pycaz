#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implements a static list of tidal potential that can be applied to a schism model.
"""

class TidalPotential(object):
    """
    The tidal potential enters the momentum equation as a body force term.
    The self.wave variable is a set of values taken from Reid 1990. 

    TODO: Implement a set of potential from the comodo-toolbox
    """
    def __init__(self):
        self.waves = {
            '2N2':{'jspc':2, 'tamp':0.006141, 'tfreq':0.0001352404964640, 'tnf':0.96720, 'tear':251.59},
            'K1':{'jspc':1, 'tamp':0.141565, 'tfreq':0.0000729211583580, 'tnf':1.10338, 'tear':324.30},
            'K2':{'jspc':2, 'tamp':0.030684, 'tfreq':0.0001458423172010, 'tnf':1.28346, 'tear':109.01},
            'L2':{'jspc':2, 'tamp':0.006931, 'tfreq':0.0001431581055310, 'tnf':0.00000, 'tear':325.06},
            'M2':{'jspc':2, 'tamp':0.242334, 'tfreq':0.0001405189025090, 'tnf':0.96720, 'tear':313.79},
            'MU2':{'jspc':2, 'tamp':0.007408, 'tfreq':0.0001355937006840, 'tnf':0.96720, 'tear':266.58},
            'N2':{'jspc':2, 'tamp':0.046397, 'tfreq':0.0001378796994870, 'tnf':0.96720, 'tear':102.69},
            'NU2':{'jspc':2, 'tamp':0.008811, 'tfreq':0.0001382329037070, 'tnf':0.96720, 'tear':117.68},
            'O1':{'jspc':1, 'tamp':0.100661, 'tfreq':0.0000675977441510, 'tnf':1.16763, 'tear':348.06},
            'P1':{'jspc':1, 'tamp':0.046848, 'tfreq':0.0000725229459750, 'tnf':1.00000, 'tear':39.25},
            'Q1':{'jspc':1, 'tamp':0.019273, 'tfreq':0.0000649585411290, 'tnf':1.16763, 'tear':136.96},
            'S2':{'jspc':2, 'tamp':0.112743, 'tfreq':0.0001454441043330, 'tnf':1.00000, 'tear':0.00},
            'T2':{'jspc':2, 'tamp':0.006608, 'tfreq':0.0001452450073530, 'tnf':1.00000, 'tear':52.32}
            }

    def value(self, wavelist='default'):
        __values = {}
        if wavelist=='default':
            __values = {wave:self.waves[wave] for wave in self.waves.keys()}
            return(__values)
        else:
            for wave in wavelist:
                if wave in self.waves.keys():
                    __values[wave] = self.waves[wave]
                else:
                    print('Wave {:s} - Not found!'.format(wave))
            return(__values)