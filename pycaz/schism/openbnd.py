# -*- encoding: utf-8 -*-

import numpy as np

class OpenBoundary(dict):
    def __init__(self, **kwargs):
        """
        A Open boundary object extended from dictonaries

        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.reset_all()
        self.update(kwargs)

    def reset_all(self):
        self.reset_properties()
        self.reset_bnd()

    def reset_properties(self):
        self.update(
            name='',
            neta=0,
            nodes=np.empty(0, dtype=int),
            xy=np.empty([0, 2], dtype=float)
        )

    def reset_bnd(self):
        self.update(
            iettype=0,  # elevation
            et={},
            ifltype=0,  # flow/current
            fl={},
            itetype=0,  # temperature
            te={},
            isatype=0,  # salinity
            sa={}
        )

    @property
    def name(self):
        return self['name']

    @name.setter
    def name(self, new_name: str):
        self['name'] = new_name

    @property
    def nodes(self):
        return self['nodes']

    @property
    def neta(self):
        return self['neta']

    @property
    def xy(self):
        return self['xy']

    @property
    def iettype(self):
        return self['iettype']

    @property
    def ifltype(self):
        return self['ifltype']

    @property
    def itetype(self):
        return self['itetype']

    @property
    def isatype(self):
        return self['isatype']

