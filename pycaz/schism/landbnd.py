# -*- encoding: utf-8 -*-
import numpy as np

class LandBoundary(dict):
    def __init__(self, **kwargs):
        """
        A Land boundary object extended from dictonaries

        Additional key-value pairs can be added using keyworded arguments.
        """
        super().__init__(self)
        self.reset()
        self.update(kwargs)

    def reset(self):
        self.update(
            name='',
            bndtype=1,
            nodes=np.empty(0, dtype=int),
            xy=np.empty([0, 2], dtype=float)
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