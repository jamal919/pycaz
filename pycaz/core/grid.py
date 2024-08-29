#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata


class Grid:
    def __init__(self, x, y, data=None):
        """
        Grid object to generate grid and provides function to find various
        values at grid points.

        x: number of rows
        y: number of columns
        data: must be of the length len(x)*len(y). It will be flatten to get
        i,j,n formation

        TODO:
            - Simplify by removing depth, it should be simple grid
            - Converts to xr.DataAarray
            - Adds functionalities to be saved to netcdf, tiff
        """
        try:
            x_ = np.asarray(x)
            y_ = np.asarray(y)
        except:
            raise TypeError(f'x, y must be coercable to array')
        else:
            self.x = x_
            self.y = y_
            self.size = (len(self.x), len(self.y))

        # Data
        if data is None:
            self.depth = 1
            self.data = np.zeros(shape=(self.size[0], self.size[1], self.depth))
        elif isinstance(data, np.ndarray):
            self.depth = int(len(data.flatten()) / self.size[0] / self.size[1])
            try:
                self.data = np.array(data).reshape((self.size[0], self.size[1], self.depth))
            except:
                raise Exception('Size mismatch')

    @property
    def meshgrid(self):
        X, Y = np.meshgrid(self.x, self.y, indexing='ij')
        return (X, Y)

    def reshape(self):
        """
        Reshape the data to conform data structure.
        """
        self.depth = int(len(self.data.flatten()) / self.size[0] / self.size[1])
        self.data = self.data.reshape((self.size[0], self.size[1], self.depth))

    def __add__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data + other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data + other
                )
            )
        elif isinstance(other, np.ndarray):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data + other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] + other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __sub__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data - other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data - other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data - other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] - other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __mul__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data * other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data * other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data * other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] * other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __truediv__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data / other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data / other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data / other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] / other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __lt__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data < other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data < other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data < other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] < other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __le__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data <= other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data <= other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data <= other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] <= other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __gt__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data > other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data > other
                )
            )
        elif isinstance(other, np.ndarray):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data > other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] > other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __ge__(self, other):
        if isinstance(other, Grid):
            try:
                assert np.all(other.size == self.size)
                assert other.depth == self.depth
            except:
                raise ValueError('Uneuqal grid object')
            else:
                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data >= other.data
                    )
                )
        elif isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data >= other
                )
            )
        elif isinstance(other, (np.ndarray)):
            try:
                assert np.all(other.shape[0:2] == self.size)
                assert len(other.shape) == 3
                assert other.shape[2] == self.depth

                return (
                    Grid(
                        x=self.x,
                        y=self.y,
                        data=self.data >= other
                    )
                )
            except:
                try:
                    assert np.all(other.shape == self.size)

                    return_data = np.copy(self.data)
                    for i in np.arange(self.depth):
                        return_data[:, :, i] = return_data[:, :, i] >= other
                    return (
                        Grid(
                            x=self.x,
                            y=self.y,
                            data=return_data
                        )
                    )
                except:
                    raise ValueError('Unequal data shape')

    def __pow__(self, other):
        if isinstance(other, (float, int)):
            return (
                Grid(
                    x=self.x,
                    y=self.y,
                    data=self.data ** other
                )
            )
        else:
            raise ValueError('Only float or int as power')

    def __repr__(self):
        return self.data.__repr__()

    def __getitem__(self, key):
        return self.data[key]

    def __getattr__(self, name):
        return getattr(self.data, name)

    def __iter__(self):
        return iter(self.data.reshape((self.size[0] * self.size[1], self.depth)))

    def apply(self, func, **kwargs):
        f = lambda x: func(x, **kwargs)

        data = np.array([f(x) for x in self])

        return (
            Grid(
                x=self.x,
                y=self.y,
                data=data
            )
        )

    def polar_coordinate(self, origin):
        """
        Calculate the polar distance from a given point of interest.

        For lon,lat values, the distance is calculated using great circle distance.
        """
        try:
            originx, originy = origin
        except:
            raise Exception('Origin must be a list of lon, lat')

        X, Y = self.meshgrid
        dfac = 60 * 1.852 * 1000
        dist_x = dfac * np.cos(np.deg2rad(Y)) * (X - originx)
        dist_y = dfac * (Y - originy)

        r = np.sqrt(dist_x ** 2 + dist_y ** 2)
        theta = np.arctan2(dist_y, dist_x)

        return (
            Grid(
                x=self.x,
                y=self.y,
                data=np.array([(rr, tt) for rr, tt in zip(r.flatten(), theta.flatten())])
            )
        )

    def interpolate(self, at, depth=0, method='linear', fill_value=np.nan, rescale=False):
        """
        Interpolate at another x,y point or grid using scipy.interpolate.griddata

        at: {list, tuple, Grid} instance
        depth: depth of grid data to interpolate
        method: {'linear', 'nearest', 'cubic'}, optional
        fill_value: value used to fill in for requested point outside of convex hull
        rescale: rescale points to unit cube before preforming interpolation

        return Grid
        """
        X, Y = self.meshgrid

        points = np.array([(x, y) for x, y in zip(X.flatten(), Y.flatten())])
        values = self[:, :, depth].flatten()

        if isinstance(at, (list, tuple)):
            # For x, y list or tuple
            return (
                griddata(
                    points, values,
                    at,
                    method=method,
                    fill_value=fill_value,
                    rescale=rescale
                )
            )

        if isinstance(at, Grid):
            return (
                Grid(
                    x=at.x,
                    y=at.y,
                    data=griddata(
                        points, values, at.meshgrid,
                        method=method,
                        fill_value=fill_value,
                        rescale=rescale
                    )
                )
            )
