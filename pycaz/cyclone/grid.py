#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Extends xarray

# def polar_coordinate(self, origin):
#     '''
#     Calculate the polar distance from a given point of interest.

#     For lon,lat values, the distance is calculated using great circle distance.
#     '''
#     try:
#         originx, originy = origin
#     except:
#         raise Exception('Origin must be a list of lon, lat')
    
#     X, Y = self.meshgrid
#     dfac = 60*1.852*1000
#     dist_x = dfac*np.cos(np.deg2rad(Y))*(X-originx)
#     dist_y = dfac*(Y-originy)
        
#     r = np.sqrt(dist_x**2 + dist_y**2)
#     theta = np.arctan2(dist_y, dist_x)

#     return(
#         Grid(
#             x=self.x,
#             y=self.y,
#             data=np.array([(rr, tt) for rr, tt in zip(r.flatten(), theta.flatten())])
#         )
#     )
