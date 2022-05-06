#!/usr/bin/env python
# -*- coding: utf-8 -*-

#%% Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.kdtree import KDTree
import cartopy.crs as ccrs

#%% Loading mesh data
fname = '/home/khan/MEGA/Models/SCHISM/Mesh_v2.x/Mesh/Mesh_v2.3_Variable_Polder_Boundary_Fixed/hgrid.gr3.orig'
with open(fname) as f:
    name = f.readline()
    nelem, nnode = [int(i) for i in f.readline().split()]

print(f'Mesh with {nelem} elements, and {nnode} nodes')
nodes = pd.read_csv(
    fname, 
    nrows=nnode, 
    skiprows=2, 
    header=None, 
    delim_whitespace=True, 
    names=['id', 'lon', 'lat', 'depth']
    ).set_index('id')
elements = pd.read_csv(
    fname, 
    nrows=nelem, 
    skiprows=2+nnode, 
    header=None,
    delim_whitespace=True, 
    names=['id', 'i34', 'node1', 'node2', 'node3']
    ).set_index('id')

# Change elements to python indexing, not used here
elements.loc[:, ['node1', 'node2', 'node3']] = elements.loc[:, ['node1', 'node2', 'node3']] - 1 

#%% Loading station file
stations = pd.read_csv(
    'station.in', 
    skiprows=2, 
    delim_whitespace=True, 
    header=None, 
    names=['id', 'lon', 'lat', 'depth', 'name']
    ).set_index('id')

#%% Create KDTree
kdtree = KDTree(nodes.loc[:, ['lon', 'lat']])
dists, inode = kdtree.query(stations.loc[:, ['lon', 'lat']])

#%% Update stations with model information
stations_model = nodes.iloc[inode]
stations['inode'] = inode
stations['model_lon'] = stations_model.lon.values
stations['model_lat'] = stations_model.lat.values
stations['dist_stn_model'] = dists
stations['depth'] = stations_model.depth.values # Replaces with local depth

#%% Save to file
stations.to_csv('stations.csv')

#%% Create a figure to see the farmost station
farmost = stations.iloc[np.argmax(dists)]

fig, ax = plt.subplots(figsize=(20, 20), facecolor='white', subplot_kw={'projection':ccrs.PlateCarree()})
ax.scatter(stations.lon, stations.lat, c='red', marker='*')
ax.scatter(stations.model_lon, stations.model_lat, c='green', marker='.')
ax.annotate(farmost['name'], farmost.loc[['model_lon', 'model_lat']] + [0.1, 0.1])
ax.coastlines()
plt.savefig('model_stations.png', dpi=150, bbox_inches='tight')
