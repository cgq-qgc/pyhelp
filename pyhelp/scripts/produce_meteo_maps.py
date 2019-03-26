# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:54:25 2018
@author: jsgosselin
"""

# ---- Standard Library Imports

from itertools import product
import os.path as osp
import os

# ---- Third Party Imports

import netCDF4
from geopandas import GeoDataFrame
import pandas as pd
from shapely.geometry import Point, Polygon
import numpy as np

dirpath_netcdf = "D:/MeteoGrilleDaily"

# %% Get lat/lon from the netCDF

filename = osp.join(dirpath_netcdf, 'GCQ_v2_2000.nc')
netcdf_dset = netCDF4.Dataset(filename, 'r+')

lat = np.array(netcdf_dset['lat'])
lon = np.array(netcdf_dset['lon'])

netcdf_dset.close()

# %% Read the weather data from the InfoClimat grid

stack_precip = []
stack_tasmax = []
stack_tasmin = []
nyear = 0
for year in range(2000, 2015):
    print("\rProcessing year %d" % year, end=' ')
    filename = osp.join(dirpath_netcdf, 'GCQ_v2_%d.nc' % year)
    netcdf_dset = netCDF4.Dataset(filename, 'r+')

    stack_precip.append(np.array(netcdf_dset['pr']))
    stack_tasmax.append(np.array(netcdf_dset['tasmax']))
    stack_tasmin.append(np.array(netcdf_dset['tasmin']))

    netcdf_dset.close()
    nyear += 1
print('')

daily_precip = np.vstack(stack_precip)
daily_tasmax = np.vstack(stack_tasmax)
daily_tasmin = np.vstack(stack_tasmin)
daily_tasavg = (daily_tasmax + daily_tasmin) / 2

yearly_avg_precip = np.sum(daily_precip, axis=0) / nyear
yearly_avg_tasavg = np.average(daily_tasavg, axis=0)
yearly_avg_tasmax = np.average(daily_tasmax, axis=0)
yearly_avg_tasmin = np.average(daily_tasmin, axis=0)

# %% Create a grid

Np = len(lat) * len(lon)

geometry = []
arr_yearly_avg_precip = np.zeros(Np)
arr_avg_yearly_tasavg = np.zeros(Np)
arr_avg_yearly_tasmax = np.zeros(Np)
arr_avg_yearly_tasmin = np.zeros(Np)
i = 0
dx = dy = 0.1/2
for j, k in product(range(len(lat)), range(len(lon))):
    print("\rProcessing cell %d of %d" % (i, Np), end=' ')
    point = Point((lon[k], lat[j]))
    # polygon = Polygon([(lon[k]-dx, lat[j]-dy),
    #                    (lon[k]-dx, lat[j]+dy),
    #                    (lon[k]+dx, lat[j]+dy),
    #                    (lon[k]+dx, lat[j]-dy)])
    geometry.append(point)

    arr_yearly_avg_precip[i] = yearly_avg_precip[j, k]
    arr_avg_yearly_tasavg[i] = yearly_avg_tasavg[j, k]
    arr_avg_yearly_tasmax[i] = yearly_avg_tasmax[j, k]
    arr_avg_yearly_tasmin[i] = yearly_avg_tasmin[j, k]
    i += 1
print("\rProcessing cell %d of %d" % (i, Np))

# %%

print('\rFormating the data in a shapefile...', end=' ')
df = pd.DataFrame(data={'precip': arr_yearly_avg_precip,
                        'tasavg': arr_avg_yearly_tasavg,
                        'tasmax': arr_avg_yearly_tasmax,
                        'tasmin': arr_avg_yearly_tasmin})
crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"
gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
print('\rFormating the data in a shapefile... done')

print('\rSaving to Shapefile...', end=' ')
path_shp_out = ("D:/MeteoGrilleDaily/grid_yearly_meteo/grid_yearly_meteo.shp")
if not osp.exists(path_shp_out):
    os.makedirs(path_shp_out)

gdf.to_file(path_shp_out)
print('\rSaving to Shapefile... done', end=' ')
