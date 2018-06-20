# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:54:25 2018
@author: jsgosselin
"""

# ---- Standard Library Imports

from itertools import product
import os.path as osp

# ---- Third Party Imports

import matplotlib.pyplot as plt
import h5py
import netCDF4
import geopandas as gpd
from geopandas import GeoDataFrame
import pandas as pd
from shapely.geometry import Point, Polygon
import numpy as np

import rasterio as rio
from rasterio import features
from affine import Affine


# ---- Local Librairies Import

from help_utils import NetCDFMeteoManager

dirpath_netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily"

# %% Get lat/lon from the netCDF

filename = osp.join(dirpath_netcdf, 'GCQ_v2_2000.nc')
netcdf_dset = netCDF4.Dataset(filename, 'r+')

lat = np.array(netcdf_dset['lat'])
lon = np.array(netcdf_dset['lon'])

netcdf_dset.close()

# %%

stack_precip = []
stack_tasmax = []
stack_tasmin = []
for year in range(2000, 2015):
    print("\rLoading year %d" % year, end=' ')
    filename = osp.join(dirpath_netcdf, 'GCQ_v2_%d.nc' % year)
    netcdf_dset = netCDF4.Dataset(filename, 'r+')

    stack_precip.append(np.array(netcdf_dset['pr']))
    stack_tasmax.append(np.array(netcdf_dset['tasmax']))
    stack_tasmin.append(np.array(netcdf_dset['tasmin']))

    netcdf_dset.close()
print('')

daily_precip = np.vstack(stack_precip)
daily_tasmax = np.vstack(stack_tasmax)
daily_tasmin = np.vstack(stack_tasmin)
daily_tasavg = (daily_tasmax + daily_tasmin)/2

yearly_avg_precip = np.sum(daily_precip, axis=0) / 15
yearly_avg_tasavg = np.average(daily_tasavg, axis=0)

# %%

Np = len(lat) * len(lon)
Ny = 15

geometry = []
arr_yearly_avg_precip = np.zeros(Np)
arr_avg_yearly_tasavg = np.zeros(Np)
i = 0
dx = dy = 0.1/2
for j, k in product(range(len(lat)), range(len(lon))):
    print("\rProcessing cell %d of %d" % (i, Np), end=' ')
    # point = Point((lon[k], lat[j]))
    polygon = Polygon([(lon[k]-dx, lat[j]-dy),
                       (lon[k]-dx, lat[j]+dy),
                       (lon[k]+dx, lat[j]+dy),
                       (lon[k]+dx, lat[j]-dy)])
    geometry.append(polygon)

    arr_yearly_avg_precip[i] = yearly_avg_precip[j, k]
    arr_avg_yearly_tasavg[i] = yearly_avg_tasavg[j, k]
    i += 1
print("\rProcessing cell %d of %d" % (i, Np))

# %%

print('\rFormating the data in a shapefile...', end=' ')
df = pd.DataFrame(data={'precip': arr_yearly_avg_precip,
                        'tasavg': arr_avg_yearly_tasavg})
crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"
gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
print('\rFormating the data in a shapefile... done')

print('\rSaving to Shapefile...', end=' ')
path_shp_out = ("C:/Users/jsgosselin/HELP/RADEAU2/grid_yearly_meteo2/"
                "grid_yearly_meteo2.shp")
gdf.to_file(path_shp_out)
print('\rSaving to Shapefile... done', end=' ')


# %% Plot of weather averages

# Read meteo grid :

print('Initializing netcdf manager...', end=' ')
path_netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily"
meteo_manager = NetCDFMeteoManager(path_netcdf)
print('done')

# Read HELP grid shapefile

print('Reading HELP grid shapefile...', end=' ')
ref_crs = ("+proj=longlat +ellps=GRS80 +datum=NAD83 "
           "+towgs84=0,0,0,0,0,0,0 +no_defs")
help_shp_fpath = ("C:/Users/jsgosselin/HELP/RADEAU2/grid_helpXYZ/"
                  "grid_helpXYZ.shp")
shp_help = gpd.read_file(help_shp_fpath)
print('done')

print('Fetching the lat/lon values from the shapefile...', end=' ')
lat_help = [p.y for p in shp_help.geometry]
lon_help = [p.x for p in shp_help.geometry]
cid_help = np.array(shp_help['id'])
print('done')

stack_idx = []
for i in range(len(cid_help)):
    progress = (i+1)/len(cid_help)*100
    msg = ("\rProcessing %d/%d (%0.1f%%)" % (i+1, len(cid_help), progress))
    print(msg, end=' ')

    lat_idx, lon_idx = meteo_manager.get_idx_from_latlon(
            lat_help[i], lon_help[i])
    if (lat_idx, lon_idx) in stack_idx:
        continue
    stack_idx.append((lat_idx, lon_idx))


stack_precip = []
stack_tasmax = []
stack_tasmin = []
stack_tasavg = []
for lat_idx, lon_idx in stack_idx:
    stack_precip.append(daily_precip[:, lat_idx, lon_idx])
    stack_tasmax.append(daily_tasmax[:, lat_idx, lon_idx])
    stack_tasmin.append(daily_tasmin[:, lat_idx, lon_idx])
    stack_tasavg.append(daily_tasavg[:, lat_idx, lon_idx])

precip_study = np.mean(np.vstack(stack_precip), axis=0)
tasmax_study = np.mean(np.vstack(stack_tasmax), axis=0)
tasmin_study = np.mean(np.vstack(stack_tasmin), axis=0)
tasavg_study = np.mean(np.vstack(stack_tasavg), axis=0) * 1000
