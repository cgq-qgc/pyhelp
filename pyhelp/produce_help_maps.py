# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 08:33:14 2018
@author: jsgosselin
"""

# ---- Third Party imports

import h5py
import os
import os.path as osp
import pandas as pd
from geopandas import GeoDataFrame
from pyhelp.maps import produce_point_geometry
from pyhelp.managers import load_grid_from_csv
from pyhelp.postprocessing import add_yrly_avg_to_shp

rootdir = "C:\\Users\\User\\pyhelp\\RADEAU2\\inputHELP_0416"

# %% Load the HELP grid

path_togrid = osp.join(rootdir, 'inputHELP_0416base.csv')
grid = load_grid_from_csv(path_togrid)

# %% Init the shapefile

print('\rInitialize the shapefile...', end=' ')
point_geo = produce_point_geometry(
    grid.as_matrix(['lat_dd']), grid.as_matrix(['lon_dd']))
cols_to_keep = ['cid', 'lat_dd', 'lon_dd', 'run', 'context', 'BT', 'Bassin']
crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"

shp = GeoDataFrame(grid[cols_to_keep], crs=crs, geometry=point_geo)
print('done')

# %% Read HELP results

for run_name in ["BTALL_inputHELP_0416base_0.35edepth",
                 "LAUALL_inputHELP_0416t3_0.15edepth"]:
    for oufile in [osp.join(rootdir, run_name, "help_%s.out" % run_name),
                   osp.join(rootdir, run_name, "surface_%s.out" % run_name)]:
        outdat = h5py.File(oufile)
        shp = add_yrly_avg_to_shp(shp, outdat)

# %% Save results in a shapefile

print('\rSaving HELP results in a shapefile...', end=' ')
path_shp_out = osp.join(
    rootdir, "help_results_0424", "help_results_0424.shp")

if not osp.exists(osp.dirname(path_shp_out)):
    os.makedirs(osp.dirname(path_shp_out))
shp.to_file(path_shp_out, driver='ESRI Shapefile')
print('done')
