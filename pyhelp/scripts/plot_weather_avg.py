# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:05:01 2019

@author: User
"""


# %% Plot of weather averages

# Read meteo grid :

print('Initializing netcdf manager...', end=' ')
path_netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily"
meteo_manager = InfoClimatGridReader(path_netcdf)
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