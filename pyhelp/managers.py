# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Standard Library Imports

import os
import os.path as osp

# ---- Third Party imports

import numpy as np
import geopandas as gpd

# ---- Local Libraries Imports

from pyhelp.preprocess import read_d10d11_file, write_d10d11_allcells
from pyhelp.help_utils import NetCDFMeteoManager
from gwhat.meteo.weather_reader import (
        save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP,
        read_cweeds_file, join_daily_cweeds_wy2_and_wy3)


class HELPInputManager(object):
    def __init__(self, path_inputdir, year_range, path_helpgrid=None):
        super(HELPInputManager, self).__init__()
        self.path_helpgrid = path_helpgrid
        self.path_inputdir = path_inputdir
        if not osp.exists(path_inputdir):
            os.makedirs(path_inputdir)
        self.year_range = year_range

        self.connect_tables = {}
        self.cellnames = []
        self.celllat = []
        self.celllon = []

        if path_helpgrid is not None:
            self.load_helpgrid_from_shapefile(path_helpgrid)

    def load_helpgrid_from_shapefile(self, path_helpgrid):
        """
        Load the shapefile that contains the grid with the location
        coordinates of the centroid of each cell.
        """
        print('Reading HELP grid shapefile...', end=' ')
        self.shp_help = shp_help = gpd.read_file(path_helpgrid)
        print('done')

        if shp_help.crs != {'init': 'epsg:4269'}:
            # Convert the geometry coordinates to lat long.
            print('Converting HELP grid to lat/lon...')
            crs = ("+proj=longlat +ellps=GRS80 +datum=NAD83 "
                   "+towgs84=0,0,0,0,0,0,0 +no_defs")
            shp_help = shp_help.to_crs(crs)
            shp_help.to_file(path_helpgrid)
            print('done')

        # Store a list of all the cell ids that are going to be run in HELP.
        self.cellnames = shp_help['cid'][shp_help['run'] == 1]

    def generate_D13_from_cweeds(self, d13fname, fpath_cweed2, fpath_cweed3):
        """
        Generate the HELP D13 input file for solar radiation from wy2 and
        wy3 CWEEDS files at a given location.
        """
        d13fpath = osp.join(self.path_inputdir, d13fname)
        print('Reading CWEEDS files...', end=' ')
        daily_wy2 = read_cweeds_file(fpath_cweed2, format_to_daily=True)
        daily_wy3 = read_cweeds_file(fpath_cweed3, format_to_daily=True)
        wy23_df = join_daily_cweeds_wy2_and_wy3(daily_wy2, daily_wy3)

        indexes = np.where((wy23_df['Years'] >= self.year_range[0]) &
                           (wy23_df['Years'] <= self.year_range[1]))[0]
        print('done')

        print('Generating HELP D13 file for solar radiation...', end=' ')
        save_solrad_to_HELP(d13fpath,
                            wy23_df['Years'][indexes],
                            wy23_df['Irradiance'][indexes],
                            'CAN_QC_MONTREAL-INTL-A_7025251',
                            wy23_df['Latitude'])
        print('done')

        if self.year_range[1] > np.max(wy23_df['Years']):
            print("Warning: there is no solar radiation data after year %d."
                  % np.max(wy23_df['Years']))
        if self.year_range[0] < np.min(wy23_df['Years']):
            print("Warning: there is no solar radiation data before year %d."
                  % np.min(wy23_df['Years']))

        # Update the connection table.
        d13_connect_table = {cid: d13fpath for cid in self.cellnames}
        self.connect_tables['D13'] = d13_connect_table

    def generate_d10d11_from_lcnp_files(self, path_d10file, path_d11file):
        """
        Prepare the D10 and D11 input datafiles for each cell.
        """
        print('Reading LCNP D10 and D11 files...', end=' ')
        d10data, d11data = read_d10d11_file(path_d10file, path_d11file)
        print('done')

        d10_conn_tbl, d11_conn_tbl = write_d10d11_allcells(
            self.path_inputdir, d10data, d11data)

        # Update the connection table.
        self.connect_tables['D10'] = d10_conn_tbl
        self.connect_tables['D11'] = d11_conn_tbl

    def generate_d4d7_from_MDELCC_grid(self, path_netcdf_dir):
        meteo_manager = NetCDFMeteoManager(path_netcdf_dir)
        d4_conn_tbl = {}
        d7_conn_tbl = {}
        N = len(self.cellnames)
        for i in range(N):
            progress = (i+1)/N * 100
            msg = ("\rGenerating HELP D4 and D7 files for cell "
                   "%d of %d (%0.1f%%)" % (i+1, N, progress))
            print(msg, end=' ')

            lat_idx, lon_idx = meteo_manager.get_idx_from_latlon(
                    self.celllat[i], self.celllon[i])
            d4fname = osp.join(self.path_inputdir,
                               '%03d_%03d.D4' % (lat_idx, lon_idx))
            d7fname = osp.join(self.path_inputdir,
                               '%03d_%03d.D7' % (lat_idx, lon_idx))

            d4_conn_tbl[self.cellnames[i]] = d4fname
            d7_conn_tbl[self.cellnames[i]] = d7fname
            if osp.exists(d4fname):
                # The help meteo files D4 and D7 already exist.
                continue
            else:
                city = 'Meteo Grid at lat/lon %0.1f ; %0.1f' % (
                        meteo_manager.lat[lat_idx], meteo_manager.lon[lon_idx])
                tasavg, precip, years = meteo_manager.get_data_from_idx(
                        lat_idx, lon_idx)
                save_precip_to_HELP(d4fname, years, precip, city)
                save_airtemp_to_HELP(d7fname, years, tasavg, city)

        # Update the connection table.
        self.connect_tables['D4'] = d4_conn_tbl
        self.connect_tables['D7'] = d7_conn_tbl