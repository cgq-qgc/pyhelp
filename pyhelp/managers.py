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
import netCDF4

# ---- Local Libraries Imports

from pyhelp.preprocessing import (write_d10d11_allcells,
                                  format_d10d11_from_excel)
from pyhelp.processing import run_help_allcells
from pyhelp.utils import savedata_to_hdf5
from gwhat.meteo.weather_reader import (
        save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP,
        read_cweeds_file, join_daily_cweeds_wy2_and_wy3)


FNAME_CONN_TABLES = 'connect_table.npy'


class HELPManager(object):
    def __init__(self, path_inputdir, year_range, path_helpgrid=None):
        super(HELPManager, self).__init__()
        self.path_helpgrid = path_helpgrid
        self.path_inputdir = path_inputdir
        if not osp.exists(path_inputdir):
            os.makedirs(path_inputdir)
        self.year_range = year_range

        self.cellnames = []
        self.celllat = []
        self.celllon = []

        if path_helpgrid is not None:
            self.load_helpgrid_from_shapefile(path_helpgrid)
        self._setup_connect_tables()

    @property
    def path_connect_tables(self):
        return osp.join(self.path_inputdir, FNAME_CONN_TABLES)

    def _setup_connect_tables(self):
        """Setup the connect tables dictionary."""
        if osp.exists(self.path_connect_tables):
            self.connect_tables = np.load(self.path_connect_tables).item()
        else:
            self.connect_tables = {}

    def _save_connect_tables(self):
        """Save the connect tables dictionary to a numpy binary file."""
        np.save(self.path_connect_tables, self.connect_tables)

    def load_helpgrid_from_shapefile(self, path_helpgrid):
        """
        Load the shapefile that contains the grid with the location
        coordinates of the centroid of each cell.
        """
        print('Reading HELP grid shapefile...', end=' ')
        self.shp_help = shp_help = gpd.read_file(path_helpgrid)
        print('done')

        if shp_help.crs != {'init': 'epsg:4269'}:
            # Convert the geometry coordinates to lat/long.
            print('Converting HELP grid to lat/lon...')
            crs = ("+proj=longlat +ellps=GRS80 +datum=NAD83 "
                   "+towgs84=0,0,0,0,0,0,0 +no_defs")
            shp_help = shp_help.to_crs(crs)
            shp_help.to_file(path_helpgrid)
            print('done')

        # Store a list of all the cell names, lat, and lon that are going to
        # be run in HELP.
        indx = np.where(shp_help['run'] == 1)[0]
        self.cellnames = shp_help['cid'][shp_help['run'] == 1].tolist()
        self.celllat = np.array([p.y for p in shp_help.geometry])[indx]
        self.celllon = np.array([p.x for p in shp_help.geometry])[indx]

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
        self._save_connect_tables()

    def generate_d10d11_from_excel_grid(self, path_excel_grid):
        """
        Prepare the D10 and D11 input datafiles for each cell from data saved
        in an Excel input grid.
        """
        path_d10d11_input = osp.join(self.path_inputdir, 'd10d11_input_files')
        if not osp.exists(path_d10d11_input):
            os.makedirs(path_d10d11_input)

        d10data, d11data = format_d10d11_from_excel(path_excel_grid)
        d10_conn_tbl, d11_conn_tbl = write_d10d11_allcells(
            path_d10d11_input, d10data, d11data)

        # Update the connection table.
        self.connect_tables['D10'] = d10_conn_tbl
        self.connect_tables['D11'] = d11_conn_tbl
        self._save_connect_tables()

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
            d4fname = osp.join(self.path_inputdir, 'd4d7_input_files',
                               '%03d_%03d.D4' % (lat_idx, lon_idx))
            d7fname = osp.join(self.path_inputdir, 'd4d7_input_files',
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
        self._save_connect_tables()

    def run_help_for(self, path_outfile, cellnames):
        """
        Run help for the cells listed in cellnames and save the result in
        an hdf5 file.
        """
        tempdir = osp.join(self.path_inputdir, ".temp")
        if not osp.exists(tempdir):
            os.makedirs(tempdir)

        if cellnames is None:
            cellnames = self.cellnames

        cellparams = {}
        for cellname in self.cellnames:
            fpath_d4 = self.connect_tables['D4'][cellname]
            fpath_d7 = self.connect_tables['D7'][cellname]
            fpath_d13 = self.connect_tables['D13'][cellname]
            fpath_d10 = self.connect_tables['D10'][cellname]
            fpath_d11 = self.connect_tables['D11'][cellname]
            fpath_out = osp.abspath(osp.join(tempdir, cellname + '.OUT'))

            daily_out = 0
            monthly_out = 1
            yearly_out = 0
            summary_out = 0

            unit_system = 2  # IP if 1 else SI
            simu_nyear = 15

            cellparams[cellname] = (fpath_d4, fpath_d7, fpath_d13, fpath_d11,
                                    fpath_d10, fpath_out, daily_out,
                                    monthly_out, yearly_out, summary_out,
                                    unit_system, simu_nyear)

        # Run HELP and save the result to an hdf5 file.

        output = run_help_allcells(cellparams)
        savedata_to_hdf5(output, path_outfile)


class NetCDFMeteoManager(object):
    def __init__(self, dirpath_netcdf):
        super(NetCDFMeteoManager, self).__init__()
        self.dirpath_netcdf = dirpath_netcdf
        self.lat = []
        self.lon = []
        self.setup_ncfile_list()
        self.setup_latlon_grid()

    def setup_ncfile_list(self):
        """Read all the available netCDF files in dirpath_netcdf."""
        self.ncfilelist = []
        for file in os.listdir(self.dirpath_netcdf):
            if file.endswith('.nc'):
                self.ncfilelist.append(osp.join(self.dirpath_netcdf, file))

    def setup_latlon_grid(self):
        if self.ncfilelist:
            netcdf_dset = netCDF4.Dataset(self.ncfilelist[0], 'r+')
            self.lat = np.array(netcdf_dset['lat'])
            self.lon = np.array(netcdf_dset['lon'])
            netcdf_dset.close()

    def get_idx_from_latlon(self, lat, lon):
        lat_idx = np.argmin(np.abs(self.lat - lat))
        lon_idx = np.argmin(np.abs(self.lon - lon))
        return lat_idx, lon_idx

    def get_data_from_idx(self, lat_idx, lon_idx):
        stack_tasmax = []
        stack_tasmin = []
        stack_precip = []
        stack_years = []
        for year in range(2000, 2015):
            filename = osp.join(self.dirpath_netcdf, 'GCQ_v2_%d.nc' % year)
            netcdf_dset = netCDF4.Dataset(filename, 'r+')
            stack_tasmax.append(
                    np.array(netcdf_dset['tasmax'])[:, lat_idx, lon_idx])
            stack_tasmin.append(
                    np.array(netcdf_dset['tasmin'])[:, lat_idx, lon_idx])
            stack_precip.append(
                    np.array(netcdf_dset['pr'])[:, lat_idx, lon_idx])
            stack_years.append(
                    np.zeros(len(stack_precip[-1])).astype(int) + year)
            netcdf_dset.close()
        tasmax = np.hstack(stack_tasmax)
        tasmin = np.hstack(stack_tasmin)
        precip = np.hstack(stack_precip)
        years = np.hstack(stack_years)
        return (tasmax + tasmin)/2, precip, years
