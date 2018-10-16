# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Standard Library Imports

import os
import os.path as osp
import csv
from datetime import datetime

# ---- Third Party imports

import numpy as np
import pandas as pd

# ---- Local Libraries Imports

from pyhelp.preprocessing import write_d10d11_allcells, format_d10d11_inputs
from pyhelp.processing import run_help_allcells
from pyhelp.utils import (savedata_to_hdf5, calc_dist_from_coord,
                          delete_folder_recursively)
from pyhelp.weather_reader import (
    save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP,
    read_cweeds_file, join_daily_cweeds_wy2_and_wy3, NetCDFMeteoManager)


FNAME_CONN_TABLES = 'connect_table.npy'


class HelpManager(object):
    """
    The :attr:`~pyhelp.HelpManager` is a class whose main purpose
    is to evaluate the component of the hydrologic water budget at the
    regional scale with the HELP model.
    """

    def __init__(self, workdir, year_range, path_togrid=None):
        super(HelpManager, self).__init__()
        self.year_range = year_range
        self.set_workdir(workdir)
        self._setup_connect_tables()

        if path_togrid is not None:
            self.load_grid(path_togrid)
        else:
            self.grid = None

    @property
    def cellnames(self):
        """Return a list with the ID numbers of all cells in the grid."""
        return [] if self.grid is None else self.grid['cid'].tolist()

    @property
    def inputdir(self):
        """
        Return the path to the folder where the HELP input files are going to
        be saved in the working directory. This folder is created in case it
        doesn't already exist in the file system.
        """
        inputdir = osp.join(self.workdir, 'help_input_files')
        if not osp.exists(inputdir):
            os.makedirs(inputdir)
        return inputdir

    @property
    def workdir(self):
        """Return the path to the current working directory."""
        return os.getcwd()

    def set_workdir(self, dirname):
        """Set the working directory of the manager."""
        if not osp.exists(dirname):
            os.makedirs(dirname)
        os.chdir(dirname)

    # ---- Connect tables

    @property
    def path_connect_tables(self):
        return osp.join(self.inputdir, FNAME_CONN_TABLES)

    def _setup_connect_tables(self):
        """Setup the connect tables dictionary."""
        if osp.exists(self.path_connect_tables):
            self.connect_tables = np.load(self.path_connect_tables).item()
        else:
            self.connect_tables = {}

    def _save_connect_tables(self):
        """Save the connect tables dictionary to a numpy binary file."""
        np.save(self.path_connect_tables, self.connect_tables)

    # ---- HELP grid

    def load_grid(self, path_togrid):
        """
        Load the grid that contains the infos required to evaluate regional
        groundwater recharge with HELP.
        """
        self.grid = load_grid_from_csv(path_togrid)
        return self.grid

    # ---- Input files creation

    def generate_d13_from_cweeds(self, d13fname, fpath_cweed2, fpath_cweed3,
                                 cellnames=None):
        """
        Generate the HELP D13 input file for solar radiation from wy2 and
        wy3 CWEEDS files at a given location.
        """
        d13fpath = osp.join(self.inputdir, d13fname)
        if cellnames is None:
            cellnames = self.cellnames
        else:
            # Keep only the cells that are in the grid.
            cellnames = self.grid['cid'][self.grid['cid'].isin(cellnames)]

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
        print("\rUpdating the connection table...", end=' ')
        d13_connect_table = {cid: d13fpath for cid in cellnames}
        self.connect_tables['D13'] = d13_connect_table
        self._save_connect_tables()
        print("done")

    def generate_d10d11_input_files(self, cellnames=None, sf_edepth=1,
                                    sf_ulai=1):
        """Prepare the D10 and D11 input datafiles for each cell."""
        d10d11_inputdir = osp.join(self.inputdir, 'd10d11_input_files')
        if not osp.exists(d10d11_inputdir):
            os.makedirs(d10d11_inputdir)

        # Only keep the cells that are going to be run in HELP because we
        # don't need the D10 or D11 input files for those that aren't.
        cellnames = self.get_run_cellnames(cellnames)

        d10data, d11data = format_d10d11_inputs(self.grid, cellnames,
                                                sf_edepth, sf_ulai)

        # Write the D10 and D11 input files.
        d10_conn_tbl, d11_conn_tbl = write_d10d11_allcells(
            d10d11_inputdir, d10data, d11data)

        # Update the connection table.
        print("\rUpdating the connection table...", end=' ')
        self.connect_tables['D10'] = d10_conn_tbl
        self.connect_tables['D11'] = d11_conn_tbl
        self._save_connect_tables()
        print("done")

    def generate_d4d7_from_MDELCC_grid(self, path_netcdf_dir, cellnames=None):
        """
        Prepare the D4 and D7 input datafiles for each cell from the
        interpolated grid of the MDDELCC.
        """
        d4d7_inputdir = osp.join(self.inputdir, 'd4d7_input_files')
        if not osp.exists(d4d7_inputdir):
            os.makedirs(d4d7_inputdir)

        cellnames = self.get_run_cellnames(cellnames)
        N = len(cellnames)

        # Get the latitudes and longitudes of the resulting cells.
        lat_dd, lon_dd = self.get_latlon_for_cellnames(cellnames)

        # Generate the connectivity table between the HELP grid and the
        # MDDELCC interpolated daily weather grid.
        print('Generating the connectivity table for each cell...', end=' ')
        meteo_manager = NetCDFMeteoManager(path_netcdf_dir)
        d4_conn_tbl = {}
        d7_conn_tbl = {}
        data = []
        for i, cellname in enumerate(cellnames):
            lat_idx, lon_idx = meteo_manager.get_idx_from_latlon(
                    lat_dd[i], lon_dd[i])

            d4fname = osp.join(
                d4d7_inputdir, '%03d_%03d.D4' % (lat_idx, lon_idx))
            d7fname = osp.join(
                d4d7_inputdir, '%03d_%03d.D7' % (lat_idx, lon_idx))

            d4_conn_tbl[cellnames[i]] = d4fname
            d7_conn_tbl[cellnames[i]] = d7fname

            data.append([lat_idx, lon_idx, d4fname, d7fname])
        print('done')

        # Fetch the daily weather data from the netCDF files.
        data = np.unique(data, axis=0)
        lat_indx = data[:, 0].astype(int)
        lon_idx = data[:, 1].astype(int)
        years = range(self.year_range[0], self.year_range[1]+1)
        tasavg, precip, years = meteo_manager.get_data_from_idx(
            lat_indx, lon_idx, years)

        # Convert and save the weather data to D4 and D7 HELP input files.
        N = len(data)
        for i in range(N):
            print(("\rGenerating HELP D4 and D7 files for location " +
                   "%d of %d (%0.1f%%)...") % (i+1, N, (i+1)/N * 100), end=' ')
            lat = meteo_manager.lat[lat_indx[i]]
            lon = meteo_manager.lon[lon_idx[i]]
            d4fname, d7fname = data[i, 2], data[i, 3]
            city = 'Meteo Grid at lat/lon %0.1f ; %0.1f' % (lat, lon)

            # Fill -999 with 0 in daily precip.
            precip_i = precip[:, i]
            precip_i[precip_i == -999] = 0

            # Fill -999 with linear interpolation in daily air temp.
            tasavg_i = tasavg[:, i]
            time_ = np.arange(len(tasavg_i))
            indx = np.where(tasavg_i != -999)[0]
            tasavg_i = np.interp(time_, time_[indx], tasavg_i[indx])

            if not osp.exists(d4fname):
                save_precip_to_HELP(d4fname, years, precip_i, city)
            if not osp.exists(d7fname):
                save_airtemp_to_HELP(d7fname, years, tasavg_i, city)
        print('done')

        # Update the connection table.
        print("\rUpdating the connection table...", end=' ')
        self.connect_tables['D4'] = d4_conn_tbl
        self.connect_tables['D7'] = d7_conn_tbl
        self._save_connect_tables()
        print('done')

    def run_help_for(self, path_outfile=None, cellnames=None, tfsoil=0):
        """
        Run help for the cells listed in cellnames and save the result in
        an hdf5 file.
        """
        # Convert from Celcius to Farenheight
        tfsoil = (tfsoil * 1.8) + 32

        tempdir = osp.join(self.inputdir, ".temp")
        if not osp.exists(tempdir):
            os.makedirs(tempdir)

        run_cellnames = self.get_run_cellnames(cellnames)
        cellparams = {}
        for cellname in run_cellnames:
            fpath_d4 = self.connect_tables['D4'][cellname]
            fpath_d7 = self.connect_tables['D7'][cellname]
            fpath_d13 = self.connect_tables['D13'][cellname]
            fpath_d10 = self.connect_tables['D10'][cellname]
            fpath_d11 = self.connect_tables['D11'][cellname]
            fpath_out = osp.abspath(osp.join(tempdir, str(cellname) + '.OUT'))

            daily_out = 0
            monthly_out = 1
            yearly_out = 0
            summary_out = 0

            unit_system = 2  # IP if 1 else SI
            simu_nyear = self.year_range[1] - self.year_range[0] + 1

            cellparams[cellname] = (fpath_d4, fpath_d7, fpath_d13, fpath_d11,
                                    fpath_d10, fpath_out, daily_out,
                                    monthly_out, yearly_out, summary_out,
                                    unit_system, simu_nyear, tfsoil)

        output = run_help_allcells(cellparams)

        if path_outfile:
            savedata_to_hdf5(output, path_outfile)

        return output

    def calc_surf_water_cells(self, evp_surf, path_netcdf_dir,
                              path_outfile=None, cellnames=None):
        cellnames = self.get_water_cellnames(cellnames)
        lat_dd, lon_dd = self.get_latlon_for_cellnames(cellnames)

        meteo_manager = NetCDFMeteoManager(path_netcdf_dir)

        N = len(cellnames)
        lat_indx = np.empty(N).astype(int)
        lon_indx = np.empty(N).astype(int)
        for i, cellname in enumerate(cellnames):
            lat_indx[i], lon_indx[i] = meteo_manager.get_idx_from_latlon(
                lat_dd[i], lon_dd[i])

        year_range = np.arange(
            self.year_range[0], self.year_range[1] + 1).astype(int)
        tasavg, precip, years = meteo_manager.get_data_from_idx(
            lat_indx, lon_indx, year_range)

        # Fill -999 with 0 in daily precip.
        precip[precip == -999] = 0

        nyr = len(year_range)
        output = {}
        for i, cellname in enumerate(cellnames):
            data = {}
            data['years'] = year_range
            data['rain'] = np.zeros(nyr)
            data['evapo'] = np.zeros(nyr) + evp_surf
            data['runoff'] = np.zeros(nyr)
            for k, year in enumerate(year_range):
                indx = np.where(years == year)[0]
                data['rain'][k] = np.sum(precip[indx, i])
                data['runoff'][k] = data['rain'][k] - evp_surf
            output[cellname] = data

        if path_outfile:
            savedata_to_hdf5(output, path_outfile)

        return output

        # # For cells for which the context is 2, convert recharge and deep
        # # subrunoff into superfical subrunoff.
        # cellnames_con_2 = cellnames[self.grid[fcon] == 2].tolist()
        # for cellname in cellnames_con_2:
        #     output[cellname]['subrun1'] += output[cellname]['subrun2']
        #     output[cellname]['subrun1'] += output[cellname]['recharge']
        #     output[cellname]['subrun2'][:] = 0
        #     output[cellname]['recharge'][:] = 0

        # # For cells for which the context is 3, convert recharge into
        # # deep runoff.
        # cellnames_con_3 = cellnames[self.grid[fcon] == 3].tolist()
        # for cellname in cellnames_con_3:
        #     output[cellname]['subrun2'] += output[cellname]['recharge']
        #     output[cellname]['recharge'][:] = 0

        # # Comput water budget for cells for which the context is 0.
        # cellnames_con_2 = cellnames[self.grid[fcon] == 0].tolist()

        # # meteo_manager = NetCDFMeteoManager(path_netcdf_dir)
        # # for cellname in cellnames_run0:

        # Save the result to an hdf5 file.

    # ---- Utilities

    def get_water_cellnames(self, cellnames):
        """
        Take a list of cellnames and return only those that are considered
        to be in a surface water area.
        """
        if cellnames is None:
            cellnames = self.cellnames
        else:
            # Keep only the cells that are in the grid.
            cellnames = self.grid['cid'][self.grid['cid'].isin(cellnames)]

        # Only keep the cells for which context is 0.
        cellnames = self.grid['cid'][cellnames][self.grid['context'] == 0]

        return cellnames.tolist()

    def get_run_cellnames(self, cellnames):
        """
        Take a list of cellnames and return only those that are in the grid
        and for which HELP can be run.
        """
        if cellnames is None:
            cellnames = self.cellnames
        else:
            # Keep only the cells that are in the grid.
            cellnames = self.grid['cid'][self.grid['cid'].isin(cellnames)]

        # Only keep the cells that are going to be run in HELP because we
        # don't need the D4 or D7 input files for those that aren't.
        cellnames = self.grid['cid'][cellnames][self.grid['run'] == 1].tolist()

        return cellnames

    def get_latlon_for_cellnames(self, cells):
        """
        Return a numpy array with latitudes and longitudes of the provided
        cells cid. Latitude and longitude for cids that are missing from
        the grid are set to nan.
        """
        lat = np.array(self.grid['lat_dd'].reindex(cells).tolist())
        lon = np.array(self.grid['lon_dd'].reindex(cells).tolist())
        return lat, lon


        """




def load_grid_from_csv(path_togrid):
    """
    Load the csv that contains the infos required to evaluate regional
    groundwater recharge with HELP.
    """
    if not osp.exists(path_togrid):
        return None

    print('Reading HELP grid from csv...', end=' ')
    grid = pd.read_csv(path_togrid)
    print('done')

    fname = osp.basename(path_togrid)
    req_keys = ['cid', 'lat_dd', 'lon_dd', 'run']
    for key in req_keys:
        if key not in grid.keys():
            raise KeyError("No attribute '%s' found in %s" % (key, fname))

    # Make sure that cid is a str.
    grid['cid'] = np.array(grid['cid']).astype(str)

    # Set 'cid' as the index of the dataframe.
    grid.set_index(['cid'], drop=False, inplace=True)

    return grid


def load_weather_from_csv(filename):
    """
    Load daily precipitation, average air temperature, or global solar
    radiation data from a correctly formatted PyHelp weather input files.
    The latitudes, longitudes, dates, and weather data values are stored in
    numpy arrays and returned as a dict with, respectively, the keys
    'lat', 'lon', 'dates', and 'data'.
    """
    if not osp.exists(filename):
        return None

    lat, lon, datestrings, data = [], [], [], []
    with open(filename, 'r') as csvfile:
        reader = list(csv.reader(csvfile, delimiter=','))

    for i, line in enumerate(reader):
        if not line or not line[0]:
            continue

        if line[0] == 'Latitude (dd)':
            lat = np.array(line[1:]).astype('float')
        elif line[0] == 'Longitude (dd)':
            lon = np.array(line[1:]).astype('float')
        elif all((len(lat), len(lon))):
            date_data = np.array(reader[i:])
            datestrings = date_data[:, 0]
            data = date_data[:, 1:].astype('float')
            break

    datetimes = [datetime.strptime(ds, "%d/%m/%Y") for ds in datestrings]
    years = [dt.year for dt in datetimes]

    if all((len(lat), len(lon), len(datestrings), len(data))):
        return {'lat': lat, 'lon': lon, 'datestrings': datestrings,
                'datetimes': datetimes, 'years': years, 'data': data}
    else:
        print("Failed to read data from {}.".format(osp.basename(filename)))
        return None

