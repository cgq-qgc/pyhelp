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
    read_cweeds_file, join_daily_cweeds_wy2_and_wy3, NetCDFMeteoManager,
    generate_input_from_cweeds)


FNAME_CONN_TABLES = 'connect_table.npy'
INPUT_PRECIP_FNAME = 'precip_input_data.csv'
INPUT_AIRTEMP_FNAME = 'airtemp_input_data.csv'
INPUT_SOLRAD_FNAME = 'solrad_input_data.csv'
INPUT_GRID_FNAME = 'input_grid.csv'


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

        self.grid = None
        self.precip_data = None
        self.airtemp_data = None
        self.solrad_data = None

        self.load_input_grid()
        self.load_weather_input_data()

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

    # ---- Grid and Input

    def load_input_grid(self):
        """
        Load the grid that contains the infos required to evaluate regional
        groundwater recharge with HELP.
        """
        grid_fname = osp.join(self.workdir, INPUT_GRID_FNAME)
        self.grid = load_grid_from_csv(grid_fname)
        return self.grid

    def load_weather_input_data(self):
        """
        Load the daily precipitation, average air temperature, and global
        solar radiation from the PyHelp weather inpy datafiles if they exists
        in the working directory.
        """
        self.precip_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_PRECIP_FNAME))
        self.airtemp_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_AIRTEMP_FNAME))
        self.solrad_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_SOLRAD_FNAME))
        return self.precip_data, self.airtemp_data, self.solrad_data

    # ---- HELP input files creation
    def clear_cache(self):
        """Delete all HELP input data files from the input folder."""
        delete_folder_recursively(self.inputdir)

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

    def _generate_d4d7d13_input_files(self, cellnames=None):
        """Generate the D4, D7, and D13 HELP input datafiles for each cell."""
        if self.grid is None:
            return

        cellnames = self.cellnames if cellnames is None else cellnames
        grid_lat, grid_lon = self.get_latlon_for_cellnames(cellnames)

        fformat = '{:3.1f}_{:3.1f}.{}'
        args = (('precip', 'D4', save_precip_to_HELP, self.precip_data),
                ('airtemp', 'D7', save_airtemp_to_HELP, self.airtemp_data),
                ('solrad', 'D13', save_solrad_to_HELP, self.solrad_data))
        for var, fext, to_help_func, data in args:
            print('Generating the connectivity table for {}...'.format(
                  var.lower()), end=' ')

            if data is None:
                print('failed')
                continue

            help_inputdir = osp.join(self.inputdir, fext + '_input_files')
            if not osp.exists(help_inputdir):
                os.makedirs(help_inputdir)

            file_conn_tbl = {}
            index_conn_tbl = {}
            for i, cellname in enumerate(cellnames):
                dist = calc_dist_from_coord(grid_lat[i], grid_lon[i],
                                            data['lat'], data['lon'])
                argmin = np.argmin(dist)

                lat, lon = data['lat'][argmin], data['lon'][argmin]
                help_input_fname = osp.join(help_inputdir,
                                            fformat.format(lat, lon, fext))
                if not osp.exists(help_input_fname):
                    city = '{} at {:3.1f} ; {:3.1f}'.format(var, lat, lon)
                    if var in ('precip', 'airtemp'):
                        to_help_func(help_input_fname, data['years'],
                                     data['data'][:, argmin], city)
                    elif var == 'solrad':
                        to_help_func(help_input_fname, data['years'],
                                     data['data'][:, argmin], city, lat)

                file_conn_tbl[cellname] = help_input_fname
                index_conn_tbl[cellname] = argmin

            self.connect_tables[fext] = file_conn_tbl
            self.connect_tables[var] = index_conn_tbl
            print('done')

        # Update the connectivity table.
        print("\rUpdating the connection table...", end=' ')
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

    # ---- Grid Utilities

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

    # ---- Input Data Utilities
    def generate_weather_inputs_from_CWEEDS(
            self, cweed2_paths, cweed3_paths, year_range=None):
        year_range = self.year_range if year_range is None else year_range
        generate_input_from_cweeds(cweed2_paths, cweed3_paths, year_range)

    def generate_weather_inputs_from_MDELCC_grid(
            self, path_to_mddelcc_grid, cellnames=None, year_range=None):
        """
        Generate weather input data files from the MDDELCC grid.

        Generate PyHelp csv data file inputs for daily precipitation and
        average air temperature using data from the MDDELCC spatially
        distributed daily precipitation and minimum and maximum air
        temperature grid for a set of lat/lon coordinates.
        """
        cellnames = self.cellnames if cellnames is None else cellnames
        year_range = self.year_range if year_range is None else year_range
        lat_dd, lon_dd = self.get_latlon_for_cellnames(cellnames)

        mddelcc_grid_mngr = NetCDFMeteoManager(path_to_mddelcc_grid)
        mddelcc_grid_mngr.generate_input_from_MDELCC_grid(
            self.workdir, lat_dd, lon_dd, year_range)


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


if __name__ == '__main__':
    workdir = "C:/Users/User/pyhelp/example"
    helpm = HelpManager(workdir, year_range=(2010, 2014))
    precip_data = helpm.precip_data
    airtemp_data = helpm.airtemp_data
    helpm.clear_cache()
    helpm._generate_d4d7d13_input_files()
