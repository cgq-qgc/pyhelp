# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Standard Library Imports

import os
import os.path as osp
import csv
from datetime import datetime
import time

# ---- Third Party imports

import numpy as np
import pandas as pd
import h5py

# ---- Local Libraries Imports

from pyhelp.preprocessing import write_d10d11_allcells, format_d10d11_inputs
from pyhelp.processing import run_help_allcells
from pyhelp.utils import (savedata_to_hdf5, calc_dist_from_coord,
                          delete_folder_recursively)
from pyhelp.weather_reader import (
    save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP,
    InfoClimatGridReader, generate_input_from_cweeds)
from pyhelp.output import HelpOutput


FNAME_CONN_TABLES = 'connect_table.npy'
INPUT_PRECIP_FNAME = 'precip_input_data.csv'
INPUT_AIRTEMP_FNAME = 'airtemp_input_data.csv'
INPUT_SOLRAD_FNAME = 'solrad_input_data.csv'
INPUT_GRID_FNAME = 'input_grid.csv'


class HelpManager(object):
    """
    The :class:`~pyhelp.HelpManager` is a class whose main purpose
    is to evaluate the component of the hydrologic water budget at the
    regional scale with the HELP model.
    """

    def __init__(self, workdir, year_range):
        super(HelpManager, self).__init__()
        self.grid = None
        self.precip_data = None
        self.airtemp_data = None
        self.solrad_data = None

        self._workdir = None
        self.set_workdir(workdir)

        self.year_range = year_range
        self._setup_connect_tables()

    @property
    def cellnames(self):
        """Return a list with the ID numbers of all cells in the grid."""
        return [] if self.grid is None else self.grid['cid'].tolist()

    # ---- Work and input dir
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
        return self._workdir

    def set_workdir(self, dirname):
        """Set the working directory of the manager."""
        if self.workdir is not None and osp.samefile(self.workdir, dirname):
            return
        if not osp.exists(dirname):
            os.makedirs(dirname)
        os.chdir(dirname)
        self._workdir = dirname
        self.load_input_grid()
        self.load_weather_input_data()

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
        Load input grid data.

        Load the grid containing the geomatic data, surface conditions, and
        soil and design data for each cell of the grid dividing the study area.
        By default, the input grid data file must be saved in the working
        directory and named :file:`input_grid.csv`.
        """
        print('Reading input data grid from csv...', end=' ')
        grid_fname = osp.join(self.workdir, INPUT_GRID_FNAME)
        self.grid = load_grid_from_csv(grid_fname)
        print('done')

    def load_weather_input_data(self):
        """
        Load input weather data.

        Load the daily precipitation, average air temperature, and
        global solar radiation from the input weather data files. By default,
        those files must be saved in the working directory and named,
        respectively, :file:`precip_input_data.csv`,
        :file:`airtemp_input_data.csv`, and :file:`solrad_input_data.csv`.
        """
        print('Reading input weather data from csv...', end=' ')
        self.precip_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_PRECIP_FNAME))
        self.airtemp_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_AIRTEMP_FNAME))
        self.solrad_data = load_weather_from_csv(
            osp.join(self.workdir, INPUT_SOLRAD_FNAME))
        print('done')

    # ---- HELP input files creation
    def clear_cache(self):
        """Delete all cached HELP input data files from the input folder."""
        print('Clearing HELP input files cache...', end=' ')
        delete_folder_recursively(self.inputdir)
        print('done')

    def build_help_input_files(self, sf_edepth=1, sf_ulai=1):
        """
        Clear all cached HELP input data files and generate new ones from the
        weather and grid input data files.
        """
        self.clear_cache()
        self._generate_d10d11_input_files(sf_edepth=sf_edepth,
                                          sf_ulai=sf_ulai)
        self._generate_d4d7d13_input_files()

    def _generate_d10d11_input_files(self, cellnames=None, sf_edepth=1,
                                     sf_ulai=1):
        """Prepare the D10 and D11 input datafiles for each cell."""
        d10d11_inputdir = osp.join(self.inputdir, 'd10d11_input_files')
        if not osp.exists(d10d11_inputdir):
            os.makedirs(d10d11_inputdir)

        # Only keep the cells that are going to be run in HELP because we
        # don't need the D10 or D11 input files for those that aren't.
        cellnames = self.get_run_cellnames(cellnames)

        # Format the data from the input grid.
        d10data, d11data = format_d10d11_inputs(
            self.grid, cellnames, sf_edepth, sf_ulai)

        # Write the D10 and D11 input files.
        d10_conn_tbl, d11_conn_tbl = write_d10d11_allcells(
            d10d11_inputdir, d10data, d11data)

        # Update the connection table.
        print("\rSaving the connectivity tables...", end=' ')
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
            print('Generating {} HELP input files for {}...'.format(
                  fext, var.lower()), end=' ')

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

    def calc_help_cells(self, path_to_hdf5=None, cellnames=None, tfsoil=0):
        """
        Calcul the water budget for all eligible cells with HELP.

        Run HELP to compute the monthly water budget for the cells listed in
        _cellnames_. Return a dict containing the resulting monthly values as
        numpy arrays. If a file name is provided in _path_outfile_, the results
        are also saved to disk in a HDF5 file.
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

        output_data = run_help_allcells(cellparams)
        output_data = self._post_process_output(output_data)
        output_grid = self.grid.loc[output_data['cid']]
        help_output = HelpOutput({'data': output_data, 'grid': output_grid})
        if path_to_hdf5:
            help_output.save_to_hdf5(path_to_hdf5)

        return help_output

    def _post_process_output(self, output):
        """
        Format and post-process HELP outputs.

        For each component of the water budget, take the data saved in a dict
        individually for each cell and generate a single 3D numpy matrix,
        where the indices i, j, and k corresponds, respectively, to the cells,
        the years, and the months.

        Also, in the process, convert groundwater recharge to subsurface
        runoff for cells that are identified as being located next to a stream.

        Return a dict that contains a 3D numpy array with the monthly values
        for precip, runoff, evapo, perco, subrun1, subrun2, and rechg. The
        dict also contains lists for the cid (cell id) and years corresponding
        to the i and j indexes of the 3D numpy arrays. A list with the
        i indexes where data is nan is also saved in the dict. Nan values
        are obtained when HELP is not able to do water budget calculations due
        to lack of data or when there is an error in the input data.
        """
        cellnames = list(output.keys())
        context = self.grid['context'][cellnames].tolist()
        years = output[cellnames[0]]['years']
        Np = len(cellnames)
        Ny = len(years)

        keys = ['precip', 'runoff', 'evapo', 'perco',
                'subrun1', 'subrun2', 'rechg']
        data = {key: np.zeros((Np, Ny, 12)) for key in keys}
        data['cid'] = cellnames
        data['years'] = years
        data['idx_nan'] = []

        for i, cellname in enumerate(cellnames):
            print("\rPost-processing cell %d of %d..." % (i+1, Np), end=' ')
            for key in keys[:-1]:
                data[key][i, :, :] = output[cellname][key]

            if np.any(np.isnan(output[cellname]['rechg'])):
                data['idx_nan'].append(i)
                data['rechg'][i, :, :] = output[cellname]['rechg']
            elif context[i] == 2:
                # Redistribute recharge as subsurface runoff if cell is
                # located next to a stream.
                if np.all(output[cellname]['subrun2'] == 0):
                    # Convert recharge as surficial subsurface runoff.
                    data['subrun1'][i, :, :] += output[cellname]['rechg']
                else:
                    # Convert recharge as deep subsurface runoff.
                    data['subrun2'][i, :, :] += output[cellname]['rechg']
            else:
                data['rechg'][i, :, :] = output[cellname]['rechg']
        print("done")

        return data

    def calc_surf_water_cells(self, evp_surf, path_outfile=None,
                              cellnames=None):
        """
        Calcul the yearly water budget for cells that are located in
        surface water bodies.
        """
        tstart = time.clock()
        print("Calculating budget for water cells...", end=' ')
        cellnames = self.get_water_cellnames(cellnames)
        lat_dd, lon_dd = self.get_latlon_for_cellnames(cellnames)

        year_range = np.arange(
            self.year_range[0], self.year_range[1] + 1).astype(int)
        nyr = len(year_range)

        output = {}
        years = self.precip_data['years']
        for i, cellname in enumerate(cellnames):
            precip_indx = self.connect_tables['precip'][cellname]
            precip = self.precip_data['data'][:, precip_indx]
            data = {}
            data['years'] = year_range
            data['rain'] = np.zeros(nyr)
            data['evapo'] = np.zeros(nyr) + evp_surf
            data['runoff'] = np.zeros(nyr)
            for k, year in enumerate(year_range):
                indx = np.where(years == year)[0]
                data['rain'][k] = np.sum(precip[indx])
                data['runoff'][k] = data['rain'][k] - evp_surf
            output[cellname] = data

        if path_outfile:
            savedata_to_hdf5(output, path_outfile)
        calcul_time = (time.clock() - tstart)
        print("done")
        print('Task completed in %0.2f sec' % calcul_time)
        return output

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
        """
        Generate global solar irradiance input data file from CWEEDS files.
        """
        year_range = self.year_range if year_range is None else year_range
        generate_input_from_cweeds(self.workdir, cweed2_paths,
                                   cweed3_paths, year_range)

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

        mddelcc_grid_mngr = InfoClimatGridReader(path_to_mddelcc_grid)
        mddelcc_grid_mngr.generate_input_from_MDELCC_grid(
            self.workdir, lat_dd, lon_dd, year_range)


def load_grid_from_csv(path_togrid):
    """
    Load the csv that contains the infos required to evaluate regional
    groundwater recharge with HELP.
    """
    if not osp.exists(path_togrid):
        return None
    grid = pd.read_csv(path_togrid)

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

    cweed2_paths = osp.join(workdir, 'CWEEDS', '94792.WY2')
    cweed3_paths = osp.join(
        workdir, 'CWEEDS', 'CAN_QC_MONTREAL-INTL-A_7025251_CWEEDS2011_T_N.WY3')
    helpm.generate_weather_inputs_from_CWEEDS(cweed2_paths, cweed3_paths)

    path_to_mddelcc_grid = "F:/MeteoGrilleDaily"
    helpm.generate_weather_inputs_from_MDELCC_grid(path_to_mddelcc_grid)

    helpm.build_help_input_files()
    path_hdf5 = osp.join(workdir, 'help_example.out')
    output_help = helpm.calc_help_cells(path_hdf5, tfsoil=-3)
    path_hdf5 = osp.join(workdir, 'surf_example.out')
    output_surf = helpm.calc_surf_water_cells(650, path_hdf5)

    # precip_data = helpm.precip_data
    # airtemp_data = helpm.airtemp_data
