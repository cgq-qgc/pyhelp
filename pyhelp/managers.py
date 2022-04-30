# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHELP Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------

from __future__ import annotations

# ---- Standard Library Imports
import calendar
import itertools
import json
import os
import os.path as osp
import csv
import time

# ---- Third Party imports
import numpy as np
import pandas as pd

# ---- Local Libraries Imports
from pyhelp.preprocessing import write_d10d11_allcells, format_d10d11_inputs
from pyhelp.processing import run_help_allcells
from pyhelp.utils import (savedata_to_hdf5, calc_dist_from_coord,
                          delete_folder_recursively)
from pyhelp.weather_reader import (
    save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP)
from pyhelp.output import HelpOutput


FNAME_CONN_TABLES = 'connect_table.json'


class HelpManager(object):
    """
    The :class:`~pyhelp.HelpManager` is a class whose main purpose
    is to evaluate the component of the hydrologic water budget at the
    regional scale with the HELP model.

    Parameters
    ----------
    workdir : str
        The path to the working directory.
    path_to_grid : str
        The path to the csv file that contains the geomatic data, surface
        conditions, and soil and design data for each cell of the grid
        dividing the study area. The default is None.
    path_to_precip : str
        The path to the csv file that contains the input data for the
        daily total precipitation.
    path_to_airtemp : str
        The path to the csv file that contains the input data for the
        daily average air temperature
    path_to_solrad : str
        The path to the csv file that contains the input data for the
        daily global solar radiation.
    """

    def __init__(self, workdir: str, path_to_grid: str,
                 path_to_precip: str, path_to_airtemp: str,
                 path_to_solrad: str):
        super().__init__()
        self.grid = None
        self.precip_data = None
        self.airtemp_data = None
        self.solrad_data = None

        self.set_workdir(workdir)
        self.load_input_grid(path_to_grid)
        self.load_weather_input_data(
            path_to_precip, path_to_airtemp, path_to_solrad)

        self._load_connect_tables()

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
        self._workdir = dirname
        if not osp.exists(dirname):
            os.makedirs(dirname)
        os.chdir(dirname)

    # ---- Connect tables
    @property
    def path_connect_tables(self):
        return osp.join(self.inputdir, FNAME_CONN_TABLES)

    def _load_connect_tables(self):
        """Load the connect tables from the json file."""
        if osp.exists(self.path_connect_tables):
            with open(self.path_connect_tables, 'r') as jsonfile:
                self.connect_tables = json.load(jsonfile)
        else:
            self.connect_tables = {}

    def _save_connect_tables(self):
        """Save the connect tables to a json text file."""
        with open(self.path_connect_tables, 'w', encoding='utf-8') as jsonfile:
            json.dump(self.connect_tables, jsonfile, indent=2,
                      separators=(",", ": "), ensure_ascii=False)

    # ---- Grid and Input
    def load_input_grid(self, path_to_grid: str):
        """
        Load input grid data.

        Parameters
        ----------
        path_to_grid : str
            The path to the csv file that contains the geomatic data, surface
            conditions, and soil and design data for each cell of the grid
            dividing the study area.
        """
        self.grid_filename = osp.abspath(path_to_grid)
        print(f'Reading grid data from {path_to_grid}...')
        if not osp.exists(path_to_grid):
            self.grid = None
            print("Grid input csv file does not exist.")
        else:
            self.grid = load_grid_from_csv(path_to_grid)
            print('Grid data read successfully from input csv file.')

    def load_weather_input_data(self, path_to_precip: str,
                                path_to_airtemp: str,
                                path_to_solrad: str):
        """
        Load input weather data.

        Parameters
        ----------
        path_to_precip : str
            The path to the csv file that contains the input data for the
            daily total precipitation.
        path_to_airtemp : str
            The path to the csv file that contains the input data for the
            daily average air temperature
        path_to_solrad : str
            The path to the csv file that contains the input data for the
            daily global solar radiation.
        """
        print(f'Reading input precip data from {path_to_precip}...')
        self.precip_data = load_weather_from_csv(path_to_precip)
        self.precip_filename = path_to_precip

        print(f'Reading input airtemp data from {path_to_airtemp}...')
        self.airtemp_data = load_weather_from_csv(path_to_airtemp)
        self.airtemp_filename = path_to_airtemp

        print(f'Reading input solrad data from {path_to_solrad}...')
        self.solrad_data = load_weather_from_csv(path_to_solrad)
        self.solrad_filename = path_to_solrad

        datasets = [self.precip_data, self.airtemp_data, self.solrad_data]
        datasets_name = {
            id(self.precip_data): 'Precipitation',
            id(self.airtemp_data): 'Air temperature',
            id(self.solrad_data): 'Solar radiation'}

        # Check that each year are complete in each dataset.
        for dataset in datasets:
            if dataset is None:
                continue

            years = dataset.index.year
            name = datasets_name[id(dataset)]
            for year in years:
                ndays = np.sum(years == year)
                if ndays != (366 if calendar.isleap(year) else 365):
                    raise ValueError((
                        "{} data for year {} only have {} daily values and "
                        "is not complete"
                        ).format(name, year, ndays))

        # Check that the datasets are synchroneous.
        for dset1, dset2 in itertools.combinations(datasets, 2):
            if dset1 is None or dset2 is None:
                continue

            x1 = dset1.index.values
            x2 = dset2.index.values
            name1 = datasets_name[id(dset1)]
            name2 = datasets_name[id(dset2)]

            # Check that the length of the datasets matches.
            if len(x1) != len(x2):
                raise ValueError((
                    "The lenght of the {} and {} data does not "
                    "match: {} != {}."
                    ).format(name1.lower(), name2.lower(),
                             len(x1), len(x2)))

            # Check that the datetimes of the datasets match.
            match = (x1 == x2)
            if not match.all():
                raise ValueError((
                    "{} and {} data does not match: {} != {}."
                    ).format(name1, name2.lower(),
                             x1[~match][0], x2[~match][0]))
        print('Input weather data read successfully.')

    # ---- HELP input files creation
    def clear_cache(self):
        """Delete all cached HELP input data files from the input folder."""
        print('Clearing HELP input files cache...', end=' ')
        delete_folder_recursively(self.inputdir)
        print('done')

    def build_help_input_files(self, cellnames: list = None,
                               sf_edepth: float = 1, sf_ulai: float = 1,
                               sf_cn: float = 1):
        """
        Clear all cached HELP input data files and generate new ones from the
        weather and grid input data files.

        Parameters
        ----------
        cellnames : list, optional
            The list of cell ids for which D10 and D11 HELP input files
            are to be generated. If None, D10 and D11 HELP input files are
            generated for each cell of the grid with a "run" value of 1.
        sf_edepth : float, optional
            Global scale factor for the Evaporative Zone Depth (applied to
            the whole grid). The default is 1.
        sf_ulai : float, optional
            Global scale factor for the Maximum Leaf Area Index (applied to
            the whole grid). The default is 1.
        sf_cn : float, optional
            Global scale factor for the Curve Number (applied to
            the whole grid). The default is 1.
        """
        self.clear_cache()
        self._generate_d10d11_input_files(
            cellnames, sf_edepth, sf_ulai, sf_cn)
        self._generate_d4d7d13_input_files(
            cellnames)

    def _generate_d10d11_input_files(self, cellnames, sf_edepth,
                                     sf_ulai, sf_cn):
        """
        Prepare the D10 and D11 input datafiles for each cell.

        D10 : Soil and Design data
        D11 : Surface condition (Evapotranspiration)

        See https://github.com/cgq-qgc/pyhelp/wiki/HELP-input-files-format-description
        """
        if sf_edepth < 0:
            raise ValueError("sf_edepth value cannot be lower than 0.")
        if sf_ulai < 0:
            raise ValueError("sf_ulai value cannot be lower than 0.")
        if sf_cn < 0:
            raise ValueError("sf_cn value cannot be lower than 0.")

        d10d11_inputdir = osp.join(self.inputdir, 'd10d11_input_files')
        if not osp.exists(d10d11_inputdir):
            os.makedirs(d10d11_inputdir)

        # Only keep the cells that are going to be run in HELP because we
        # don't need the D10 or D11 input files for those that aren't.
        cellnames = self.get_run_cellnames(cellnames)

        # Apply scaling factors to the grid.
        grid = self.grid.copy()
        grid['EZD'] = grid['EZD'] * sf_edepth
        grid['LAI'] = grid['LAI'] * sf_ulai
        grid['CN'] = grid['CN'] * sf_cn

        # Format the data from the input grid.
        d10data, d11data = format_d10d11_inputs(grid, cellnames)

        # Write the D10 and D11 input files.
        d10_conn_tbl, d11_conn_tbl = write_d10d11_allcells(
            d10d11_inputdir, d10data, d11data)

        # Update the connection table.
        print("\rSaving the connectivity tables...", end=' ')
        self.connect_tables['D10'] = d10_conn_tbl
        self.connect_tables['D11'] = d11_conn_tbl
        self._save_connect_tables()
        print("done")

    def _generate_d4d7d13_input_files(self, cellnames):
        """
        Generate the D4, D7, and D13 HELP input datafiles for each cell.

        D4: Total precipitation.
        D7: Mean air temperature.
        D13: Solar radiation.
        """
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
            data_lat = data.columns.get_level_values('lat_dd').values
            data_lon = data.columns.get_level_values('lon_dd').values
            for i, cellname in enumerate(cellnames):
                dist = calc_dist_from_coord(grid_lat[i], grid_lon[i],
                                            data_lat, data_lon)
                argmin = int(np.argmin(dist))

                lat = data_lat[argmin]
                lon = data_lon[argmin]
                help_input_fname = osp.join(
                    help_inputdir, fformat.format(lat, lon, fext))
                if not osp.exists(help_input_fname):
                    city = '{} at {:3.1f} ; {:3.1f}'.format(var, lat, lon)
                    if var in ('precip', 'airtemp'):
                        to_help_func(help_input_fname,
                                     data.index.year.values,
                                     data.values[:, argmin],
                                     city)
                    elif var == 'solrad':
                        to_help_func(help_input_fname,
                                     data.index.year.values,
                                     data.values[:, argmin],
                                     city,
                                     lat)

                file_conn_tbl[cellname] = help_input_fname
                index_conn_tbl[cellname] = argmin

            self.connect_tables[fext] = file_conn_tbl
            self.connect_tables[var] = index_conn_tbl
            print('done')

        # Update the connectivity table.
        print("\rUpdating the connection table...", end=' ')
        self._save_connect_tables()
        print('done')

    def calc_help_cells(self, path_to_hdf5=None, cellnames=None, tfsoil=0,
                        sf_edepth: float = 1, sf_ulai: float = 1,
                        sf_cn: float = 1,
                        build_help_input_files: bool = True) -> HelpOutput:
        """
        Calcul the water budget for all eligible cells with HELP.

        Run HELP to compute the monthly water budget for the cells listed in
        "cellnames". If a file name is provided in _path_outfile_, the results
        are also saved to disk in a HDF5 file.

        Parameters
        ----------
        path_to_hdf5: str, optional
            File path where to save the results to disk in a HDF5 file.
        tfsoil: float, optional
            The average air temperature, in Celcius degrees,
            below which the soil is assumed to be freezing. The default is 0.
        sf_edepth : float, optional
            Global scale factor for the Evaporative Zone Depth (applied to
            the whole grid). The default is 1.
        sf_ulai : float, optional
            Global scale factor for the Maximum Leaf Area Index (applied to
            the whole grid). The default is 1.
        sf_cn : float, optional
            Global scale factor for the Curve Number (applied to
            the whole grid). The default is 1.
        build_help_input_files: bool
            A flag to indicate whether to generate the basic HELP input
            files before running the simulation.
        """
        if build_help_input_files:
            self.build_help_input_files(cellnames, sf_edepth, sf_ulai, sf_cn)

        # Convert from Celcius to Farenheight
        tfsoil = (tfsoil * 1.8) + 32

        tempdir = osp.join(self.inputdir, ".temp")
        if not osp.exists(tempdir):
            os.makedirs(tempdir)

        run_cellnames = self.get_run_cellnames(cellnames)
        cellparams = {}
        skipped_cells = []
        for cellname in run_cellnames:
            fpath_d4 = self.connect_tables['D4'][cellname]
            fpath_d7 = self.connect_tables['D7'][cellname]
            fpath_d13 = self.connect_tables['D13'][cellname]
            fpath_d10 = self.connect_tables['D10'][cellname]
            fpath_d11 = self.connect_tables['D11'][cellname]
            fpath_out = osp.abspath(osp.join(tempdir, str(cellname) + '.OUT'))

            if fpath_d10 is None or fpath_d11 is None:
                skipped_cells.append(cellname)
                continue

            daily_out = 0
            monthly_out = 1
            yearly_out = 0
            summary_out = 0

            unit_system = 2  # IP if 1 else SI

            year_start = self.precip_data.index.year.min()
            year_end = self.precip_data.index.year.max()
            simu_nyear = year_end - year_start + 1

            cellparams[cellname] = (fpath_d4, fpath_d7, fpath_d13, fpath_d11,
                                    fpath_d10, fpath_out, daily_out,
                                    monthly_out, yearly_out, summary_out,
                                    unit_system, simu_nyear, tfsoil)

        skipped_cells = list(set(skipped_cells))
        if skipped_cells:
            print('-' * 25)
            msg = "Warning: calcul for "
            msg += "cell " if len(skipped_cells) == 1 else "cells "
            msg += ", ".join(skipped_cells)
            msg += " will be skipped due to problems with the input data."
            print(msg)
            print('-' * 25)

        output_data = run_help_allcells(cellparams)
        output_data = self._post_process_output(output_data)
        help_output = HelpOutput(output_data)
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
        to the i and j indexes of the 3D numpy arrays.

        A list with the indexes of the cells where data is nan is also saved
        in the dict. Nan values are obtained when HELP is not able to do water
        budget calculations due to lack of data or when there is an error in
        the input data.
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
        data['lat_dd'] = self.grid.loc[cellnames]['lat_dd'].values
        data['lon_dd'] = self.grid.loc[cellnames]['lon_dd'].values

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
        tstart = time.perf_counter()
        print("Calculating budget for water cells...", end=' ')
        cellnames = self.get_water_cellnames(cellnames)
        lat_dd, lon_dd = self.get_latlon_for_cellnames(cellnames)

        year_start = self.precip_data.index.year.min()
        year_end = self.precip_data.index.year.max()
        year_range = np.arange(year_start, year_end + 1).astype(int)
        nyr = len(year_range)

        output = {}
        years = self.precip_data.index.year.values
        for i, cellname in enumerate(cellnames):
            precip_indx = self.connect_tables['precip'][cellname]
            precip = self.precip_data.values[:, precip_indx]
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
        calcul_time = (time.perf_counter() - tstart)
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

    def get_run_cellnames(self, cellnames=None):
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

        # Only keep the cells that have  a nlayer > 0 because we cannot run
        # those in HELP anyway.
        cellnames = (
            self.grid['cid'][cellnames][self.grid['nlayer'] > 0].tolist())

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


def load_grid_from_csv(path_togrid):
    """
    Load the csv that contains the infos required to evaluate regional
    groundwater recharge with HELP.
    """
    grid = pd.read_csv(path_togrid, dtype={'cid': 'str'})

    fname = osp.basename(path_togrid)
    req_keys = ['cid', 'lat_dd', 'lon_dd', 'run', 'context']
    for key in req_keys:
        if key not in grid.keys():
            raise KeyError("No attribute '{}' found in {}".format(key, fname))

    # Set 'cid' as the index of the dataframe.
    grid.set_index(['cid'], drop=False, inplace=True)

    return grid


def load_weather_from_csv(filename: str) -> pd.DataFrame:
    """
    Load daily precipitation, average air temperature, or global solar
    radiation data from a correctly formatted PyHELP weather input files.

    The data are saved in a pandas dataframe where the index corresponds
    to the dates and the columns to latitudes and longitudes of each
    data series.
    """
    with open(filename, 'r') as csvfile:
        reader = list(csv.reader(csvfile, delimiter=','))

    # Get the latitudes and longitudes and the number of lines that
    # are located above the data block.
    latitudes = None
    longitudes = None
    for i, line in enumerate(reader):
        if not line or not line[0]:
            continue

        if line[0].lower().strip().startswith('lat'):
            latitudes = np.array(line[1:]).astype('float')
        elif line[0].lower().strip().startswith('lon'):
            longitudes = np.array(line[1:]).astype('float')
        elif latitudes is not None and longitudes is not None:
            break

    dataf = pd.read_csv(
        filename,
        header=None,
        index_col=[0],
        skiprows=i,
        skip_blank_lines=False,
        parse_dates=True,
        infer_datetime_format=True,
        dayfirst=True)
    dataf.index.name = 'date'
    dataf.columns = pd.MultiIndex.from_tuples(
        [(lat, lon) for lat, lon in zip(latitudes, longitudes)],
        names=['lat_dd', 'lon_dd'])

    return dataf


if __name__ == '__main__':
    workdir = "C:/Users/User/pyhelp/example"
    helpm = HelpManager(workdir)

    helpm.build_help_input_files()
    path_hdf5 = osp.join(workdir, 'help_example.out')
    output_help = helpm.calc_help_cells(path_hdf5, tfsoil=-3)
    path_hdf5 = osp.join(workdir, 'surf_example.out')
    output_surf = helpm.calc_surf_water_cells(650, path_hdf5)
