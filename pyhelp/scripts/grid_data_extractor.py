# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 14:29:22 2021
@author: Jean-Sébastien Gosselin
"""

# ---- Standard Library imports
from __future__ import annotations
import functools
import time
import os
import os.path as osp
from pathlib import Path
import datetime
import re

# ---- Third Party imports
import numpy as np
import pandas as pd
import netCDF4


def timethis(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print('runtime: {} -> {:0.5} sec'.format(func.__name__, end-start))
        return result
    return wrapper


def calc_dist_from_coord(lat1, lon1, lat2, lon2):
    """
    Compute the  horizontal distance in km between a location given in
    decimal degrees and a set of locations also given in decimal degrees.

    https://en.wikipedia.org/wiki/Haversine_formula
    https://www.nhc.noaa.gov/gccalc.shtml
    """
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)

    r = 6373  # Earth radius in km.

    # Note that the units used for r determine the return value units.

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    # Note that arctan2(sqrt(a), sqrt(1-a)) is the same as arcsin(sqrt(a)) in
    # this case.

    return r * c


class ConnectTable(pd.DataFrame):
    # https://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas

    @property
    def _constructor(self):
        return ConnectTable

    def save_to_csv(self, filepath):
        """
        Parameters
        ----------
        filepath : str | Path
            The path of the csv file where to save this connect table.
        """
        filepath = osp.abspath(filepath)
        if not osp.exists(osp.dirname(filepath)):
            os.makedirs(osp.dirname(filepath))

        self.to_csv(filepath, sep=',', line_terminator='\n', encoding='utf-8')


class ExtractedDataFrame(pd.DataFrame):
    # https://pandas.pydata.org/pandas-docs/stable/development/extending.html#extending-subclassing-pandas
    _metadata = ["varname", "source"]

    @property
    def _constructor(self):
        return ExtractedDataFrame

    def save_to_csv(self, filepath: str | Path):
        """
        Parameters
        ----------
        filepath : str | Path
            The path of the csv file where to save this climate dataframe.
        """
        filepath = osp.abspath(filepath)
        if not osp.exists(osp.dirname(filepath)):
            os.makedirs(osp.dirname(filepath))

        print("Saving '{}' data to '{}'...".format(self.varname, filepath))

        now_str = datetime.datetime.now().strftime('%Y-%m-%d')
        index_strlist = self.columns.get_level_values(0).astype(str).tolist()
        lat_strlist = self.columns.get_level_values(1).astype(str).tolist()
        lon_strlist = self.columns.get_level_values(2).astype(str).tolist()

        # Define the content of the file header.
        fheader = [
            [self.varname],
            [''],
            ['Created on:', now_str],
            ['Source:', self.source],
            [''],
            ['Grid index:'] + index_strlist,
            ['Latitude (dd):'] + lat_strlist,
            ['Longitude (dd):'] + lon_strlist,
            [''],
            ['']]
        fheader = '\n'.join([','.join(row) for row in fheader])

        # Define the content of the file data.
        fdata = self.to_csv(
            header=None, sep=',', line_terminator='\n', encoding='utf-8')

        # Save the content of the file buffer to disk.
        with open(filepath, mode='w', encoding='utf-8') as f:
            f.write(fheader + fdata)


class GridExtractor(object):
    """
    A base class to read and extract daily time series from a grid.

    Parameters
    ----------
    gridpath : str | Path
        The path of the directory where the grid netcdf files are stored.
    filename_pattern : str
        The filename pattern of the yearly netcdf files.
        Example : solar_radiation_flux_{year}.nc
    varname : str
        The name of the variable for which data is stored in the grid.
    """

    def __init__(self, gridpath: str | Path, filename_pattern: str,
                 varname: str):
        super().__init__()
        if not osp.exists(gridpath):
            raise ValueError(f'{gridpath} does not exist.')
        self.gridpath = gridpath
        self.filename_pattern = filename_pattern
        self.varname = varname

        self._grid_latdd = np.array([])
        self._grid_londd = np.array([])
        self._nan_grid_mask = np.array([])

        self.setup_grid()

    def setup_grid(self):
        """
        Setup the latitude and longitude coordinates of the grid cells.
        """
        # The grid need to be the same for all netcdf files, so we
        # simply load the grid lat and lon from the first netcdf file
        # found 'gridpath'.
        for file in os.listdir(self.gridpath):
            pattern = self.filename_pattern.format(year='(.*)')
            if re.match(pattern, file):
                break
        else:
            raise ValueError(f"There is no netcdf file in {self.gridpath}.")

        with netCDF4.Dataset(osp.join(self.gridpath, file), 'r+') as ncdset:
            lat = np.array(ncdset['lat'])
            lon = np.array(ncdset['lon'])

            # Define a mask to identify the cells of the grid with no data.
            array = np.array(ncdset[self.varname])
            array = array.reshape(array.shape[0], -1)
            array[array <= -999] = np.nan
            self._nan_grid_mask = np.isnan(array).all(axis=0)

        if len(np.shape(lat)) == 1 and len(np.shape(lon)) == 1:
            self._grid_latdd = np.repeat(lat, len(lon)).flatten()
            self._grid_londd = np.tile(lon, len(lat)).flatten()
        elif len(np.shape(lat)) == 2 and len(np.shape(lon)) == 2:
            self._grid_latdd = lat.flatten()
            self._grid_londd = lon.flatten()
        else:
            raise ValueError(
                "The lat and lon data is not formatted correctly. ")

    def get_idx_from_latlon(self, latitudes: list[float],
                            longitudes: list[float]) -> list[int]:
        """
        Return the logical indexes of the cells of the flatened grid
        containing each pair of latitude and longitude coordinates.

        Parameters
        ----------
        latitudes : array-like, Iterable, or scalar value
            Contains the latitude values, in decimal degrees, of the
            coordinates for which we want to find the logical indexes of
            the corresponding grid cells.
        longitudes : array-like, Iterable, or scalar value
            Contains the longitude values, in decimal degrees, of the
            coordinates for which we want to find the logical indexes of
            the corresponding grid cells.

        Returns
        -------
        list[int]
            A list containing the flattened grid indexes of the cells
            containint the location at which climate data is to be extracted.
        """
        lat_max = np.max(latitudes) + 1
        lat_min = np.min(latitudes) - 1
        lon_max = np.max(longitudes) + 1
        lon_min = np.min(longitudes) - 1

        mask = (
            (self._grid_latdd <= lat_max) &
            (self._grid_latdd >= lat_min) &
            (self._grid_londd <= lon_max) &
            (self._grid_londd >= lon_min) &
            (~self._nan_grid_mask)
            )

        masked_grid_latdd = self._grid_latdd[mask]
        masked_grid_londd = self._grid_londd[mask]

        grid_idx = np.arange(len(self._grid_latdd))
        masked_grid_idx = [
            np.argmin(calc_dist_from_coord(
                lat, lon, masked_grid_latdd, masked_grid_londd)) for
            lat, lon in zip(latitudes, longitudes)
            ]

        return grid_idx[mask][masked_grid_idx]

    def create_connect_table(self, lat_dd: list[float], lon_dd: list[float],
                             loc_id: list = None) -> ConnectTable:
        """
        Create a connection table that contains the relation between a set
        of location coordinates and the cell of the grids.

        Parameters
        ----------
        latitudes : list[float]
            Contains the latitude coordinates, in decimal degrees, of the
            locations for which you want to extract climate data from the grid.
        longitudes : list[float]
            Contains the longitude coordinates, in decimal degrees, of the
            locations for which you want to extract climate data from the grid
        loc_id : list
            An optional list of custom ID to identify each location for which
            you want to extract climate data from the grid. If not 'loc_id'
            are specified, a list of incremental integer values will be
            generated and used by default.

        Returns
        -------
        ConnectTable
            A pandas dataframe containing the following columns:

            * loc_id: The identifiers of the locations for which climate data
              is to be extracted from the grid.
            * loc_lat_dd: The latitude coordinates of the locations for which
              climate data is to be extracted from the grid.
            * loc_lon_dd : The latitude coordinates of the locations for which
              climate data is to be extracted from the grid.
            * grid_idx: The cell indexes of the flattened grid containing
              the locations for which climate data is to be extracted.
            * grid_lat_dd: The latitude coordinates of the cells containing
              the location for which climate data is to be extracted.
            * grid_lon_dd: The longitude coordinates of the cells containing
              the locations for which climate data is to be extracted.
            * dist_km: The distance in km between the locations for which
              climate data is to be extracted and the location of the
              corresponding cells of the grid.
        """
        print("Creating a connection table...")
        connect_table = ConnectTable([])

        if loc_id is None:
            connect_table['loc_id'] = list(range(len(connect_table)))
        else:
            connect_table['loc_id'] = loc_id
        connect_table['loc_lat_dd'] = lat_dd
        connect_table['loc_lon_dd'] = lon_dd

        indexes = self.get_idx_from_latlon(lat_dd, lon_dd)

        connect_table['grid_idx'] = indexes
        connect_table['grid_lat_dd'] = self._grid_latdd[indexes]
        connect_table['grid_lon_dd'] = self._grid_londd[indexes]
        connect_table['dist_km'] = calc_dist_from_coord(
            connect_table['loc_lat_dd'].values,
            connect_table['loc_lon_dd'].values,
            connect_table['grid_lat_dd'].values,
            connect_table['grid_lon_dd'].values
            ).round(2)

        return connect_table

    def get_data(self, connect_table: pd.DataFrame,
                 first_year: int, last_year: int) -> ExtractedDataFrame:
        """
        Extract data from the grid.

        Parameters
        ----------
        varname : str
            The variable for which values are to be extracted from the grid.
        connect_table: ConnectTable
            The connection table containing the information on the
            locations for which data is to be extracted and their
            relation with the nodes of the grid.
        first_year: int
            The first year of the period for which data is to be
            extracted from the grid.
        last_year: int
            The last year of the period for which data is to be
            extracted from the grid.

        Returns
        -------
        extracted_data: ExtractedDataFrame
            A pandas dataframe containing the data extracted from the grid.
        """
        years = np.arange(first_year, last_year + 1)

        grid_extract_info = (
            connect_table.copy()
            .drop_duplicates(subset='grid_idx')
            .sort_values('grid_idx'))

        grid_lat_dd = grid_extract_info['grid_lat_dd'].values
        grid_lon_dd = grid_extract_info['grid_lon_dd'].values
        grid_idx = grid_extract_info['grid_idx'].values.tolist()

        data_stack = []
        index_stack = pd.Index([], dtype='datetime64[ns]')
        source = ''
        for year in years:
            print('Fetching daily {} data for year {}...'.format(
                self.varname, year))

            ncfilename = self.filename_pattern.format(year=year)
            ncfilepath = osp.join(self.gridpath, ncfilename)
            if not osp.exists(ncfilepath):
                print("'{}' does not exist: skipping.".format(ncfilename))
                continue

            with netCDF4.Dataset(ncfilepath, 'r+') as ncdset:
                array = np.array(ncdset[self.varname])
                source = ncdset.source

            data_stack.append(
                array.reshape(array.shape[0], -1)[:, grid_idx])
            index_stack = index_stack.append(
                pd.date_range(start=datetime.datetime(year, 1, 1),
                              end=datetime.datetime(year, 12, 31)))

        extracted_data = ExtractedDataFrame(
            data=np.vstack(data_stack),
            index=index_stack,
            columns=pd.MultiIndex.from_tuples(
                zip(grid_idx, grid_lat_dd, grid_lon_dd),
                names=['Grid index', 'Latitude(dd)', 'Longitude (dd)'])
            )
        extracted_data = extracted_data.sort_index()
        extracted_data.varname = self.varname
        extracted_data.source = source

        # Apply a mask to set missing values to NaN.
        # Missing values can be attributed a value of
        # -1.7976931348623157e+308 in the netcdf files.
        extracted_data = extracted_data[extracted_data > -999].copy()

        return extracted_data
