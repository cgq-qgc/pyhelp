# -*- coding: utf-8 -*-
"""
A Script to download solar radiation data from CDS.

https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-agrometeorological-indicators?tab=overview
"""

# https://pypi.org/project/cdsapi/
# https://github.com/ecmwf/cdsapi

from __future__ import annotations

import os.path as osp
import cdsapi
import netCDF4
import numpy as np
import pandas as pd
import zipfile
import datetime
import re
import tempfile


def read_zipped_solrad_netcdf(path_to_zip: str, bbox: list):
    """
    Read solar radiation data from a zip archive.
    """
    with zipfile.ZipFile(path_to_zip) as thezip:
        namelist = thezip.namelist()
        namelist.sort()

        time = []
        datastack = []
        for filename in namelist:
            time.append(
                datetime.datetime.strptime(
                    re.findall('AgERA5_(.*?)_final', filename)[0], '%Y%m%d'))
            with thezip.open(filename, mode='r') as thefile:
                with netCDF4.Dataset(
                        'dummy', mode='r', memory=thefile.read()) as nc:
                    dataf = pd.DataFrame(
                        data=np.array(nc['Solar_Radiation_Flux'])[0, :, :],
                        index=np.array(nc['lat']),
                        columns=np.array(nc['lon']))
                mask_index = (
                    (dataf.index <= bbox[1]) &
                    (dataf.index >= bbox[3]))
                mask_columns = (
                    (dataf.columns >= bbox[0]) &
                    (dataf.columns <= bbox[2]))
                dataf = dataf.loc[mask_index, mask_columns]
                datastack.append(dataf.values)

    data = np.stack(datastack, axis=0)
    time = np.array(time)
    lat = dataf.index.values
    lon = dataf.columns.values

    return data, time, lat, lon


def save_solrad_to_nc(filename: str, data: np.ndarray, time: np.ndarray,
                      lat: np.ndarray, lon: np.ndarray):
    """
    Save a year of solar radiation data to a netcdf file.
    """
    with netCDF4.Dataset(filename, 'w', format="NETCDF4") as ncfile:
        ncfile.year = time[0].year
        ncfile.source = "https://doi.org/10.24381/cds.6c68c9bb"

        ncfile.createDimension('time', len(time))
        ncfile.createDimension('lat', len(lat))
        ncfile.createDimension('lon', len(lon))

        var_time = ncfile.createVariable('time', 'f4', ('time',))
        var_lat = ncfile.createVariable('lat', 'f8', ('lat',))
        var_lon = ncfile.createVariable('lon', 'f8', ('lon',))
        var_solrad = ncfile.createVariable(
            'solrad', 'f8', ('time', 'lat', 'lon',))

        var_time.unit = 'day of year'
        var_lat.unit = 'decimal degrees'
        var_lon.unit = 'decimal degrees'
        var_solrad.unit = 'MJ m-2 day-1'

        var_solrad.description = (
            "Total amount of energy provided by solar radiation at "
            "the surface over the period 00-24h local time per unit "
            "area and time.")

        var_lat[:] = lat
        var_lon[:] = lon
        var_time[:] = np.arange(1, len(time) + 1)
        var_solrad[:] = data


def get_solrad_data(year: int, bbox: list[float, float, float, float],
                    filename: str = None):
    """
    Get a year of solar radiation data from CDS.

    Parameters
    ----------
    year : int
        The year for which to download the solar radiation data.
    bbox : list[float, float, float, float]
        The bounding box within which to extract the data. The first element
        corresponds to the westernmost longitude of the bbox, the second to
        the northernmost latitude, the third to the easternmost longitude, and
        the fourth to the southernmost latitude.
    filename: str
        A string corresponding to a netCDF file path where to save the data.

    Returns
    -------
    data : np.ndarray
        A 3D numpy matrix containing the daily solar radiation data, where
        axes 0 corresponds to the time, axes 1 to the latitude and axes 2 to
        the longitude.
    time : np.ndarray
        A 1D array containing the time of the data.
    lat : np.ndarray
        A 1D array containing the latitudes of the data.
    lon : np.ndarray
        A 1D array containing the longitudes of the data.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        zippath = osp.join(tmpdir, 'download.zip')

        timestack = []
        datastack = []
        cds = cdsapi.Client()
        for month in range(1, 13):
            cds.retrieve(
                'sis-agrometeorological-indicators',
                {"variable": "solar_radiation_flux",
                 "year": ['{:0.0f}'.format(year)],
                 'month': ['{0:02.0f}'.format(month)],
                 "day": ['{0:02.0f}'.format(day) for day in range(1, 32)],
                 "format": "zip"},
                zippath)
            data, time, lat, lon = read_zipped_solrad_netcdf(zippath, bbox)
            timestack.append(time)
            datastack.append(data)

    data = np.vstack(datastack)
    time = np.hstack(timestack)

    # Convert solar radiation data to MJ/m2/day (from J/m2/day).
    data = data / 10**6

    if filename:
        save_solrad_to_nc(filename, data, time, lat, lon)

    return data, time, lat, lon


if __name__ == '__main__':
    year = 2018
    filename_pattern = 'solar_radiation_flux_{year}.nc'
    savedir = 'C:/Users/jean-/Documents/Data/CDS_Solar_Radiation_Flux'
    savefile = osp.join(savedir, filename_pattern.format(year=year))

    data, time, lat, lon = get_solrad_data(
        year=year,
        bbox=[-81, 63, -55, 44],
        filename=savefile)

    from grid_data_extractor import GridExtractor
    grid_extractor = GridExtractor(
        gridpath=savedir,
        filename_pattern=filename_pattern,
        varname='solrad'
        )

    loc_id = ['loc1', 'loc2', 'loc3']
    lat_dd = [45.42571, 49.1564, 45.43753]
    lon_dd = [-73.0764, -68.24755, -73.0813]

    connect_table = grid_extractor.create_connect_table(
        lat_dd, lon_dd, loc_id)
    solrad_data = grid_extractor.get_data(
        connect_table, first_year=2019, last_year=2019)

    print(solrad_data)
    print()
    print(connect_table)

    # connect_table.save_to_csv('connect_table.csv')
    # solrad_data.save_to_csv('solrad_data_2019.csv')
