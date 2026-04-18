# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHELP Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------


# ---- Standard Library imports
import os.path as osp
import calendar

# ---- Third Party imports
import numpy as np
import pandas as pd

# ---- Local imports
from pyhelp.utils import save_content_to_csv


def save_precip_to_HELP(filename: str, precip: pd.Series, city: str):
    """
    Formats and saves a daily precipitation time series in mm
    to the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D4' else filename + '.D4'

    fheader = format_weather_header_for_HELP(3, 2, city)

    # Calculate and add monthly normals to header.
    monthly_normals = precip.groupby(precip.index.month).mean().values
    fheader.append([''.join([f'{x:>6.2f}' for x in monthly_normals])])

    fdata = format_timeseries_for_HELP(
        precip.index.year.values,
        precip.values,
        '{0:>10}', '{0:>5.1f}')
    save_content_to_csv(filename, fheader + fdata)


def save_airtemp_to_HELP(filename: str, airtemp: pd.Series, city: str):
    """
    Formats and saves a daily average air temperature time series in
    Celcius to the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D7' else filename + '.D7'

    fheader = format_weather_header_for_HELP(3, 2, city)

    # Calculate and add monthly normals to header.
    monthly_normals = airtemp.groupby(airtemp.index.month).mean().values
    fheader.append([''.join([f'{x:>6.2f}' for x in monthly_normals])])

    fdata = format_timeseries_for_HELP(
        airtemp.index.year.values,
        airtemp.values,
        '{0:>5}', '{0:>6.1f}')
    save_content_to_csv(filename, fheader + fdata)


def save_solrad_to_HELP(
        filename: str, solrad: pd.Series, city: str, lat: float
        ):
    """
    Formats and saves a daily global solar radiation time series in MJ/m2/day
    to the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D13' else filename + '.D13'

    fheader = format_weather_header_for_HELP(3, 2, city, lat)
    fdata = format_timeseries_for_HELP(
        solrad.index.year.values,
        solrad.values,
        '{0:>5}', '{0:>6.2f}')
    save_content_to_csv(filename, fheader + fdata)


def format_weather_header_for_HELP(itype, iunits, city, lat=None):
    """
    Prepare the header for the precipitation, air temperature and
    global solar radiation input weather datafile for HELP. The format of the
    header is defined in the subroutine READIN of the HELP Fortran source code.
    """
    fheader = [['{0:>2}'.format(itype)],  # 3: data was entered by the user.
               ['{0:>2}'.format(iunits)],  # 1 for IP and 2 for SI
               ['{0:<40}'.format(city[:40])],
               ]
    if lat is not None:
        # Append the latitude if the data are solar radiation.
        fheader.append(['{0:>6.2f}'.format(lat)])
    return fheader


def format_timeseries_for_HELP(years, data, year_format, data_format):
    fdata = []
    for year in np.unique(years):
        # Selects the data and asserts that the data are complete for
        # that year :

        indexes = np.where(years == year)[0]
        days_in_year = 366 if calendar.isleap(year) else 365
        assert len(indexes) == days_in_year

        # Adds zeros to complete de last row and reshape the data
        # in a 37 x 10 grid:

        year_data = data[indexes]
        year_data = np.hstack(
            [year_data, np.zeros(10 - len(year_data) % 10)])
        year_data = year_data.reshape(37, 10).tolist()

        # Save the data in a format compatible with HELP :

        for line_data in year_data:
            formated_line = year_format.format(year)
            for i in range(10):
                formated_line += data_format.format(line_data[i])
            fdata.append([formated_line])
    return fdata
