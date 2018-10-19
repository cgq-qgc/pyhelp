# -*- coding: utf-8 -*-

# Copyright © 2014-2018 GWHAT Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of GWHAT (Ground-Water Hydrograph Analysis Toolbox).
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library imports
import os
import os.path as osp
import csv
import calendar
from calendar import monthrange
import datetime
from time import strftime

# ---- Third Party imports
import numpy as np
import netCDF4
from xlrd.xldate import xldate_from_datetime_tuple

# ---- Local imports
from pyhelp import __namever__
from pyhelp.utils import save_content_to_csv, nan_as_text_tolist


class InfoClimatGridReader(object):
    """
    The :attr:`~pyhelp.weather_reader.InfoClimatGridReader` is a class
    to read and format precipitation and air temperature data from the
    interpolated grid produced by the `Info-climat service`_ of the MDDELCC.

    .. _Info-climat service:
       http://www.mddelcc.gouv.qc.ca/climat/surveillance/produits.htm
    """

    def __init__(self, dirpath_netcdf):
        super(InfoClimatGridReader, self).__init__()
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

    def get_idx_from_latlon(self, latitudes, longitudes, unique=False):
        """
        Get the i and j indexes of the grid meshes from a list of latitude
        and longitude coordinates. If unique is True, only the unique pairs of
        i and j indexes will be returned.
        """
        try:
            lat_idx = [np.argmin(np.abs(self.lat - lat)) for lat in latitudes]
            lon_idx = [np.argmin(np.abs(self.lon - lon)) for lon in longitudes]
            if unique:
                ijdx = np.vstack({(i, j) for i, j in zip(lat_idx, lon_idx)})
                lat_idx = ijdx[:, 0].tolist()
                lon_idx = ijdx[:, 1].tolist()
        except TypeError:
            lat_idx = np.argmin(np.abs(self.lat - latitudes))
            lon_idx = np.argmin(np.abs(self.lon - longitudes))

        return lat_idx, lon_idx

    def get_data_from_latlon(self, latitudes, longitudes, years):
        """
        Return the daily minimum, maximum and average air temperature and daily
        precipitation
        """
        lat_idx, lon_idx = self.get_idx_from_latlon(latitudes, longitudes)
        return self.get_data_from_idx(lat_idx, lon_idx, years)

    def get_data_from_idx(self, lat_idx, lon_idx, years):
        try:
            len(lat_idx)
        except TypeError:
            lat_idx, lon_idx = [lat_idx], [lon_idx]

        tasmax_stacks = []
        tasmin_stacks = []
        precip_stacks = []
        years_stack = []
        for year in years:
            print('\rFetching daily weather data for year %d...' % year,
                  end=' ')
            filename = osp.join(self.dirpath_netcdf, 'GCQ_v2_%d.nc' % year)
            netcdf_dset = netCDF4.Dataset(filename, 'r+')

            tasmax_stacks.append(
                np.array(netcdf_dset['tasmax'])[:, lat_idx, lon_idx])
            tasmin_stacks.append(
                np.array(netcdf_dset['tasmin'])[:, lat_idx, lon_idx])
            precip_stacks.append(
                np.array(netcdf_dset['pr'])[:, lat_idx, lon_idx])
            years_stack.append(
                np.zeros(len(precip_stacks[-1][:])).astype(int) + year)

            netcdf_dset.close()
        print('done')

        tasmax = np.vstack(tasmax_stacks)
        tasmin = np.vstack(tasmin_stacks)
        precip = np.vstack(precip_stacks)
        years = np.hstack(years_stack)

        return (tasmax + tasmin)/2, precip, years

    def generate_input_from_MDELCC_grid(self, outdir, lat_dd, lon_dd,
                                        year_range):
        """
        Generate input data files from the MDDELCC grid.

        Generate PyHelp csv data file inputs for daily precipitation and
        average air temperature  using data from the MDDELCC spatially
        distributed daily precipitation and minimum and maximum air
        temperature grid for a set of lat/lon coordinates.
        """
        if not osp.exists(outdir):
            os.makedirs(outdir)

        lat_idx, lon_idx = self.get_idx_from_latlon(
            lat_dd, lon_dd, unique=True)
        lat_dd = [self.lat[i] for i in lat_idx]
        lon_dd = [self.lon[i] for i in lon_idx]

        # Fetch the daily weather data from the netCDF files.
        years = range(year_range[0], year_range[1] + 1)
        tasavg, precip, years = self.get_data_from_idx(lat_idx, lon_idx, years)

        # Create an array of datestring and lat/lon
        Ndt, Ndset = np.shape(tasavg)
        start = datetime.datetime(years[0], 1, 1)
        datetimes = [start + datetime.timedelta(days=i) for i in range(Ndt)]
        datestrings = [dt.strftime("%d/%m/%Y") for dt in datetimes]

        # Fill -999 with 0 in daily precip.
        precip[:, :][precip[:, :] == -999] = 0

        # Fill -999 with linear interpolation in daily air temp.
        time_ = np.arange(Ndt)
        for i in range(Ndset):
            indx = np.where(tasavg[:, i] != -999)[0]
            tasavg[:, i] = np.interp(time_, time_[indx], tasavg[:, i][indx])

        #  Convert and save the weather data to PyHelp csv input files.
        for var in ['precip', 'airtemp']:
            if var == 'precip':
                varname = 'Precipitation in mm'
                data = nan_as_text_tolist(precip)
            elif var == 'airtemp':
                varname = 'Average daily air temperature in \u00B0C'
                data = nan_as_text_tolist(tasavg)
            fname = osp.join(outdir, var + '_input_data.csv')

            print('Saving {} data to {}...'.format(var, fname), end=' ')
            fheader = [
                [varname],
                ['', ''],
                ['Created by ' + __namever__],
                ['Created on ' + strftime("%d/%m/%Y")],
                ['Created from MDDELCC grid'],
                ['', ''],
                ['Latitude (dd)'] + lat_dd,
                ['Longitude (dd)'] + lon_dd,
                ['', '']]
            fdata = [[datestrings[i]] + data[i] for i in range(Ndt)]
            fcontent = fheader + fdata
            save_content_to_csv(fname, fcontent)
            print('done')


# ---- Read CWEEDS Files
def generate_input_from_cweeds(outdir, cweed2_paths, cweed3_paths, year_range):
    """Generate an input PyHelp data file from CWEED files."""
    if not isinstance(cweed2_paths, (list, tuple)):
        cweed2_paths = [cweed2_paths]
    if not isinstance(cweed3_paths, (list, tuple)):
        cweed3_paths = [cweed3_paths]

    print('Reading CWEEDS files...', end=' ')
    lat_dd = []
    lon_dd = []
    stations = []
    data = []
    for cweed2, cweed3 in zip(cweed2_paths, cweed3_paths):
        daily_wy2 = read_cweeds_file(cweed2, format_to_daily=True)
        daily_wy3 = read_cweeds_file(cweed3, format_to_daily=True)
        wy23_df = join_daily_cweeds_wy2_and_wy3(daily_wy2, daily_wy3)

        lat_dd.append(wy23_df['Latitude'])
        lon_dd.append(wy23_df['Longitude'])
        stations.append(wy23_df['Location'])

        indexes = np.where((wy23_df['Years'] >= year_range[0]) &
                           (wy23_df['Years'] <= year_range[1]))[0]
        data.append(wy23_df['Irradiance'][indexes])
    data = nan_as_text_tolist(np.array(data).astype(float).transpose())
    print('done')

    fname = osp.join(outdir, 'solrad_input_data.csv')
    print('Saving {} data to {}...'.format('solrad', fname), end=' ')

    # Create an array of datestring and lat/lon
    Ndt = len(wy23_df['Years'][indexes])
    start = datetime.datetime(year_range[0], 1, 1)
    datetimes = [start + datetime.timedelta(days=i) for i in range(Ndt)]
    datestrings = [dt.strftime("%d/%m/%Y") for dt in datetimes]

    # Save the data to file.
    fheader = [['Global solar irradiance in MJ/m²'],
               ['', ''],
               ['Created by ' + __namever__],
               ['Created on ' + strftime("%d/%m/%Y")],
               ['Created from CWEED files'],
               ['', ''],
               ['Stations'] + stations,
               ['Latitude (dd)'] + lat_dd,
               ['Longitude (dd)'] + lon_dd,
               ['', '']]
    fdata = [[datestrings[i]] + data[i] for i in range(Ndt)]
    fcontent = fheader + fdata
    save_content_to_csv(fname, fcontent)
    print('done')


def read_cweeds_file(filename, format_to_daily=True):
    """
    Reads and formats data from a CWEEDS file, either version WY2 or WY3.
    Returns a dictionary, which includes a numpy array of the global
    solar irradiance in MJ/m², as well as corresponding arrays of the years,
    months, days, and hours. By default, the hourly data from the CWEEDS file
    are formated to daily values. The data are kept in a hourly format if
    format_to_daily is set to False.
    """
    # Determine if the CWEEDS file is in the WY2 or WY3 format.
    root, ext = osp.splitext(filename)
    ext = ext.replace('.', '')
    if ext not in ['WY2', 'WY3']:
        raise ValueError("%s is not a valid file extension. CWEEHDS files must"
                         " have either a WY2 or WY3 extension" % ext)

    # Open and format the data from the CWEEDS file.
    with open(filename, 'r') as f:
        reader = list(csv.reader(f))

    header_df = {}
    if ext == 'WY3':
        # We remove the header line from the data if the format is WY3.
        header_list = reader.pop(0)
        header_df['HORZ version'] = header_list[0]
        header_df['Location'] = header_list[1]
        header_df['Province'] = header_list[2]
        header_df['Country'] = header_list[3]
        header_df['Station ID'] = header_list[4]
        header_df['Latitude'] = float(header_list[5])
        header_df['Longitude'] = float(header_list[6])
        header_df['Time Zone'] = float(header_list[7])
        header_df['Elevation'] = float(header_list[8])

    char_offset = 0 if ext == 'WY2' else 2
    hourly_df = {}
    hourly_df['Years'] = np.empty(len(reader)).astype(int)
    hourly_df['Months'] = np.empty(len(reader)).astype(int)
    hourly_df['Days'] = np.empty(len(reader)).astype(int)
    hourly_df['Hours'] = np.empty(len(reader)).astype(int)
    hourly_df['Time'] = np.empty(len(reader)).astype('float64')
    # Global horizontal irradiance, kJ/m²
    hourly_df['Irradiance'] = np.empty(len(reader)).astype('float64')

    for i, line in enumerate(reader):
        hourly_df['Years'][i] = year = int(line[0][char_offset:][6:10])
        hourly_df['Months'][i] = month = int(line[0][char_offset:][10:12])
        hourly_df['Days'][i] = day = int(line[0][char_offset:][12:14])
        hourly_df['Hours'][i] = hour = int(line[0][char_offset:][14:16]) - 1
        # The global horizontal irradiance is converted from kJ/m² to MJ/m².
        hourly_df['Irradiance'][i] = float(line[0][char_offset:][20:24])/1000

        # Compute time in Excel numeric format :
        hourly_df['Time'][i] = xldate_from_datetime_tuple(
                (year, month, day, hour, 0, 0), 0)

    if format_to_daily:
        # Convert the hourly data to daily format.
        assert len(hourly_df['Irradiance']) % 24 == 0
        new_shape = (len(hourly_df['Irradiance'])//24, 24)

        daily_df = {}
        daily_df['Irradiance'] = np.sum(
                hourly_df['Irradiance'].reshape(new_shape), axis=1)
        for key in ['Years', 'Months', 'Days', 'Time']:
            daily_df[key] = hourly_df[key].reshape(new_shape)[:, 0]
        daily_df['Hours'] = np.zeros(len(daily_df['Irradiance']))

        daily_df.update(header_df)
        daily_df['Time Format'] = 'daily'
        daily_df['CWEEDS Format'] = ext
        return daily_df
    else:
        hourly_df.update(header_df)
        hourly_df['Time Format'] = 'hourly'
        hourly_df['CWEEDS Format'] = ext
        return hourly_df


def join_daily_cweeds_wy2_and_wy3(wy2_df, wy3_df):
    """
    Join a CWEEDS dataset in the wy2 format to another cweeds dataset in the
    wy3 format.
    """
    assert wy2_df['CWEEDS Format'] == 'WY2'
    assert wy3_df['CWEEDS Format'] == 'WY3'
    assert wy2_df['Time Format'] == wy3_df['Time Format']

    time_wy23 = np.hstack([wy2_df['Time'], wy3_df['Time']])
    time_wy23 = np.unique(time_wy23)
    time_wy23 = np.sort(time_wy23)

    wy23_df = {}
    wy23_df['Time Format'] = wy3_df['Time Format']
    wy23_df['CWEEDS Format'] = 'WY2+WY3'

    # Copy the header info from WY3 dataset :

    for key in ['HORZ version', 'Location', 'Province', 'Country',
                'Station ID', 'Latitude', 'Longitude', 'Time Zone',
                'Elevation']:
        wy23_df[key] = wy3_df[key]

    # Merge the two datasets :

    wy23_df['Time'] = time_wy23
    wy23_df['Years'] = np.empty(len(time_wy23)).astype(int)
    wy23_df['Months'] = np.empty(len(time_wy23)).astype(int)
    wy23_df['Days'] = np.empty(len(time_wy23)).astype(int)
    wy23_df['Hours'] = np.empty(len(time_wy23)).astype(int)
    wy23_df['Irradiance'] = np.empty(len(time_wy23)).astype('float64')

    for dataset in [wy2_df, wy3_df]:
        indexes = np.digitize(dataset['Time'], time_wy23, right=True)
        for key in ['Years', 'Months', 'Days', 'Hours', 'Irradiance']:
            wy23_df[key][indexes] = dataset[key]

    return wy23_df


# ---- Export to HELP format

def save_precip_to_HELP(filename, years, precip, city):
    """
    Formats and saves a daily precipitation time series in mm
    to the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D4' else filename + '.D4'

    fheader = format_weather_header_for_HELP(3, 2, city)
    fdata = format_timeseries_for_HELP(years, precip, '{0:>10}', '{0:>5.1f}')
    save_content_to_csv(filename, fheader + fdata)


def save_airtemp_to_HELP(filename, years, precip, city):
    """
    Formats and saves a daily average air temperature time series in Celcius to
    the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D7' else filename + '.D7'

    fheader = format_weather_header_for_HELP(3, 2, city)
    fdata = format_timeseries_for_HELP(years, precip, '{0:>5}', '{0:>6.1f}')
    save_content_to_csv(filename, fheader + fdata)


def save_solrad_to_HELP(filename, years, precip, city, lat):
    """
    Formats and saves a daily global solar radiation time series in MJ/m2/day
    to the HELP format.
    """
    root, ext = osp.splitext(filename)
    filename = filename if ext == '.D13' else filename + '.D13'

    fheader = format_weather_header_for_HELP(3, 2, city, lat)
    fdata = format_timeseries_for_HELP(years, precip, '{0:>5}', '{0:>6.2f}')
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
    else:
        fheader.append([])
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


def save_data_to_HELP_format(filename, years, data, city, lat=None):
    """Formats a time series to the HELP format."""
    root, ext = osp.splitext(filename)
    ext = ext[1:]
    if ext == 'D4':  # precipitation
        year_format = '{0:>10}'
        data_format = '{0:>5.1f}'
    elif ext == 'D7':  # air temperature
        year_format = '{0:>5}'
        data_format = '{0:>6.1f}'
    elif ext == 'D13':  # global solar radiation
        year_format = '{0:>5}'
        data_format = '{0:>6.2f}'
        if lat is None:
            raise ValueError("A value must be specified for lat.")
    else:
        raise ValueError("%s is not a valid file extension." % ext)

    # ---- Format Header

    itype = 3   # Precipitation data for {city} was entered by the user.
    iunits = 2  # 1 for IP and 2 for SI
    fcontent = [['{0:>2}'.format(itype)],
                ['{0:>2}'.format(iunits)],
                ['{0:<40}'.format(city[:40])],
                ]
    if ext == 'D13':
        # Append the latitude if the data are solar radiation.
        fcontent.append(['{0:>6.2f}'.format(lat)])
    else:
        fcontent.append([])

    # ---- Format Data

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
            fcontent.append([formated_line])

    save_content_to_csv(filename, fcontent)


# ---- Base functions: monthly downscaling

def calc_monthly_sum(yy_dly, mm_dly, x_dly):
    """
    Calcul monthly cumulative values from daily values, where yy_dly are the
    years, mm_dly are the months (1 to 12), and x_dly are the daily values.
    """
    return calc_monthly(yy_dly, mm_dly, x_dly, np.sum)


def calc_monthly_mean(yy_dly, mm_dly, x_dly):
    """
    Calcul monthly mean values from daily values, where yy_dly are the
    years, mm_dly are the months (1 to 12), and x_dly are the daily values.
    """
    return calc_monthly(yy_dly, mm_dly, x_dly, np.mean)


def calc_monthly(yy_dly, mm_dly, x_dly, func):
    yy = np.unique(yy_dly)
    mm = range(1, 13)

    yy_mly = np.repeat(yy, len(mm))
    mm_mly = np.tile(mm, len(yy))
    x_mly = np.zeros(len(mm)*len(yy))

    for i in range(len(mm)*len(yy)):
        indx = np.where((yy_dly == yy_mly[i]) & (mm_dly == mm_mly[i]))[0]
        if len(indx) < monthrange(yy_mly[i], mm_mly[i])[1]:
            x_mly[i] = np.nan  # incomplete dataset for this month
        else:
            x_mly[i] = func(x_dly[indx])

    return yy_mly, mm_mly, x_mly


def calcul_monthly_normals(years, months, x_mly, yearmin=None, yearmax=None):
    """Calcul the monthly normals from monthly values."""
    if len(years) != len(months) != len(x_mly):
        raise ValueError("The dimension of the years, months, and x_mly array"
                         " must match exactly.")
    if np.min(months) < 1 or np.max(months) > 12:
        raise ValueError("Months values must be between 1 and 12.")

    # Mark as nan monthly values that are outside the year range that is
    # defined by yearmin and yearmax :
    x_mly = np.copy(x_mly)
    if yearmin is not None:
        x_mly[years < yearmin] = np.nan
    if yearmax is not None:
        x_mly[years > yearmax] = np.nan

    # Calcul the monthly normals :
    x_norm = np.zeros(12)
    for i, mm in enumerate(range(1, 13)):
        indx = np.where((months == mm) & (~np.isnan(x_mly)))[0]
        if len(indx) > 0:
            x_norm[i] = np.mean(x_mly[indx])
        else:
            x_norm[i] = np.nan

    return x_norm


# ----- Base functions: yearly downscaling

def calc_yearly_sum(yy_dly, x_dly):
    """
    Calcul yearly cumulative values from daily values, where yy_dly are the
    years and x_dly are the daily values.
    """
    return calc_yearly(yy_dly, x_dly, np.sum)


def calc_yearly_mean(yy_dly, x_dly):
    """
    Calcul yearly mean values from daily values, where yy_dly are the years
    and x_dly are the daily values.
    """
    return calc_yearly(yy_dly, x_dly, np.mean)


def calc_yearly(yy_dly, x_dly, func):
    yy_yrly = np.unique(yy_dly)
    x_yrly = np.zeros(len(yy_yrly))
    for i in range(len(yy_yrly)):
        indx = np.where(yy_dly == yy_yrly[i])[0]
        x_yrly[i] = func(x_dly[indx])

    return yy_yrly, x_yrly


# ----- Base functions: secondary variables

def calcul_rain_from_ptot(Tavg, Ptot, Tcrit=0):
    rain = np.copy(Ptot)
    rain[np.where(Tavg < Tcrit)[0]] = 0
    return rain


if __name__ == '__main__':
    outdir = "C:\\Users\\User\\pyhelp\\example"
    cweed2_paths = osp.join(outdir, 'CWEEDS', '94792.WY2')
    cweed3_paths = osp.join(
        outdir, 'CWEEDS', 'CAN_QC_MONTREAL-INTL-A_7025251_CWEEDS2011_T_N.WY3')
    year_range = [2010, 2014]
    generate_input_from_cweeds(outdir, cweed2_paths, cweed3_paths, year_range)
