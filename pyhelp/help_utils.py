# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:37:59 2018
@author: jsgosselin
"""

# ---- Standard Library Imports

import sys
import os
import os.path as osp
import time
import calendar
import multiprocessing
import shutil
import subprocess
import csv
from collections import OrderedDict


# ---- Third Party imports

import netCDF4
import geopandas as gpd

import numpy as np
import h5py
from PyQt5.QtCore import pyqtSlot as QSlot
from PyQt5.QtCore import pyqtSignal as QSignal
from PyQt5.QtCore import Qt, QObject, QThread
from PyQt5.QtWidgets import QApplication


# ---- Local Libraries Imports

from pyhelp import HELP3O
from pyhelp.meteo.weather_reader import (
        save_precip_to_HELP, save_airtemp_to_HELP, save_solrad_to_HELP,
        read_cweeds_file, join_daily_cweeds_wy2_and_wy3)


DIRNAME = '.help_threads'
HELPNAME = 'H3OFS32F.exe'


def run_help(path_to_exe):
    if osp.exists(path_to_exe):
        subprocess.call(path_to_exe)


class HelpThreadPoolManager(QObject):
    """
    An object that parallelize the HELP calculation and post-processing
    over multiple threads.
    """

    def __init__(self, nthread=None, path_hdf5=None):
        super(QObject, self).__init__()
        self._thread_pool = []
        self._worker_pool = []
        # self._hdf5file = None
        self._output = {}

        self.d10data = None
        self.d11data = None
        self.connect_tables = None
        self.cellnames = []

        # Output from HELP :

        self._daily_out = False
        self._monthly_out = True
        self._yearly_out = False
        self._summary_out = False

        # Number of threads to use for the calculation :
        self.nthread = (multiprocessing.cpu_count() if nthread is None
                        else nthread)
        self.setup_threadpool()

        self.path_hdf5 = (osp.abspath('HELP.OUT') if path_hdf5 is None
                          else path_hdf5)
        self.setup_hdf5_output_file()

    def setup_threadpool(self):
        """Setup the threads in which HELP calculations will be run."""
        for i in range(self.nthread):
            clonedir = osp.join(DIRNAME, 'thread%d' % i)
            if not os.path.exists(clonedir):
                os.makedirs(clonedir)
            clonename = osp.join(clonedir, HELPNAME)
            if not osp.exists(clonename):
                shutil.copy2(HELPNAME, clonename)

            help_thread = QThread()
            help_worker = HelpWorker(clonename, self)
            self._thread_pool.append(help_thread)
            self._worker_pool.append(help_worker)

            help_worker.moveToThread(help_thread)
            help_thread.started.connect(help_worker.run_help)
            help_worker.sig_cellchunk_finished.connect(
                    self._end_thread_calculation)
            help_worker.sig_singlecell_finished.connect(
                    self._handle_singlecell_result)

    def setup_hdf5_output_file(self):
        """Setup the hdf5 file where HELP simulation output are saved."""
        if not osp.exists(osp.dirname(self.path_hdf5)):
            os.makedirs(osp.dirname(self.path_hdf5))
        # self._hdf5file = h5py.File(self.path_hdf5, mode='w')

    def load_help_D10D11_inputs(self, path_d10file, path_d11file):
        """
        Read the D10 and D11 simulation inputs from files that were produced
        with the LCNP Help application and return the data for all cells
        in an OrderedDict, where the keys are the names of the cells and
        the value the content of the D10 and D11 csv files for each cell.
        """
        self.d10data, self.d11data = read_d10d11_file(
                path_d10file, path_d11file)
        self.cellnames = list(self.d11data.keys())

    def load_meteo_connect_tables(self, filename):
        """
        Load the table that connects the D4, D7, and D13 input weather HELP
        files to each cell.
        """
        self.connect_tables = np.load(filename).item()

    def is_calcul_running(self):
        """Return whether calculations are still being run."""
        for thread in self._thread_pool:
            if thread.isRunning():
                return True
        else:
            return False

    def start_calculation(self):
        """
        Divide the cells in chunks and assign them to a thread and start the
        calculation process for each thread.
        """
        # self.cellnames = self.cellnames[:10000]
        self._output = {}
        self.__calcul_progress = 0
        self.__start_calcul_time = time.clock()
        print("Total number of cells: %d" % len(self.cellnames))
        print("Number of threads: %d" % self.nthread)
        cellchunks = np.array_split(self.cellnames, self.nthread)
        self._nbr_of_thread_running = 0
        for i in range(self.nthread):
            self._nbr_of_thread_running += 1
            self._worker_pool[i].cellchunk = cellchunks[i]
            self._thread_pool[i].start()

    @QSlot(QObject, float)
    def _end_thread_calculation(self, worker, calcul_time):
        self._nbr_of_thread_running += -1
        worker.thread().quit()

    def _store_monthly_values(self, filename, cellname):
        self._output[cellname] =  read_monthly_help_output(filename)
        # data = read_monthly_help_output(filename)
        # if 'years' not in list(self._hdf5file.keys()):
        #     self._hdf5file.create_dataset('years', data=data['years'])
        # cellgrp = self._hdf5file.create_group(cellname)
        # cellgrp.create_dataset('rain', data=data['rain'])
        # cellgrp.create_dataset('runoff', data=data['runoff'])
        # cellgrp.create_dataset('evapo', data=data['evapo'])
        # cellgrp.create_dataset('sub-runoff', data=data['sub-runoff'])
        # cellgrp.create_dataset('percolation', data=data['percolation'])
        # cellgrp.create_dataset('recharge', data=data['recharge'])
        # self._hdf5file.flush()

    def _store_daily_values(self, filename, cellname):
        data = read_daily_help_output(filename)
        # if 'years' not in list(self._hdf5file.keys()):
        #     self._hdf5file.create_dataset('years', data=data['years'])
        # if 'days' not in list(self._hdf5file.keys()):
        #     self._hdf5file.create_dataset('days', data=data['days'])
        # cellgrp = self._hdf5file.create_group(cellname)
        # cellgrp.create_dataset('rain', data=data['rain'])
        # cellgrp.create_dataset('runoff', data=data['runoff'])
        # cellgrp.create_dataset('et', data=data['et'])
        # cellgrp.create_dataset('ezone', data=data['ezone'])
        # cellgrp.create_dataset('head first', data=data['head first'])
        # cellgrp.create_dataset('drain first', data=data['drain first'])
        # cellgrp.create_dataset('leak first', data=data['leak first'])
        # cellgrp.create_dataset('leak last', data=data['leak last'])
        # self._hdf5file.flush()

    @QSlot(str, str)
    def _handle_singlecell_result(self, filename, cellname):
        self._store_monthly_values(filename, cellname)

        os.remove(filename)

        self.__calcul_progress += 1
        progress_pct = self.__calcul_progress/len(self.cellnames)*100
        tpassed = time.clock() - self.__start_calcul_time
        tremain = (100-progress_pct)*tpassed/progress_pct/60
        print('\r%0.1f%% (%d min remaining)' % (progress_pct, tremain), end='')

        if self.__calcul_progress == len(self.cellnames):
            calcul_time = (time.clock() - self.__start_calcul_time)
            print('\nCalculation time: %0.2fs\n' % calcul_time)
            # self._hdf5file.close()
            savedata_to_hdf5(self._output, 'monthly_help_all.out')
            app.quit()


class HelpWorker(QObject):

    sig_cellchunk_finished = QSignal(QObject, float)
    sig_singlecell_finished = QSignal(str, str)

    def __init__(self, help_exe_path, manager):
        super(HelpWorker, self).__init__()
        self.manager = manager
        self.help_exe_path = osp.abspath(help_exe_path)
        self.thread_dir = osp.dirname(self.help_exe_path)
        self.path_outparam = osp.join(self.thread_dir, 'OUTPARAM.DAT')
        self.path_datad10 = osp.join(self.thread_dir, 'DATA10.D10')
        self.path_datad11 = osp.join(self.thread_dir, 'DATA11.D11')
        self.cellchunk = []

    def run_help(self):
        """Run HELP for all the cells in cellchunk."""
        tstart = time.clock()
        for cellname in self.cellchunk:
            self._update_d10d11(cellname)
            self._update_outparam(cellname)
            subprocess.call(self.help_exe_path,
                            cwd=osp.dirname(self.help_exe_path))
            outdir = osp.join(self.thread_dir, cellname+'.OUT')
            self.sig_singlecell_finished.emit(outdir, cellname)
        tend = time.clock()
        self.sig_cellchunk_finished.emit(self, tend-tstart)

    def _update_d10d11(self, cellname):
        """Update the data in the D10 and D11 input files for cellname."""
        with open(self.path_datad10, 'w') as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerows(self.manager.d10data[cellname])
        with open(self.path_datad11, 'w') as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerows(self.manager.d11data[cellname])

    def _update_outparam(self, cellname):
        d4name = self.manager.connect_tables['D4'][cellname]
        d7name = self.manager.connect_tables['D7'][cellname]
        d13name = self.manager.connect_tables['D13'][cellname]
        shutil.copyfile(d4name, osp.join(self.thread_dir, 'DATA4.D4'))
        shutil.copyfile(d7name, osp.join(self.thread_dir, 'DATA7.D7'))
        shutil.copyfile(d13name, osp.join(self.thread_dir, 'DATA13.D13'))
        outputparam = [['DATA4.D4'],
                       ['DATA7.D7'],
                       ['DATA13.D13'],
                       ['DATA11.D11'],
                       ['DATA10.D10'],
                       [cellname + '.OUT'],
                       [2],
                       [15],
                       [int(self.manager._daily_out)],
                       [int(self.manager._monthly_out)],
                       [int(self.manager._yearly_out)],
                       [int(self.manager._summary_out)]
                       ]

        with open(self.path_outparam, 'w') as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerows(outputparam)


class HelpResultReader(QObject):
    def __init__(self, path_hdf5):
        super(HelpResultReader, self).__init__()
        self.path_hdf5 = path_hdf5
        self._hdf5 = h5py.File(self.path_hdf5, mode='r+')


def read_monthly_help_output(filename):
    """
    Read the monthly output from .OUT HELP file and return the data as
    numpy arrays stored in a dictionary.
    """
    with open(filename, 'r') as csvfile:
        csvread = list(csv.reader(csvfile))

    arr_years = []
    vstack_precip = []
    vstack_runoff = []
    vstack_evapo = []
    vstack_subrunoff = []
    vstack_percol = []
    vstack_rechg = []

    year = None
    i = 0
    while True:
        if i+1 >= len(csvread):
            break
        if len(csvread[i]) == 0:
            i += 1
            continue

        line = csvread[i][0]
        if 'MONTHLY TOTALS' in line:
            year = int(line.split()[-1])
            arr_years.append(year)
            subrunoff = None
            percol = None
            while True:
                i += 1
                if len(csvread[i]) == 0:
                    continue
                line = csvread[i][0]
                if '**********' in line:
                    break
                if len(csvread[i+1]) == 0:
                    continue

                nline = csvread[i+1][0]
                if 'PRECIPITATION' in line:
                    precip = line.split()[-6:] + nline.split()[-6:]
                elif 'RUNOFF' in line:
                    runoff = line.split()[-6:] + nline.split()[-6:]
                elif 'EVAPOTRANSPIRATION' in line:
                    evapo = line.split()[-6:] + nline.split()[-6:]
                elif 'LATERAL DRAINAGE' in line and subrunoff is None:
                    subrunoff = line.split()[-6:] + nline.split()[-6:]
                elif 'PERCOLATION' in line:
                    if percol is None:
                        percol = line.split()[-6:] + nline.split()[-6:]
                    rechg = line.split()[-6:] + nline.split()[-6:]
            vstack_precip.append(np.array(precip).astype('float32'))
            vstack_runoff.append(np.array(runoff).astype('float32'))
            vstack_evapo.append(np.array(evapo).astype('float32'))
            vstack_rechg.append(np.array(rechg).astype('float32'))
            vstack_percol.append(np.array(percol).astype('float32'))
            if subrunoff is None:
                vstack_subrunoff.append(np.zeros(12).astype('float32'))
            else:
                vstack_subrunoff.append(np.array(subrunoff).astype('float32'))
        elif 'FINAL WATER STORAGE' in line:
            break

        i += 1

    data = {'years': np.array(arr_years).astype('uint16'),
            'rain': np.vstack(vstack_precip),
            'runoff': np.vstack(vstack_runoff),
            'evapo': np.vstack(vstack_evapo),
            'sub-runoff': np.vstack(vstack_subrunoff),
            'percolation': np.vstack(vstack_percol),
            'recharge': np.vstack(vstack_rechg)}
    return data


def read_daily_help_output(filename):
    """
    Read the daily output from .OUT HELP file and return the data as
    numpy arrays stored in a dictionary.
    """
    with open(filename, 'r') as csvfile:
        csvread = list(csv.reader(csvfile))

    nlay = None
    arr_years = []
    arr_days = []
    arr_rain = []
    arr_ru = []
    arr_et = []
    arr_ezone = []
    arr_headfirst = []
    arr_drainfirst = []
    arr_leakfirst = []
    arr_leaklast = []

    year = None
    nlay = nsub = None
    for i, line in enumerate(csvread):
        if line:
            line = line[0]
            if 'TOTAL NUMBER OF LAYERS' in line:
                nlay = int(line.split()[-1])
            elif 'TOTAL NUMBER OF SUBPROFILES' in line:
                nsub = int(line.split()[-1])
            elif 'DAILY OUTPUT FOR YEAR' in line:
                year = int(line.split()[-1])
                days_in_year = 366 if calendar.isleap(year) else 365
            elif year is not None:
                try:
                    day = int(line[2:5])
                    rain = float(line[13:19])
                    ru = float(line[19:26])
                    et = float(line[26:33])
                    ezone = float(line[33:41])
                    headfirst = float(line[41:51])
                    drainfirst = float(line[51:61])
                    leakfirst = float(line[61:71])
                    leaklast = float(line[-10:])
                except ValueError:
                    pass
                else:
                    arr_years.append(year)
                    arr_days.append(day)
                    arr_rain.append(rain)
                    arr_ru.append(ru)
                    arr_et.append(et)
                    arr_ezone.append(ezone)
                    arr_headfirst.append(headfirst)
                    arr_drainfirst.append(drainfirst)
                    arr_leakfirst.append(leakfirst)
                    arr_leaklast.append(leaklast)
                    if day == days_in_year:
                        year = None

    dataf = {'years': np.array(arr_years).astype('uint16'),
             'days': np.array(arr_days).astype('uint16'),
             'rain': np.array(arr_rain).astype('float32'),
             'runoff': np.array(arr_ru).astype('float32'),
             'et': np.array(arr_et).astype('float32'),
             'ezone': np.array(arr_ezone).astype('float32'),
             'head first': np.array(arr_headfirst).astype('float32'),
             'drain first': np.array(arr_drainfirst).astype('float32'),
             'leak first': np.array(arr_leakfirst).astype('float32'),
             'leak last': np.array(arr_leaklast).astype('float32')
             }
    return dataf


def read_d10d11_file(path_d10file, path_d11file):
    """
    Read the concatenated D10 and D11 files that contain the D10 and D11 inputs
    for all the cells.
    """

    enc = 'iso-8859-1'  # 'utf-8-sig', 'iso-8859-1', 'utf-8', 'utf-16'

    # ---- Read and Format D11 File

    with open(path_d11file, 'r', encoding=enc) as csvfile:
        d11reader = list(csv.reader(csvfile))

    d11dat = OrderedDict()
    cellnames = []
    for i in range(0, len(d11reader), 3):
        cell_name, zone_name = d11reader[i+1][0].split()
        # cell_name = cell_name[3:]
        d11dat[cell_name] = [d11reader[i],
                             d11reader[i+1],
                             d11reader[i+2]]
        cellnames.append(cell_name)

    # ---- Read and Format D10 File

    with open(path_d10file, 'r', encoding=enc) as csvfile:
        d10reader = list(csv.reader(csvfile))

    d10dat = OrderedDict()
    i = 0
    curcell = None
    nextcell = cellnames[i]
    for line in d10reader:
        if nextcell in line[0]:
            d10dat[nextcell] = [line]
            curcell = nextcell
            nextcell = 'None' if i >= len(cellnames)-1 else cellnames[i+1]
            i += 1
        else:
            d10dat[curcell].append(line)

    return d10dat, d11dat


def write_d10d11_input_files(dirpath, d10data, d11data):
    """
    Write the content of each cell in individual D10 and D11 files.
    """
    if not osp.exists(dirpath):
        os.makedirs(dirpath)

    d10_connect_table = {}
    d11_connect_table = {}

    # Write evapotranspiration input file :

    N = len(d10data)
    for i, key in enumerate(d10data.keys()):
        print("\rProcessing D10 input file for cell %d of %d (%0.1f%%)." %
              (i+1, N, (i+1)/N*100), end=' ')
        d10fname = osp.join(dirpath, key+'.D10')
        d10_connect_table[key] = d10fname
        if not osp.exists(d10fname):
            with open(d10fname, 'w') as csvfile:
                writer = csv.writer(csvfile, lineterminator='\n')
                writer.writerows(d10data[key])
    print("\rProcessing D10 input file for cell %d of %d (%0.1f%%)." %
          (i+1, N, (i+1)/N*100))

    # Write the soil and design input file :

    N = len(d11data)
    for i, key in enumerate(d11data.keys()):
        print("\rProcessing D11 input file for cell %d of %d (%0.1f%%)." %
              (i+1, N, (i+1)/N*100), end=' ')
        d11fname = osp.join(dirpath, key+'.D11')
        d11_connect_table[key] = d11fname
        if not osp.exists(d11fname):
            with open(d11fname, 'w') as csvfile:
                writer = csv.writer(csvfile, lineterminator='\n')
                writer.writerows(d11data[key])
    print("\rProcessing D11 input file for cell %d of %d (%0.1f%%)." %
          (i+1, N, (i+1)/N*100))

    filename = osp.join(dirpath, 'connect_tables.npy')
    conn_tbl = {'D10': d10_connect_table, 'D11': d11_connect_table}
    np.save(filename, conn_tbl)


# ---- Meteo

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


def generate_helpfile_from_meteogrid(dirname_out, suffix='cnb'):
    if not osp.exists(dirname_out):
        os.makedirs(dirname_out)

    # ---- Generate HELP D13 for Solar Rad

    d13fname = osp.join(dirname_out, 'MTLINTL.D13')
    if not osp.exists(d13fname):
        print('Generating HELP D13 file for Solar Radiation...', end=' ')
        filename = "C:/Users/jsgosselin/HELP/RADEAU2/CWEEDS/94792.WY2"
        daily_wy2 = read_cweeds_file(filename, format_to_daily=True)

        filename = ("C:/Users/jsgosselin/HELP/RADEAU2/CWEEDS/"
                    "CAN_QC_MONTREAL-INTL-A_7025251_CWEEDS2011_T_N.WY3")
        daily_wy3 = read_cweeds_file(filename, format_to_daily=True)
        wy23_df = join_daily_cweeds_wy2_and_wy3(daily_wy2, daily_wy3)

        indexes = np.where((wy23_df['Years'] >= 2000) &
                           (wy23_df['Years'] <= 2015))[0]
        save_solrad_to_HELP(d13fname, wy23_df['Years'][indexes],
                            wy23_df['Irradiance'][indexes],
                            'CAN_QC_MONTREAL-INTL-A_7025251',
                            wy23_df['Latitude'])

    # ---- Read meteo grid

    print('Initializing netcdf manager...', end=' ')
    path_netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily"
    meteo_manager = NetCDFMeteoManager(path_netcdf)
    print('done')

    # ---- Read HELP grid shapefile

    print('Reading HELP grid shapefile...', end=' ')
    ref_crs = ("+proj=longlat +ellps=GRS80 +datum=NAD83 "
               "+towgs84=0,0,0,0,0,0,0 +no_defs")
    help_shp_fpath = ("C:/Users/jsgosselin/HELP/RADEAU2/grid_helpXYZ/"
                      "grid_helpXYZ.shp")
    shp_help = gpd.read_file(help_shp_fpath)
    print('done')

    # Convert the geometry coordinates to lat long from UTM18 :
    if shp_help.crs != {'init': 'epsg:4269'}:
        print('Converting HELP grid to lat/lon...')
        shp_help = shp_help.to_crs(ref_crs)
        shp_help.to_file(help_shp_fpath)
        print('done')

    print('Fetching the lat/lon values from the shapefile...', end=' ')
    lat_help = [p.y for p in shp_help.geometry]
    lon_help = [p.x for p in shp_help.geometry]
    cid_help = np.array(shp_help['id'])
    print('done')

    d4_linktable = {}
    d7_linktable = {}
    d13_linktable = {}
    for i in range(len(cid_help)):
        progress = (i+1)/len(cid_help)*100
        msg = ("\rGenerating HELP D4 and D7 files for Meteo..."
               "%d/%d (%0.1f%%)" % (i, len(cid_help), progress))
        print(msg, end=' ')

        lat_idx, lon_idx = meteo_manager.get_idx_from_latlon(
                lat_help[i], lon_help[i])
        d4fname = osp.join(dirname_out, '%03d_%03d.D4' % (lat_idx, lon_idx))
        d7fname = osp.join(dirname_out, '%03d_%03d.D7' % (lat_idx, lon_idx))

        d4_linktable[suffix+str(cid_help[i])] = d4fname
        d7_linktable[suffix+str(cid_help[i])] = d7fname
        d13_linktable[suffix+str(cid_help[i])] = d13fname
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
    print('\rGenerating HELP D4 and D7 files for Meteo... done')

    print('\rSaving the connectivity tables to file... ', end=' ')
    filename = osp.join(dirname_out, 'connect_tables.npy')
    conn_tbl = {'D13': d13_linktable, 'D4': d4_linktable, 'D7': d7_linktable}
    np.save(filename, conn_tbl)
    print('\rSaving the connectivity tables to file... done' + ' '*25)

    return d4_linktable, d7_linktable, d13_linktable


def savedata_to_hdf5(data, hdf5_filename):
    hdf5file = h5py.File(hdf5_filename, mode='w')
    for cid in data.keys():
        cellgrp = hdf5file.create_group(cid)
        for key in data[cid].keys():
            cellgrp.create_dataset(key, data=data[cid][key])
    hdf5file.close()


# %% if __name__ == '__main__'

if __name__ == '__main__':
    dirname_out = ("C:/Users/jsgosselin/HELP/help/.help_threads/"
                   "help_meteofiles")

    # filename = "C:/Users/jsgosselin/HELP/help/.help_threads/thread0/"
    # data164 = read_monthly_help_output(filename+"164.OUT")
    # data641 = read_monthly_help_output(filename+"641.OUT")

    if False:
        generate_helpfile_from_meteogrid(dirname_out, suffix='rad')
        d10data, d11data = read_d10d11_file("RADEAU2.D10", "RADEAU2.D11")
        write_d10d11_input_files("help_input_d10d11", d10data, d11data)

    if False:
        app = QApplication(sys.argv)
        help_thread_pool_mngr = HelpThreadPoolManager(nthread=7)
        help_thread_pool_mngr.load_help_D10D11_inputs(
                'RADEAU2.D10', 'RADEAU2.D11')
        help_thread_pool_mngr.load_meteo_connect_tables(
                osp.join(dirname_out, 'connect_tables.npy'))
        help_thread_pool_mngr.start_calculation()
        sys.exit(app.exec_())

    if False:
        outpath = "C:/Users/jsgosselin/HELP/help/HELP.OUT"
        help_reader = HelpResultReader(outpath)

    # dirpath_netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily"
    # netcdf_meteo_reader = NetCDFMeteoReader(dirpath_netcdf)

#    dataf, nlay, nsub = read_daily_help_output(
#            'C:/Users/jsgosselin/HELP/help/RCRA_TEST.OUT')
#

    # http://proj4.org/parameters.html
    # crs1 = "+proj=lcc +lon_0=-68.5 +lat_1=46 +lat_2=60 +lat_0=53 +ellps=GRS80 +datum=NAD83 +no_defs"
    # crs2 = "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

    # from geopandas import GeoDataFrame
    # import geopandas as gpd
    # import pandas as pd
    # from shapely.geometry import Point, Polygon
    # from itertools import product
    # import rasterio as rio
    # from rasterio import features
    # from affine import Affine

    # # print('\rReading netcdf dataset...', end=' ')
    # path__netcdf = "C:/Users/jsgosselin/MeteoGrilleDaily/GCQ_v2_2000.nc"
    # netcdf_dset = netCDF4.Dataset(path__netcdf, 'r+')

    # lat = np.array(netcdf_dset['lat'])
    # lon = np.array(netcdf_dset['lon'])
    # # tasmax = np.array(netcdf_dset['tasmax'])

    # # nt = netcdf_dset.dimensions['time'].size
    # # ny = netcdf_dset.dimensions['lat'].size
    # # nx = netcdf_dset.dimensions['lon'].size
    # # print('done')

    # # print('\rFormating the data in a shapefile...', end=' ')
    # # tasmax_yavg = np.empty(ny*nx).astype('float64')
    # # geometry = []
    # # i = 0
    # # dx = dy = 0.1/2
    # # for j, k in product(range(ny), range(nx)):
    # #     # point = Point((lon[k], lat[j]))
    # #     polygon = Polygon([(lon[k]-dx, lat[j]-dy),
    # #                        (lon[k]-dx, lat[j]+dy),
    # #                        (lon[k]+dx, lat[j]+dy),
    # #                        (lon[k]+dx, lat[j]-dy)])
    # #     geometry.append(polygon)
    # #     tasmax_yavg[i] = np.average(tasmax[:, j, k])
    # #     i += 1

    # # df = pd.DataFrame(data={'tasmax_yavg': tasmax_yavg})
    # crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"
    # # gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
    # # print('done')

    # # print('\rSaving to Shapefile...', end=' ')
    # # out = "C:/Users/jsgosselin/MeteoGrilleDaily/" + 'GCQ_v2_2000.shp'
    # # gdf.to_file(out)
    # # print('done')

    # # shp_fpath = "C:/Users/jsgosselin/MeteoGrilleDaily/" + 'GCQ_v2_2000.shp'
    # # shp = gpd.read_file(shp_fpath)
    
    # rst_fpath = "C:/Users/jsgosselin/MeteoGrilleDaily/" + 'test_meteogrille.tif'
    # rst = rio.open(rst_fpath)
    # meta = rst.meta.copy()
    # meta.update(compress='lzw')
    
    # meta = {'driver': 'GTiff', 
    #         'dtype': 'float64',
    #         'nodata': None,
    #         'width': 3000,
    #         'height': 3000,
    #         'count': 1,
    #         'crs': crs,
    #         'compress': 'lzw',
    #         'transform': Affine(0.00883333333333333, 0.0, -81.55,
    #                             0.0, -0.0066666666666666645, 62.949999999999996)}
    
    
    # rst_fpath = "C:/Users/jsgosselin/MeteoGrilleDaily/" + 'test3_raster.tif'
    # with rio.open(rst_fpath, 'w', **meta) as out:
    #     out_arr = out.read(1)
    #     shapes = ((geom, value) for geom, value in zip(shp.geometry, shp['tasmax_yav']))

    # # # # # # this is where we create a generator of geom, value pairs to use in rasterizing

    #     burned = features.rasterize(shapes=shapes, fill=0, out=out_arr,
    #                                 transform=out.transform)
    #     out.write_band(1, burned)
    
    # meta = rst.meta.copy()
    # meta.update(compress='lzw')
    # out_fn = "C:/Users/jsgosselin/MeteoGrilleDaily/" + 'template_raster.tif'
    # with rasterio.open(out_fn, 'w', **meta) as out:
    #     out_arr = out.read(1)

    #     # this is where we create a generator of geom, value pairs to use in rasterizing
    #     shapes = ((geom,value) for geom, value in zip(gdf.geometry, gdf.tasmax_yavg))

    #     burned = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    #     out.write_band(1, burned)
    
