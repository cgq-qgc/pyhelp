# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:37:59 2018
@author: jsgosselin
"""

# ---- Standard library imports

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

import numpy as np
import h5py
from PyQt5.QtCore import pyqtSlot as QSlot
from PyQt5.QtCore import pyqtSignal as QSignal
from PyQt5.QtCore import Qt, QObject, QThread
from PyQt5.QtWidgets import QApplication

DIRNAME = '.help_threads'
HELPNAME = 'H3OFS32F.exe'


def run_help():
    if osp.exists(HELPNAME):
        subprocess.call(HELPNAME)


class HelpThreadPoolManager(QObject):
    """
    An object that parallelize the HELP calculation and post-processing
    over multiple threads.
    """

    def __init__(self, nthread=None, path_hdf5=None):
        super(QObject, self).__init__()
        self._thread_pool = []
        self._worker_pool = []
        self._hdf5file = None

        self.d10data = None
        self.d11data = None
        self.cellnames = []

        # Output from HELP :

        self._daily_out = True
        self._monthly_out = False
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
        self._hdf5file = h5py.File(self.path_hdf5, mode='w')

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
        self.cellnames = self.cellnames[:16]
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

    @QSlot(str, str)
    def _handle_singlecell_result(self, filename, cellname):
        data = read_daily_help_output(filename)

        if 'years' not in list(self._hdf5file.keys()):
            self._hdf5file.create_dataset('years', data=data['years'])
        if 'days' not in list(self._hdf5file.keys()):
            self._hdf5file.create_dataset('days', data=data['days'])
        cellgrp = self._hdf5file.create_group(cellname)
        cellgrp.create_dataset('rain', data=data['rain'])
        cellgrp.create_dataset('runoff', data=data['runoff'])
        cellgrp.create_dataset('et', data=data['et'])
        cellgrp.create_dataset('ezone', data=data['ezone'])
        cellgrp.create_dataset('head1', data=data['head1'])
        cellgrp.create_dataset('drain1', data=data['drain1'])
        cellgrp.create_dataset('leak1', data=data['leak1'])
        self._hdf5file.flush()

        os.remove(filename)

        self.__calcul_progress += 1
        progress_pct = self.__calcul_progress/len(self.cellnames)*100
        print('\r%0.1f%%' % progress_pct, end='')
        if self.__calcul_progress == len(self.cellnames):
            calcul_time = (time.clock() - self.__start_calcul_time)
            print('\nCalculation time: %0.2fs\n' % calcul_time)
            self._hdf5file.close()
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
        outputparam = [['DATA4.D4'],
                       ['DATA7.D7'],
                       ['DATA13.D13'],
                       ['DATA11.D11'],
                       ['DATA10.D10'],
                       [cellname+'.OUT'],
                       [2],
                       [5],
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


def read_daily_help_output(filename):
    """
    Read the daily output from .OUT HELP file and return the data as
    numpy arrays stored in a dictionary.
    """
    with open(filename, 'r') as csvfile:
        csvread = list(csv.reader(csvfile))

    arr_years = []
    arr_days = []
    arr_rain = []
    arr_ru = []
    arr_et = []
    arr_ezone = []
    arr_head1 = []
    arr_drain1 = []
    arr_leak1 = []

    year = None
    for i, line in enumerate(csvread):
        if line:
            line = line[0]
            if 'DAILY OUTPUT FOR YEAR' in line:
                year = int(line.split()[-1])
                days_in_year = 366 if calendar.isleap(year) else 365
            elif len(line) == 70 and year is not None:
                try:
                    day = int(line[2:5])
                    rain = float(line[13:19])
                    ru = float(line[19:26])
                    et = float(line[26:33])
                    ezone = float(line[33:41])
                    head1 = float(line[41:51])
                    drain1 = float(line[51:61])
                    leak1 = float(line[61:71])
                except ValueError:
                    pass
                else:
                    arr_years.append(year)
                    arr_days.append(day)
                    arr_rain.append(rain)
                    arr_ru.append(ru)
                    arr_et.append(et)
                    arr_ezone.append(ezone)
                    arr_head1.append(head1)
                    arr_drain1.append(drain1)
                    arr_leak1.append(leak1)
                    if day == days_in_year:
                        year = None

    dataf = {'years': np.array(arr_years).astype('uint8'),
             'days': np.array(arr_days).astype('uint8'),
             'rain': np.array(arr_rain).astype('float32'),
             'runoff': np.array(arr_ru).astype('float32'),
             'et': np.array(arr_et).astype('float32'),
             'ezone': np.array(arr_ezone).astype('float32'),
             'head1': np.array(arr_head1).astype('float32'),
             'drain1': np.array(arr_drain1).astype('float32'),
             'leak1': np.array(arr_leak1).astype('float32')
             }
    return dataf


def read_d10d11_file(path_d10file, path_d11file):
    """
    Read the concatenated D10 and D11 files that contain the D10 and D11 inputs
    for all the cells.
    """

    # ---- Read and Format D11 File

    with open(path_d11file, 'r') as csvfile:
        d11reader = list(csv.reader(csvfile))

    d11dat = OrderedDict()
    cellnames = []
    for i in range(0, len(d11reader), 3):
        cell_name, zone_name = d11reader[i+1][0].split()
        d11dat[cell_name] = [d11reader[i],
                             d11reader[i+1],
                             d11reader[i+2]]
        cellnames.append(cell_name)

    # ---- Read and Format D10 File

    with open(path_d10file, 'r') as csvfile:
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


# %% if __name__ == '__main__'

if __name__ == '__main__':
    app = QApplication(sys.argv)
    help_thread_pool_mngr = HelpThreadPoolManager(nthread=8)
    help_thread_pool_mngr.load_help_D10D11_inputs('H3b.D10', 'H3b.D11')
    help_thread_pool_mngr.start_calculation()
    sys.exit(app.exec_())

    # dataf = read_daily_help_output('C:/Users/jsgosselin/HELP/HELP'
                                   # '/.help_threads/thread2/H3b4.OUT')
