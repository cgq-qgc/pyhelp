# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library Imports
import os
import os.path as osp
from multiprocessing import Pool
import multiprocessing as mp
import time
import csv
import calendar

# ---- Third Party imports
import numpy as np

# ---- Local Libraries Imports
from pyhelp import HELP3O

DEL_TEMPFILES = True


# ---- Run HELP

def run_help_singlecell(item):
    """Run HELP for a single cell."""
    cellname, outparam = item
    HELP3O.run_simulation(*outparam)
    results = read_monthly_help_output(outparam[5])
    if DEL_TEMPFILES:
        os.remove(outparam[5])
    return (cellname, results)


def run_help_allcells(cellparams, ncore=None):
    """Run HELP in batch for multiple cells."""
    output = {}
    ncore = max(mp.cpu_count() if ncore is None else ncore, 1)
    tstart = time.perf_counter()
    calcul_progress = 0
    N = len(cellparams)
    pool = Pool(ncore)
    for cell in pool.imap_unordered(run_help_singlecell, cellparams.items()):
        output[cell[0]] = cell[1]
        calcul_progress += 1
        progress_pct = calcul_progress/N*100
        tpassed = time.perf_counter() - tstart
        tremain = (100-progress_pct)*tpassed/progress_pct/60
        print(('\rHELP simulation in progress: %3.1f%% (%0.1f min remaining)'
               "     ") % (progress_pct, tremain), end='')
    calcul_time = (time.perf_counter() - tstart)
    print('\nTask completed in %0.2f sec' % calcul_time)

    return output


# ---- Read HELP output

def read_monthly_help_output(filename):
    """
    Read the monthly output from .OUT HELP file and return the data as
    numpy arrays stored in a dictionary. Support the output format that was
    modified from HELP 3.07 (see PR#2).
    """
    with open(filename, 'r') as csvfile:
        csvread = list(csv.reader(csvfile))

    arr_years = []
    vstack_precip = []
    vstack_runoff = []
    vstack_evapo = []
    vstack_subrun1 = []
    vstack_subrun2 = []
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
            subrun1 = None
            subrun2 = np.zeros(12).astype('float32')
            percol = None
            while True:
                i += 1
                if len(csvread[i]) == 0:
                    continue
                line = csvread[i][0]

                if '**********' in line:
                    break
                elif 'PRECIPITATION' in line:
                    precip = np.array(line.split()[-12:]).astype('float32')
                elif 'RUNOFF' in line:
                    runoff = np.array(line.split()[-12:]).astype('float32')
                elif 'EVAPOTRANSPIRATION' in line:
                    evapo = np.array(line.split()[-12:]).astype('float32')
                elif 'LAT. DRAINAGE' in line:
                    if subrun1 is None:
                        subrun1 = np.array(
                            line.split()[-12:]).astype('float32')
                    else:
                        subrun2 += np.array(
                            line.split()[-12:]).astype('float32')
                elif 'PERCOLATION' in line:
                    if percol is None:
                        percol = np.array(line.split()[-12:]).astype('float32')
                    rechg = np.array(line.split()[-12:]).astype('float32')

            vstack_precip.append(precip)
            vstack_runoff.append(runoff)
            vstack_evapo.append(np.array(evapo).astype('float32'))
            vstack_rechg.append(np.array(rechg).astype('float32'))
            vstack_percol.append(np.array(percol).astype('float32'))
            if subrun1 is None:
                vstack_subrun1.append(np.zeros(12).astype('float32'))
            else:
                vstack_subrun1.append(subrun1)
            vstack_subrun2.append(subrun2)
        elif 'FINAL WATER STORAGE' in line:
            break

        i += 1

    data = {'years': np.array(arr_years).astype('uint16'),
            'precip': np.vstack(vstack_precip),
            'runoff': np.vstack(vstack_runoff),
            'evapo': np.vstack(vstack_evapo),
            'subrun1': np.vstack(vstack_subrun1),
            'subrun2': np.vstack(vstack_subrun2),
            'perco': np.vstack(vstack_percol),
            'rechg': np.vstack(vstack_rechg)}
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
            if 'DAILY OUTPUT FOR YEAR' in line:
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
