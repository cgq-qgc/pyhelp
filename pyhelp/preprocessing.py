# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library Imports

import csv
import time
import os
import os.path as osp
from collections import OrderedDict
import multiprocessing as mp
from multiprocessing import Pool


# ---- Third Party imports

import numpy as np
import geopandas as gpd
import netCDF4
import xlrd


# ---- Evapotranspiration and Soil and Design data (D10 and D11)

def _read_data_from_excel(filename):
    """
    Read the evapotranspiration and soil and design data from an excel sheet.
    """
    with xlrd.open_workbook(filename, on_demand=True) as wb:
        sheet = wb.sheet_by_index(0)

        data = [sheet.row_values(rowx, start_colx=0, end_colx=None) for
                rowx in range(3, sheet.nrows)]
    return data


def _format_d11_singlecell(row, sf_edepth, sf_ulai):
    """
    Format the D11 input data for a single cell (one row in the excel file).
    """
    iu11 = 2
    city = row[0]
    ulat = float(row[3])
    ipl, ihv = int(row[9]), int(row[10])
    ulai = float(row[12]) * sf_ulai
    edepth = max(float(row[13]) * sf_edepth, 10)
    wind = float(row[4])
    hum1 = float(row[5])
    hum2 = float(row[6])
    hum3 = float(row[7])
    hum4 = float(row[8])

    d11dat = []

    # READ (11, 5050) IU11, CITY11
    # 5050 FORMAT (I2/A40)

    d11dat.append(['{0:>2}'.format(iu11)])
    d11dat.append(['{0:<40}'.format(city)])

    # READ (11, 5060) ULAT, IPL, IHV, ULAI, EDEPTH, WIND, HUM1, HUM2,
    #                 HUM3, HUM4
    # 5060 FORMAT (F10.2,I4,I4,F7.0,F8.0,5F5.0)

    d11dat.append(['{0:<10.2f}'.format(ulat) +
                   '{0:>4}'.format(ipl) +
                   '{0:>4}'.format(ihv) +
                   '{0:>7.2f}'.format(ulai) +
                   '{0:>8.1f}'.format(edepth) +
                   '{0:>5.1f}'.format(wind) +
                   '{0:>5.1f}'.format(hum1) +
                   '{0:>5.1f}'.format(hum2) +
                   '{0:>5.1f}'.format(hum3) +
                   '{0:>5.1f}'.format(hum4)])

    return d11dat


def _format_d10_singlecell(row):
    """
    Format the D10 input data for a single cell (one row in the excel file).
    """
    title = row[0]
    iu10 = 2
    ipre = 0
    irun = 1
    osno = 0     # initial snow water
    area = 6.25  # area projected on horizontal plane
    frunof = 100
    irun = 1
    runof = float(row[14])

    d10dat = []

    # READ (10, 5070) TITLE
    # 5070 FORMAT(A60)

    d10dat.append(['{0:<60}'.format(title)])

    # READ (10, 5080) IU10, IPRE, OSNO, AREA, FRUNOF, IRUN
    # 5080 FORMAT(I2,I2,2F10.0,F6.0,I2)

    d10dat.append(['{0:>2}'.format(iu10) +
                   '{0:>2}'.format(ipre) +
                   '{0:>10.0f}'.format(osno) +
                   '{0:>10.0f}'.format(area) +
                   '{0:>6.0f}'.format(frunof) +
                   '{0:>2}'.format(irun)])

    # IF (IRUN .EQ. 1) READ (10, 5090) CN2
    # 5090 FORMAT(F7.0)

    d10dat.append(['{0:>7.0f}'.format(runof)])

    # Format the layer properties.

    nlayers = int(row[11])
    layers = row[15:]
    for lay in range(nlayers):
        layer = int(layers.pop(0))
        thick = float(layers.pop(0))
        isoil = 0
        poro = float(layers.pop(0))
        fc = float(layers.pop(0))
        wp = float(layers.pop(0))
        sw = ''
        rc = float(layers.pop(0))
        xleng = float(layers.pop(0))
        slope = float(layers.pop(0))

        # READ (10, 5120) LAYER (J), THICK (J), ISOIL (J),
        #   PORO (J), FC (J), WP (J), SW (J), RC (J)
        # 5120 FORMAT(I2,F7.0,I4,4F6.0,F16.0)

        d10dat.append(['{0:>2}'.format(layer) +
                       '{0:>7.0f}'.format(thick) +
                       '{0:>4}'.format(isoil) +
                       '{0:>6.3f}'.format(poro) +
                       '{0:>6.3f}'.format(fc) +
                       '{0:>6.3f}'.format(wp) +
                       '{0:>6}'.format(sw) +
                       '{0:>16.14f}'.format(rc)])
        recir = subin = phole = defec = ipq = trans = ''
        layr = 0
        line_layr = ''

        # READ (10, 5130) XLENG (J), SLOPE (J), RECIR (J), LAYR (J),
        #   SUBIN (J), PHOLE (J), DEFEC (J), IPQ (J), TRANS (J)
        # 5130 FORMAT(F7.0,2F6.0,I3,F13.0,2F7.0,I2,G14.6)

        d10dat.append(['{0:>7.0f}'.format(xleng) +
                       '{0:>6.2f}'.format(slope) +
                       '{0:>6}'.format(recir) +
                       '{0:>3}'.format(layr) +
                       '{0:>13}'.format(subin) +
                       '{0:>7}'.format(phole) +
                       '{0:>7}'.format(defec) +
                       '{0:>2}'.format(ipq) +
                       '{0:>14}'.format(trans)])

    return d10dat


def format_d10d11_from_excel(filename, sf_edepth=1, sf_ulai=1):
    """
    Format the evapotranspiration (D11) and soil and design data (D11) in a
    format that is compatible by HELP.
    """
    print('\rReading D10 an D11 data from Excel file...', end=' ')
    data = _read_data_from_excel(filename)
    print('done')

    d11dat = {}
    d10dat = {}
    N = len(data)
    for i, row in enumerate(data):
        print("\rFormatting D10 and D11 data for cell %d of %d (%0.1f%%)" %
              (i+1, N, (i+1)/N*100), end=' ')

        cellname = row[0]
        d11dat[cellname] = _format_d11_singlecell(row, sf_edepth, sf_ulai)
        d10dat[cellname] = _format_d10_singlecell(row)

    print("\rFormatting D10 and D11 data for cell %d of %d (%0.1f%%)" %
          (i+1, N, (i+1)/N*100))

    return d10dat, d11dat


def read_concatenated_d10d11_file(path_d10file, path_d11file):
    """
    Read the concatenated D10 and D11 files that contain the D10 and D11 inputs
    for all the cells. Return a formatted dictionary where the data
    relative to each cell is saved as a list at the key that correspond to the
    unique id of the cell.
    """

    enc = 'iso-8859-1'

    # ---- Read and Format D11 File

    with open(path_d11file, 'r', encoding=enc) as csvfile:
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


def write_d10d11_singlecell(packed_data):
    fname, cid, d10data = packed_data
    if not osp.exists(fname):
        with open(fname, 'w') as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerows(d10data)
    return {cid: fname}


def write_d10d11_allcells(dirpath, d10data, d11data, ncore=None):
    """
    Write the content of each cell in individual D10 and D11 files.
    """
    ncore = max(mp.cpu_count() if ncore is None else ncore, 1)
    pool = Pool(ncore)

    # Prepare soil and design input files :

    tic = time.clock()
    iterable = [(osp.join(dirpath, str(cid) + '.D10'), cid, d10data[cid]) for
                cid in d10data.keys()]
    d10_connect_table = {}
    calcul_progress = 0
    N = len(iterable)
    for i in pool.imap_unordered(write_d10d11_singlecell, iterable):
        d10_connect_table.update(i)
        calcul_progress += 1
        progress_pct = calcul_progress/N*100
        print("\rCreating D10 input file for cell %d of %d (%0.1f%%)" %
              (calcul_progress, N, progress_pct), end=' ')
    tac = time.clock()
    print('\nTask completed in %0.2f sec' % (tac-tic))

    # Prepare evapotranspiration input files :

    tic = time.clock()
    iterable = [(osp.join(dirpath, str(cid) + '.D11'), cid, d11data[cid]) for
                cid in d10data.keys()]
    d11_connect_table = {}
    calcul_progress = 0
    N = len(iterable)
    for i in pool.imap_unordered(write_d10d11_singlecell, iterable):
        d11_connect_table.update(i)
        calcul_progress += 1
        progress_pct = calcul_progress/N*100
        print("\rCreating D11 input file for cell %d of %d (%0.1f%%)" %
              (calcul_progress, N, progress_pct), end=' ')
    tac = time.clock()
    print('\nTask completed in %0.2f sec' % (tac-tic))

    return d10_connect_table, d11_connect_table
