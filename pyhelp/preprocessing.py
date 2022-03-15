# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------


# ---- Standard Library Imports

import csv
import time
import os.path as osp
import multiprocessing as mp
from multiprocessing import Pool


MINEDEPTH = 3
MAXEDEPTH = 80
MINTHICK = 10


# ---- Evapotranspiration and Soil and Design data (D10 and D11)

def _format_d11_singlecell(row, sf_edepth, sf_ulai):
    """
    Format the D11 input data for a single cell (one row in the excel file).
    """
    nlayers = int(row['nlayer'])
    if nlayers == 0:
        # This means this cell cannot be run in HELP.
        return None

    iu11 = 2
    try:
        city = str(int(row['cid']))
    except ValueError:
        city = str(row['cid'])
    ulat = float(row['lat_dd'])
    ipl = int(row['growth_start'])
    ihv = int(row['growth_end'])
    ulai = float(row['LAI']) * sf_ulai
    edepth = float(row['EZD']) * sf_edepth
    edepth = min(max(edepth, MINEDEPTH), MAXEDEPTH)
    wind = float(row['wind'])
    hum1 = float(row['hum1'])
    hum2 = float(row['hum2'])
    hum3 = float(row['hum3'])
    hum4 = float(row['hum4'])

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
    Format the D10 input data for a single cell (corresponds to a single row
    in the input csv file).
    """
    nlayers = int(row['nlayer'])
    if nlayers == 0:
        # This means this cell cannot be run in HELP.
        return None
    try:
        title = str(int(row['cid']))
    except ValueError:
        title = str(row['cid'])
    iu10 = 2
    ipre = 0
    irun = 1
    osno = 0     # initial snow water
    area = 6.25  # area projected on horizontal plane
    frunof = 100
    runof = float(row['CN'])

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
    for i in range(nlayers):
        lay = str(i+1)
        layer = int(row['lay_type'+lay])
        thick = max(float(row['thick'+lay]), MINTHICK)
        isoil = 0
        poro = float(row['poro'+lay])
        fc = float(row['fc'+lay])
        wp = float(row['wp'+lay])
        sw = ''
        rc = float(row['ksat'+lay])
        xleng = float(row['dist_dr'+lay])
        slope = float(row['slope'+lay])

        # Check that all values are valid for the layer.
        check = [val == -9999 for val in
                 (thick, poro, fc, wp, rc, xleng, slope)]
        if any(check):
            return None

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


def format_d10d11_inputs(grid, cellnames, sf_edepth=1, sf_ulai=1):
    """
    Format the evapotranspiration (D11) and soil and design data (D10) in a
    format that is compatible with HELP.
    """
    tic = time.perf_counter()
    d11dat = {}
    d10dat = {}
    N = len(cellnames)
    for i, cid in enumerate(cellnames):
        print("\rFormatting D10 and D11 data for cell %d of %d (%0.1f%%)" %
              (i+1, N, (i+1)/N*100), end=' ')

        row = grid.loc[cid]
        d11dat[cid] = _format_d11_singlecell(row, sf_edepth, sf_ulai)
        d10dat[cid] = _format_d10_singlecell(row)

    print("\rFormatting D10 and D11 data for cell %d of %d (%0.1f%%)" %
          (i+1, N, (i+1)/N*100))
    tac = time.perf_counter()
    print('Task completed in %0.2f sec' % (tac-tic))

    warnings = [cid for cid, val in d10dat.items() if val is None]
    if warnings:
        print('-' * 25)
        msg = "Warning: the data for "
        msg += "cell " if len(warnings) == 1 else "cells "
        msg += ", ".join(warnings)
        msg += " are not formatted correctly."
        print(msg)
        print('-' * 25)

    return d10dat, d11dat


def write_d10d11_singlecell(packed_data):
    """Write the content of cell in a D10 and D11 file."""
    fname, cid, d10data = packed_data
    if d10data is None:
        fname = None
    else:
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

    # Prepare soil and design input files (D10).

    tic = time.perf_counter()
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
    tac = time.perf_counter()
    print('\nTask completed in %0.2f sec' % (tac-tic))

    # Prepare evapotranspiration input files (D11).

    tic = time.perf_counter()
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
    tac = time.perf_counter()
    print('\nTask completed in %0.2f sec' % (tac-tic))

    return d10_connect_table, d11_connect_table
