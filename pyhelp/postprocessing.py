# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Third Party imports

import numpy as np
import time


def add_yrly_avg_to_shp(shp, outdat):
    tic = time.clock()

    cellnames = list(outdat.keys())
    years = outdat[cellnames[0]]['years'].value
    nyears = len(years)
    ncells = len(cellnames)

    avg_rain = np.zeros(ncells)
    avg_runoff = np.zeros(ncells)
    avg_evapo = np.zeros(ncells)
    avg_perco = np.zeros(ncells)
    avg_subrun1 = np.zeros(ncells)
    avg_subrun2 = np.zeros(ncells)
    avg_recharge = np.zeros(ncells)

    for i, cellname in enumerate(cellnames):
        print("\rProcessing cell %d of %d..." % (i+1, ncells), end=' ')
        data = outdat[cellname]
        avg_rain[i] = np.sum(data['rain'].value) / nyears
        avg_runoff[i] = np.sum(data['runoff'].value) / nyears
        avg_evapo[i] = np.sum(data['evapo'].value) / nyears
        if shp['context'][cellname] == 0:
            # Handle surface water results.
            avg_perco[i] = 0
            avg_subrun1[i] = 0
            avg_subrun2[i] = 0
            avg_recharge[i] = 0
        else:
            # Handle HELP results.
            avg_perco[i] = np.sum(data['percolation'].value) / nyears
            avg_subrun1[i] = np.sum(data['subrun1'].value) / nyears
            avg_subrun2[i] = np.sum(data['subrun2'].value) / nyears
            avg_recharge[i] = np.sum(data['recharge'].value) / nyears
            if shp['context'][cellname] == 2:
                # Convert recharge to runoff when cells are close to a stream.
                if avg_subrun2[-1] == 0:
                    # Convert recharge as surficial subrunoff.
                    avg_subrun1[-1] = avg_subrun1[-1] + avg_recharge[-1]
                else:
                    # This means there is a layer of sand above the clay layer.
                    # Convert recharge as deep runoff.
                    avg_subrun2[-1] = avg_subrun2[-1] + avg_recharge[-1]
                avg_recharge[-1] = 0

    # Insert yearly average values in the grid.
    shp.loc[cellnames, 'rain'] = avg_rain
    shp.loc[cellnames, 'runoff'] = avg_runoff
    shp.loc[cellnames, 'evapo'] = avg_evapo

    shp.loc[cellnames, 'percolation'] = avg_perco
    shp.loc[cellnames, 'subrun1'] = avg_subrun1
    shp.loc[cellnames, 'subrun2'] = avg_subrun2
    shp.loc[cellnames, 'recharge'] = avg_recharge
    print("\rProcessing cell %d of %d... done in %0.1fs" %
          (i+1, ncells, (time.clock() - tic)))
    return shp
