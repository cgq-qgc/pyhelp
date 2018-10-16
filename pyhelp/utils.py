# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Third party imports
import h5py
import numpy as np


def savedata_to_hdf5(data, hdf5_filename):
    """
    Save the data contains in a dictionary into a hdf5 file.
    """
    hdf5file = h5py.File(hdf5_filename, mode='w')
    for cid in data.keys():
        cellgrp = hdf5file.create_group(str(cid))
        for key in data[cid].keys():
            cellgrp.create_dataset(key, data=data[cid][key])
    hdf5file.close()


def nan_as_text_tolist(arr):
    """
    Convert the float nan to text while converting a numpy 2d array to a
    list, so that it is possible to save to an Excel file.
    """
    if np.isnan(arr).any():
        m, n = np.shape(arr)
        list_ = []
        for i in range(m):
            list_.append(['nan' if np.isnan(x) else x for x in arr[i, :]])
    else:
        list_ = arr.tolist()
    return list_


