# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports
import os
import os.path as osp
import csv
from shutil import rmtree

# ---- Third party imports
import h5py
import numpy as np


def delete_folder_recursively(dirpath):
    """Try to delete all files and sub-folders below the given dirpath."""
    for filename in os.listdir(dirpath):
        filepath = os.path.join(dirpath, filename)
        try:
            rmtree(filepath)
        except OSError:
            os.remove(filepath)


def savedata_to_hdf5(data, hdf5_filename, grid=None):
    """
    Save the data contains in a dictionary into a hdf5 file.
    """
    hdf5file = h5py.File(hdf5_filename, mode='w')
    for cid in data.keys():
        cellgrp = hdf5file.create_group(str(cid))
        for key in data[cid].keys():
            cellgrp.create_dataset(key, data=data[cid][key])

        if grid is not None:
            # Add the grid data to the group.
            row = grid.loc[cid]
            for key in list(row.keys()):
                cellgrp.attrs[key] = row[key]
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


def save_content_to_csv(fname, fcontent, mode='w', delimiter=',',
                        encoding='utf8'):
    """
    Save fcontent in a csv file with the specifications provided
    in arguments.
    """
    with open(fname, mode, encoding='utf8') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter, lineterminator='\n')
        writer.writerows(fcontent)


def calc_dist_from_coord(lat1, lon1, lat2, lon2):
    """
    Compute the  horizontal distance in km between a location given in
    decimal degrees and a set of locations also given in decimal degrees.
    """
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)

    r = 6373  # r is the Earth radius in km

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    return r * c
