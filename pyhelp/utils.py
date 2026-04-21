# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from matplotlib.axis import Axis

# ---- Standard imports
import os
import os.path as osp
import csv
from shutil import rmtree

# ---- Third party imports
import h5py
from matplotlib.transforms import Bbox
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


def get_ticklabel_extents(axis: Axis, renderer) -> tuple[Bbox, Bbox]:
    """
    Get the extents of the tick labels on either side of the axes.

    This function replaces the deprecated `axis.get_ticklabel_extents()`
    method that was removed in matplotlib 3.6.

    See: https://github.com/matplotlib/matplotlib/pull/23190

    Parameters
    ----------
    axis : matplotlib.axis.Axis
        The axis (xaxis or yaxis) to get tick label extents for
    renderer : matplotlib.backend_bases.RendererBase
        The renderer to use for calculating bounding boxes

    Returns
    -------
    bbox_primary : matplotlib.transforms.Bbox
        Bounding box for labels on the primary side
        (bottom for x-axis, left for y-axis)
    bbox_secondary : matplotlib.transforms.Bbox
        Bounding box for labels on the secondary side
        (top for x-axis, right for y-axis)

    Notes
    -----
    Returns empty bboxes (0, 0, 0, 0) instead of None when no labels exist.
    This maintains compatibility with code expecting bbox objects.
    """
    if not hasattr(axis, '_get_ticklabel_bboxes'):
        # This is required for backward compatiblity with matplotlib < 3.6.
        return axis.get_ticklabel_extents(renderer)

    # Use matplotlib's internal method (available in 3.6+).
    ticks_to_draw = axis._update_ticks()
    tlb1, tlb2 = axis._get_ticklabel_bboxes(ticks_to_draw, renderer)

    # Create bounding boxes (empty if no labels)
    bbox1 = Bbox.union(tlb1) if tlb1 else Bbox.from_extents(0, 0, 0, 0)
    bbox2 = Bbox.union(tlb2) if tlb2 else Bbox.from_extents(0, 0, 0, 0)

    return bbox1, bbox2
