# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Third Party imports

import h5py


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
