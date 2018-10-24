# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library imports
import os.path as osp
from collections.abc import Mapping


# ---- Third party imports
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import pandas as pd
import numpy as np
import h5py


class HelpOutput(Mapping):
    def __init__(self, path_or_dict):
        super(HelpOutput, self).__init__()
        if isinstance(path_or_dict, dict):
            self.data = path_or_dict['data']
            self.grid = path_or_dict['grid']
        elif isinstance(path_or_dict, str) and osp.exists(path_or_dict):
            # Load the data from an HDF5 file saved on disk.
            hdf5 = h5py.File(path_or_dict, mode='r+')
            self.data = {}
            for key in list(hdf5['data'].keys()):
                if key == 'cid':
                    self.data[key] = hdf5['data'][key].value.tolist()
                else:
                    self.data[key] = hdf5['data'][key].value
            hdf5.close()

            # Load the grid from an HDF5 file saved on disk.
            self.grid = pd.read_hdf(path_or_dict, 'grid')
        else:
            self.data = None
            self.grid = None

    def __getitem__(self):
        pass

    def __iter__(self):
        pass

    def __len__(self):
        return len(self.data['cid'])

    def save_to_hdf5(self, path_to_hdf5):
        """Save the data and grid to a HDF5 file at the specified location."""
        print("Saving data to {}...".format(osp.basename(path_to_hdf5)),
              end=" ")

        # Save the data.
        hdf5file = h5py.File(path_to_hdf5, mode='w')
        datagrp = hdf5file.create_group('data')
        for key in list(self.data.keys()):
            if key == 'cid':
                # This is required to avoid a "TypeError: No conversion path
                # for dtype: dtype('<U5')".
                # See https://github.com/h5py/h5py/issues/289
                datagrp.create_dataset(
                    key, data=[np.string_(i) for i in self.data['cid']])
            else:
                datagrp.create_dataset(key, data=self.data[key])
        hdf5file.close()

        # Save the grid.
        self.grid.to_hdf(path_to_hdf5, key='grid', mode='a')

        print('done')

