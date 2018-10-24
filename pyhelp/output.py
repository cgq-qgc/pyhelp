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


