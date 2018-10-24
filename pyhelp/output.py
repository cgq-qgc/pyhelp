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
import numpy as np
import h5py


class HelpOutput(Mapping):
    def __init__(self, help_output):
        super(HelpOutput, self).__init__()
        if isinstance(help_output, dict):
            self.store = help_output
        elif isinstance(help_output, str) and osp.exists(help_output):
            self.store = h5py.File(help_output, mode='r+')

        self.avg_monthly = self.calc_area_monthly_avg()

    def __getitem__(self):
        pass

    def __iter__(self):
        pass

    def __len__(self):
        pass

