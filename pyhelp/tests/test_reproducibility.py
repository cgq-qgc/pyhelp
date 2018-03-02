# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

# ---- Standard Library Imports

import os
import os.path as osp

# ---- Third Party Imports

import numpy as np
import pytest

# ---- Local Library Imports

from pyhelp import HELP3O
from pyhelp.help_utils import read_monthly_help_output

# ---- Standard library imports


def test_help_reproducibility():
    dirname = "C:\\users\\jsgosselin\\pyhelp\\pyhelp\\tests"
    fpath_precip = osp.join(dirname, "RCRA.D4")
    fpath_tasavg = osp.join(dirname, "RCRA.D7")
    fpath_solrad = osp.join(dirname, "RCRA.D13")
    fpath_evapo = osp.join(dirname, "RCRA.D11")
    fpath_soil = osp.join(dirname, "RCRA.D10")
    fpath_output = osp.join(dirname, "rcra_results.out")

    daily_out = 0
    monthly_out = 1
    yearly_out = 0
    summary_out = 0

    unit_system = 1  # IP if 1 else SI
    simu_nyear = 3

    if osp.exists(fpath_output):
        os.remove(fpath_output)

    HELP3O.run_simulation(fpath_precip, fpath_tasavg, fpath_solrad,
                          fpath_evapo, fpath_soil, fpath_output, daily_out,
                          monthly_out, yearly_out, summary_out, unit_system,
                          simu_nyear)

    # Assert that the output file was correctly generated :

    assert(osp.exists(fpath_output))

    # Assert that the monthly results are the same :

    fpath_expected_output = osp.join(dirname, "rcra_results.out")
    expected_results = read_monthly_help_output(fpath_expected_output)
    results = read_monthly_help_output(fpath_output)

    assert np.array_equal(expected_results['years'], results['years'])

    for var in ['rain', 'runoff', 'evapo', 'sub-runoff', 'percolation',
                'recharge']:
        max_err = np.max(np.abs(expected_results[var] - results[var]))
        assert max_err <= 0.01


if __name__ == '__main__':
    pytest.main([os.path.basename(__file__)])
