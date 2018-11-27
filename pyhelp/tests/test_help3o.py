# -*- coding: utf-8 -*-

# Copyright © 2018 PyHelp Project Contributors
# https://github.com/jnsebgosselin/gwhat
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library Imports
import os
import os.path as osp


# ---- Third party imports
import numpy as np
import pytest


# ---- Local library Imports
from pyhelp import HELP3O
from pyhelp import __rootdir__
from pyhelp.processing import run_help_singlecell


# ---- Fixtures
@pytest.fixture(scope="module")
def test_folder():
    return osp.join(osp.dirname(__rootdir__), 'pyhelp', 'tests')


@pytest.fixture(scope="module")
def rca_folder(test_folder):
    return osp.join(test_folder, 'rca_original_testcase_1997')


@pytest.fixture(scope="module")
def rca_params(rca_folder):
    daily_out = 0
    monthly_out = 1
    yearly_out = 0
    summary_out = 0
    tfsoil = 32.0  # Must be in Fahrenheit
    unit_system = 2  # IP if 1 else SI
    simu_nyear = 3
    return (osp.join(rca_folder, 'RCRA.D4'),
            osp.join(rca_folder, 'RCRA.D7'),
            osp.join(rca_folder, 'RCRA.D13'),
            osp.join(rca_folder, 'RCRA.D11'),
            osp.join(rca_folder, 'RCRA.D10'),
            osp.join(rca_folder, 'NEW_RCA.OUT'),
            daily_out,
            monthly_out,
            yearly_out,
            summary_out,
            unit_system,
            simu_nyear,
            tfsoil)


# ---- Test HelpManager
def test_run_help3o(rca_params):
    """
    Test that the HELP3O extension run and create an output file as expected.
    """
    if osp.exists(rca_params[5]):
        os.remove(rca_params[5])
    assert not osp.exists(rca_params[5])

    HELP3O.run_simulation(*rca_params)
    assert osp.exists(rca_params[5])


