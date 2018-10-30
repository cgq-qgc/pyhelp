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
from pyhelp.managers import HelpManager


# ---- Fixtures

@pytest.fixture(scope="module")
def example_folder():
    return osp.join(osp.dirname(__rootdir__), 'example')


@pytest.fixture(scope="module")
def input_files(example_folder):
    return {'airtemp': osp.join(example_folder, 'airtemp_input_data.csv'),
            'precip': osp.join(example_folder, 'precip_input_data.csv'),
            'solrad': osp.join(example_folder, 'solrad_input_data.csv'),
            'grid': osp.join(example_folder, 'input_grid.csv')
            }


@pytest.fixture
def helpm(example_folder):
    manager = HelpManager(example_folder, year_range=(2000, 2010))
    return manager


# ---- Tests
def test_autoread_input(helpm):
    """
    Test that the input files are read automatically when instantiating
    the HelpManager.
    """
    assert helpm.precip_data is not None
    assert helpm.airtemp_data is not None
    assert helpm.solrad_data is not None
    assert helpm.grid is not None


def test_calc_help_cells(helpm):
    """
    Test that the HelpManager is able to run calculation.
    """
    cellnames = helpm.cellnames[:100]
    help_output_hdf5 = osp.join(helpm.workdir, 'help_example.out')
    helpm.calc_help_cells(help_output_hdf5, cellnames)

    assert osp.exists(help_output_hdf5)


if __name__ == '__main__':
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
