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
from pyhelp.output import HelpOutput


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


@pytest.fixture(scope="module")
def output_file(example_folder):
    return osp.join(example_folder, 'help_example.out')


@pytest.fixture
def helpm(example_folder):
    manager = HelpManager(example_folder, year_range=(2000, 2010))
    return manager


# ---- Test HelpManager
def test_autoread_input(helpm):
    """
    Test that the input files are read automatically when instantiating
    the HelpManager.
    """
    assert helpm.precip_data is not None
    assert helpm.airtemp_data is not None
    assert helpm.solrad_data is not None
    assert helpm.grid is not None


def test_calc_help_cells(helpm, output_file):
    """Test that the HelpManager is able to run water budget calculation."""
    cellnames = helpm.cellnames[:100]
    helpm.build_help_input_files()
    helpm.calc_help_cells(output_file, cellnames, tfsoil=-3)
    assert osp.exists(output_file)


# ---- Test HelpOutput
def test_validate_results(output_file):
    """Test that the water budget results are as expected. """
    output = HelpOutput(output_file)
    area_yrly_avg = output.calc_area_yearly_avg()
    expected_results = {'precip': 11614.46,
                        'perco': 2767.51,
                        'evapo': 6034.42,
                        'rechg': 1432.13,
                        'runoff': 2334.89,
                        'subrun1': 509.02,
                        'subrun2': 1243.24}
    for key in list(expected_results.keys()):
        assert abs(np.sum(area_yrly_avg[key]) - expected_results[key]) < 1, key


def test_plot_water_budget(output_file, example_folder):
    """
    Test that the water budget plots are created and saved as expected.

    Regression test for Issue #29.
    """
    output = HelpOutput(output_file)

    figfilename = osp.join(example_folder, 'area_monthly_avg.pdf')
    output.plot_area_monthly_avg(figfilename)
    assert osp.exists(figfilename)

    figfilename = osp.join(example_folder, 'area_yearly_avg.pdf')
    output.plot_area_yearly_avg(figfilename)
    assert osp.exists(figfilename)

    figfilename = osp.join(example_folder, 'area_yearly_series.pdf')
    output.plot_area_yearly_series(figfilename)
    assert osp.exists(figfilename)


if __name__ == '__main__':
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
