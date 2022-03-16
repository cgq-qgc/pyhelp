# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHELP Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------

# ---- Standard library imports
import os
import os.path as osp


# ---- Third party imports
import numpy as np
import pytest


# ---- Local library imports
from pyhelp import __rootdir__
from pyhelp.managers import HelpManager
from pyhelp.output import HelpOutput

EXAMPLE_FOLDER = osp.join(osp.dirname(__rootdir__), 'example')
INPUT_FILES = {
    'airtemp': osp.join(EXAMPLE_FOLDER, 'airtemp_input_data.csv'),
    'precip': osp.join(EXAMPLE_FOLDER, 'precip_input_data.csv'),
    'solrad': osp.join(EXAMPLE_FOLDER, 'solrad_input_data.csv'),
    'grid': osp.join(EXAMPLE_FOLDER, 'input_grid.csv')}


# =============================================================================
# ---- Fixtures
# =============================================================================
@pytest.fixture(scope="module")
def output_dir(tmp_path_factory):
    return tmp_path_factory.getbasetemp()


@pytest.fixture(scope="module")
def output_file(output_dir):
    return osp.join(output_dir, 'help_example.out')


@pytest.fixture
def helpm():
    manager = HelpManager(EXAMPLE_FOLDER)
    return manager


# =============================================================================
# ---- Tests
# =============================================================================
def test_calc_help_cells(helpm, output_file):
    """Test that the HelpManager is able to run water budget calculation."""
    cellnames = helpm.cellnames[:100]
    helpm.build_help_input_files()
    helpm.calc_help_cells(output_file, cellnames, tfsoil=-3)
    assert osp.exists(output_file)


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


def test_plot_water_budget(output_dir, output_file):
    """
    Test that the water budget plots are created and saved as expected.

    Regression test for Issue #29.
    """
    output = HelpOutput(output_file)

    figfilename = osp.join(output_dir, 'area_monthly_avg.pdf')
    output.plot_area_monthly_avg(figfilename)
    assert osp.exists(figfilename)

    figfilename = osp.join(output_dir, 'area_yearly_avg.pdf')
    output.plot_area_yearly_avg(figfilename)
    assert osp.exists(figfilename)

    figfilename = osp.join(output_dir, 'area_yearly_series.pdf')
    output.plot_area_yearly_series(figfilename)
    assert osp.exists(figfilename)


if __name__ == '__main__':
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
