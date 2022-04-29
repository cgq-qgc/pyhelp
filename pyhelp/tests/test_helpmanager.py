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
import pandas as pd
import pytest
from pandas.api.types import is_string_dtype

# ---- Local library imports
from pyhelp import __rootdir__
from pyhelp.managers import HelpManager
from pyhelp.output import HelpOutput, VARNAMES

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
def helpm(output_dir):
    manager = HelpManager(
        workdir=output_dir,
        path_to_grid=INPUT_FILES['grid'],
        path_to_precip=INPUT_FILES['precip'],
        path_to_airtemp=INPUT_FILES['airtemp'],
        path_to_solrad=INPUT_FILES['solrad'])
    return manager


# =============================================================================
# ---- Tests
# =============================================================================
def test_read_input_grid(helpm, output_file):
    """Test that the input grid is read as expected."""
    assert is_string_dtype(helpm.grid.index)
    assert is_string_dtype(helpm.grid['cid'])


def test_calc_help_cells(helpm, output_file):
    """Test that the HelpManager is able to run water budget calculation."""
    cellnames = helpm.cellnames[:100]

    helpm.calc_help_cells(output_file, cellnames, tfsoil=-3)
    assert osp.exists(output_file)

    inputdir = osp.join(helpm.inputdir, 'd10d11_input_files')
    assert len(os.listdir(inputdir)) == 98 * 2
    inputdir = osp.join(helpm.inputdir, 'D4_input_files')
    assert len(os.listdir(inputdir)) == 2
    inputdir = osp.join(helpm.inputdir, 'D7_input_files')
    assert len(os.listdir(inputdir)) == 2
    inputdir = osp.join(helpm.inputdir, 'D13_input_files')
    assert len(os.listdir(inputdir)) == 1

    # Assert that the results are as expected.
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


@pytest.mark.parametrize('fig_title', [None, 'Exemple figure title'])
def test_plot_area_monthly_avg(output_dir, output_file, fig_title):
    """
    Test that the water budget plots are created and saved as expected.
    """
    output = HelpOutput(output_file)

    figfilename = osp.join(output_dir, 'area_monthly_avg.pdf')
    fig = output.plot_area_monthly_avg(
        figfilename, year_from=2003, year_to=2009, fig_title=fig_title)

    assert fig is not None
    assert osp.exists(figfilename)

    children = fig.axes[0].get_children()
    assert (children[0].get_ydata().sum() - 1086.8448950125246) < 0.01  # precip
    assert (children[1].get_ydata().sum() - 136.12068516893297) < 0.01  # rechg
    assert (children[2].get_ydata().sum() - 226.9635476845988) < 0.01   # runoff
    assert (children[3].get_ydata().sum() - 550.11140519815) < 0.01     # evapo
    assert (children[4].get_ydata().sum() - 47.97935577404126) < 0.01   # subrun1
    assert (children[5].get_ydata().sum() - 121.66539490800443) < 0.01  # subrun2

    if fig_title is None:
        assert fig._suptitle is None
    else:
        assert fig._suptitle.get_text() == fig_title


@pytest.mark.parametrize('fig_title', [None, 'Exemple figure title'])
def test_plot_area_yearly_series(output_dir, output_file, fig_title):
    """
    Test that plotting the yearly values is working expected.
    """
    output = HelpOutput(output_file)

    figfilename = osp.join(output_dir, 'area_yearly_series.pdf')
    fig = output.plot_area_yearly_series(
        figfilename, year_from=2003, year_to=2009, fig_title=fig_title)

    assert fig is not None
    assert osp.exists(figfilename)

    expected_xdata = [2003, 2004, 2005, 2006, 2007, 2008, 2009]

    children = fig.axes[0].get_children()
    for i in range(12):
        assert list(children[i].get_xdata()) == expected_xdata
    assert (children[0].get_ydata().mean() - 1086.8448950125246) < 0.01   # precip
    assert (children[2].get_ydata().mean() - 136.12068516893297) < 0.01   # rechg
    assert (children[4].get_ydata().mean() - 226.9635476845988) < 0.01    # runoff
    assert (children[6].get_ydata().mean() - 550.11140519815) < 0.01      # evapo
    assert (children[8].get_ydata().mean() - 47.97935577404126) < 0.01    # subrun1
    assert (children[10].get_ydata().mean() - 121.66539490800443) < 0.01  # subrun2

    if fig_title is None:
        assert fig._suptitle is None
    else:
        assert fig._suptitle.get_text() == fig_title


@pytest.mark.parametrize('fig_title', [None, 'Exemple figure title'])
def test_plot_area_yearly_avg(output_dir, output_file, fig_title):
    """
    Test that plotting the yearly averages is working expected.
    """
    output = HelpOutput(output_file)

    figfilename = osp.join(output_dir, 'area_yearly_avg.pdf')
    fig = output.plot_area_yearly_avg(
        figfilename, year_from=2003, year_to=2009, fig_title=fig_title)

    assert fig is not None
    assert osp.exists(figfilename)

    children = fig.axes[0].get_children()
    assert (children[0].get_height() - 1086.8448950125246) < 0.01   # precip
    assert children[1].get_text() == '1087\nmm/an'
    assert (children[2].get_height() - 136.12068516893297) < 0.01   # rechg
    assert children[3].get_text() == '136\nmm/an'
    assert (children[4].get_height() - 226.9635476845988) < 0.01    # runoff
    assert children[5].get_text() == '227\nmm/an'
    assert (children[6].get_height() - 550.11140519815) < 0.01      # evapo
    assert children[7].get_text() == '550\nmm/an'
    assert (children[8].get_height() - 47.97935577404126) < 0.01    # subrun1
    assert children[9].get_text() == '48\nmm/an'
    assert (children[10].get_height() - 121.66539490800443) < 0.01  # subrun2
    assert children[11].get_text() == '122\nmm/an'

    if fig_title is None:
        assert fig._suptitle is None
    else:
        assert fig._suptitle.get_text() == fig_title


def test_calc_cells_yearly_avg(output_file):
    """
    Test that the method to calculate cells yearly average values is
    working as expected.
    """
    output = HelpOutput(output_file)

    # Test calc_cells_yearly_avg without providing any value for the
    # year_from and year_to argument.
    yearly_avg = output.calc_cells_yearly_avg()
    expected_results = {
        'precip': 1055.8597373343132,
        'perco': 251.5843470406275,
        'evapo': 548.5954285729573,
        'rechg': 130.18263914385366,
        'runoff': 212.26214164981408,
        'subrun1': 46.26946108344004,
        'subrun2': 113.02721698440918}
    for varname in list(expected_results.keys()):
        result = np.sum(yearly_avg[varname]) / len(yearly_avg[varname])
        assert abs(result - expected_results[varname]) < 1, varname

    # Test calc_cells_yearly_avg with non null year_from and year_to argument.
    yearly_avg = output.calc_cells_yearly_avg(year_from=2003, year_to=2009)
    expected_results = {
        'precip': 1086.8448950125246,
        'perco': 259.85654385728577,
        'evapo': 550.11140519815,
        'rechg': 136.12068516893297,
        'runoff': 226.9635476845988,
        'subrun1': 47.97935577404126,
        'subrun2': 121.66539490800443}
    for varname in list(expected_results.keys()):
        result = np.sum(yearly_avg[varname]) / len(yearly_avg[varname])
        assert abs(result - expected_results[varname]) < 1, varname

    # Test calc_cells_yearly_avg with year_from == year_to.
    yearly_avg = output.calc_cells_yearly_avg(year_from=2003, year_to=2003)
    expected_results = {
        'precip': 1144.4142919267927,
        'perco': 324.15252048559057,
        'evapo': 492.4243657442988,
        'rechg': 148.72946963740077,
        'runoff': 164.32582637467374,
        'subrun1': 56.33706154407843,
        'subrun2': 140.96849990912418}
    for varname in list(expected_results.keys()):
        result = np.sum(yearly_avg[varname]) / len(yearly_avg[varname])
        assert abs(result - expected_results[varname]) < 1, varname


def test_save_output_to_csv(output_dir, output_file):
    """
    Test that saving yearly results to csv is working as expected.
    """
    output = HelpOutput(output_file)

    # Save yearly results to csv.
    csvfilename = osp.join(output_dir, 'test_help_yearly_results.csv')
    assert not osp.exists(csvfilename)
    output.save_to_csv(csvfilename)
    assert osp.exists(csvfilename)

    # Assert that the content of the csv is as expected.
    df = pd.read_csv(csvfilename, dtype={'cid': 'str'})
    df = df.set_index('cid', drop=True)
    assert list(df.columns) == ['lat_dd', 'lon_dd'] + VARNAMES
    assert df.index.name == 'cid'
    assert len(df) == 98

    assert df.index[0] == output.data['cid'][0]
    assert df.iloc[0]['lat_dd'] == output.data['lat_dd'][0]
    assert df.iloc[0]['lon_dd'] == output.data['lon_dd'][0]

    expected_results = {
        'precip': 1055.86,
        'perco': 251.59,
        'evapo': 548.58,
        'rechg': 130.19,
        'runoff': 212.26,
        'subrun1': 46.27,
        'subrun2': 113.02}
    for key in list(expected_results.keys()):
        result = df[key].sum() / len(df)
        expected_result = expected_results[key]
        assert abs(result - expected_result) < 1, key


if __name__ == '__main__':
    pytest.main(['-x', __file__, '-v', '-rw'])
