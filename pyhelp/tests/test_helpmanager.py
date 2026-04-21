# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHELP Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------

# ---- Standard library imports
import datetime
import os
import os.path as osp

# ---- Third party imports
import matplotlib
matplotlib.use("Agg")

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
    'grid': osp.join(EXAMPLE_FOLDER, 'input_grid.csv')
    }


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


@pytest.mark.parametrize(
    'write_help_input_files, write_help_output_files',
    [(True, True), (False, False)]
    )
def test_calc_help_cells(
        helpm, output_file, write_help_input_files, write_help_output_files
        ):
    """Test that the HelpManager is able to run water budget calculation."""
    cellnames = helpm.cellnames[:100]

    helpm.calc_help_cells(
        output_file, cellnames, tfsoil=-3,
        write_help_input_files=write_help_input_files,
        write_help_output_files=write_help_output_files
        )

    assert osp.exists(output_file)

    inputdir = osp.join(helpm.inputdir, 'd10d11_input_files')
    assert osp.exists(inputdir) == write_help_input_files
    if write_help_input_files:
        assert len(os.listdir(inputdir)) == 98 * 2

    inputdir = osp.join(helpm.inputdir, 'D4_input_files')
    assert len(os.listdir(inputdir)) == 2
    inputdir = osp.join(helpm.inputdir, 'D7_input_files')
    assert len(os.listdir(inputdir)) == 2
    inputdir = osp.join(helpm.inputdir, 'D13_input_files')
    assert len(os.listdir(inputdir)) == 1

    outputdir = osp.join(helpm.inputdir, 'HELP30_output_files')
    assert osp.exists(outputdir) == write_help_output_files
    if write_help_output_files:
        assert len(os.listdir(outputdir)) == 98

    daily_outdir = osp.join(helpm.inputdir, 'Daily_output_files')
    assert not osp.exists(daily_outdir)

    # Assert that the results are as expected.
    output = HelpOutput(output_file)

    area_yrly_avg = output.calc_area_yearly_avg()
    expected_results = {'precip': 11614.457110677447,
                        'perco': 2968.4506514697555,
                        'evapo': 5769.752614678167,
                        'rechg': 1531.6787247646182,
                        'runoff': 2348.9837706497565,
                        'subrun1': 561.5879183655345,
                        'subrun2': 1350.018408411377}
    for name, value in expected_results.items():
        result = np.sum(area_yrly_avg[name])
        assert abs(result - value) < 1, f'{name}: {result} vs {value}'


def test_write_daily_output(helpm, output_file):
    """
    Test that daily output are saved as expected when
    `write_daily_output` is set to `True`.
    """
    cellnames = helpm.cellnames[:100]

    helpm.calc_help_cells(
        output_file, cellnames, tfsoil=-3,
        write_daily_output=True,
        )

    daily_outdir = osp.join(helpm.inputdir, 'Daily_output_files')
    assert osp.exists(daily_outdir)
    assert len(os.listdir(daily_outdir)) == 98

    daily_data = pd.read_csv(
        osp.join(daily_outdir, '37728.csv'),
        index_col=0,
        parse_dates=True
        )

    assert list(daily_data.columns) == [
        'RAIN', 'RUNOFF', 'ET', 'E_ZONE_WATER', 'SNOW_SURF',
        'TAIR', 'TSOIL_SURF', 'TSOIL_EDEPTH', 'FROZEN_SOIL',
        'HEAD1_ON_LAY3', 'DRAIN1_FROM_LAY2', 'LEAK1_THROUGH_LAY3',
        'HEAD2_ON_LAY5', 'DRAIN2_FROM_LAY4', 'LEAK2_THROUGH_LAY5'
        ]

    assert daily_data.index[0] == datetime.datetime(2000, 1, 1)
    assert daily_data.index[-1] == datetime.datetime(2010, 12, 31)

    expected_sums = {
        'RAIN': 11567.8,
        'RUNOFF': 2808.92,
        'ET': 5357.29,
        'DRAIN1_FROM_LAY2': 0.0238081869138438,
        'LEAK1_THROUGH_LAY3': 3438.1293521999996,
        'DRAIN2_FROM_LAY4': 0.0074378913938849996,
        'LEAK2_THROUGH_LAY5': 3396.4980834,
        }
    for col, expected in expected_sums.items():
        actual = daily_data[col].values.sum()
        err = abs(expected - actual)
        assert err < 1, f"Mismatch in col '{col}': {actual}"

    expected_sums = {
        'DRAIN1_FROM_LAY2': 0.0238081869138438,
        'DRAIN2_FROM_LAY4': 0.0074378913938849996,
        }
    for col, expected in expected_sums.items():
        actual = daily_data[col].values.sum()
        err = abs(expected - actual)
        assert err < 0.01, f"Mismatch in col '{col}': {actual}"

    assert daily_data.FROZEN_SOIL.sum() == 1248

    expected_means = {
        'E_ZONE_WATER': 0.20224736187157788,
        'SNOW_SURF': 30.35667496266799,
        'TAIR': 5.921279243404678,
        'TSOIL_SURF': 1.6576630164260828,
        'TSOIL_EDEPTH': 2.5955550024888003
        }
    for col, expected in expected_means.items():
        actual = daily_data[col].values.sum()
        err = abs(expected == actual)
        assert err < 0.1, f"Mismatch in col '{col}': {actual}"

    expected_means = {
        'HEAD1_ON_LAY3': 0.31645661523145846,
        'HEAD2_ON_LAY5': 1.0583286809357888,
        }
    for col, expected in expected_means.items():
        actual = daily_data[col].values.sum()
        err = abs(expected == actual)
        assert err < 0.01, f"Mismatch in col '{col}': {actual}"


def test_monthly_normals_in_weather_headers(helpm, output_file):
    """
    Test that monthly climate normals are correctly calculated and injected
    into the 4th line of the D4 and D7 input files headers.

    See cgq-qgc/pyhelp#127
    """
    # Generate the input files for the first cell
    cellnames = helpm.cellnames[:1]
    cellname = cellnames[0]

    # -------------------------------------------------------------------------
    # Test Precipitation (D4)
    # -------------------------------------------------------------------------
    d4_file = helpm.connect_tables['D4'][cellname]
    d4_col_idx = helpm.connect_tables['precip'][cellname]

    # Calculate the expected monthly normals directly from the pandas Series
    precip = helpm.precip_data.iloc[:, d4_col_idx]

    expected_d4_total = precip.resample("ME").sum()
    expected_d4_normals = (
        expected_d4_total.groupby(expected_d4_total.index.month)
        .mean()
        .sort_index()
        .values
        .round(2)
        )

    with open(d4_file, 'r') as f:
        d4_lines = f.readlines()

    assert len(d4_lines) == (37 * 11) + 4

    # The 4th line (index 3) should contain the 12 normals
    d4_normals_line = d4_lines[3].rstrip('\n')

    # Format F6.2 for 12 values means the line must be exactly
    # 72 characters long
    assert len(d4_normals_line) == 72

    # Parse the strings back into floats
    d4_parsed_normals = [
        float(d4_normals_line[i:i+6]) for i in range(0, 72, 6)]
    assert len(d4_parsed_normals) == 12

    # Assert they match our expected values.
    assert np.max(np.abs(d4_parsed_normals - expected_d4_normals)) < 0.001

    # -------------------------------------------------------------------------
    # Test Air Temperature (D7)
    # -------------------------------------------------------------------------
    d7_file = helpm.connect_tables['D7'][cellname]
    d7_col_idx = helpm.connect_tables['airtemp'][cellname]

    # Calculate the expected monthly normals directly from the pandas Series
    airtemp = helpm.airtemp_data.iloc[:, d7_col_idx]
    expected_d7_normals = (
        airtemp.groupby(airtemp.index.month)
        .mean()
        .round(2)
        .values
        )

    with open(d7_file, 'r') as f:
        d7_lines = f.readlines()

    assert len(d7_lines) == (37 * 11) + 4

    d7_normals_line = d7_lines[3].rstrip('\n')
    assert len(d7_normals_line) == 72

    d7_parsed_normals = [
        float(d7_normals_line[i:i+6]) for i in range(0, 72, 6)]
    assert len(d7_parsed_normals) == 12

    # Assert they match our expected values.
    max_abs_err = np.max(np.abs(d7_parsed_normals - expected_d7_normals))
    assert max_abs_err < 0.001, max_abs_err


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
    expected_values = {
        'precip': 1086.8448950125246,
        'rechg': 145.8440860111447,
        'runoff': 228.41669613037263,
        'evapo': 523.7231069076861,
        'subrun1': 53.17822738799963,
        'subrun2': 132.172476817029
        }
    for i, (name, value) in enumerate(expected_values.items()):
        result = children[i].get_ydata().sum()
        assert abs(result - value) < 0.1, f'{name}: {result} vs {value}'

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

    expected_values = {
        'precip': 1086.8448950125246,
        'rechg': 145.8440860111447,
        'runoff': 228.41669613037263,
        'evapo': 523.7231069076861,
        'subrun1': 53.17822738799963,
        'subrun2': 132.172476817029
        }
    for i, (name, value) in enumerate(expected_values.items()):
        result = children[i * 2].get_ydata().mean()
        assert abs(result - value) < 0.1, f'{name}: {result} vs {value}'

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
    expected_values = {
        'precip': 1086.8448950125246,
        'rechg': 145.8440860111447,
        'runoff': 228.41669613037263,
        'evapo': 523.7231069076861,
        'subrun1': 53.17822738799963,
        'subrun2': 132.172476817029
        }
    for i, (name, value) in enumerate(expected_values.items()):
        height = children[i * 2].get_height()
        assert abs(height - value) < 0.1, f'{name}: {height} vs {value}'
        assert children[i * 2 + 1].get_text() == f'{round(value)}\nmm/an'

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
        'perco': 269.85915013361415,
        'evapo': 524.5229649707427,
        'rechg': 139.2435204331471,
        'runoff': 213.5439791499779,
        'subrun1': 51.0534471241395,
        'subrun2': 122.72894621921611}
    for name, value in expected_results.items():
        result = np.mean(yearly_avg[name])
        assert abs(result - value) < 1, f'{name}: {result} vs {value}'

    # Test calc_cells_yearly_avg with non null year_from and year_to argument.
    yearly_avg = output.calc_cells_yearly_avg(year_from=2003, year_to=2009)
    expected_results = {
        'precip': 1086.8448950125246,
        'perco': 279.6540188589911,
        'evapo': 523.7231069076861,
        'rechg': 145.84408601114473,
        'runoff': 228.41669613037257,
        'subrun1': 53.178227387999634,
        'subrun2': 132.17247681702898}
    for name, value in expected_results.items():
        result = np.mean(yearly_avg[name])
        assert abs(result - value) < 1, f'{name}: {result} vs {value}'

    # Test calc_cells_yearly_avg with year_from == year_to.
    yearly_avg = output.calc_cells_yearly_avg(year_from=2003, year_to=2003)
    expected_results = {
        'precip': 1144.4142919267927,
        'perco': 338.01189818333636,
        'evapo': 474.3865899455791,
        'rechg': 157.384908103213,
        'runoff': 165.76093843767458,
        'subrun1': 60.32943865514778,
        'subrun2': 150.2097857308441}
    for name, value in expected_results.items():
        result = np.mean(yearly_avg[name])
        assert abs(result - value) < 1, f'{name}: {result} vs {value}'


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
        'precip': 1055.859737334313,
        'perco': 267.6719784568864,
        'evapo': 527.1941984146556,
        'rechg': 137.92593783036364,
        'runoff': 213.45404069293568,
        'subrun1': 50.6600028195145,
        'subrun2': 121.90106029782127}
    for key in list(expected_results.keys()):
        result = df[key].sum() / len(df)
        expected_result = expected_results[key]
        assert abs(result - expected_result) < 1, key


if __name__ == '__main__':
    pytest.main(['-x', __file__, '-vv', '-rw'])
