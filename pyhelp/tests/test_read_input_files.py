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
import shutil

# ---- Third party imports
import numpy as np
import pytest

# ---- Local library imports
from pyhelp import __rootdir__
from pyhelp.managers import HelpManager

DATAFOLDER = osp.join(__rootdir__, 'tests', 'data')
VARNAMES = ['precip', 'airtemp', 'solrad']


# =============================================================================
# ---- Fixtures
# =============================================================================
@pytest.fixture
def workdir(tmp_path):
    return tmp_path


# =============================================================================
# ---- Tests
# =============================================================================
def test_read_input(workdir):
    """
    Test that the input files are read as expected by the HelpManager.
    """
    helpm = HelpManager(
        workdir,
        path_to_grid=osp.join(
            DATAFOLDER, 'input_grid.csv'),
        path_to_precip=osp.join(
            DATAFOLDER, 'precip_input_data_2001-2002.csv'),
        path_to_airtemp=osp.join(
            DATAFOLDER, 'airtemp_input_data_2001-2002.csv'),
        path_to_solrad=osp.join(
            DATAFOLDER, 'solrad_input_data_2001-2002.csv')
        )

    assert len(helpm.airtemp_data) == 365 * 2
    assert list(helpm.airtemp_data.index.year.unique()) == [2001, 2002]
    assert list(helpm.airtemp_data.columns.get_level_values('lat_dd')) == [
        46.1, 45.9, 45.8, 45.9, 46, 46.1]
    assert list(helpm.airtemp_data.columns.get_level_values('lon_dd')) == [
        -74, -74.4, -74.1, -74, -74.1, -74.1]

    assert len(helpm.precip_data) == 365 * 2
    assert list(helpm.precip_data.index.year.unique()) == [2001, 2002]
    assert list(helpm.precip_data.columns.get_level_values('lat_dd')) == [
        46.1, 45.9, 45.8]
    assert list(helpm.precip_data.columns.get_level_values('lon_dd')) == [
        -74, -74.4, -74.1]

    assert len(helpm.solrad_data) == 365 * 2
    assert list(helpm.solrad_data.index.year.unique()) == [2001, 2002]
    assert list(helpm.solrad_data.columns.get_level_values('lat_dd')) == [
        45.47]
    assert list(helpm.solrad_data.columns.get_level_values('lon_dd')) == [
        -73.74]


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_len_nomatch(testvar, workdir):
    """
    Test that the HelpManager throws an error if the len of the weather
    data do not match.
    """
    kwargs = {'workdir': workdir,
              'path_to_grid': osp.join(DATAFOLDER, 'input_grid.csv')}
    for var in VARNAMES:
        if var == testvar:
            kwargs[f'path_to_{var}'] = osp.join(
                DATAFOLDER, f'{var}_input_data_2001.csv')
        else:
            kwargs[f'path_to_{var}'] = osp.join(
                DATAFOLDER, f'{var}_input_data_2001-2002.csv')

    with pytest.raises(ValueError):
        HelpManager(**kwargs)


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_dates_nomatch(testvar, workdir):
    """
    Test that the HelpManager throws an error if the len of the weather
    data do not match.
    """
    kwargs = {'workdir': workdir,
              'path_to_grid': osp.join(DATAFOLDER, 'input_grid.csv')}
    for var in VARNAMES:
        if var == testvar:
            kwargs[f'path_to_{var}'] = osp.join(
                DATAFOLDER, f'{var}_input_data_2002-2003.csv')
        else:
            kwargs[f'path_to_{var}'] = osp.join(
                DATAFOLDER, f'{var}_input_data_2001-2002.csv')

    with pytest.raises(ValueError):
        HelpManager(**kwargs)


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_incomplete(testvar, workdir):
    """
    Test that the HelpManager throws an error if the weather are not complete
    for one or more years in the dataset.
    """
    kwargs = {'workdir': workdir,
              'path_to_grid': osp.join(DATAFOLDER, 'input_grid.csv')}
    for var in ['precip', 'airtemp', 'solrad']:
        if var == testvar:
            filename = f'{var}_input_data_incomplete.csv'
        else:
            filename = f'{var}_input_data_2001-2002.csv'
        kwargs[f'path_to_{var}'] = osp.join(DATAFOLDER, filename)

    with pytest.raises(ValueError):
        HelpManager(**kwargs)


if __name__ == '__main__':
    pytest.main(['-x', __file__, '-v', '-rw'])
