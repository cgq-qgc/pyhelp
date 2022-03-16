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
def helpm(tmp_path):
    manager = HelpManager(tmp_path, year_range=(2000, 2010))

    assert manager.precip_data is None
    assert manager.airtemp_data is None
    assert manager.solrad_data is None
    assert manager.grid is None

    return manager


# =============================================================================
# ---- Tests
# =============================================================================
def test_read_input(helpm):
    """
    Test that the input files are read as expected by the HelpManager.
    """
    # Copy input weather data in the temp folder.
    for var in VARNAMES:
        shutil.copyfile(
            osp.join(DATAFOLDER, f'{var}_input_data_2001-2002.csv'),
            osp.join(helpm.workdir, f'{var}_input_data.csv'))

    helpm.load_weather_input_data()

    assert len(helpm.airtemp_data['data']) == 365 * 2
    assert np.unique(helpm.airtemp_data['years']).tolist() == [2001, 2002]
    assert helpm.airtemp_data['lat'].tolist() == [
        46.1, 45.9, 45.8, 45.9, 46, 46.1]
    assert helpm.airtemp_data['lon'].tolist() == [
        -74, -74.4, -74.1, -74, -74.1, -74.1]

    assert len(helpm.precip_data['data']) == 365 * 2
    assert np.unique(helpm.precip_data['years']).tolist() == [2001, 2002]
    assert helpm.precip_data['lat'].tolist() == [46.1, 45.9, 45.8]
    assert helpm.precip_data['lon'].tolist() == [-74, -74.4, -74.1]

    assert len(helpm.solrad_data['data']) == 365 * 2
    assert np.unique(helpm.solrad_data['years']).tolist() == [2001, 2002]
    assert helpm.solrad_data['lat'].tolist() == [45.47]
    assert helpm.solrad_data['lon'].tolist() == [-73.74]


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_len_nomatch(helpm, testvar):
    """
    Test that the HelpManager throws an error if the len of the weather
    data do not match.
    """
    for var in VARNAMES:
        if var == testvar:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_2001.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
        else:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_2001-2002.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
    with pytest.raises(ValueError):
        helpm.load_weather_input_data()


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_dates_nomatch(helpm, testvar):
    """
    Test that the HelpManager throws an error if the len of the weather
    data do not match.
    """
    for var in VARNAMES:
        if var == testvar:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_2002-2003.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
        else:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_2001-2002.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
    with pytest.raises(ValueError):
        helpm.load_weather_input_data()

    assert (len(helpm.airtemp_data['data']) ==
            len(helpm.precip_data['data']) ==
            len(helpm.solrad_data['data']))


@pytest.mark.parametrize('testvar', VARNAMES)
def test_read_weather_incomplete(helpm, testvar):
    """
    Test that the HelpManager throws an error if the weather are not complete
    for one or more years in the dataset.
    """
    for var in ['precip', 'airtemp', 'solrad']:
        if var == testvar:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_incomplete.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
        else:
            shutil.copyfile(
                osp.join(DATAFOLDER, f'{var}_input_data_2001-2002.csv'),
                osp.join(helpm.workdir, f'{var}_input_data.csv'))
    with pytest.raises(ValueError):
        helpm.load_weather_input_data()


if __name__ == '__main__':
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
