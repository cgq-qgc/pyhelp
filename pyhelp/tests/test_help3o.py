# -*- coding: utf-8 -*-
# =============================================================================
# Copyright © PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the MIT License.
# =============================================================================


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
def rca_folder():
    return osp.join(__rootdir__, 'tests', 'rca_original_testcase_1997')


@pytest.fixture
def rca_params(rca_folder, tmp_path):
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
            osp.join(tmp_path, 'NEW_RCA.OUT'),
            daily_out,
            monthly_out,
            yearly_out,
            summary_out,
            unit_system,
            simu_nyear,
            tfsoil)


# ---- Tests
def test_run_help3o(rca_params):
    """
    Test that the HELP3O extension run and create an output file as expected.
    """
    HELP3O.run_simulation(*rca_params)
    assert osp.exists(rca_params[5])


def test_run_help_singlecell(rca_params):
    """
    Run HELP for a single cell using the RCA test case and assert that the
    results are as expected.
    """
    cellname, results = run_help_singlecell(('rca', rca_params))
    assert not osp.exists(rca_params[5])
    assert cellname == 'rca'

    # Precipitations.
    precip = np.sum(results['precip'] * 0.0393701, axis=1)
    for i, expected_result in enumerate([48.53, 58.32, 56.71]):
        assert abs(precip[i] - expected_result) < 0.1, 'precip year %i' % i

    # Runoff.
    runoff = np.sum(results['runoff'] * 0.0393701, axis=1)
    for i, expected_result in enumerate([0.596, 4.130, 2.397]):
        assert abs(runoff[i] - expected_result) < 0.1, 'runoff year %i' % i

    # Recharge.
    rechg = np.sum(results['rechg'] * 0.0393701, axis=1)
    for i, expected_result in enumerate([0.000004, 0.000006, 0.000007]):
        assert abs(rechg[i] - expected_result) < 0.1, 'rechg year %i' % i

    # Evapotranspiration.
    evapo = np.sum(results['evapo'] * 0.0393701, axis=1)
    for i, expected_result in enumerate([32.178, 33.420, 36.540]):
        assert abs(evapo[i] - expected_result) < 0.1, 'evapo year %i' % i

    # Superficial subsurface runoff.
    # (DRAINAGE COLLECTED FROM LAYER  2)
    subrun1 = np.sum(results['subrun1'] * 0.0393701, axis=1)
    for i, expected_result in enumerate([15.4971, 22.8698, 19.1360]):
        assert abs(subrun1[i] - expected_result) < 0.1, 'subrun1 year %i' % i

    # Deep subsurface runoff.
    # (DRAINAGE COLLECTED FROM LAYER  7 + DRAINAGE COLLECTED FROM LAYER  9)
    subrun2 = np.sum(results['subrun2'] * 0.0393701, axis=1)
    expected_results = (0.0543 + 0.0833, 0.1481 + 0.1658, 0.1835 + 0.1935)
    for i, expected_result in enumerate(expected_results):
        assert abs(subrun2[i] - expected_result) < 0.1, 'subrun2 year %i' % i


if __name__ == '__main__':
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
