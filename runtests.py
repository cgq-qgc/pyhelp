# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

"""
File for running tests programmatically.
"""
import sys
import pytest


def main():
    """
    Run pytest tests.
    """
    errno = pytest.main(['-x', 'pyhelp',  '-v', '-rw', '--durations=10',
                         '--cov=pyhelp'])
    if errno != 0:
        raise SystemExit(errno)


if __name__ == '__main__':
    main()
