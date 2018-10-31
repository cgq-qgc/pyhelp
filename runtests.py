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
    sys.path.insert(
        0, "C:\\mingw-w64\\x86_64-7.2.0-posix-seh-rt_v5-rev1\\mingw64\\bin")

    errno = pytest.main(['-x', 'pyhelp',  '-v', '-rw', '--durations=10',
                         '--cov=pyhelp'])
    if errno != 0:
        raise SystemExit(errno)


if __name__ == '__main__':
    main()
