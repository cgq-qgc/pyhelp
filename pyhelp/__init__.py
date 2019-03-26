# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

import os

version_info = (0, 1, 2)
__version__ = '.'.join(map(str, version_info))
__appname__ = 'PyHELP'
__namever__ = __appname__ + " " + __version__
__date__ = '23/11/2018'
__project_url__ = "https://github.com/cgq-qgc/pyhelp"
__releases_url__ = __project_url__ + "/releases"
__releases_api__ = "https://api.github.com/repos/cgq-qgc/pyhelp/releases"

__rootdir__ = os.path.dirname(os.path.realpath(__file__))

try:
    from pyhelp.managers import HelpManager
except ImportError as e:
    # We need to do this to avoid an error when building the
    # help extension with setup.py
    print('ImportError:', e)
