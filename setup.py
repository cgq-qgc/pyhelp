# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

from numpy.distutils.core import setup, Extension
import numpy


setup(ext_modules=[Extension(name='pyhelp.HELP3O',
                             sources=['pyhelp/HELP3O.FOR'],
                             extra_link_args=['-static'])],
      include_dirs=[numpy.get_include()])
