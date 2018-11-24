# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.

"""Installation script """

import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from distutils.command.sdist import sdist
import numpy
from pyhelp import __version__, __project_url__


LONG_DESCRIPTION = ("PyHELP is an object oriented Python library that "
                    "provides a set of tools to evaluate spatially "
                    "distributed groundwater recharge at the regional "
                    "scale using the HELP model. This work is based on the "
                    "method that was originally developed by "
                    "Croteau et al. (2011) to assess spatially distributed "
                    "groundwater recharge in the Chateauguay River "
                    "Watershed, Quebec, Canada.")

HELPEXT = Extension(name='pyhelp.HELP3O',
                    sources=['pyhelp/HELP3O.FOR'],
                    include_dirs=[numpy.get_include()],
                    extra_link_args=["-static",
                                     "-static-libgfortran",
                                     "-static-libgcc"]
                    )

setup(name='pyhelp',
      version=__version__,
      description=('A Python library for groundwater recharge assessment in '
                   'regional studies with the HELP  model'),
      long_description=LONG_DESCRIPTION,
      license='MIT',
      author='PyHELP Project Contributors',
      author_email='jean-sebastien.gosselin@ete.inrs.ca',
      url=__project_url__,
      ext_modules=[HELPEXT],
      packages=setuptools.find_packages(),
      package_data={'pyhelp': ['*.pyd']},
      include_package_data=True,
      cmdclass={"sdist": sdist},
      classifiers=["Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: Microsoft :: Windows",
                   "Programming Language :: Python :: 3.6",
                   "Intended Audience :: Science/Research",
                   "Intended Audience :: Education",
                   "Topic :: Scientific/Engineering"],
      )
