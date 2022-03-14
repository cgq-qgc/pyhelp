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
import numpy
from pyhelp import __version__, __project_url__

URL_HELP = ("https://www.epa.gov/land-research/"
            "hydrologic-evaluation-landfill-performance-help-model")
URL_CROTEAU = "https://www.tandfonline.com/doi/abs/10.4296/cwrj3504451"
LONG_DESCRIPTION = (
    "PyHELP is an object oriented Python library providing a set of tools "
    "to estimate spatially distributed groundwater recharge and other "
    "hydrological components (runoff and evapotranspiration) using "
    "the HELP ([Hydrologic Evaluation of Landfill Performance]({})) model."
    "PyHELP integrates weather data (from grids or stations), "
    "land conditions defined by a series of GIS maps as well as soil "
    "and geological material properties into HELP input files. "
    "PyHELP also processes HELP simulation results and outputs them as "
    "maps and graphs, including comparisons of simulation results with "
    "stream hydrographs. PyHELP thus accompanies users through the entire "
    "workflow from input file assembly to model calibration and to the "
    "documentation of results. This workflow is based on the method "
    "originally developed by [Croteau et al. (2011)]({}) to assess spatially "
    "distributed groundwater recharge at the regional scale."
    ).format(URL_HELP, URL_CROTEAU)

HELPEXT = Extension(name='pyhelp.HELP3O',
                    sources=['pyhelp/HELP3O.FOR'],
                    include_dirs=[numpy.get_include()],
                    extra_link_args=["-static",
                                     "-static-libgfortran",
                                     "-static-libgcc"]
                    )

INSTALL_REQUIRES = [
    'numpy',
    'scipy'
    'pandas'
    'h5py>=3'
    'matplotlib']

setup(name='pyhelp',
      version=__version__,
      description=('A Python library for the assessment of spatially '
                   'distributed groundwater recharge and hydrological '
                   'components with HELP.'),
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      license='MIT',
      author='PyHELP Project Contributors',
      author_email='jean-sebastien.gosselin@ete.inrs.ca',
      url=__project_url__,
      ext_modules=[HELPEXT],
      packages=setuptools.find_packages(),
      package_data={'pyhelp': ['*.pyd']},
      include_package_data=True,
      install_requires=INSTALL_REQUIRES,
      classifiers=["License :: OSI Approved :: MIT License",
                   "Operating System :: Microsoft :: Windows",
                   "Programming Language :: Python :: 3.7",
                   "Programming Language :: Python :: 3.8",
                   "Intended Audience :: Science/Research",
                   "Intended Audience :: Education",
                   "Topic :: Scientific/Engineering"],
      )
