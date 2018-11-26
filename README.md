![PyHELP](./images/pyhelp_banner.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Latest release](https://img.shields.io/github/release/cgq-qgc/pyhelp.svg)](https://github.com/cgq-qgc/pyhelp/releases)
[![Documentation Status](https://readthedocs.org/projects/pyhelp/badge/?version=latest)](http://pyhelp.readthedocs.io)
[![Build status](https://ci.appveyor.com/api/projects/status/ns6s8x0hkd31ffb3/branch/master?svg=true)](https://ci.appveyor.com/project/jnsebgosselin/pyhelp-rd625/branch/master)
[![codecov](https://codecov.io/gh/cgq-qgc/pyhelp/branch/master/graph/badge.svg)](https://codecov.io/gh/cgq-qgc/pyhelp)


PyHELP is a [Python](https://www.python.org/) library for the assessment of groundwater recharge at the regional scale with the [Hydrologic Evaluation of Landfill Performance](https://www.epa.gov/land-research/hydrologic-evaluation-landfill-performance-help-model) (HELP) model.

Please consult the [documentation](http://pyhelp.readthedocs.io) for more details or contact us at [jean-sebastien.gosselin@ete.inrs.ca](mailto:jean-sebastien.gosselin@ete.inrs.ca).

## Install from built distributions

Pip wheels and conda packages are both available for Python 3.6 on the Windows 64bits plateform. If you need to use PyHELP with Python 3.7 or are working on Linux or MacOS, you will have to build and install PyHELP from source.

[![Anaconda-Server Badge](https://anaconda.org/cgq-qgc/pyhelp/badges/installer/conda.svg)](https://anaconda.org/cgq-qgc/pyhelp)

The easiest method to install the released version of PyHELP on Windows is with [Conda](https://conda.io/docs/index.html). To do so, you will need first to download and install the [Anaconda distribution](https://www.anaconda.com/distribution/) on your computer. Anaconda comes with the most important Python scientific libraries (i.e. Numpy, Pandas, Matplotlib, IPython, etc), including all PyHELP dependencies, in a single, easy to use environment.

PyHELP can then be installed, along with all its dependencies, by executing the following `conda` command in a terminal:

`conda install -c cgq-qgc pyhelp scipy geopandas xlrd netcdf4 h5py pytables matplotlib`

[![Anaconda-Server Badge](https://anaconda.org/cgq-qgc/pyhelp/badges/installer/pypi.svg)](https://pypi.org/project/pyhelp/)

It is also possible to install PyHELP with [pip](https://pypi.org/project/pip/), but be aware that pip installations are for advanced users. PyHELP depends on several low-level libraries for geospatial analysis, and this may cause dependency conflicts if you are not careful.

First, you will need to download and install [Python 3.6](https://www.python.org/downloads/release/python-367/) on your computer. Then you will need to download and install with `pip`, in that specific order, `numpy`, `matplotlib`, `scipy`, `pandas`, `shapely`, `fiona`, `pyproj`, `geopandas`, `h5py`, `pytables`, `xlrd`, `netcdf4`from Christopher Gohlke's [Unofficial Windows Binaries for Python Extension Packages](https://www.lfd.uci.edu/~gohlke/pythonlibs/). Be carefull to install the packages that were built for Python 3.6 on Windows 64bits.

Finally, you can install PyHELP by executing the following Python command in a terminal:

`python -m pip install `







