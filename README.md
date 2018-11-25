![PyHELP](./images/pyhelp_banner.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Latest release](https://img.shields.io/github/release/cgq-qgc/pyhelp.svg)](https://github.com/cgq-qgc/pyhelp/releases)
[![Documentation Status](https://readthedocs.org/projects/pyhelp/badge/?version=latest)](http://pyhelp.readthedocs.io)
[![Build status](https://ci.appveyor.com/api/projects/status/ns6s8x0hkd31ffb3/branch/master?svg=true)](https://ci.appveyor.com/project/jnsebgosselin/pyhelp-rd625/branch/master)
[![codecov](https://codecov.io/gh/cgq-qgc/pyhelp/branch/master/graph/badge.svg)](https://codecov.io/gh/cgq-qgc/pyhelp)


PyHELP is a [Python](https://www.python.org/) library for the assessment of groundwater recharge at the regional scale with the [Hydrologic Evaluation of Landfill Performance](https://www.epa.gov/land-research/hydrologic-evaluation-landfill-performance-help-model) (HELP) model.

Please consult the [documentation](http://pyhelp.readthedocs.io) for more details or contact us at [jean-sebastien.gosselin@ete.inrs.ca](mailto:jean-sebastien.gosselin@ete.inrs.ca).

# Install
[![Anaconda-Server Badge](https://anaconda.org/cgq-qgc/pyhelp/badges/installer/conda.svg)](https://anaconda.org/cgq-qgc/pyhelp)
[![Anaconda-Server Badge](https://anaconda.org/cgq-qgc/pyhelp/badges/installer/pypi.svg)](https://pypi.org/project/pyhelp/)

If you are new to Python or the Scientific Python ecosystem, we strongly recommend you to to install the released version of PyHELP with [Conda](https://conda.io/docs/index.html). To do so, you will need first to download and install the [Anaconda distribution](https://www.anaconda.com/distribution/) on your computer. Anaconda comes with the most important Python scientific libraries (i.e. Numpy, Pandas, Matplotlib, IPython, etc), including all PyHELP dependencies, in a single, easy to use environment. Then, the easiest way to install PyHELP and all its depencies is by executing the following `conda` command in a terminal:

`conda install -c cgq-qgc pyhelp scipy geopandas xlrd netcdf4 h5py pytables matplotlib`

We also support [pip](https://pypi.org/project/pip/), but be aware that pip installations are for advanced users. PyHELP depends on several low-level libraries for geospatial analysis, and this may cause dependency conflicts if you are not careful.



