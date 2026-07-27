![PyHELP](https://raw.githubusercontent.com/cgq-qgc/pyhelp/master/images/pyhelp_banner_v2.png)

[![license](https://img.shields.io/pypi/l/pyhelp.svg)](./LICENSE)
[![Latest release](https://img.shields.io/github/release/cgq-qgc/pyhelp.svg)](https://github.com/cgq-qgc/pyhelp/releases)
[![Documentation Status](https://readthedocs.org/projects/pyhelp/badge/?version=latest)](http://pyhelp.readthedocs.io)
[![tests](https://github.com/cgq-qgc/pyhelp/actions/workflows/python-test.yml/badge.svg?branch=master)](https://github.com/cgq-qgc/pyhelp/actions/workflows/python-test.yml)
[![codecov](https://codecov.io/gh/cgq-qgc/pyhelp/branch/master/graph/badge.svg)](https://codecov.io/gh/cgq-qgc/pyhelp)

PyHELP is an object oriented Python library providing a set of tools to
estimate spatially distributed groundwater recharge and other hydrological
components (runoff and evapotranspiration) using the HELP
([Hydrologic Evaluation of Landfill Performance](https://www.epa.gov/land-research/hydrologic-evaluation-landfill-performance-help-model))
model.
PyHELP integrates weather data (from grids or stations), land conditions
defined by a series of GIS maps as well as soil and geological material
properties into HELP input files.
PyHELP also processes HELP simulation results and outputs them as
maps and graphs, including comparisons of simulation results with
stream hydrographs.
PyHELP thus accompanies users through the entire workflow from input file
assembly to model calibration and to the documentation of results.
This workflow is based on the method originally developed by
[Croteau et al. (2011)](https://www.tandfonline.com/doi/abs/10.4296/cwrj3504451)
to assess spatially distributed groundwater recharge at the regional scale.

Please consult the [documentation](http://pyhelp.readthedocs.io) for more
details.

## Installation
[![pypi version](https://img.shields.io/pypi/v/pyhelp.svg)](https://pypi.org/project/pyhelp/)

Pip Wheels are available for Python 3.11 to 3.13 on the Windows 64bits plateform. If you need to use PyHELP with a version of Python older than 3.11 or more recent than 3.13, or if you are working on Linux or macOS, you will have to build and install PyHELP from source.

To install PyHELP, along with all its dependencies, simply run the following command in a terminal:

`python -m pip install pyhelp`

### Building from source (Linux / macOS / development)

On Linux and macOS, and when developing PyHELP, you need to build the compiled
HELP3O extension from source. This requires a C and a **Fortran** compiler:

- Debian/Ubuntu: `sudo apt install gfortran`
- macOS: `brew install gcc`

PyHELP supports Python 3.11 to 3.13 (Python 3.14 and newer are not yet
supported), so make sure you work inside a virtual environment using one of
these versions:

```
python3.12 -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
```

Then install the build toolchain and PyHELP in editable mode:

```
pip install -r requirements-dev.txt
pip install --no-build-isolation --editable .
```

`--no-build-isolation` makes meson-python install an import hook that
**recompiles the HELP3O extension automatically** whenever the package is
imported, so your changes to both the Python and the Fortran sources are picked
up without reinstalling. To see the rebuild output, add
`--config-settings=editable-verbose=true` to the editable install command.

Run the test suite with:

```
python runtests.py
```
