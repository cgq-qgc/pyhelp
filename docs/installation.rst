Installation
=================================

Pip Wheels and Conda packages are both available for Python 3.6 and 3.7 on the
Windows 64bits plateform.
If you need to use PyHELP with a version of Python older than 3.6 or if you
are working on Linux or macOS, you will have to build and install PyHELP from
source.

Install with Conda
---------------------------------

.. warning:: This installation method is currently only supported for the
             Windows 64bits plateform.

The easiest method to install a released version of PyHELP on Windows is
with Conda. To do so, you will need first to download and install the
`Anaconda distribution`_ on your computer.
Anaconda comes with the most important Python scientific libraries
(i.e. Numpy, Pandas, Matplotlib, IPython, etc), including all PyHELP
dependencies, in a single, easy to use environment. It also includes the
`Anaconda Navigator`_, which is a graphical user interface to the Conda
package and environment manager.

First, to avoid installation problems and dependency conflicts, we advise you
to install PyHELP in a fresh new conda environment.
This can be done with the Anaconda Navigator or with Conda by executing the
following command in a terminal:

.. code-block:: bash

   conda create -n my_new_env_name python=3.7.*

If you want to use Python 3.6 instead, simply replace the ``python=3.7.*``
argument by ``python=3.6.*``.

Then, PyHELP can be installed, along with all its dependencies, by executing
the following command in a terminal:

.. code-block:: bash

   conda install -c cgq-qgc pyhelp
   
When a new released version of PyHELP is made available on the
`Anaconda cgq-qgc channel`_, PyHELP can be updated by executing the following
command in a terminal:

.. code-block:: bash

   conda update -c cgq-qgc pyhelp

   
If you need more guidance on how to install packages or manage Conda
environments with Conda or the Anaconda Navigator, please consult the 
`Getting started with conda`_ or `Getting started with Navigator`_ guide.
            
Install with Pip
---------------------------------

.. warning:: This installation method is currently only supported for the
             Windows 64bits plateform.

It is also possible to install PyHELP with `pip`_, but be aware that pip
installations are for advanced users.
PyHELP depends on several low-level libraries for geospatial analysis, and
this may cause dependency conflicts if you are not careful.

First, you will need to download and install `Python 3.6 or 3.7`_ on your
computer.
Then you will need to install all the dependencies that are listed in
the section :ref:`sec_requirements` below.
Unless you really know what you are doing, we strongly recommand against
installing these dependencies directly from the `The Python Package Index (PyPI)`_
with pip, because you will most likely run into installation problems and
dependency conflicts.
The easiest and safest way to install PyHELP's depencies on Windows is to
download Wheels from Christopher Gohlke's
`Unofficial Windows Binaries for Python Extension Packages`_ and
`install them with pip`_.
Be carefull to install the packages that were built for Windows 64bits and
the version of Python that you downloaded and installed on your computer.

Then, you can install PyHELP with pip by executing the following command
in a terminal:

.. code-block:: bash
   
   python -m pip install pyhelp
   
.. _sec_install_from_source:

Install PyHELP from source
---------------------------------

If you need to use PyHELP with a version of Python older than 3.6 or
if you are working on Linux or macOS, you will have to build and install
PyHELP from source.
Below is a step-by-step guide that describe how to achieve this.

#. Install PyHELP's `requirements`_.

   The recommended and easiest way to do this is with `Anaconda`_, a free
   and open source distribution of Python. If you decide to do so,
   PyHELP's `requirements`_ can be installed one by one with the Anaconda
   Navigator or with Conda by executing the following command in a terminal:
   
   .. code-block:: bash

      conda install scipy geopandas xlrd netcdf4 h5py pytables matplotlib

#. Install `Git`_, a powerful source control management tool, or install one
   of the numerous `GUI client`_ that exists for it .

#. Clone the PyHELP source code repository by running this command :

   .. code-block:: bash

      git clone https://github.com/cgq-qgc/pyhelp.git <path-to-target-dir>
    
   or do it with your GUI client if you are using one.

#. Build and install PyHELP by executing the following commands
   in a terminal from inside your cloned directory:
   
   .. code-block:: bash

      python setup.py build_ext
      python setup.py install
      
   To do the above, you will need to have a Fortran and C++ compiler installed
   on your computer. If you are using Anaconda, you can achieve that simply by
   installing the conda package named `m2w64-toolchain`.
   If you do not use Anaconda and are working on Linux or macOS, you can
   install the free and open source `GNU Compiler Collection (GCC)`_ and
   the `GNU Fortran compiler (gfortran)`_ with the package manager of your
   operating system.
   If you are on Windows, you can download and install `mingw-w64`_, which is
   a complete runtime environment for gcc to support binaries native to
   Windows 64-bit and 32-bit operating systems.

#. Open Python and start using PyHELP.

#. To keep your PyHELP repository up-to-date, run ``git pull`` inside the
   cloned directory or do it with your GUI client.
   You then need to re-build and re-install PyHELP, so that the pulled
   changes are applied to the PyHELP installation used by your Python
   installation. 

.. _sec_requirements:

Requirements
-----------------------------------------------

- `Python <https://www.python.org/>`_ == 3.6
- `Matplotlib <https://matplotlib.org/>`_
- `Numpy <https://www.numpy.org/>`_
- `Pandas <https://pandas.pydata.org/>`_
- `Scipy <https://www.scipy.org/>`_
- `xlrd <https://github.com/python-excel/xlrd/>`_
- `netCDF4 <http://unidata.github.io/netcdf4-python/>`_
- `H5py <https://www.h5py.org/>`_
- `GeoPandas <http://geopandas.org/>`_
- `PyTables <https://www.pytables.org/>`_

.. _Anaconda cgq-qgc channel: https://anaconda.org/cgq-qgc/pyhelp
.. _Anaconda: https://www.anaconda.com/download/
.. _Anaconda distribution: <https://www.anaconda.com/download/
.. _Anaconda Navigator: https://docs.anaconda.com/anaconda/navigator/
.. _mingw-w64: https://sourceforge.net/projects/mingw-w64/
.. _Getting started with conda: https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html
.. _Getting started with Navigator: https://docs.anaconda.com/anaconda/navigator/getting-started/
.. _Git: https://git-scm.com/downloads
.. _GitHub repository: https://github.com/jnsebgosselin/pyhelp
.. _GNU Fortran compiler (gfortran): https://gcc.gnu.org/wiki/GFortran
.. _GNU Compiler Collection (GCC): https://gcc.gnu.org/
.. _GUI client: https://git-scm.com/download/gui/windows
.. _install them with pip:  https://pip.pypa.io/en/stable/user_guide/#installing-from-wheels
.. _pip: https://pypi.org/project/pip/
.. _Python 3.6 or 3.7: https://www.python.org/downloads/
.. _The Python Package Index (PyPI): https://pypi.org/
.. _Unofficial Windows Binaries for Python Extension Packages: https://www.lfd.uci.edu/~gohlke/pythonlibs/
