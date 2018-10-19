Installation
=================================

PyHelp is in an early stage of development and no version has been published
yet.
Therefore, the only way to use PyHelp currently is to clone it from our
`GitHub repository`_ and run it from the source files after installing the
`requirements`_ that are listed below.

.. note::  PyHelp is currently available only on the Windows plateform and on
           Python 3.6 64-bits.

In summary:

#. Install the PyHelp `requirements`_.

   The recommended and easiest way to do this is with `Anaconda`_, a free
   and open source distribution of Python.

#. Install `Git`_, a powerful source control management tool, or install one
   of the numerous`GUI client`_ that exists for it .

#. Clone the PyHelp source code repository by running this command :

   .. code-block:: bash

      git clone https://github.com/jnsebgosselin/pyhelp.git <path-to-target-dir>
    
   or do it from your GUI client if you are using one.

#. Add your cloned PyHelp directory to your Windoes `PYTHONPATH`.

#. To keep your repository up-to-date, run ``git pull`` inside the cloned
   directory or do it directly from your GUI client.
   
#. Open Python and starts using the latest version of PyHelp.


Requirements
=================================

- `Python <https://www.python.org/>`_ == 3.6
- `Matplotlib <https://matplotlib.org/>`_
- `Numpy <https://www.numpy.org/>`_
- `Pandas <https://pandas.pydata.org/>`_
- `Scipy <https://www.scipy.org/>`_
- `xlrd <https://github.com/python-excel/xlrd/>`_
- `netCDF4 <http://unidata.github.io/netcdf4-python/>`_
- `H5py  <https://www.h5py.org/>`_


.. _Anaconda: https://www.anaconda.com/download/
.. _Git: https://git-scm.com/downloads
.. _GitHub repository: https://github.com/jnsebgosselin/pyhelp
.. _GUI client: https://git-scm.com/download/gui/windows