Data Input
=================================

This section presents how to prepare the input data required to run the `HELP`_
infiltration model with the :class:`~pyhelp.HelpManager` class of the PyHelp
module.

Weather data needed by PyHelp include precipitation, average air temperature,
and global solar radiation at a daily time step.
Geomatics data needed by PyHelp include a uniform grid that must include, for
each cell of the grid, the location coordinates, the surface conditions used
to compute evapotranspiration and runoff, and finally the soil and design data
needed to compute infiltration, subrunoff, and groundwater recharge.

Format of the weather input data
---------------------------------

Daily precipitation data, average air temperature and global solar irradiance
are needed by PyHelp to calculate the water balance for each cell of the grid
with the HELP model.
Data from any sources can be used, as long as the data are correctly formatted
and saved in the working directory as coma-separated text files with the
following names:

- 1. Precipitation: :file:`precip_input_data.csv`
- 2. Average air temperature: :file:`airtemp_input_data.csv`
- 3. Global solar radiation: :file:`solrad_input_data.csv`

An example of correctly formatted input weather data file is shown in
:numref:`weather_datafile_example`.

The first column of the data is used to store the dates in a `dd/mm/yyyy`
format.
Additional columns are used to store the data of the daily meteorological
time series for one or more geographic locations.
At least one daily time series is required.
When more than one time series is available, PyHelp takes the time series
closest to each cell of the grid to do the calculations.
The latitude and longitude coordinates, in decimal degree, must be provided for
each time series in the header of the file.
All other information in the header is optional and does not need to be
included in the input file.

The data must be saved as text using the coma (,) as separator and the dot (.)
as decimal symbol.
Saving the data using the utf8 encoding is strongly recommended.
The input weather data must be saved in the International System of Units (SI).
More specifically, daily precipitation must be saved in mm, average air
temperature in °C, and global solar irradiance in MJ/m².

By default, the weather data are automatically loaded when instantiating the
:class:`~pyhelp.HelpManager` class or when setting the working directory
with :meth:`~pyhelp.HelpManager.set_workdir`. It is also possible to
explicitly force a reload of the weather data from the input data files
with the :meth:`~pyhelp.HelpManager.load_weather_input_data` method of the
:class:`~pyhelp.HelpManager` class.

Finally, tools are available to generate automatically the weather input data
files from various sources of data. This is covered in more details in the
section `Tools to generate input data files`_.

.. _weather_datafile_example:
.. figure:: img/weather_input_data.*
    :align: center
    :width: 85%
    :alt: weather_input_datafile_example.png
    :figclass: align-center

    Example of a correctly formatted input weather data file.
    

Format of the grid input data
---------------------------------

.. _sec_utils_data:

Tools to generate input data files
-----------------------------------

Tools are available in PyHelp to generate automatically the weather input data
files from various sources of data.

The :meth:`~pyhelp.HelpManager.generate_weather_inputs_from_MDELCC_grid` method
of the :class:`~pyhelp.HelpManager` class can be used to generate automatically
the precipitation and average air temperature input data files using data from
the MDDELCC spatially distributed daily meteo grid.


Similarly, the :meth:`~pyhelp.HelpManager.generate_weather_inputs_from_CWEEDS`
method of the :class:`~pyhelp.HelpManager` class can be used to generate
automatically the global solar irradiance input data file from a set of
CWEEDS files.

Please consult the documentation of each method for more details.


Example
---------------------------------

Import and instantiate the :class:`~pyhelp.HelpManager` class ::

    >>> from pyhelp import HelpManager
    >>> helpm = HelpManager("C:/path_to_pyhelp_project")

Generate precipitation and air temperature input files from the MDDELCC
weather grid ::

    >>> helpm.generate_weather_inputs_from_MDELCC_grid("C:/path_to_mddelcc_grid") 

Generate global solar irradiance input file from CWEEDS files ::

     >>> cweed2_paths = "C:/path_to_cweed2_file"
     >>> cweed3_paths = "C:/path_to_cweed3_file"
     >>> helpm.generate_weather_inputs_from_CWEEDS(cweed2_paths, cweed3_paths) 
     
     
.. _HELP: https://www.epa.gov/land-research/hydrologic-evaluation-landfill-performance-help-model