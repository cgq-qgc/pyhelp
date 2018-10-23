.. _sec_data_input:

Data Input
=================================
The PyHelp module requires multiple weather and geomatics data to perform
spatially distributed water balance calculations at the regional scale with
the `HELP`_ infiltration model.

More specifically, weather data include precipitation, average air temperature,
and global solar irradiance at a daily time step.
Geomatics data consists in a uniformly spaced grid where the location
coordinates must be defined for each cell of the grid, as well as the surface
conditions required to compute evapotranspiration and runoff, and the soil and
design data needed to compute infiltration, subrunoff, and groundwater
recharge.

The input data must be formatted and saved as coma-separated text files
in the working directory of PyHelp with the following names:

- 1. :file:`precip_input_data.csv` for daily precipitation in mm
- 2. :file:`airtemp_input_data.csv` for average air temperature in °C
- 3. :file:`solrad_input_data.csv` for global solar irradiance in MJ/m²
- 4. :file:`input_grid.csv` for the grid and geomatics data.

This section presents how to prepare these input data files, either manually
or automatically with the tools that are available in PyHelp.

Format of the weather input data
---------------------------------

Daily precipitation data, average air temperature and global solar irradiance
are needed by PyHelp to calculate the water balance for each cell of the grid
with the HELP model.
Data from any sources can be used, as long as the data are correctly formatted
and saved in the working directory as coma-separated text files named as
:file:`precip_input_data.csv` for daily precipitation, 
:file:`airtemp_input_data.csv` for average air temperature, and
:file:`solrad_input_data.csv` for global solar irradiance.

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

Format of the grid input data
---------------------------------

The geomatics data required to run HELP calculations for each cell of
the grid must be formatted and saved in the working directory as a
coma-separated text file named :file:`grid_input.csv`.
An example of correctly formatted input grid data file is shown in
:numref:`grid_datafile_example`.
:numref:`table_grid_field_desc` presents the required information that must be
provided for each cell of the grid in the input grid data file.
Note that the name of the fields must be respected faithfully in the data
header of the file, as well as the units of the data.

The field `run` is used to identify cells that must be run with HELP. All 
cells with a `run` value of 0 are skipped when executing
:meth:`pyhelp.HelpManager.calc_help_cells`. The method 
:meth:`pyhelp.HelpManager.get_run_cellnames` can be used to get a list of cell
ids for wich the `run` value is 1.
The field `context` is used to identify cells that are consisered to be
located in surface water bodies. It is also used to identify cells that are
located near a stream, in urban areas, and cells for which data are
incomplete.

In addition, any field can be added to the grid for cell selection purpose.
For example, a field could be added to faciliate the selection of cells
by watershed or region. These selection fields are particularly useful for
the calibration of the model.

.. _grid_datafile_example:
.. figure:: img/grid_input_data.*
    :align: center
    :width: 85%
    :alt: grid_input_datafile_example.png
    :figclass: align-center

    Example of a correctly formatted grid input data file.

.. _table_grid_field_desc:
.. table:: Field description of the :file:`grid_input.csv`
   :widths: auto

   +--------------+-----------------+----------------------------------------+
   | Field Name   | Units           | Description                            |
   +==============+=================+========================================+
   | cid          |                 | Unique cell ID                         |
   +--------------+-----------------+----------------------------------------+
   | lat_dd       | Decimal degrees | Latitude of the cell centroid          |
   +--------------+-----------------+----------------------------------------+
   | lon_dd       | Decimal degrees | Longitude of the cell centroid         |
   +--------------+-----------------+----------------------------------------+
   | wind         | km/h            | Average annual wind speed              |
   +--------------+-----------------+----------------------------------------+
   | hum1         | %               | Average quaterly relative humidity     |
   |              |                 | (jan to mar)                           |
   +--------------+-----------------+----------------------------------------+
   | hum2         | %               | Average quaterly relative humidity ()  |
   |              |                 | (apr to jun)                           |
   +--------------+-----------------+----------------------------------------+
   | hum3         | %               | Average quaterly relative humidity     |
   |              |                 | (jul to sep)                           |
   +--------------+-----------------+----------------------------------------+
   | hum4         | %               | Average quaterly relative humidity     |
   |              |                 | (oct to dec)                           |
   +--------------+-----------------+----------------------------------------+
   | growth_start | Julian day      | First day of the growing season        |
   +--------------+-----------------+----------------------------------------+
   | growth_start | Julian day      | Last day of the growing season         |
   +--------------+-----------------+----------------------------------------+
   | LAI          |                 | Maximum leaf area index                |
   +--------------+-----------------+----------------------------------------+
   | EZD          | cm              | Evaporative zone depth                 |
   +--------------+-----------------+----------------------------------------+
   | CN           |                 | Curve Number                           |
   +--------------+-----------------+----------------------------------------+
   | nlayer       |                 | Number of hydrostratigraphic layers in |
   |              |                 | the soil profile at cell cid           |
   +--------------+-----------------+----------------------------------------+
   | lay_type{i}  |                 | Type of HELP layer of the ith soil     |
   |              |                 | layer                                  |
   +--------------+-----------------+----------------------------------------+
   | thick{i}     |cm               | Thickness of the ith soil layer        |
   +--------------+-----------------+----------------------------------------+
   | poro{i}      | m³/m³           | Total porosity of the ith soil layer   |
   +--------------+-----------------+----------------------------------------+
   | fc{i}        | m³/m³           | Field capacity of the ith soil layer   |
   +--------------+-----------------+----------------------------------------+
   | wp{i}        | m³/m³           | Wilting point of the ith soil layer    |
   +--------------+-----------------+----------------------------------------+
   | ksat         | cm/s            | Saturated hydraulic conductivity of    |
   |              |                 | the ith soil layer                     |
   +--------------+-----------------+----------------------------------------+
   | dist_dr      | m               | Distance to discharge                  |
   +--------------+-----------------+----------------------------------------+
   | slope        | %               | Average slope                          |
   +--------------+-----------------+----------------------------------------+
   | run          |                 | Identify cells that need to be run with|
   |              |                 | the HELP model                         |
   +--------------+-----------------+----------------------------------------+
   | context      |                 | Identify cells by context:             |
   |              |                 |                                        |
   |              |                 | 0. Water cell                          |
   |              |                 | 1. Normal cell                         |
   |              |                 | 2. Stream edge cell with superficial   |
   |              |                 |    hypodermic runoff                   |
   |              |                 | 3. River edge cell with deep           |
   |              |                 |    hypodermic runoff                   |
   |              |                 | 4. Urban cell                          |
   |              |                 | 5. Cell not mapped                     |
   +--------------+-----------------+----------------------------------------+

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