Example
=================================

The first step is to import and create an instance of the
:class:`~pyhelp.HelpManager` class, where we need to pass as argument a path
to a working directory.
The working directory is where the input, output and temporary files are read
and saved by default by the :class:`~pyhelp.HelpManager`.
It must be a location where you have `Write and Read` permissions.
The path to the working directory can be changed at any time with the
:meth:`~pyhelp.HelpManager.set_workdir` method.
Here, we will use the path to the folder `example` that is distributed along
the PyHelp module.

    >>> import os.path as osp
    >>> from pyhelp import HelpManager
    >>> workdir = "C:/Users/User/pyhelp/example"
    >>> helpm = HelpManager(workdir, year_range=(2000, 2010))
    Reading input data grid data from csv... done
    Reading input weather data files... done

During initialization or when setting a new working directory with
:meth:`~pyhelp.HelpManager.set_workdir`, :class:`~pyhelp.HelpManager`
automatically looks in the specified directory and loads the input
data grid and weather data from any valid existing input files.
Input data files with a formatted grid and weather data are provided in the
folder `example`.
Please read the :ref:`sec_data_input` section for more details on how
to prepare the input data manually or with the available tools to generate
these files automatically from other existing sources of data.

Once :class:`~pyhelp.HelpManager` has been instantiated and the input
data grid and weather data loaded successfully, the D4, D7, D10, D11, and D13
input data files that are required by the HELP model to run need to be
generated for each cell of the grid.
These files can be automatically generated from the input grid and weather
data by doing ::

    >>> helpm.build_help_input_files()
    Clearing HELP input files cache... done
    Formatting D10 and D11 data for cell 10 of 10 (100.0%) 
    Task completed in 0.01 sec
    Creating D10 input file for cell 10 of 10 (100.0%) 
    Task completed in 0.99 sec
    Creating D11 input file for cell 10 of 10 (100.0%) 
    Task completed in 0.01 sec
    Saving the connectivity tables... done
    Generating D4 HELP input files for precip... done
    Generating D7 HELP input files for airtemp... done
    Generating D13 HELP input files for solrad... done
    Updating the connection table... done

Note that by default, these files are saved in the folder `help_input_files`
of the working directory.

We can now tell :class:`~pyhelp.HelpManager` to calcul the monthly water budget
for each cells of the grid by doing ::

    >>> help_output_hdf5 = osp.join(workdir, 'help_example.out')
    >>> output_help = helpm.calc_help_cells(help_output_hdf5, tfsoil=-3)
    HELP simulation in progress: 100.0% (0.0 min remaining)     
    Task completed in 1.29 sec
    
The :meth:`~pyhelp.HelpManager.calc_help_cells` method returns a dictionary
with monthtly values of the water budget components for every cell of the
grid. Also, since we passed a hdf5 file name as argument, the data are also
saved in a HDF5 file, with the same structure as the `output_help` dict.
    
Similarly, yearly values of the surface water budget can be calculated for 
the cells that are assumed to be located in surface water bodies by doing ::

    >>> evp_surf = 650
    >>> surf_output_hdf5 = osp.join(workdir, 'surf_example.out')
    >>> output_surf = helpm.calc_surf_water_cells(evp_surf, surf_output_hdf5)
    Calculating budget for water cells... done
    Task completed in 0.02 sec
    
The :meth:`~pyhelp.HelpManager.calc_surf_water_cells` method returns a
dictionary with yearly values of the water budget components for every cell
of the grid that is assumed to be located in surface water bodies.

Various scripts are avaible in `postprocessing.py`, `produce_help_maps.py`, 
`produce_meteo_maps.py`, and `produce_water_budget.py`
to produce a shapefile and various graphs from the results. Note that
the code in these files are in an early stage of development and are subject
to change without notice in the near futur.
