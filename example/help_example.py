# -*- coding: utf-8 -*-
"""
This example shows how to use PyHELP to calculate the monthly water balance
for a section of the Rivière du Nord watershed in the Laurentians, Quebec, Can.
"""

import os.path as osp
from pyhelp.managers import HelpManager


if __name__ == '__main__':
    # For an explanation of why (on Windows) the if __name__ == '__main__'
    # part is necessary, please see :
    #    https://docs.python.org/3.6/library/
    #    multiprocessing.html#programming-guidelines

    # Define the directory where the weather and grid input files are saved.
    workdir = osp.dirname(__file__)

    # Instantiate the HelpManager.
    helpm = HelpManager(workdir)

    # Note that you can access the weather input data through
    # the 'precip_data', 'airtemp_data', and 'solrad_data'  attributes.

    # Generates the input files required by the HELP model
    # for each cell of the grid.
    helpm.build_help_input_files(sf_edepth=0.15)

    # We want to run HELP only for the cells that are located within
    # a jauged subsection of the Rivière du Nord watershed.

    # The field "Bassin" was added to the grid input data to identify the
    # cells that are located within this jauged subsection of the watershed.
    cellnames = helpm.grid.index[helpm.grid['Bassin'] == 1]

    # Run HELP simulation for all the cells in cellnames.

    # Note that the monthly output data will be automatically saved to
    # the HDF5 file define in filepath.
    output = helpm.calc_help_cells(
        path_to_hdf5=osp.join(workdir, 'help_example.out'),
        cellnames=cellnames,
        tfsoil=-3)

    # Export and save annual averages of HELP output values to a csv file.
    output.save_to_csv(osp.join(workdir, 'help_example_yearly.csv'))

    # Plot some results.
    output.plot_area_monthly_avg()
    output.plot_area_yearly_avg()
    output.plot_area_yearly_series()

    # Calculate the yearly water budget for surface water cells.
    output_surf = helpm.calc_surf_water_cells(
        evp_surf=650,
        path_outfile=osp.join(workdir, 'surf_example.out'))
