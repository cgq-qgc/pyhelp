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

    # Generates the input files required by the HELP model for each cell
    # of the grid.
    helpm.build_help_input_files()

    # It is possible to run HELP only for a subset of the study area
    # grid cells. To do this, we need to pass a list of cell ids as
    # an argument to the calc_help_cells method.
    # For example, we can run here the results for the first 100 cells
    # of the grid.
    cellnames = helpm.cellnames[:100]

    # If we pass a filename in argument to the calc_help_cells method, the
    # monthly output data will be automatically saved to disk as an HDF5 file.
    filepath = osp.join(workdir, 'help_example.out')
    output = helpm.calc_help_cells(filepath, cellnames, tfsoil=-3)

    # Export and save the data to an ESRI shapefile.
    filepath = osp.join(workdir, 'help_example_yearly.csv')
    output.save_to_csv(filepath)

    # Plot some results.
    output.plot_area_monthly_avg()
    output.plot_area_yearly_avg()
    output.plot_area_yearly_series()

    # Calculate the yearly water budget for surface water cells.
    evp_surf = 650
    filepath = osp.join(workdir, 'surf_example.out')
    output_surf = helpm.calc_surf_water_cells(evp_surf, filepath)
