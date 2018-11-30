# -*- coding: utf-8 -*-
"""
This example shows how to use PyHELP to calculate the monthly water balance
for a section of the Rivière du Nord watershed in the Laurentians, Quebec,
Canada.
"""

import os.path as osp
from pyhelp.managers import HelpManager


def generate_solrad_input_data_from_CWEEDS(helpm):
    """
    Generate input global solar radiation from the two CWEEDS files that are
    distributed with the example.
    CWEEDS files for more location across Canada are available here:
    http://climate.weather.gc.ca/prods_servs/engineering_e.html
    """
    cweed2_paths = osp.join(workdir, 'CWEEDS', '94792.WY2')
    cweed3_paths = osp.join(workdir, 'CWEEDS', '94792.WY3')
    helpm.generate_weather_inputs_from_CWEEDS(cweed2_paths, cweed3_paths)
    helpm.load_weather_input_data()


def generate_weather_input_data_from_InfoClimat(helpm):
    """
    Generate input daily precipitation and air temperature data from a
    spatially distributed grid produced by the Info-Climat service of the
    Governement of Quebec.
    http://www.environnement.gouv.qc.ca/climat/surveillance/produits.htm
    """
    path_to_mddelcc_grid = "path/to/folder/with/infoclimat/netcdf/files"
    helpm.generate_weather_inputs_from_MDELCC_grid(path_to_mddelcc_grid)
    helpm.load_weather_input_data()


def run_help(helpm):
    """
    Calculate the monthly water budget with HELP.
    """
    # It is possible to run HELP only for a subset of the study area
    # grid cells. To do this, we need to pass a list of cell ids as
    # an argument to the calc_help_cells method.
    # For example, we can run here the results for the first 100 cells
    # of the grid.
    cellnames = helpm.cellnames[:100]

    # If we pass a filename in argument to the calc_help_cells method, the
    # monthly output data will be automatically saved to disk as an HDF5
    # file.
    help_output_hdf5 = osp.join(helpm.workdir, 'help_example.out')

    helpm.build_help_input_files()
    output = helpm.calc_help_cells(help_output_hdf5, cellnames, tfsoil=-3)
    return output


if __name__ == '__main__':
    # For an explanation of why (on Windows) the if __name__ == '__main__'
    # part is necessary, please see :
    #    https://docs.python.org/3.6/library/
    #    multiprocessing.html#programming-guidelines

    # Instantiate the HelpManager and calculate the monthrly water budget
    # for the 2000-2010 period.
    workdir = "C:/Users/User/pyhelp/example"
    helpm = HelpManager(workdir, year_range=(2000, 2010))
    output = run_help(helpm)

    # Export and save the data to an ESRI shapefile.
    help_output_shp = osp.join(workdir, 'help_example.shp')
    output.save_to_shp(help_output_shp)

    # Plot some results.
    output.plot_area_monthly_avg()
    output.plot_area_yearly_avg()
    output.plot_area_yearly_series()

    # Calculate the yearly water budget for surface water cells.
    evp_surf = 650
    surf_output_hdf5 = osp.join(workdir, 'surf_example.out')
    output_surf = helpm.calc_surf_water_cells(evp_surf, surf_output_hdf5)
