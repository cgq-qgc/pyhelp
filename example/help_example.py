# -*- coding: utf-8 -*-
"""
This example shows how to use PyHELP to calculate the monthly water balance
for a section of the Rivière du Nord watershed in the Laurentians, Quebec, Can.

Updated for PyHELP version 0.3.1
"""

import os.path as osp
import pandas as pd
from pyhelp.managers import HelpManager
import pyhelp.bilan as HelpBilan

if __name__ == '__main__':
    # For an explanation of why (on Windows) the if __name__ == '__main__'
    # part is necessary, please see :
    #    https://docs.python.org/3.6/library/
    #    multiprocessing.html#programming-guidelines

    # Define the working directory.
    workdir = "D:/Projets/pyhelp/example/"

    # Instantiate the HelpManager and provide the paths to the grid and
    # weather input data files so that they are loaded automatically.
    helpm = HelpManager(
        workdir,
        path_to_grid=workdir + 'input_grid.csv',
        path_to_precip=workdir + 'precip_input_data.csv',
        path_to_airtemp=workdir + 'airtemp_input_data.csv',
        path_to_solrad=workdir + 'solrad_input_data.csv')

    # Note that you can access the grid input data through
    # the 'grid' attribute of the HelpManager.

    # Note that you can access the weather input data through the
    # 'precip_data', 'airtemp_data', and 'solrad_data' attributes
    # of the HelpManager.

    # =========================================================================
    # Run HELP simulation for all the cells in cellnames.
    # =========================================================================

    # We want to run HELP only for the cells that are located within
    # a jauged subsection of the Rivière du Nord watershed.

    # The field "Bassin" was added to the grid input data to identify the
    # cells that are located within this jauged subsection of the watershed.
    cellnames = helpm.grid.index[helpm.grid['Bassin'] == 1]

    # Note that the monthly output data will be automatically saved to
    # the HDF5 file define in filepath.
    output_help = helpm.calc_help_cells(
        path_to_hdf5=workdir + 'help_example.out',
        cellnames=cellnames,
        tfsoil=-3,
        sf_edepth=0.15,
        sf_ulai=1)

    # Export and save annual averages of HELP output values to a csv file.
    output_help.save_to_csv(osp.join(workdir, 'help_example_yearly.csv'))

    # Plot some results.
    output_help.plot_area_monthly_avg(fig_title="PyHELP Example")
    output_help.plot_area_yearly_avg(fig_title="PyHELP Example")
    output_help.plot_area_yearly_series(fig_title="PyHELP Example")

    # =========================================================================
    # Compare with river total and base streamflow
    # =========================================================================

    # Calculate the yearly water budget for surface water cells.
    output_surf = helpm.calc_surf_water_cells(
        evp_surf=650,
        path_outfile=osp.join(workdir, 'surf_example.out'))

    # Read observed yearly total and base streamflow (in mm/year).
    obs_qflow = pd.read_csv(
        osp.join(workdir, 'obs_yearly_river_flow.csv'),
        index_col=0)

    # Calcul simulated early total and base streamflow (in mm/year).
    sim_qflow = HelpBilan.calc_yearly_streamflow(output_help, output_surf)

    # Plot results.
    HelpBilan.plot_sim_vs_obs_yearly_streamflow(
        sim_qflow, obs_qflow, fig_title="PyHELP Example")
    HelpBilan.plot_streamflow_scatter(
        sim_qflow, obs_qflow, fig_title="PyHELP Example")
