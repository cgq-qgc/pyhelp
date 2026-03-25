# -*- coding: utf-8 -*-
"""
This example shows how to use PyHELP to calculate the monthly water balance
for a section of the Rivière du Nord watershed in the Laurentians, Quebec, Can.

Updated for PyHELP version 0.5.0
"""

from pathlib import Path
import pandas as pd
from pyhelp import __rootdir__, HelpManager
import pyhelp.bilan as HelpBilan

if __name__ == '__main__':
    # For an explanation of why (on Windows) the if __name__ == '__main__'
    # part is necessary, please see :
    #    https://docs.python.org/3.6/library/
    #    multiprocessing.html#programming-guidelines

    # Define the working directory here.
    workdir = Path(__rootdir__) / 'example'
    print(workdir)

    # Instantiate the HelpManager and provide the paths to the grid and
    # weather input data files so that they are loaded automatically.
    helpm = HelpManager(
        str(workdir),
        path_to_grid=str(workdir / 'input_grid.csv'),
        path_to_precip=str(workdir / 'precip_input_data.csv'),
        path_to_airtemp=str(workdir / 'airtemp_input_data.csv'),
        path_to_solrad=str(workdir / 'solrad_input_data.csv')
        )

    # Note that you can access the grid input data through
    # the 'grid' attribute of the HelpManager.

    # Note that you can access the weather input data through the
    # 'precip_data', 'airtemp_data', and 'solrad_data' attributes
    # of the HelpManager.

    # %%

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
        path_to_hdf5=str(workdir / 'help_example.out'),
        cellnames=cellnames,
        tfsoil=-3,
        sf_edepth=0.15,
        sf_ulai=1,
        sf_cn=1.15,
        # Don't write D10/D11 HELP30 input files. See the documentation.
        write_help_input_files=False,
        # Don' write HELP30 .OUT files. See the documentation.
        write_help_output_files=False,
        )

    # Export and save annual averages of HELP output values to a csv file.
    output_help.save_to_csv(workdir / 'help_example_yearly.csv')

    # Plot some results.
    output_help.plot_area_monthly_avg(fig_title="PyHELP Example")
    output_help.plot_area_yearly_avg(fig_title="PyHELP Example")
    output_help.plot_area_yearly_series(fig_title="PyHELP Example")

    # %%

    # =========================================================================
    # Compare with river total and base streamflow
    # =========================================================================

    # Calculate the yearly water budget for surface water cells.
    output_surf = helpm.calc_surf_water_cells(
        cellnames=cellnames,
        evp_surf=650,
        path_outfile=workdir / 'surf_example.out'
        )

    # Read observed yearly total and base streamflow (in mm/year).
    obs_qflow = pd.read_csv(
        workdir / 'obs_yearly_river_flow.csv',
        index_col=0)

    # Calcul simulated early total and base streamflow (in mm/year).
    sim_qflow = HelpBilan.calc_yearly_streamflow(output_help, output_surf)

    # Plot results.
    HelpBilan.plot_sim_vs_obs_yearly_streamflow(
        sim_qflow, obs_qflow, fig_title="PyHELP Example")
    HelpBilan.plot_streamflow_scatter(
        sim_qflow, obs_qflow, fig_title="PyHELP Example")
