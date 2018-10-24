# -*- coding: utf-8 -*-
"""
A simple example for PyHelp.
"""

import os.path as osp
from pyhelp.managers import HelpManager

# %% Instantiate HelpManager

workdir = "C:/Users/User/pyhelp/example"
helpm = HelpManager(workdir, year_range=(2000, 2010))

# %% Generate input global solar radiation from CWEEDS

cweed2_paths = osp.join(workdir, 'CWEEDS', '94792.WY2')
cweed3_paths = osp.join(workdir, 'CWEEDS', '94792.WY3')
helpm.generate_weather_inputs_from_CWEEDS(cweed2_paths, cweed3_paths)
helpm.load_weather_input_data()

# %% Generate precip and air temp from MDDELCC meteo grid

path_to_mddelcc_grid = "F:/MeteoGrilleDaily"
helpm.generate_weather_inputs_from_MDELCC_grid(path_to_mddelcc_grid)
helpm.load_weather_input_data()

# %% Calculate the monthly water budget with HELP

helpm.build_help_input_files()

help_output_hdf5 = osp.join(workdir, 'help_example.out')
output_help = helpm.calc_help_cells(help_output_hdf5, tfsoil=-3)

# %% Plot some results

output_help.plot_area_monthly_avg()
output_help.plot_area_yearly_avg

# %% Calculate the yearly water budget for surface water cells.

evp_surf = 650
surf_output_hdf5 = osp.join(workdir, 'surf_example.out')
output_surf = helpm.calc_surf_water_cells(evp_surf, surf_output_hdf5)
