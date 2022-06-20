# -*- coding: utf-8 -*-
"""
Script permettant de faire le post-processing des données de la température
moyenne journalière pour une simulation pyhelp donnée.
"""

from pyhelp.managers import (
    load_weather_from_csv, calc_dist_from_coord, HelpOutput)
import numpy as np
import pandas as pd


# Chargement de la grille des données journalières de la température
# moyenne de l'air.
airtemp_data = load_weather_from_csv(
    'D:/Projets/pyhelp/example/airtemp_input_data.csv')
airtemp_lat = airtemp_data.columns.get_level_values('lat_dd').values
airtemp_lon = airtemp_data.columns.get_level_values('lon_dd').values

# Chargement des résultats de la simulation pyhelp pour laquelle nous voulons
# faire le post-processing de la température moyenne de l'air.
output_help = HelpOutput(
    "D:/Projets/pyhelp/example/help_example.out")

grid_cells = output_help.data['cid']
grid_lat = output_help.data['lat_dd']
grid_lon = output_help.data['lon_dd']

# Association des cellules de la grille de la zone d'étude aux cellules de
# la grille de la température moyenne de l'air.
connect_table = {}
for grid_index in range(len(grid_cells)):
    dist = calc_dist_from_coord(
        grid_lat[grid_index], grid_lon[grid_index],
        airtemp_lat, airtemp_lon)
    airtemp_index = int(np.argmin(dist))

    if airtemp_index not in connect_table:
        connect_table[airtemp_index] = [grid_index]
    else:
        connect_table[airtemp_index].append(grid_index)

# Mise en forme des données journalières de la température moyenne de l'air
# dans un bloc 3D de données où l'axe 0 correspond aux cellules de la
# grille pyhelp, l'axe 1 aux années et l'axe 2 aux mois.
airtemp_monthly = np.zeros(
    (len(grid_cells), len(airtemp_data.index.year.unique()), 12)
    ) * np.nan
for airtemp_index, grid_indexes in connect_table.items():
    subset_airtemp_daily = pd.DataFrame(
        data=airtemp_data.values[:, airtemp_index],
        index=airtemp_data.index
        )
    subset_monthly_airtemp = subset_airtemp_daily.groupby(
        [subset_airtemp_daily.index.year, subset_airtemp_daily.index.month]
        ).mean().unstack().values

    airtemp_monthly[grid_indexes, :, :] = subset_monthly_airtemp


# Calcul des valeurs annuelles moyennes aux cellules de la
# grille pyhelp pour chaque année.
airtemp_yearly = pd.DataFrame(
    data=np.mean(airtemp_monthly, axis=2),
    index=grid_cells,
    columns=airtemp_data.index.year.unique()
    )
airtemp_yearly.index.name = 'cid'
airtemp_yearly.insert(0, 'lon_dd', grid_lat)
airtemp_yearly.insert(0, 'lat_dd', grid_lon)

airtemp_yearly.to_csv(
    "D:/Projets/pyhelp/example/airtemp_yearly.csv")

# Calcul des valeurs annuelles moyennes aux cellules de la grille
# pyhelp sur une période donnée.
year_from = 2003
year_to = 2006
airtemp_yearly_period = airtemp_yearly[['lat_dd', 'lon_dd']].copy()
airtemp_yearly_period[f'{year_from}-{year_to}'] = airtemp_yearly[
    list(range(year_from, year_to + 1))].mean(axis=1)

airtemp_yearly_period.to_csv(
    f"D:/Projets/pyhelp/example/airtemp_yearly_{year_from}_{year_to}.csv")

# Calcul et sauvegarde des valeurs annuelles moyennes de la température de
# l'air pour l'ensemble la région d'étude.
airtemp_yearly_avg = pd.DataFrame(
    data=np.mean(airtemp_monthly, axis=(0, 2)),
    index=airtemp_data.index.year.unique(),
    columns=['airtemp_avg']
    )
airtemp_yearly_avg.index.name = 'year'

airtemp_yearly_avg.to_csv(
    "D:/Projets/pyhelp/example/airtemp_yearly_mean.csv")

# Calcul et sauvegarde des valeurs mensuelles moyennes de la température
# de l'air pour toute la région d'étude.
airtemp_monthly_avg = pd.DataFrame(
    data=np.mean(airtemp_monthly, axis=0),
    index=airtemp_data.index.year.unique(),
    columns=range(1, 13)
    )
airtemp_monthly_avg.index.name = 'year'
airtemp_monthly_avg.to_csv(
    "D:/Projets/pyhelp/example/airtemp_monthly_mean.csv")
