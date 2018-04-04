# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 08:33:14 2018
@author: jsgosselin
"""

# ---- Third Party imports

import matplotlib.pyplot as plt
import h5py
import netCDF4
import geopandas as gpd
import numpy as np
import os
import os.path as osp

# %% Read HELP results

path_hdf5 = "E:/pyhelp/RADEAU2/help_output_grid0320_BT.out"
hdf5 = h5py.File(path_hdf5, mode='r+')
cellnames = list(hdf5.keys())

# %% Read HELP grid

print('\rReading HELP grid shapefile...', end=' ')
path_shp_help = "E:/pyhelp/RADEAU2/grid_helpXYZ/grid_helpXYZ_bassin.shp"
shp_help = gpd.read_file(path_shp_help)
print('\rReading HELP grid shapefile... done')

# %% Add the results to the dataframe

# Fields : years, rain, runoff, evapo, sub-runoff, percolation, recharge

avg_yearly_precip = np.zeros(len(shp_help['cid']))
avg_yearly_rechg = np.zeros(len(shp_help['cid']))
avg_yearly_runoff = np.zeros(len(shp_help['cid']))
avg_yearly_evapo = np.zeros(len(shp_help['cid']))
avg_year_perco = np.zeros(len(shp_help['cid']))
avg_year_subrun1 = np.zeros(len(shp_help['cid']))
avg_year_subrun2 = np.zeros(len(shp_help['cid']))

N = 15
for i, cellname in enumerate(cellnames):
    print("\rProcessing cell %d of %d" % (i, len(cellnames)), end=' ')
    cid = int(cellname[3:])

    data = hdf5[cellname]
    rain = data['rain'].value
    recharge = data['recharge'].value
    runoff = data['runoff'].value
    evapo = data['evapo'].value
    perco = data['percolation'].value
    subrun1 = data['subrun1'].value
    subrun2 = data['subrun2'].value

    avg_yearly_precip[cid] = np.sum(rain)/N
    avg_yearly_rechg[cid] = np.sum(recharge)/N
    avg_yearly_runoff[cid] = np.sum(runoff)/N
    avg_yearly_evapo[cid] = np.sum(evapo)/N
    avg_year_perco[cid] = np.sum(perco)/N
    avg_year_subrun1[cid] = np.sum(subrun1)/N
    avg_year_subrun2[cid] = np.sum(subrun2)/N

shp_help['precip'] = avg_yearly_precip
shp_help['rechg'] = avg_yearly_rechg
shp_help['runoff'] = avg_yearly_runoff
shp_help['evapo'] = avg_yearly_evapo
shp_help['perco'] = avg_year_perco
shp_help['subrun1'] = avg_year_subrun1
shp_help['subrun2'] = avg_year_subrun2

# %% Save results in a shapefile

crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"
print('\rSaving HELP results in a shapefile...', end=' ')
path_shp_out = ("E:/pyhelp/RADEAU2/"
                "grid_help_result_grid0320/"
                "grid_help_result_grid0320_80reduc_etp.shp")
if not osp.exists(osp.dirname(path_shp_out)):
    os.makedirs(osp.dirname(path_shp_out))
shp_help.to_file(path_shp_out, driver='ESRI Shapefile', crs_wkt=crs)
print('\rSaving HELP results in a shapefile... done')

# %% Produce yearly time series

# cellnames = shp_help['cid'][shp_help['Bassin'] == 2][shp_help['run'] == 1].tolist()
# cellnames = shp_help['cid'][shp_help['BT'] == 0][shp_help['run'] == 1].tolist()
Np = len(cellnames)
avg_monthly_precip = np.zeros((15, 12))
avg_monthly_rechg = np.zeros((15, 12))
avg_monthly_runoff = np.zeros((15, 12))
avg_monthly_evapo = np.zeros((15, 12))
avg_monthly_subrun1 = np.zeros((15, 12))
avg_monthly_perco = np.zeros((15, 12))
nan_cells = []

for i, cellname in enumerate(cellnames):
    print("\rProcessing cell %d of %d" % (i+1, Np), end=' ')
    data = hdf5[cellname]
    if np.any(np.isnan(data['recharge'].value)):
        nan_cells.append(cellname)
        continue

    avg_monthly_precip += data['rain'].value
    avg_monthly_rechg += data['recharge'].value
    avg_monthly_runoff += data['runoff'].value
    avg_monthly_evapo += data['evapo'].value
    avg_monthly_subrun1 += data['sub-runoff'].value
    avg_monthly_perco += data['percolation'].value
print("\rProcessing cell %d of %d" % (i, Np))
print(nan_cells)

avg_monthly_subrun = (avg_monthly_perco - avg_monthly_rechg +
                      avg_monthly_subrun1)


years = range(2000, 2015)
avg_yearly_precip = np.sum(avg_monthly_precip, axis=1) / (Np - len(nan_cells))
avg_yearly_rechg = np.sum(avg_monthly_rechg, axis=1) / (Np - len(nan_cells))
avg_yearly_runoff = np.sum(avg_monthly_runoff, axis=1) / (Np - len(nan_cells))
avg_yearly_evapo = np.sum(avg_monthly_evapo, axis=1) / (Np - len(nan_cells))
avg_yearly_subrun = np.sum(avg_monthly_subrun, axis=1) / (Np - len(nan_cells))


# %%

plt.close('all')
fwidth, fheight = 9, 6.5
fig, ax = plt.subplots()
fig.set_size_inches(fwidth, fheight)

# Setup axe margins :

left_margin = 1.5/fwidth
right_margin = 0.25/fwidth
top_margin = 1./fheight
bot_margin = 0.7/fheight
ax.set_position([left_margin, bot_margin,
                 1 - left_margin - right_margin, 1 - top_margin - bot_margin])

l1, = ax.plot(years, avg_yearly_precip, marker='o', mec='white', clip_on=False,
              lw=2)
l2, = ax.plot(years, avg_yearly_rechg, marker='o', mec='white', clip_on=False,
              lw=2)
l3, = ax.plot(years, avg_yearly_runoff, marker='o', mec='white', clip_on=False,
              lw=2)
l4, = ax.plot(years, avg_yearly_evapo, marker='o', mec='white', clip_on=False,
              lw=2)
l5, = ax.plot(years, avg_yearly_subrun, marker='o', mec='white', clip_on=False,
              lw=2)

ax.tick_params(axis='both', direction='out', labelsize=12)
ax.set_ylabel('Composantes du bilan hydrologique\n(mm/an)',
              fontsize=16, labelpad=10)
ax.set_xlabel('Années', fontsize=16, labelpad=10)
ax.axis(ymin=0, ymax=1400)
ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)

lines = [l1, l2, l3, l4, l5]
labels = ["Précipitations totales", "Recharge au roc",
          "Ruissellement de surface", "Évapotranspiration",
          "Ruissellement hypodermique"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                   borderaxespad=0, loc='lower left', borderpad=0.5,
                   bbox_to_anchor=(0, 1), ncol=2)
legend.draw_frame(False)
# fig.savefig("bilan_hydro_annuel_LAU.pdf")

# %%

plt.close('all')
fwidth, fheight = 8, 6.5
fig, ax = plt.subplots()
fig.set_size_inches(fwidth, fheight)

# Setup axe margins :

left_margin = 1.5/fwidth
right_margin = 0.25/fwidth
top_margin = 0.25/fheight
bot_margin = 0.25/fheight
ax.set_position([left_margin, bot_margin,
                 1 - left_margin - right_margin, 1 - top_margin - bot_margin])

l1 = ax.bar(1, np.mean(avg_yearly_precip), 0.85, align='center')
l2 = ax.bar(2, np.mean(avg_yearly_rechg), 0.85, align='center')
l3 = ax.bar(3, np.mean(avg_yearly_runoff), 0.85, align='center')
l4 = ax.bar(4, np.mean(avg_yearly_evapo), 0.85, align='center')
l5 = ax.bar(5, np.mean(avg_yearly_subrun), 0.85, align='center')
ax.axis(ymin=0, ymax=1400, xmin=0, xmax=6)
ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
ax.set_axisbelow(True)

ax.text(1, np.mean(avg_yearly_precip)+10,
        "%d\nmm/an" % np.mean(avg_yearly_precip), ha='center', va='bottom')
ax.text(2, np.mean(avg_yearly_rechg)+10,
        "%d\nmm/an" % np.mean(avg_yearly_rechg), ha='center', va='bottom')
ax.text(3, np.mean(avg_yearly_runoff)+10,
        "%d\nmm/an" % np.mean(avg_yearly_runoff), ha='center', va='bottom')
ax.text(4, np.mean(avg_yearly_evapo)+10,
        "%d\nmm/an" % np.mean(avg_yearly_evapo), ha='center', va='bottom')
ax.text(5, np.mean(avg_yearly_subrun)+10,
        "%d\nmm/an" % np.mean(avg_yearly_subrun), ha='center', va='bottom')

ax.tick_params(axis='y', direction='out', labelsize=12)
ax.tick_params(axis='x', direction='out', length=0)
ax.set_ylabel('Composantes du bilan hydrologique\n(mm/an)',
              fontsize=16, labelpad=10)

ax.set_xticklabels([])

lines = [l1, l2, l3, l4, l5]
labels = ["Précipitations totales", "Recharge au roc",
          "Ruissellement de surface", "Évapotranspiration",
          "Ruissellement hypodermique"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                   borderaxespad=0, loc='upper right', borderpad=0.5,
                   bbox_to_anchor=(1, 1), ncol=2)
legend.draw_frame(False)
# fig.savefig("bilan_hydro_moyen_annuel_LAU.pdf")


# %%

plt.close('all')
fwidth, fheight = 9, 6.5
fig, ax = plt.subplots()
fig.set_size_inches(fwidth, fheight)

# Setup axe margins :

left_margin = 1.5/fwidth
right_margin = 0.25/fwidth
top_margin = 1./fheight
bot_margin = 0.7/fheight
ax.set_position([left_margin, bot_margin,
                 1 - left_margin - right_margin, 1 - top_margin - bot_margin])

months = range(1, 13)
l1, = ax.plot(months, np.mean(avg_monthly_precip, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l2, = ax.plot(months, np.mean(avg_monthly_rechg, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l3, = ax.plot(months, np.mean(avg_monthly_runoff, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l4, = ax.plot(months, np.mean(avg_monthly_evapo, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l5, = ax.plot(months, np.mean(avg_monthly_subrun, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l6, = ax.plot(months, np.mean(avg_monthly_perco, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)

ax.set_ylabel('Composantes du bilan hydrologique\n(mm/mois)',
              fontsize=16, labelpad=10)
ax.set_xlabel('Mois', fontsize=16, labelpad=10)
ax.axis(ymin=-20, ymax=140)
ax.grid(axis='both', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
ax.set_xticks(months)
ax.set_xticklabels(['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Jun', 'Jul', 'Aoû',
                    'Sep', 'Oct', 'Nov', 'Déc'])
ax.tick_params(axis='both', direction='out', labelsize=12)

lines = [l1, l2, l3, l4, l5, l6]
labels = ["Précipitations totales", "Recharge au roc",
          "Ruissellement de surface", "Évapotranspiration",
          "Ruissellement hypodermique", "Percolation"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                    borderaxespad=0, loc='lower left', borderpad=0.5,
                    bbox_to_anchor=(0, 1), ncol=2)
legend.draw_frame(False)
# fig.savefig('bilan_hydro_mensuel_LAU.pdf')
