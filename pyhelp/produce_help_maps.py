# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 08:33:14 2018
@author: jsgosselin
"""

# ---- Third Party imports

import matplotlib.pyplot as plt
import h5py
import geopandas as gpd
import numpy as np
import os
import os.path as osp

rname = "BT-0406_edepth60"

os.chdir("C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018")
workdir = "C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018\\%s\\" % rname
path_help_output = workdir + "help_%s.out" % rname
path_surf_output = workdir + "surface_%s.out" % rname

# %% Read HELP results

help_output = h5py.File(path_help_output, mode='r+')
helpcells = list(help_output.keys())
years = help_output[helpcells[0]]['years'].value

surf_output = h5py.File(path_surf_output, mode='r+')
surfcells = list(surf_output.keys())

# %% Read HELP grid

print('\rReading HELP grid shapefile...', end=' ')
path_shp = ("C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018" +
            "\\XYgrille0405\\XYgrille0405.shp")
shp_help = gpd.read_file(path_shp)
print('\rReading HELP grid shapefile... done')

# %% Add the results to the dataframe

# Fields : years, rain, runoff, evapo, sub-runoff, percolation, recharge

ncells = len(shp_help['maille'])
nyears = len(years)

avg_yearly_precip = np.zeros(ncells) * np.nan
avg_yearly_rechg = np.zeros(ncells) * np.nan
avg_yearly_runoff = np.zeros(ncells) * np.nan
avg_yearly_evapo = np.zeros(ncells) * np.nan
avg_year_perco = np.zeros(ncells) * np.nan
avg_year_subrun1 = np.zeros(ncells) * np.nan
avg_year_subrun2 = np.zeros(ncells) * np.nan

# Handle HELP results.

for i, cellname in enumerate(helpcells):
    print("\rProcessing cell %d of %d..." % (i, len(helpcells)), end=' ')
    cid = int(cellname)

    data = help_output[cellname]
    rain = np.sum(data['rain'].value) / nyears
    runoff = np.sum(data['runoff'].value) / nyears
    evapo = np.sum(data['evapo'].value) / nyears
    perco = np.sum(data['percolation'].value) / nyears
    subrun1 = np.sum(data['subrun1'].value) / nyears
    subrun2 = np.sum(data['subrun2'].value) / nyears
    recharge = np.sum(data['recharge'].value) / nyears

    avg_yearly_precip[cid] = rain
    avg_yearly_runoff[cid] = runoff
    avg_yearly_evapo[cid] = evapo
    avg_year_perco[cid] = perco

    if shp_help['Contexte'][cid] == 2:
        # Convert recharge to runoff when cells are close to a stream.
        if subrun2 == 0:
            # Convert recharge as surficial subrunoff.
            subrun1 = subrun1 + recharge
        else:
            # This means there is a layer of sand above the clay layer.
            # Convert recharge as deep runoff.
            subrun2 = subrun2 + recharge
        recharge = 0

    avg_year_subrun1[cid] = subrun1
    avg_year_subrun2[cid] = subrun2
    avg_yearly_rechg[cid] = recharge
print("done")

# Handle Surface Water results.

for i, cellname in enumerate(surfcells):
    print("\rProcessing surface cell %d of %d..." % (i, len(surfcells)),
          end=' ')
    cid = int(cellname)

    data = surf_output[cellname]
    avg_yearly_precip[cid] = np.sum(data['rain'].value) / nyears
    avg_yearly_runoff[cid] = np.sum(data['runoff'].value) / nyears
    avg_yearly_evapo[cid] = np.sum(data['evapo'].value) / nyears

    avg_year_perco[cid] = 0
    avg_year_subrun1[cid] = 0
    avg_year_subrun2[cid] = 0
    avg_yearly_rechg[cid] = 0
print("done")

# Add the results to the shapefile.

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
path_shp_out = ("C:/Users/User/pyhelp/RADEAU2/"
                "help_result_BT-0406_edepth60/"
                "help_result_BT-0406_edepth60.shp")
if not osp.exists(osp.dirname(path_shp_out)):
    os.makedirs(osp.dirname(path_shp_out))
shp_help.to_file(path_shp_out, driver='ESRI Shapefile')
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
