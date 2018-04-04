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
import os.path as osp
import os

rname = 'RVC-0403'
figname_sufix = rname
CONTEXT = True

os.chdir("C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018")
workdir = "C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018\\%s\\" % rname
path_help_output = workdir + "help_%s.out" % rname
path_surf_output = workdir + "surface_%s.out" % rname

# %% Load the HELP grid

grid = gpd.read_file("C:\\Users\\User\\pyhelp\\RADEAU2\\Calage_27mars2018" +
                     "\\XYgrille0327\\XYgrille0327.shp")
grid.set_index(['cid'], drop=False, inplace=True)

# %% Compute yearly values

# Compute monthly values for cells that were run in HELP.

help_output = h5py.File(path_help_output, mode='r+')
cellnames = list(help_output.keys())
years = help_output[cellnames[0]]['years'].value

Np = len(cellnames)
Ny = len(years)

avg_monthly_precip = np.zeros((Ny, 12))
avg_monthly_runoff = np.zeros((Ny, 12))
avg_monthly_evapo = np.zeros((Ny, 12))
avg_monthly_perco = np.zeros((Ny, 12))
avg_monthly_subrun2 = np.zeros((Ny, 12))
avg_monthly_rechg = np.zeros((Ny, 12))
nan_cells = []

for i, cellname in enumerate(cellnames):
    print("\rProcessing cell %d of %d..." % (i+1, Np), end=' ')
    data = help_output[cellname]
    if np.any(np.isnan(data['recharge'].value)):
        nan_cells.append(cellname)
        continue

    avg_monthly_precip += data['rain'].value
    avg_monthly_evapo += data['evapo'].value
    avg_monthly_perco += data['percolation'].value
    avg_monthly_runoff += data['runoff'].value + data['subrun1'].value
    avg_monthly_subrun2 += data['subrun2'].value
    if grid['context'][int(cellname)] == 2 and CONTEXT:
        # Convert recharge to runoff.
        if np.sum(data['subrun2'].value) == 0:
            # Convert recharge as surficial runoff.
            avg_monthly_runoff += data['recharge'].value
        else:
            # This means there is a layer of sand above the clay layer.
            # Convert recharge as deep runoff.
            avg_monthly_subrun2 += data['recharge'].value
    else:
        avg_monthly_rechg += data['recharge'].value
print("done")

avg_yearly_precip = np.sum(avg_monthly_precip, axis=1)
avg_yearly_perco = np.sum(avg_monthly_perco, axis=1)
avg_yearly_rechg = np.sum(avg_monthly_rechg, axis=1)
avg_yearly_runoff = np.sum(avg_monthly_runoff, axis=1)
avg_yearly_evapo = np.sum(avg_monthly_evapo, axis=1)
avg_yearly_subrun2 = np.sum(avg_monthly_subrun2, axis=1)

# Add the yearly values for surface water cells.

surf_output = h5py.File(path_surf_output, mode='r+')
surf_cellnames = list(surf_output.keys())

Nsc = len(surf_cellnames)

for i, cellname in enumerate(surf_cellnames):
    print("\rProcessing surf. cell %d of %d..." % (i+1, Nsc), end=' ')
    data = surf_output[cellname]
    avg_yearly_precip += data['rain'].value
    avg_yearly_runoff += data['runoff'].value
    avg_yearly_evapo += data['evapo'].value
print("done")

# Scale the results with the total numer of cells so we get mm/year.

Ntot = Nsc + Np - len(nan_cells)

avg_yearly_precip /= Ntot
avg_yearly_runoff /= Ntot
avg_yearly_evapo /= Ntot
avg_yearly_perco /= Ntot
avg_yearly_subrun2 /= Ntot
avg_yearly_rechg /= Ntot

summary = [np.mean(avg_yearly_precip),
           np.mean(avg_yearly_runoff),
           np.mean(avg_yearly_evapo),
           np.mean(avg_yearly_perco),
           np.mean(avg_yearly_subrun2),
           np.mean(avg_yearly_rechg)]


# %% Yearly averages time series


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
l5, = ax.plot(years, avg_yearly_subrun2, marker='o', mec='white',
              clip_on=False, lw=2)

# Plot the observations

# Riv. du Chêne
base_years = np.array([1980, 1981, 1982, 1983, 1984, 2011, 2012, 2013, 2014])
base_evapo = [642.6, 498.5, 509.2, 545.8, 331.0, 551.2, 593.9, 490.5, 519.4]

# Riv. du Nord
# base_years = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
#               2010, 2011, 2012, 2013, 2014]
# base_evapo = [424.2, 400.8, 452.1, 517.4, 398.7, 381.4, 380.8, 521.0, 453.5,
#               262.5, 503.6, 459.4, 572.9, 411.8, 417.6]

l6, = ax.plot(base_years, base_evapo, marker='^', mec='white',
              clip_on=False, lw=2, linestyle='--')

ax.tick_params(axis='both', direction='out', labelsize=12)
ax.set_ylabel('Composantes du bilan hydrologique\n(mm/an)',
              fontsize=16, labelpad=10)
ax.set_xlabel('Années', fontsize=16, labelpad=10)
ax.axis(ymin=0, ymax=1600)
ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)

lines = [l1, l2, l3, l4, l5, l6]
labels = ["Précipitations totales", "Recharge au roc",
          "Ruissellement de surface", "Évapotranspiration",
          "Ruissellement hypodermique", "Évapotranspiration observée"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                   borderaxespad=0, loc='lower left', borderpad=0.5,
                   bbox_to_anchor=(0, 1), ncol=2)
legend.draw_frame(False)
fig.savefig("bilan_hydro_annuel_%s.pdf" % figname_sufix)

# %% Yearly Average Barplot

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
l5 = ax.bar(5, np.mean(avg_yearly_subrun2), 0.85, align='center')
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
ax.text(5, np.mean(avg_yearly_subrun2)+10,
        "%d\nmm/an" % np.mean(avg_yearly_subrun2), ha='center', va='bottom')

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
fig.savefig("bilan_hydro_moyen_annuel_%s.pdf" % figname_sufix)

# %% Comparison with stream hydrograph

plt.close('all')

base_years = [2011, 2012, 2013, 2014]
base_evapo = [551.2, 593.9, 490.5, 519.4]
base_runoff = [320.6, 182.4, 315.1, 313.9]
base_baseflow = [316.0, 189.9, 315.7, 304.0] 

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

l0, = ax.plot(years, avg_yearly_precip, marker='o', mec='white', clip_on=False,
              lw=2)
    
l1, = ax.plot(years, avg_yearly_evapo, marker='o', mec='white', clip_on=False,
              lw=2, color='green')
l1b, = ax.plot(base_years, base_evapo, marker='s', mec='white', clip_on=False,
               lw=2, linestyle='--', color='green')


# l1, = ax.plot(years, avg_yearly_evapo, marker='o', mec='white', clip_on=False,
#               lw=2)
# l1b, = ax.plot(base_years, base_evapo, marker='o', mec='white', clip_on=False,
#                lw=2, linestyle='--')

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
l5, = ax.plot(months, np.mean(avg_monthly_subrun2, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)
l6, = ax.plot(months, np.mean(avg_monthly_perco, axis=0)/(Np-len(nan_cells)),
              marker='o', mec='white', clip_on=False, lw=2)

ax.set_ylabel('Composantes du bilan hydrologique\n(mm/mois)',
              fontsize=16, labelpad=10)
ax.set_xlabel('Mois', fontsize=16, labelpad=10)
ax.axis(ymin=-5, ymax=140)
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
fig.savefig('bilan_hydro_mensuel_%s.pdf' % figname_sufix)
