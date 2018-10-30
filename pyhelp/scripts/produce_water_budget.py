# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 08:33:14 2018
@author: jsgosselin
"""

# ---- Third Party imports

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import h5py
import geopandas as gpd
import numpy as np
import os.path as osp
import os
from pyhelp.managers import HelpManager
from scipy.stats import linregress

# rname = "LAUALL_inputHELP_0416t3_0.15edepth"
rname = "BTALL_inputHELP_0416base_0.35edepth"
figname_sufix = rname
riv = 2

rootdir = "C:\\Users\\User\\pyhelp\\RADEAU2\\calage_inputHELP_0416"
workdir = osp.join(rootdir, rname)
os.chdir(workdir)

path_help_output = osp.join(workdir, "help_%s.out" % rname)
path_surf_output = osp.join(workdir, "surface_%s.out" % rname)

# %% Load the HELP grid

path_togrid = osp.join(rootdir, 'inputHELP_0416t3.csv')
helpm = HelpManager(workdir, year_range=(1965, 2014))
grid = helpm.load_grid(path_togrid)

# %% Compute yearly values

# Compute monthly values for cells that were run in HELP.

help_output = h5py.File(path_help_output, mode='r+')
cellnames = list(help_output.keys())
years = help_output[cellnames[0]]['years'].value

# urbancells = grid['maille'][grid['Urbain'] == 1][grid['run'] == 1][grid['BT'] == 1].tolist()
# # ruralcells = grid['maille'][grid['Urbain'] == 0][grid['BT'] == 1]
# urbancells = grid['maille'][grid['Urbain'] == 1][grid['Contexte'] == 0][grid['BT'] == 1].tolist()

Np = len(cellnames)
Ny = len(years)

avg_monthly_precip = np.zeros((Ny, 12))
avg_monthly_runoff = np.zeros((Ny, 12))
avg_monthly_evapo = np.zeros((Ny, 12))
avg_monthly_perco = np.zeros((Ny, 12))
avg_monthly_subrun = np.zeros((Ny, 12))
avg_monthly_rechg = np.zeros((Ny, 12))
nan_cells = []


for i, cellname in enumerate(cellnames):
    cellname = str(int(cellname))
    print("\rProcessing cell %d of %d..." % (i+1, Np), end=' ')
    data = help_output[cellname]
    if np.any(np.isnan(data['recharge'].value)):
        nan_cells.append(cellname)
        continue

    avg_monthly_precip += data['rain'].value
    avg_monthly_evapo += data['evapo'].value
    avg_monthly_perco += data['percolation'].value
    avg_monthly_runoff += data['runoff'].value
    avg_monthly_subrun += data['subrun2'].value + data['subrun1'].value
    if grid['context'][int(cellname)] == 2:
        # Convert recharge to subsurface runoff.
        avg_monthly_subrun += data['recharge'].value
    else:
        avg_monthly_rechg += data['recharge'].value
print("done")

avg_yearly_precip = np.sum(avg_monthly_precip, axis=1)
avg_yearly_perco = np.sum(avg_monthly_perco, axis=1)
avg_yearly_rechg = np.sum(avg_monthly_rechg, axis=1)
avg_yearly_runoff = np.sum(avg_monthly_runoff, axis=1)
avg_yearly_evapo = np.sum(avg_monthly_evapo, axis=1)
avg_yearly_subrun = np.sum(avg_monthly_subrun, axis=1)

# Add the yearly values for surface water cells.

surf_output = h5py.File(path_surf_output, mode='r+')
surf_cellnames = list(surf_output.keys())
# surf_cellnames = grid['maille'][grid['Urbain'] == 1][grid['Contexte'] == 0][grid['BT'] == 1].tolist()

Nsc = len(surf_cellnames)

for i, cellname in enumerate(surf_cellnames):
    cellname = str(int(cellname))
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
avg_yearly_subrun /= Ntot
avg_yearly_rechg /= Ntot

summary = [np.mean(avg_yearly_precip),
           np.mean(avg_yearly_runoff),
           np.mean(avg_yearly_evapo),
           np.mean(avg_yearly_perco),
           np.mean(avg_yearly_subrun),
           np.mean(avg_yearly_rechg),
           '']

# %% Plot debits

if riv == 1:
    years_debits = np.array([
        1965, 1966, 1967, 1968, 1969,
        1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979,
        1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989,
        1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999,
        2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
        2010, 2011, 2012, 2013, 2014])
    qcehq_tot = np.array([
        562, 583, 701, 477, 676,
        579, 533, 913, 827, 769, 676, 782, 642, 587, 884,
        657, 827, 533, 799, 672, 553, 622, 514, 706, 583,
        691, 569, 633, 736, 689, 537, 858, 616, 604, 668,
        696, 497, 522, 715, 556, 727, 1004, 538, 825, 848,
        641, 739, 505, 760, 812])
    qcehq_base = np.array([
        282, 290, 350, 241, 337,
        291, 265, 456, 412, 386, 335, 395, 319, 295, 437,
        333, 414, 261, 404, 335, 279, 311, 256, 353, 292,
        343, 287, 316, 367, 346, 269, 426, 312, 300, 333,
        348, 248, 263, 353, 279, 366, 499, 271, 408, 429,
        319, 369, 254, 381, 404])
elif riv == 2:
    years_debits = np.array([
        1973, 1974, 1975, 1976, 1977, 1978, 1979,
        1980, 1981, 1982, 1983, 1984, np.nan,
        2011, 2012, 2013, 2014])
    qcehq_tot = np.array([
        646, 572, 584, 615, 514, 482, 631,
        398, 639, 417, 640, 551, np.nan,
        637, 372, 631, 617])
    qcehq_base = np.array([
        321, 286, 293, 309, 256, 242, 308,
        205, 320, 203, 322, 277, np.nan,
        316, 190, 316, 303])

qhelp_total = avg_yearly_runoff + avg_yearly_subrun + avg_yearly_rechg
qhelp_base = avg_yearly_subrun + avg_yearly_rechg

fwidth, fheight = 9, 5.5
fig, ax = plt.subplots()
fig.set_size_inches(fwidth, fheight)

left_margin = 1.5/fwidth
right_margin = 0.25/fwidth
top_margin = 0.5/fheight
bot_margin = 0.7/fheight
ax.set_position([left_margin, bot_margin,
                 1 - left_margin - right_margin, 1 - top_margin - bot_margin])

# Streamflow total

l1, = ax.plot(years_debits, qcehq_tot, marker=None, mec='white',
              lw=2, linestyle='--',  clip_on=True, color='#3690c0')

l2, = ax.plot(years, qhelp_total, marker=None, mec='white', clip_on=True,
              lw=2, linestyle='-', color='#034e7b')

# Streamflow base

l3, = ax.plot(years_debits, qcehq_base, marker=None, mec='white',
              lw=2, linestyle='--',  clip_on=True, color='#ef6548')

l4, = ax.plot(years, qhelp_base, marker=None, mec='white', clip_on=True,
              lw=2, linestyle='-', color='#990000')


ax.tick_params(axis='both', direction='out', labelsize=12)
ax.set_ylabel('Débit par unité de surface\n(mm/an)',
              fontsize=16, labelpad=10)
ax.set_xlabel('Années', fontsize=16, labelpad=10)
ax.axis(ymin=0, ymax=1600)
ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
ax.axis([1960, 2017, 0, 1200])

lines = [l1, l2, l3, l4]
labels = ["CEHQ débit total", "HELP débit total",
          "CEHQ débit base", "HELP débit base"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                   borderaxespad=0, loc='upper left', borderpad=0.5,
                   bbox_to_anchor=(0, 1), ncol=2)
legend.draw_frame(False)

# Add a graph title.
offset = transforms.ScaledTranslation(0/72, 12/72, fig.dpi_scale_trans)
ax.text(0.5, 1, figname_sufix, fontsize=16, ha='center', va='bottom',
        transform=ax.transAxes+offset)

fig.savefig('debit_temps_%s.pdf' % figname_sufix)

# %% PLot debits scatter

# Calcul the RMSE

years_debits = years_debits[~np.isnan(years_debits)]
qcehq_tot = qcehq_tot[~np.isnan(qcehq_tot)]
qcehq_base = qcehq_base[~np.isnan(qcehq_base)]

indx = np.digitize(years_debits, years, right=True)
rmse_qtot = np.mean((qcehq_tot - qhelp_total[indx])**2)**0.5
me_qtot = np.mean(qhelp_total[indx] - qcehq_tot)
rmse_qbase = np.mean((qcehq_base - qhelp_base[indx])**2)**0.5
me_qbase = np.mean(qhelp_base[indx] - qcehq_base)

summary.append(rmse_qbase)
summary.append(rmse_qtot)
summary.append(me_qbase)
summary.append(me_qtot)

# Plot the results

fwidth, fheight = 5, 5
fig, ax = plt.subplots()
fig.set_size_inches(fwidth, fheight)

left_margin = 1/fwidth
right_margin = 0.5/fwidth
top_margin = 0.5/fheight
bot_margin = 1/fheight
ax.set_position([left_margin, bot_margin,
                 1 - left_margin - right_margin, 1 - top_margin - bot_margin])

xymin, xymax = 100, 1100
ax.axis([xymin, xymax, xymin, xymax])
ax.set_ylabel('Débits HELP (mm/an)', fontsize=16, labelpad=15)
ax.set_xlabel('Débits CEHQ (mm/an)', fontsize=16, labelpad=15)

ax.tick_params(axis='both', direction='out', labelsize=12)

l1, = ax.plot([xymin, xymax], [xymin, xymax], '--', color='black', lw=1)
l2, = ax.plot(qcehq_tot, qhelp_total[indx], '.', color='#3690c0')
l3, = ax.plot(qcehq_base, qhelp_base[indx], '.', color='#990000')

# Plot the model fit stats.

dx, dy = 3, -3
offset = transforms.ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
ax.text(0, 1, "RMSE débit total = %0.1f mm/an" % rmse_qtot,
        transform=ax.transAxes+offset, ha='left', va='top')
dy += -12
offset = transforms.ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
ax.text(0, 1, "ME débit total = %0.1f mm/an" % me_qtot,
        transform=ax.transAxes+offset, ha='left', va='top')

dy += -18
offset = transforms.ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
ax.text(0, 1, "RMSE débit base = %0.1f mm/an" % rmse_qbase,
        transform=ax.transAxes+offset, ha='left', va='top')
dy += -12
offset = transforms.ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
ax.text(0, 1, "ME débit base = %0.1f mm/an" % me_qbase,
        transform=ax.transAxes+offset, ha='left', va='top')

# Add a graph title.
offset = transforms.ScaledTranslation(0/72, 12/72, fig.dpi_scale_trans)
ax.text(0.5, 1, figname_sufix, fontsize=16, ha='center', va='bottom',
        transform=ax.transAxes+offset)

# Add a legend.
lines = [l1, l2, l3]
labels = ["1:1", "Débit total", "Débit base"]
legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                   borderaxespad=0, loc='lower right', borderpad=0.5,
                   bbox_to_anchor=(1, 0), ncol=1)
legend.draw_frame(False)


fig.savefig('calage_'+rname+'.pdf')
