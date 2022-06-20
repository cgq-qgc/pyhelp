# -*- coding: utf-8 -*-
"""
A script to compare simulated with observed river yearly total and base flow.
"""

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd


def calc_yearly_streamflow(help_output, surf_output):
    """
    Calcul the simulated yearly total and base streamflow in mm/year.
    """
    # Calculate the cumulative yearly values for all the cells of the
    # study area that were un with HELP.
    help_cellnames = help_output.data['cid']
    nan_cells = help_output.data['idx_nan']
    years = help_output.data['years']

    zone_yearly_precip = np.sum(help_output.data['precip'], axis=(0, 2))
    zone_yearly_perco = np.sum(help_output.data['perco'], axis=(0, 2))
    zone_yearly_rechg = np.sum(help_output.data['rechg'], axis=(0, 2))
    zone_yearly_runoff = np.sum(help_output.data['runoff'], axis=(0, 2))
    zone_yearly_evapo = np.sum(help_output.data['evapo'], axis=(0, 2))
    zone_yearly_subrun1 = np.sum(help_output.data['subrun1'], axis=(0, 2))
    zone_yearly_subrun2 = np.sum(help_output.data['subrun2'], axis=(0, 2))

    # Calculate the yearly values for the surface water cells.
    surf_cellnames = list(surf_output.keys())
    for i, cellname in enumerate(surf_cellnames):
        cellname = str(int(cellname))
        data = surf_output[cellname]
        zone_yearly_precip += np.array(data['rain'])
        zone_yearly_runoff += np.array(data['runoff'])
        zone_yearly_evapo += np.array(data['evapo'])

    # Scale the results with the total numer of cells so we get mm/year.
    Ntot = len(help_cellnames) + len(surf_cellnames) - len(nan_cells)

    zone_yearly_precip /= Ntot
    zone_yearly_perco /= Ntot
    zone_yearly_rechg /= Ntot
    zone_yearly_runoff /= Ntot
    zone_yearly_evapo /= Ntot
    zone_yearly_subrun1 /= Ntot
    zone_yearly_subrun2 /= Ntot

    # Calcul simulated yearly total and base river flow.
    yearly_qflow = pd.DataFrame(index=years, columns=['qtot_sim'])
    yearly_qflow.index.name = 'years'

    yearly_qflow['qtot_sim'] = (
        zone_yearly_runoff +
        zone_yearly_subrun1 +
        zone_yearly_subrun2 +
        zone_yearly_rechg)

    yearly_qflow['qbase_sim'] = (
        zone_yearly_rechg +
        zone_yearly_subrun1 +
        zone_yearly_subrun2)

    return yearly_qflow


def plot_sim_vs_obs_yearly_streamflow(sim_qflow, obs_qflow, fig_title,
                                      figname=None):
    """
    Plot simulated vs observed yearly total and base streamflow.

    Parameters
    ----------
    figname : str, optional
        The abolute path of the file where to save the figure to disk.
        Note that the format of the file is inferred from the extension of
        "figname".
    """
    fwidth, fheight = 9, 5.5
    fig, ax = plt.subplots()
    fig.set_size_inches(fwidth, fheight)

    left_margin = 1.5/fwidth
    right_margin = 0.25/fwidth
    top_margin = 0.5/fheight
    bot_margin = 0.7/fheight
    ax.set_position([left_margin, bot_margin,
                     1 - left_margin - right_margin,
                     1 - top_margin - bot_margin])

    # Streamflow total
    l1, = ax.plot(obs_qflow['qtot_obs'], marker=None, mec='white',
                  lw=2, linestyle='--', clip_on=True, color='#3690c0')

    l2, = ax.plot(sim_qflow['qtot_sim'], marker=None, mec='white',
                  clip_on=True, lw=2, linestyle='-', color='#034e7b')

    # Streamflow base
    l3, = ax.plot(obs_qflow['qbase_obs'], marker=None, mec='white',
                  lw=2, linestyle='--', clip_on=True, color='#ef6548')

    l4, = ax.plot(sim_qflow['qbase_sim'], marker=None, mec='white',
                  clip_on=True, lw=2, linestyle='-', color='#990000')

    ax.tick_params(axis='both', direction='out', labelsize=12)
    ax.set_ylabel('Débit par unité de surface\n(mm/an)',
                  fontsize=16, labelpad=10)
    ax.set_xlabel('Années', fontsize=16, labelpad=10)
    ax.axis(ymin=0, ymax=1600)
    ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
    ax.axis([1960, 2017, 0, 1200])

    lines = [l1, l2, l3, l4]
    labels = ["Débit total observé", "Débit total simulé",
              "Débit base observé", "Débit base simulé"]
    legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                       borderaxespad=0, loc='upper left', borderpad=0.5,
                       bbox_to_anchor=(0, 1), ncol=2)
    legend.draw_frame(False)

    # Add a graph title.
    offset = transforms.ScaledTranslation(0/72, 12/72, fig.dpi_scale_trans)
    ax.text(0.5, 1, fig_title, fontsize=16, ha='center', va='bottom',
            transform=ax.transAxes+offset)

    # Save figure to file.
    if figname is not None:
        fig.savefig(figname)

    return fig


def plot_streamflow_scatter(sim_qflow, obs_qflow, fig_title, figname=None):
    """
    Create a scatter plot comparing simulated total and base
    streamflow yearly values with observed values.

    Parameters
    ----------
    figname : str, optional
        The abolute path of the file where to save the figure to disk.
        Note that the format of the file is inferred from the extension of
        "figname".
    """
    # Join both dataframe in a single dataframe.
    yearly_qflow = sim_qflow.copy()
    yearly_qflow = yearly_qflow.reindex(index=obs_qflow.index)
    yearly_qflow.loc[obs_qflow.index, 'qtot_obs'] = obs_qflow['qtot_obs']
    yearly_qflow.loc[obs_qflow.index, 'qbase_obs'] = obs_qflow['qbase_obs']

    fwidth, fheight = 5, 5
    fig, ax = plt.subplots()
    fig.set_size_inches(fwidth, fheight)

    left_margin = 1/fwidth
    right_margin = 0.5/fwidth
    top_margin = 0.5/fheight
    bot_margin = 1/fheight
    ax.set_position([left_margin, bot_margin,
                     1 - left_margin - right_margin,
                     1 - top_margin - bot_margin])

    xymin, xymax = 100, 1100
    ax.axis([xymin, xymax, xymin, xymax])
    ax.set_ylabel('Débits Simulés (mm/an)', fontsize=16, labelpad=15)
    ax.set_xlabel('Débits Observés (mm/an)', fontsize=16, labelpad=15)

    ax.tick_params(axis='both', direction='out', labelsize=12)

    l1, = ax.plot([xymin, xymax], [xymin, xymax], '--', color='black', lw=1)
    l2, = ax.plot(yearly_qflow['qtot_obs'],
                  yearly_qflow['qtot_sim'],
                  '.', color='#3690c0')
    l3, = ax.plot(yearly_qflow['qbase_obs'],
                  yearly_qflow['qbase_sim'],
                  '.', color='#990000')

    # Plot the model fit stats.
    rmse_qtot = np.nanmean(
        (yearly_qflow['qtot_sim'] - yearly_qflow['qtot_obs'])**2)**0.5
    me_qtot = np.nanmean(
        yearly_qflow['qtot_sim'] - yearly_qflow['qtot_obs'])

    rmse_qbase = np.nanmean(
        (yearly_qflow['qbase_sim'] - yearly_qflow['qbase_obs'])**2)**0.5
    me_qbase = np.nanmean(
        yearly_qflow['qbase_sim'] - yearly_qflow['qbase_obs'])

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
    ax.text(0.5, 1, fig_title, fontsize=16, ha='center', va='bottom',
            transform=ax.transAxes+offset)

    # Add a legend.
    lines = [l1, l2, l3]
    labels = ["1:1", "Débit total", "Débit base"]
    legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                       borderaxespad=0, loc='lower right', borderpad=0.5,
                       bbox_to_anchor=(1, 0), ncol=1)
    legend.draw_frame(False)

    # Save figure to file.
    if figname is not None:
        fig.savefig(figname)

    return fig
