# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/jnsebgosselin/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library imports
import os.path as osp
from collections.abc import Mapping


# ---- Third party imports
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import pandas as pd
import numpy as np
import h5py


class HelpOutput(Mapping):
    def __init__(self, path_or_dict):
        super(HelpOutput, self).__init__()
        if isinstance(path_or_dict, dict):
            self.data = path_or_dict['data']
            self.grid = path_or_dict['grid']
        elif isinstance(path_or_dict, str) and osp.exists(path_or_dict):
            # Load the data from an HDF5 file saved on disk.
            hdf5 = h5py.File(path_or_dict, mode='r+')
            self.data = {}
            for key in list(hdf5['data'].keys()):
                if key == 'cid':
                    self.data[key] = hdf5['data'][key].value.tolist()
                else:
                    self.data[key] = hdf5['data'][key].value
            hdf5.close()

            # Load the grid from an HDF5 file saved on disk.
            self.grid = pd.read_hdf(path_or_dict, 'grid')
        else:
            self.data = None
            self.grid = None

    def __getitem__(self):
        pass

    def __iter__(self):
        pass

    def __len__(self):
        return len(self.data['cid'])

    def save_to_hdf5(self, path_to_hdf5):
        """Save the data and grid to a HDF5 file at the specified location."""
        print("Saving data to {}...".format(osp.basename(path_to_hdf5)),
              end=" ")

        # Save the data.
        hdf5file = h5py.File(path_to_hdf5, mode='w')
        datagrp = hdf5file.create_group('data')
        for key in list(self.data.keys()):
            if key == 'cid':
                # This is required to avoid a "TypeError: No conversion path
                # for dtype: dtype('<U5')".
                # See https://github.com/h5py/h5py/issues/289
                datagrp.create_dataset(
                    key, data=[np.string_(i) for i in self.data['cid']])
            else:
                datagrp.create_dataset(key, data=self.data[key])
        hdf5file.close()

        # Save the grid.
        self.grid.to_hdf(path_to_hdf5, key='grid', mode='a')

        print('done')

    def calc_area_monthly_avg(self):
        """
        Calcul water budget average monthly values for the study area.

        Return a dictionary that contains a 2D numpy array for each
        component of the water budget with average values calculated over
        the study area for each month of every year for which data is
        available.
        """
        Np = len(self.data['cid']) - len(self.data['idx_nan'])
        keys = ['precip', 'runoff', 'evapo', 'perco',
                'subrun1', 'subrun2', 'rechg']
        if self.data is None:
            return {key: np.zeros((1, 12)) for key in keys}
        else:
            return {key: np.nansum(self.data[key], axis=0)/Np for key in keys}

    def calc_area_yearly_avg(self):
        """
        Calcul water budget average yearly values for the study area.

        Return a dictionary that contains a numpy array for each
        component of the water budget with average values calculated over
        the study area for every year for which data is available.
        """
        monthly_avg = self.calc_area_monthly_avg()
        keys = list(monthly_avg.keys())
        return {key: np.sum(monthly_avg[key], axis=1) for key in keys}

    def calc_cells_yearly_avg(self):
        """
        Calcul water budget average yearly values for each cell.

        Return a dictionary that contains a numpy array for each
        component of the water budget with average values calculated for
        each cell for which data are available.
        """
        keys = ['precip', 'runoff', 'evapo', 'perco',
                'subrun1', 'subrun2', 'rechg']
        return {key: np.mean(np.sum(self.data[key], axis=2), axis=1)
                for key in keys}

    def plot_area_monthly_avg(self, figname=None):
        """
        Plot water budget average monthly values for the study area.

        Plot the monthly average values calculated for each component of the
        water budget for the entire study area.
        """
        fwidth, fheight = 9, 6.5
        fig, ax = plt.subplots()
        fig.set_size_inches(fwidth, fheight)

        # Setup axe margins :
        left_margin = 1.5/fwidth
        right_margin = 0.25/fwidth
        top_margin = 1/fheight
        bot_margin = 0.7/fheight
        ax.set_position([left_margin, bot_margin,
                         1 - left_margin - right_margin,
                         1 - top_margin - bot_margin])

        avg_monthly = self.calc_area_monthly_avg()
        months = range(1, 13)
        l1, = ax.plot(months, np.mean(avg_monthly['precip'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)
        l2, = ax.plot(months, np.mean(avg_monthly['rechg'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)
        l3, = ax.plot(months, np.mean(avg_monthly['runoff'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)
        l4, = ax.plot(months, np.mean(avg_monthly['evapo'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)
        l5, = ax.plot(months, np.mean(avg_monthly['subrun1'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)
        l6, = ax.plot(months, np.mean(avg_monthly['subrun2'], axis=0),
                      marker='o', mec='white', clip_on=False, lw=2)

        ax.set_ylabel('Composantes du bilan hydrologique\n(mm/mois)',
                      fontsize=16, labelpad=10)
        ax.set_xlabel('Mois', fontsize=16, labelpad=10)
        ax.axis(ymin=-5, ymax=140)
        ax.grid(axis='both', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
        ax.set_xticks(months)
        ax.set_xticklabels(['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Jun', 'Jul',
                            'Aoû', 'Sep', 'Oct', 'Nov', 'Déc'])
        ax.tick_params(axis='both', direction='out', labelsize=12)

        lines = [l1, l2, l3, l4, l5, l6]
        labels = ["Précipitations totales", "Recharge au roc",
                  "Ruissellement de surface", "Évapotranspiration",
                  "Ruissellement hypodermique superficiel",
                  "Ruissellement hypodermique profond"]
        legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                           borderaxespad=0, loc='lower left', borderpad=0.5,
                           bbox_to_anchor=(0, 1), ncol=2)
        legend.draw_frame(False)

        # if figname is not None:
        #     fig.savefig('bilan_hydro_mensuel_%s.pdf' % figname_sufix)

    def plot_area_yearly_avg(self):
        """
        Plot water budget average yearly values for the study area.

        Plot the yearly average values calculated for each component of the
        water budget for the entire study area.
        """
        fwidth, fheight = 8, 6.5
        fig, ax = plt.subplots()
        fig.set_size_inches(fwidth, fheight)

        # Setup axe margins.
        left_margin = 1.5/fwidth
        right_margin = 0.25/fwidth
        top_margin = 0.5/fheight
        bot_margin = 0.25/fheight
        ax.set_position([left_margin, bot_margin,
                         1 - left_margin - right_margin,
                         1 - top_margin - bot_margin])

        avg_yearly = self.calc_area_yearly_avg()
        avg_yearly_precip = np.mean(avg_yearly['precip'])
        avg_yearly_rechg = np.mean(avg_yearly['rechg'])
        avg_yearly_runoff = np.mean(avg_yearly['runoff'])
        avg_yearly_evapo = np.mean(avg_yearly['evapo'])
        avg_yearly_subrun1 = np.mean(avg_yearly['subrun1'])
        avg_yearly_subrun2 = np.mean(avg_yearly['subrun2'])

        l1 = ax.bar(1, avg_yearly_precip, 0.85, align='center')
        l2 = ax.bar(2, avg_yearly_rechg, 0.85, align='center')
        l3 = ax.bar(3, avg_yearly_runoff, 0.85, align='center')
        l4 = ax.bar(4, avg_yearly_evapo, 0.85, align='center')
        l5 = ax.bar(5, avg_yearly_subrun1, 0.85, align='center')
        l6 = ax.bar(6, avg_yearly_subrun2, 0.85, align='center')

        ax.axis(ymin=0, ymax=1200, xmin=0, xmax=7)
        ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
        ax.set_axisbelow(True)

        ax.text(1, avg_yearly_precip + 10, "%d\nmm/an" % avg_yearly_precip,
                ha='center', va='bottom')
        ax.text(2, avg_yearly_rechg + 10, "%d\nmm/an" % avg_yearly_rechg,
                ha='center', va='bottom')
        ax.text(3, avg_yearly_runoff + 10, "%d\nmm/an" % avg_yearly_runoff,
                ha='center', va='bottom')
        ax.text(4, avg_yearly_evapo + 10, "%d\nmm/an" % avg_yearly_evapo,
                ha='center', va='bottom')
        ax.text(5, avg_yearly_subrun1 + 10, "%d\nmm/an" % avg_yearly_subrun1,
                ha='center', va='bottom')
        ax.text(6, avg_yearly_subrun2 + 10, "%d\nmm/an" % avg_yearly_subrun2,
                ha='center', va='bottom')

        ax.tick_params(axis='y', direction='out', labelsize=12)
        ax.tick_params(axis='x', direction='out', length=0)
        ax.set_ylabel('Composantes du bilan hydrologique\n(mm/an)',
                      fontsize=16, labelpad=10)
        ax.set_xticklabels([])

        lines = [l1, l2, l3, l4, l5, l6]
        labels = ["Précipitations totales", "Recharge au roc",
                  "Ruissellement de surface", "Évapotranspiration",
                  "Ruissellement hypodermique superficiel",
                  "Ruissellement hypodermique profond"]
        legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                           borderaxespad=0, loc='upper right', borderpad=0.5,
                           bbox_to_anchor=(1, 1), ncol=1)
        legend.draw_frame(False)

        # Add a graph title.
        offset = transforms.ScaledTranslation(0/72, 12/72, fig.dpi_scale_trans)
        # ax.text(0.5, 1, figname_sufix, fontsize=16, ha='center', va='bottom',
                # transform=ax.transAxes+offset)
        # fig.savefig("hist_bilan_hydro_moyen_annuel_%s.pdf" % figname_sufix)


if __name__ == "__main__":
    output_fpath = "C:/Users/User/pyhelp/example/help_example.out"
    hout = HelpOutput(output_fpath)
    data = hout.data
    grid = hout.grid
    hout.plot_area_monthly_avg()
    hout.plot_area_yearly_avg()
