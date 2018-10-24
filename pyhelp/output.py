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
        the study area for each month of the year.
        """
        Np = len(self.data['cid']) - len(self.data['idx_nan'])
        keys = ['precip', 'runoff', 'evapo', 'perco',
                'subrun1', 'subrun2', 'rechg']
        if self.data is None:
            return {key: np.zeros((1, 12)) for key in keys}
        else:
            return {key: np.nansum(self.data[key], axis=0)/Np for key in keys}

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


if __name__ == "__main__":
    output_fpath = "C:/Users/User/pyhelp/example/help_example.out"
    hout = HelpOutput(output_fpath)
    data = hout.data
    grid = hout.grid
    hout.plot_avg_monthly_budget()
