# -*- coding: utf-8 -*-

# Copyright © PyHelp Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHelp.
# Licensed under the terms of the GNU General Public License.


# ---- Standard Library imports
import os
import os.path as osp
from collections.abc import Mapping


# ---- Third party imports
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import pandas as pd
from geopandas import GeoDataFrame
import numpy as np
import h5py
from scipy.stats import linregress

# ---- Local imports
from pyhelp.maps import produce_point_geometry


class HelpOutput(Mapping):
    """
    A container to read and post-process monthly water budget results produced
    with the :class:`~pyhelp.HelpManager` class.
    """

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
                    self.data[key] = hdf5['data'][key].value.astype(str)
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
        """Save the data and grid to an HDF5 file at the specified location."""
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

    def save_to_shp(self, path_to_shp):
        """
        Save the grid data and cell yearly average values for each component
        of the water budget to a shapefile.
        """
        print('\rInitialize the shapefile...', end=' ')
        point_geo = produce_point_geometry(
            self.grid['lat_dd'].values, self.grid['lon_dd'].values)
        crs = ("+proj=longlat +ellps=GRS80 +datum=NAD83 "
               "+towgs84=0,0,0,0,0,0,0 +no_defs")
        shp = GeoDataFrame(self.grid, crs=crs, geometry=point_geo)
        print('done')
        print('\rAdding results to the shapefile...', end=' ')
        yearly_avg = self.calc_cells_yearly_avg()
        for key, value in yearly_avg.items():
            shp.loc[self.data['cid'], key] = value
        print('done')
        print('\rSaving data to the shapefile...', end=' ')
        if not osp.exists(osp.dirname(path_to_shp)):
            os.makedirs(osp.dirname(path_to_shp))
        shp.to_file(path_to_shp, driver='ESRI Shapefile')
        print('done')

    # ---- Calcul

    def calc_area_monthly_avg(self):
        """
        Calcul the monthly values of the water budget in mm/month for the
        whole study area.

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
        Calcul the average yearly values of the water budget in mm/year
        for the whole study area.

        Return a dictionary that contains a numpy array for each
        component of the water budget with average values calculated over
        the study area for every year for which data is available.
        """
        monthly_avg = self.calc_area_monthly_avg()
        keys = list(monthly_avg.keys())
        return {key: np.sum(monthly_avg[key], axis=1) for key in keys}

    def calc_cells_yearly_avg(self):
        """
        Plot the yearly values of the water budget in mm/year for the whole
        study area.
        Calcul water budget average yearly values for each cell.

        Return a dictionary that contains a numpy array for each
        component of the water budget with average values calculated for
        each cell for which data are available.
        """
        keys = ['precip', 'runoff', 'evapo', 'perco',
                'subrun1', 'subrun2', 'rechg']
        return {key: np.mean(np.sum(self.data[key], axis=2), axis=1)
                for key in keys}

    # ---- Plots

    def _create_figure(self, fsize=None, margins=None):
        """
        Create and return a figure and an axe mpl object using the
        specified settings.
        """
        fig, ax = plt.subplots()

        # Setup figure size.
        if fsize is not None:
            fwidth, fheight = fsize
        else:
            fwidth, fheight = fig.get_size_inches()
        fig.set_size_inches(*fsize)

        # Setup axe margins.
        if margins is not None:
            left_margin = margins[0]/fwidth
            top_margin = margins[1]/fheight
            right_margin = margins[2]/fwidth
            bot_margin = margins[3]/fheight
        ax.set_position([left_margin, bot_margin,
                         1 - left_margin - right_margin,
                         1 - top_margin - bot_margin])

        return fig, ax

    def plot_area_monthly_avg(self, figname=None):
        """
        Plot the monthly average values of the water budget in mm/month
        for the whole study area.
        """
        fig, ax = self._create_figure(
            fsize=(9, 6.5), margins=(1.5, 1, 0.25, 0.7))

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

        if figname is not None:
            fig.savefig(figname)

    def plot_area_yearly_avg(self, figname=None):
        """
        Plot the average yearly values of the water budget in mm/year
        for the whole study area.
        """
        fig, ax = self._create_figure(
            fsize=(8, 6.5), margins=(1.5, 0.5, 0.25, 0.25))

        area_yearly_avg = self.calc_area_yearly_avg()
        avg_yearly_precip = np.mean(area_yearly_avg['precip'])
        avg_yearly_rechg = np.mean(area_yearly_avg['rechg'])
        avg_yearly_runoff = np.mean(area_yearly_avg['runoff'])
        avg_yearly_evapo = np.mean(area_yearly_avg['evapo'])
        avg_yearly_subrun1 = np.mean(area_yearly_avg['subrun1'])
        avg_yearly_subrun2 = np.mean(area_yearly_avg['subrun2'])

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

        if figname is not None:
            fig.savefig(figname)

        # Add a graph title.
        # offset = transforms.ScaledTranslation(0/72, 12/72, fig.dpi_scale_trans)
        # ax.text(0.5, 1, figname_sufix, fontsize=16, ha='center', va='bottom',
                # transform=ax.transAxes+offset)

    def plot_area_yearly_series(self, figname=None):
        """
        Plot the yearly values of the water budget in mm/year for the whole
        study area.
        """
        fig, ax = self._create_figure(
            fsize=(9, 6.5), margins=(1.5, 1, 0.25, 0.7))

        years = self.data['years']
        yearly_avg = self.calc_area_yearly_avg()

        # Precipitation
        precip = yearly_avg['precip']
        l1, = ax.plot(years, precip, marker='o', mec='white',
                      clip_on=False, lw=2, color='#1f77b4')
        slope, intercept, r_val, p_val, std_err = linregress(years, precip)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#1f77b4')

        # Recharge
        rechg = yearly_avg['rechg']
        l2, = ax.plot(years, rechg, marker='o', mec='white', clip_on=False,
                      lw=2, color='#ff7f0e')
        slope, intercept, r_val, p_val, std_err = linregress(years, rechg)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#ff7f0e')

        # Runoff
        runoff = yearly_avg['runoff']
        l3, = ax.plot(years, runoff, marker='o', mec='white', clip_on=False,
                      lw=2, color='#2ca02c')
        slope, intercept, r_val, p_val, std_err = linregress(years, runoff)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#2ca02c')

        # Evapotranspiration
        evapo = yearly_avg['evapo']
        l4, = ax.plot(years, evapo, marker='o', mec='white', clip_on=False,
                      lw=2, color='#d62728')
        slope, intercept, r_val, p_val, std_err = linregress(years, evapo)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#d62728')

        # Superficial subsurface runoff
        subrun1 = yearly_avg['subrun1']
        l5, = ax.plot(years, subrun1, marker='o', mec='white',
                      clip_on=False, lw=2, color='#9467bd')
        slope, intercept, r_val, p_val, std_err = linregress(years, subrun1)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#9467bd')

        # Superficial subsurface runoff
        subrun2 = yearly_avg['subrun2']
        l6, = ax.plot(years, subrun2, marker='o', mec='white',
                      clip_on=False, lw=2, color='#8c564b')
        slope, intercept, r_val, p_val, std_err = linregress(years, subrun2)
        ax.plot(years, years * slope + intercept, marker=None, mec='white',
                clip_on=False, lw=1, dashes=[5, 3], color='#8c564b')

        ax.tick_params(axis='both', direction='out', labelsize=12)
        ax.set_ylabel('Composantes du bilan hydrologique\n(mm/an)',
                      fontsize=16, labelpad=10)
        ax.set_xlabel('Années', fontsize=16, labelpad=10)
        ax.axis(ymin=0, ymax=1600)
        ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)

        lines = [l1, l2, l3, l4, l5, l6]
        labels = ["Précipitations totales", "Recharge au roc",
                  "Ruissellement de surface", "Évapotranspiration",
                  "Ruissellement hypodermique superficiel",
                  "Ruissellement hypodermique profond"]

        legend = ax.legend(lines, labels, numpoints=1, fontsize=12,
                           borderaxespad=0, loc='lower left', borderpad=0.5,
                           bbox_to_anchor=(0, 1), ncol=2)
        legend.draw_frame(False)

        if figname is not None:
            fig.savefig(figname)


if __name__ == "__main__":
    output_fpath = "C:/Users/User/pyhelp/example/help_example.out"
    hout = HelpOutput(output_fpath)
    data = hout.data
    grid = hout.grid
    hout.plot_area_monthly_avg()
    hout.plot_area_yearly_avg()
    hout.plot_area_yearly_series()
    cells_yearly_avg = hout.calc_cells_yearly_avg()

    shp_fname = "C:/Users/User/pyhelp/example/help_example.shp"
    hout.save_to_shp(shp_fname)
