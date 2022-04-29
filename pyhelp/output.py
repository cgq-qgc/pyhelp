# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright © PyHELP Project Contributors
# https://github.com/cgq-qgc/pyhelp
#
# This file is part of PyHELP.
# Licensed under the terms of the MIT License.
# -----------------------------------------------------------------------------

from __future__ import annotations

# ---- Standard Library imports
import os.path as osp

# ---- Third party imports
from matplotlib.figure import Figure
from matplotlib import transforms
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import h5py
from scipy.stats import linregress


VARNAMES = ['precip', 'rechg', 'runoff', 'evapo',
            'subrun1', 'subrun2', 'perco']
LABELS = {'precip': "Précipitations totales",
          'rechg': "Recharge au roc",
          'runoff': "Ruissellement de surface",
          'evapo': "Évapotranspiration",
          'subrun1': "Ruissellement hypodermique superficiel",
          'subrun2': "Ruissellement hypodermique profond"}
COLORS = {'precip': '#1f77b4',
          'rechg': '#ff7f0e',
          'runoff': '#2ca02c',
          'evapo': '#d62728',
          'subrun1': '#9467bd',
          'subrun2': '#8c564b'}


class HelpOutput(object):
    """
    A container to read and post-process monthly water budget results produced
    with the :class:`~pyhelp.HelpManager` class.
    """

    def __init__(self, path_or_data: str | dict):
        if isinstance(path_or_data, dict):
            self.data = path_or_data
        elif isinstance(path_or_data, str) and osp.exists(path_or_data):
            self.load_from_hdf5(path_or_data)
        else:
            self.data = None

    def load_from_hdf5(self, path_to_hdf5: str):
        """Read data and grid from an HDF5 file at the specified location."""
        print(f"Loading data and grid from {path_to_hdf5}")
        hdf5 = h5py.File(path_to_hdf5, mode='r+')
        try:
            # Load the data.
            self.data = {}
            for key in list(hdf5['data'].keys()):
                values = np.array(hdf5['data'][key])
                if key == 'cid':
                    values = values.astype(str)
                self.data[key] = values
        finally:
            hdf5.close()
        print("Data and grid loaded successfully.")

    def save_to_hdf5(self, path_to_hdf5: str):
        """Save the data and grid to an HDF5 file at the specified location."""
        print("Saving data to {}...".format(osp.basename(path_to_hdf5)))
        hdf5file = h5py.File(path_to_hdf5, mode='w')
        try:
            # Save the data.
            group = hdf5file.create_group('data')
            for key in list(self.data.keys()):
                if key == 'cid':
                    # See http://docs.h5py.org/en/latest/strings.html as to
                    # why this is necessary to do this in order to save a list
                    # of strings in a dataset with h5py.
                    group.create_dataset(
                        key,
                        data=self.data[key],
                        dtype=h5py.string_dtype())
                else:
                    group.create_dataset(key, data=self.data[key])
        finally:
            hdf5file.close()
        print("Data saved successfully.")

    def save_to_csv(self, path_to_csv: str,
                    year_from: int = -np.inf,
                    year_to: int = np.inf) -> None:
        """
        Save in a csv file the annual average values of the components of the
        water budget calculated for each cell of the grid.

        Parameters
        ----------
        year_from : int, optional
            Minimum year of the period over which the average annual values
            are calculated. The default is -np.inf.
        year_to : int, optional
            Maximum year of the period over which the average annual values
            are calculated. The default is np.inf.
        """
        print("Saving data to {}...".format(osp.basename(path_to_csv)))
        df = pd.DataFrame(index=self.data['cid'])
        df.index.name = 'cid'

        df['lat_dd'] = self.data['lat_dd']
        df['lon_dd'] = self.data['lon_dd']

        yearly_avg = self.calc_cells_yearly_avg(year_from, year_to)
        for key, value in yearly_avg.items():
            df[key] = value

        df.to_csv(path_to_csv, encoding='utf8')
        print("Data saved successfully.")

    # ---- Calcul
    def calc_area_monthly_avg(self):
        """
        Calcul the water budget average monthly values in mm/month for the
        whole study area.

        Returns
        -------
        dict
            A dictionary that contains a pandas dataframe for each
            component of the water budget with monthly average values
            calculated over the whole study area for each month (columns) and
            every year (index) for which data is available.
        """
        Np = len(self.data['cid']) - len(self.data['idx_nan'])

        monthly_avg = {}
        for varname in VARNAMES:
            if self.data is None:
                df = pd.DataFrame(
                    columns=list(range(1, 13)))
            else:
                df = pd.DataFrame(
                    data=np.nansum(self.data[varname], axis=0) / Np,
                    index=self.data['years'],
                    columns=list(range(1, 13)))
            df.columns.name = 'month'
            df.index.name = 'year'
            monthly_avg[varname] = df
        return monthly_avg

    def calc_area_yearly_avg(self):
        """
        Calcul the water budget average yearly values in mm/year for the
        whole study area.

        Returns
        -------
        dict
            A dictionary that contains pandas dataframe for each
            component of the water budget with with yearly average values
            calculated over the whole study area for each year (index) for
            which data is available.
        """
        monthly_avg = self.calc_area_monthly_avg()
        return {varname: monthly_avg[varname].sum(axis=1) for
                varname in VARNAMES}

    def calc_cells_yearly_avg(self, year_from: int = -np.inf,
                              year_to: int = np.inf) -> dict:
        """
        Calcul the water budget average yearly values for each cell.

        Parameters
        ----------
        year_from : int, optional
            Minimum year of the period over which the average annual values
            are calculated. The default is -np.inf.
        year_to : int, optional
            Maximum year of the period over which the average annual values
            are calculated. The default is np.inf.

        Returns
        -------
        dict
           A dictionary that contains, for each component of the water budget,
           a numpy array of the average yearly values calculated for each cell
           of the grid.
        """
        years_mask = (
            (self.data['years'] >= year_from) &
            (self.data['years'] <= year_to))

        yearly_avg = {}
        for varname in VARNAMES:
            vardata = self.data[varname][:, years_mask, :]
            yearly_avg[varname] = np.mean(np.sum(vardata, axis=2), axis=1)
        return yearly_avg

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
            fig.set_size_inches(*fsize)
        else:
            fwidth, fheight = fig.get_size_inches()

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

    def plot_area_monthly_avg(self, figname: str = None,
                              year_from: int = -np.inf,
                              year_to: int = np.inf,
                              fig_title: str = None) -> Figure:
        """
        Plot the monthly average values of the water budget in mm/month
        for the whole study area.

        Parameters
        ----------
        figname : str, optional
            The abolute path of the file where to save the figure to disk.
            Note that the format of the file is inferred from the extension of
            "figname".
        year_from : int, optional
            Minimum year of the period over which the average monthly values
            are calculated. The default is -np.inf.
        year_to : int, optional
            Maximum year of the period over which the average monthly values
            are calculated. The default is np.inf.
        fig_title : str, optional
            A title to add at the top of the figure.

        Returns
        -------
        Figure
            The matplotlib figure instance created by this method.
        """
        avg_monthly = self.calc_area_monthly_avg()

        fig, ax = self._create_figure(fsize=(9, 6.5))

        mask_years = (
            (self.data['years'] >= year_from) &
            (self.data['years'] <= year_to))
        months = list(range(1, 13))
        for varname in VARNAMES[:-1]:
            vardataf = avg_monthly[varname]
            ax.plot(months, vardataf.loc[mask_years, :].mean(axis=0),
                    marker='o', ms=5, mec='white', clip_on=False, lw=2,
                    label=LABELS[varname], color=COLORS[varname])

        ax.set_ylabel(
            'Composantes mensuelles moyennes\ndu bilan hydrologique (mm/mois)',
            fontsize=16, labelpad=10)
        ax.axis(ymin=-5)
        ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
        ax.set_xticks(months)

        # http://bdl.oqlf.gouv.qc.ca/bdl/gabarit_bdl.asp?id=3619
        ax.set_xticklabels(['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Jun', 'Jul',
                            'Aoû', 'Sep', 'Oct', 'Nov', 'Déc'])
        ax.tick_params(axis='both', direction='out', labelsize=12)

        ax.legend(
            numpoints=1, fontsize=12, frameon=False, borderaxespad=0,
            loc='lower left', borderpad=0.5, bbox_to_anchor=(0, 1), ncol=2)

        # Add the figure title if provided in argument.
        if fig_title is not None:
            fig.suptitle(fig_title, fontsize=16)

        # Add note about year considered for the hydrologic budget.
        masked_years = self.data['years'][mask_years]
        year_min = masked_years.min()
        year_max = masked_years.max()
        if year_min == year_max:
            text = f"Année considérée pour le bilan : {year_min:0.0f}"
        else:
            text = "Années considérées pour le bilan : "
            text += f"{year_min:0.0f} - {year_max:0.0f}"

        fig.canvas.draw()
        bbox_bottom, _ = ax.xaxis.get_ticklabel_extents(
            fig.canvas.get_renderer())
        y0 = ax.transAxes.inverted().transform(bbox_bottom)[0][1]
        offset = transforms.ScaledTranslation(0, -12/72, fig.dpi_scale_trans)
        ax.text(0, y0, text, transform=ax.transAxes + offset,
                va='top', ha='left')

        # We call tight_layout two times to make sure the layout is
        # adjusted correctly.
        fig.tight_layout()
        fig.tight_layout()

        # Save figure to file.
        if figname is not None:
            fig.savefig(figname)

        return fig

    def plot_area_yearly_avg(self, figname: str = None,
                             year_from: int = -np.inf,
                             year_to: int = np.inf,
                             fig_title: str = None) -> Figure:
        """
        Plot the yearly average values of the water budget in mm/year
        for the whole study area.

        Parameters
        ----------
        figname : str, optional
            The abolute path of the file where to save the figure to disk.
            Note that the format of the file is inferred from the extension of
            "figname".
        year_from : int, optional
            Minimum year of the period over which the average annual values
            are calculated. The default is -np.inf.
        year_to : int, optional
            Maximum year of the period over which the average annual values
            are calculated. The default is np.inf.
        fig_title : str, optional
            A title to add at the top of the figure.

        Returns
        -------
        Figure
            The matplotlib figure instance created by this method.
        """
        fig, ax = self._create_figure(fsize=(9, 5.5))

        text_offset = transforms.ScaledTranslation(
            0, 3/72, fig.dpi_scale_trans)

        area_yearly_avg = self.calc_area_yearly_avg()
        mask_years = (
            (self.data['years'] >= year_from) &
            (self.data['years'] <= year_to))

        x = 0
        text_handles = []
        for varname in VARNAMES[:-1]:
            x += 1

            vardataf = area_yearly_avg[varname]
            var_avg_yearly = vardataf.loc[mask_years].mean()

            ax.bar(x, var_avg_yearly, 0.85, align='center',
                   label=LABELS[varname], color=COLORS[varname])
            text_handles.append(
                ax.text(x, var_avg_yearly, "%d\nmm/an" % round(var_avg_yearly),
                        ha='center', va='bottom',
                        transform=ax.transData + text_offset))
        fig.canvas.draw()

        # Setup axis limits.
        ymax = 0
        for handle in text_handles:
            bbox = handle.get_window_extent(fig.canvas.get_renderer())
            bbox = bbox.transformed(ax.transData.inverted())
            ymax = max(ymax, bbox.y1)
        ymax = np.ceil(ymax * 1.05)
        ax.axis(ymin=0, ymax=ymax, xmin=0.25, xmax=6.75)

        ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)
        ax.set_axisbelow(True)

        ax.tick_params(axis='y', direction='out', labelsize=12)
        ax.tick_params(axis='x', direction='out', length=0)
        ax.set_ylabel(
            'Composantes annuelles moyennes\ndu bilan hydrologique (mm/an)',
            fontsize=16, labelpad=10)
        ax.set_xticklabels([])

        ax.legend(
            numpoints=1, fontsize=12, frameon=False, borderaxespad=0,
            loc='lower left', borderpad=0.5, bbox_to_anchor=(0, 1, 1, 1),
            ncol=2)

        # Add the figure title if provided in argument.
        if fig_title is not None:
            fig.suptitle(fig_title, fontsize=16)

        # Add note about year considered for the hydrologic budget.
        masked_years = self.data['years'][mask_years]
        year_min = masked_years.min()
        year_max = masked_years.max()
        if year_min == year_max:
            text = f"Année considérée pour le bilan : {year_min:0.0f}"
        else:
            text = "Années considérées pour le bilan : "
            text += f"{year_min:0.0f} - {year_max:0.0f}"

        fig.canvas.draw()
        bbox_bottom, _ = ax.xaxis.get_ticklabel_extents(
            fig.canvas.get_renderer())
        y0 = ax.transAxes.inverted().transform(bbox_bottom)[0][1]
        offset = transforms.ScaledTranslation(0, -6/72, fig.dpi_scale_trans)
        ax.text(0, y0, text, transform=ax.transAxes + offset,
                va='top', ha='left')

        # We call tight_layout two times to make sure the layout is
        # adjusted correctly.
        fig.tight_layout()
        fig.tight_layout()

        if figname is not None:
            fig.savefig(figname)

        return fig

    def plot_area_yearly_series(self, figname: str = None,
                                year_from: int = -np.inf,
                                year_to: int = np.inf,
                                fig_title: str = None) -> Figure:
        """
        Plot the yearly values of the water budget in mm/year calculated
        over the whole study area.

        Parameters
        ----------
        figname : str, optional
            The abolute path of the file where to save the figure to disk.
            Note that the format of the file is inferred from the extension of
            "figname".
        year_from : int, optional
            Minimum year of the period over which the annual values
            are plotted. The default is -np.inf.
        year_to : int, optional
            Maximum year of the period over which the annual values
            are calculated. The default is np.inf.
        fig_title : str, optional
            A title to add at the top of the figure.

        Returns
        -------
        Figure
            The matplotlib figure instance created by this method.
        """
        fig, ax = self._create_figure(
            fsize=(9, 6.5), margins=(1.5, 1, 0.25, 0.7))

        years = self.data['years']
        yearly_avg = self.calc_area_yearly_avg()
        mask_years = (years >= year_from) & (years <= year_to)

        for varname in VARNAMES[:-1]:
            masked_data = yearly_avg[varname].loc[mask_years]
            masked_years = masked_data.index.values.astype('int')

            ax.plot(masked_years, masked_data, marker='o', mec='white', ms=5,
                    clip_on=False, lw=2, color=COLORS[varname],
                    label=LABELS[varname])

            slope, intercept, r_val, p_val, std_err = linregress(
                masked_years, masked_data.values)

            ax.plot(masked_years, masked_years * slope + intercept,
                    marker=None, mec='white', clip_on=False, lw=1,
                    dashes=[5, 3], color=COLORS[varname])

        ax.tick_params(axis='both', direction='out', labelsize=12)
        ax.set_ylabel('Composantes annuelles\ndu bilan hydrologique (mm/an)',
                      fontsize=16, labelpad=10)
        ax.set_xlabel('Années', fontsize=16, labelpad=10)
        ax.axis(ymin=0)
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.grid(axis='y', color=[0.35, 0.35, 0.35], ls='-', lw=0.5)

        ax.legend(numpoints=1, fontsize=12, frameon=False,
                  borderaxespad=0, loc='lower left', borderpad=0.5,
                  bbox_to_anchor=(0, 1), ncol=2)

        # Add the figure title if provided in argument.
        if fig_title is not None:
            fig.suptitle(fig_title, fontsize=16)

        # We call tight_layout two times to make sure the layout is
        # adjusted correctly.
        fig.tight_layout()
        fig.tight_layout()

        if figname is not None:
            fig.savefig(figname)

        return fig


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
