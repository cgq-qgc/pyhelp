# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:21:49 2022
@author: Jean-Sébastien Gosselin
"""

from pyhelp.output import HelpOutput
import matplotlib.pyplot as plt

plt.close('all')
dirname = "D:/Projets/pyhelp/example/"

# On charge en mémoire les résultats PyHELP d'une simulation antérieure
# dont nous avions sauvegarder les résultats sur le disque.

# À noter que cela revient au même que d'utiliser l'output de la commande
# "output_help = helpm.calc_help_cells()" directement.
output_help = HelpOutput(dirname + "help_example.out")

# On peut ensuite exporter, dans un fichier CSV, les moyennes annuelles du
# bilan hydrologique.

# On peut également spécifier (disponible dans pyhelp 0.3) les années pour
# lesquelles ont souhaite calculer les moyennes annuelles du bilan. Par
# exemple, pour calculer les moyennes entre 2003 et 2009 inclusivement, on
# écrirait :
output_help.save_to_csv(
    dirname + "help_example_yearly_2003-2009.csv",
    year_from=2003,
    year_to=2007)

# On peut également produires des graphiques des moyennes annuelles et
# ensuelles du bilan hydrologique de même que les valeurs annuelles du
# bilan hydrologique à l'aide des 3 fonctions ci-dessous. Ici encore,
# Il est possible de spécifier la période pour laquelle ont veut calculer
# les moyennes ou afficher les valeurs annuelles.
output_help.plot_area_monthly_avg(year_from=2003, year_to=2007)
output_help.plot_area_yearly_avg(year_from=2005, year_to=2005)
output_help.plot_area_yearly_series(year_from=2004, year_to=2009)

# Les valeurs numériques derrière les graphiques produits ci-dessus
# peuvent être accédées avec les fonctions suivantes:
yearly_avg = output_help.calc_area_yearly_avg()
monthly_avg = output_help.calc_area_monthly_avg()

# Il est très facile d'exporter ces résultats dans des csv. Par exemple,
# pour exporter les valeurs mensuelles moyennes de recharge, nous
# ferions tout simplement:
monthly_avg['rechg'].to_csv(dirname + "monthly_rechg.csv")

# Similairement, pour exporter les valeurs moyennes annuelles
# d'évapotranspiration, on pourrait faire:
yearly_avg['evapo'].to_csv(dirname + "yearly_evapo.csv")
