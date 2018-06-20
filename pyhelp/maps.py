# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 11:18:39 2018
@author: User
"""

from itertools import product
from shapely.geometry import Point
import numpy as np


def produce_point_geometry(lat, lon):
    """Return a list of shapely points from lat, lon coordinates."""
    return [Point(x, y) for x, y in zip(lon, lat)]
