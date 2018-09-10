#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
compute voronoi diagram for a set of 2D points
"""

import sys
sys.path.insert(0, '../lib/')
import numpy as np
import matplotlib.pyplot as plt
from triangulation_DS import Triangulation
# checking correctness
from scipy.spatial import Voronoi, voronoi_plot_2d


def voronoi_diagram():
    """create voronoi diagram of points in file"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    # 1. our voronoi
    tri = Triangulation(points)
    tri.compute_voronoi()
    figure = plt.figure(figsize=(8, 8))
    plt.subplot(1, 1, 1)
    tri.plot_voronoi()
    # 2. scipy's voronoi
    # vor = Voronoi(tri.vertices)
    # voronoi_plot_2d(vor)
    plt.show()


if __name__ == '__main__':
    voronoi_diagram()
