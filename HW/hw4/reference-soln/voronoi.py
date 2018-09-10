#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# create voronoi diagram from 
# a set of points
##################################

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def voronoi_diagram():
    """create voronoi diagram of points in file"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    filename = sys.argv[1]
    points = np.loadtxt(filename)
    name = filename[:filename.rfind('.')]
    # 1. scipy's voronoi
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes()
    vor = Voronoi(points)
    voronoi_plot_2d(vor=vor, ax=ax, line_colors='red')
    ax.plot(vor.points[:, 0], vor.points[:, 1], 'bo')
    plt.savefig(name+"_voronoi.png") #, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    voronoi_diagram()
