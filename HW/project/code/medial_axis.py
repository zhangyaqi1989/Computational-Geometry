#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
compute medial axis of a simple polygon
"""

# standard library
import sys

# third party library
import numpy as np
import matplotlib.pyplot as plt

# local modules
sys.path.insert(0, '../../lib/')
from polygon_medial_axis import compute_polygon_medial_axis,\
        plot_polygon_medial_axis


def main():
    if len(sys.argv) != 2:
        print('Usage: >> python {} <polygon_file>'.format(sys.argv[0]))
        sys.exit(1)
    polygon = np.loadtxt(sys.argv[1])
    fig, ax = plt.subplots(figsize=(8, 8))
    medial_axis = compute_polygon_medial_axis(polygon, h=0.1)
    plot_polygon_medial_axis(polygon, medial_axis, ax=ax)
    ax.axis('equal')
    ax.set_title('Delaunay Triangulation')
    plt.show()

if __name__ == "__main__":
    main()
