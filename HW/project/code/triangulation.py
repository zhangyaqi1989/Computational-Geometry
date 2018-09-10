#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
implement constrained delaunay triangulation (CST)
of a simple polygon
"""

# standard library
import sys

# third party library
import numpy as np
import matplotlib.pyplot as plt

# local modules
sys.path.insert(0, '../../lib/')
from polygon_delaunay_triangulation import plot_triangulation, polygon_dt

def main():
    if len(sys.argv) != 2:
        print('Usage: >> python {} <polygon_file>'.format(sys.argv[0]))
        sys.exit(1)
    polygon = np.loadtxt(sys.argv[1])
    tri = polygon_dt(polygon)
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_triangulation(tri, ax=ax)
    ax.axis('equal')
    ax.set_title('Delaunay Triangulation')
    plt.show()


if __name__ == "__main__":
    main()
