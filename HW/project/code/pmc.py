#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
read in polygon and points, do point membership check
and compute minimum distance
"""

# standard library
import sys
sys.path.insert(0, '../../lib/')

# third party library
import numpy as np
import matplotlib.pyplot as plt

# local modules
from marching_square import compute_point_to_polygon, point_membership_check

def main():
    if len(sys.argv) != 3:
        print("Usage: python {} <polygon_file> <points_file>".format(sys.argv[0]))
        sys.exit(1)
    polygon = np.loadtxt(sys.argv[1])
    points_query = np.loadtxt(sys.argv[2])
    for i in range(points_query.shape[0]):
        P = points_query[i]
        min_distance = compute_point_to_polygon(polygon, P)
        print('{}: {}, {:0.6f}'.format(P, point_membership_check(polygon, P), min_distance))
    fig, ax = plt.subplots(figsize=(8, 8))
    temp_polygon = np.vstack((polygon, polygon[0:1, :]))
    ax.plot(temp_polygon[:, 0], temp_polygon[:, 1], 'b-o')
    ax.plot(points_query[:, 0], points_query[:, 1], 'ro')
    ax.axis('equal')
    for i, (x, y) in enumerate(points_query):
        ax.annotate(i, (x, y))
    plt.show()


if __name__ == "__main__":
    main()
