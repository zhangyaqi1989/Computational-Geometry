#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
# crust algorithm
##################################

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, Delaunay


def crust():
    """apply crust algorithm to a set of points defined in file"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    filename = sys.argv[1]
    points = np.loadtxt(filename)
    name = filename[:filename.rfind('.')]
    points_set = set([tuple(row) for row in points])
    tri1 = Delaunay(points)
    vor = Voronoi(points)
    new_points = np.vstack((points, vor.vertices))
    tri2 = Delaunay(new_points)
    # plot the crust of the shape
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes()
    n_triangles = tri2.simplices.shape[0]
    for i in range(n_triangles):
        triangle = tri2.points[tri2.simplices[i], :]
        for j in range(3):
            k = (j + 1) % 3
            A, B = triangle[j], triangle[k]
            if tuple(A) in points_set and tuple(B) in points_set:
                ax.plot([A[0], B[0]], [A[1], B[1]], 'b-')
    ax.plot(points[:, 0], points[:, 1], 'ro')
    # plt.savefig(name + "_crust.png") #, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    crust()
