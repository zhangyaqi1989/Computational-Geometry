#!/usr/bin/env python3
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
compute delaunay triangulation of a set of 2D points
"""

import sys
sys.path.insert(0, '../lib/')
import numpy as np
from convex_hull import convex_hull
import matplotlib.pyplot as plt
from triangulation_DS import Triangulation


def _is_left(A, B, point, tol=1e-16):
    """check if a point is on the left of AB or not"""
    x0, y0 = A
    x1, y1 = B
    x, y = point
    v = (x1 - x0) * (y - y0) - (y1 - y0) * (x - x0)
    if abs(v) < tol:
        return 0
    elif v > 0:
        return 1
    else:
        return -1


def hw3_2():
    """script for hw3-2()"""
    if len(sys.argv) != 2:
        print("Usage: python {} <filename>".format(sys.argv[0]))
        sys.exit(1)
    points = np.loadtxt(sys.argv[1])
    tri = Triangulation(points)
    tri.create_delaunay()
    # tri.plot_naive_triangulation()
    tri.plot_delaunay()
    hull = convex_hull(points, alg=2)
    M = len(hull)
    n = points.shape[0]
    k = 0 # compute how many points are on CH
    for p in range(n):
        point = points[p, :]
        for i in range(M):
            j = (i + 1) % M
            if _is_left(hull[i], hull[j], (point[0], point[1])) == 0:
                k += 1
                break
    t = tri.get_nfaces()
    e = tri.get_nedges()
    print("number of triangles = {}".format(t))
    print("number of edges = {}".format(e))
    print("number of vertices = {}".format(n))
    print("Checking topology")
    print("(n: number of vertices; t: number of triangles; e: number of edges)")
    print("(t = 2n - 2 - k; e = 3n - 3 - k)")
    print("{} {} 2 * {} - 2 - {}".format(t, "==" if t == 2*n-2-k else "!=", n, k))
    print("{} {} 3 * {} - 3 - {}".format(e, "==" if e == 3*n-3-k else "!=", n, k))
    plt.axis("equal")
    plt.show()


if __name__ == '__main__':
    hw3_2()
